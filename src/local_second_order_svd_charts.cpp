#include <geosmooth/local_second_order_svd_charts.h>
#include <geosmooth/local_pca_charts.h>

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace geosmooth {

namespace {

Eigen::RowVectorXd second_order_chart_center(
    const Eigen::MatrixXd& local,
    const Eigen::RowVectorXd& anchor,
    const std::string& center_mode
) {
    if (center_mode == "anchor") {
        return anchor;
    }
    if (center_mode == "mean") {
        return local.colwise().mean();
    }
    Rcpp::stop("'center.mode' must be either 'anchor' or 'mean'.");
}

void orient_second_order_basis_columns(Eigen::MatrixXd& basis) {
    for (int j = 0; j < basis.cols(); ++j) {
        Eigen::Index max_index = 0;
        basis.col(j).cwiseAbs().maxCoeff(&max_index);
        if (basis(max_index, j) < 0.0) {
            basis.col(j) *= -1.0;
        }
    }
}

bool matrix_is_finite(const Eigen::MatrixXd& X) {
    return X.array().isFinite().all();
}

bool vector_is_finite(const Eigen::VectorXd& x) {
    return x.array().isFinite().all();
}

bool row_vector_is_finite(const Eigen::RowVectorXd& x) {
    return x.array().isFinite().all();
}

Eigen::MatrixXd nan_matrix(int rows, int cols) {
    Eigen::MatrixXd out(rows, cols);
    out.setConstant(std::numeric_limits<double>::quiet_NaN());
    return out;
}

void fill_quadratic_monomials(
    int m,
    std::vector<int>& monomial_a,
    std::vector<int>& monomial_b,
    std::vector<double>& monomial_multiplier
) {
    monomial_a.clear();
    monomial_b.clear();
    monomial_multiplier.clear();
    monomial_a.reserve(m * (m + 1) / 2);
    monomial_b.reserve(m * (m + 1) / 2);
    monomial_multiplier.reserve(m * (m + 1) / 2);

    for (int a = 0; a < m; ++a) {
        monomial_a.push_back(a);
        monomial_b.push_back(a);
        monomial_multiplier.push_back(1.0);
    }
    for (int a = 0; a < m; ++a) {
        for (int b = a + 1; b < m; ++b) {
            monomial_a.push_back(a);
            monomial_b.push_back(b);
            monomial_multiplier.push_back(2.0);
        }
    }
}

Eigen::VectorXd sanitize_weights(const Eigen::VectorXd* weights, int k) {
    Eigen::VectorXd out(k);
    if (weights == nullptr) {
        out.setOnes();
        return out;
    }
    for (int i = 0; i < k; ++i) {
        const double wi = (*weights)(i);
        out(i) = (std::isfinite(wi) && wi > 0.0) ? wi : 0.0;
    }
    return out;
}

int positive_weight_count(const Eigen::VectorXd& weights) {
    int count = 0;
    for (int i = 0; i < weights.size(); ++i) {
        if (weights(i) > 0.0) {
            ++count;
        }
    }
    return count;
}

Eigen::MatrixXd apply_sqrt_weights(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& weights
) {
    Eigen::MatrixXd out = X;
    for (int i = 0; i < out.rows(); ++i) {
        out.row(i) *= std::sqrt(weights(i));
    }
    return out;
}

struct rank_result_t {
    int rank = 0;
    double cutoff = 0.0;
    double sigma_max = 0.0;
    double sigma_min_kept = 0.0;
    double condition = std::numeric_limits<double>::infinity();
};

rank_result_t relative_svd_rank(
    const Eigen::VectorXd& singular_values,
    int nrow,
    int ncol,
    double rank_tolerance,
    double rank_absolute_tolerance
) {
    rank_result_t out;
    if (singular_values.size() < 1) {
        return out;
    }

    out.sigma_max = singular_values.maxCoeff();
    if (!std::isfinite(out.sigma_max) ||
        out.sigma_max <= rank_absolute_tolerance) {
        out.rank = 0;
        out.cutoff = rank_absolute_tolerance;
        return out;
    }

    out.cutoff = rank_tolerance *
        static_cast<double>(std::max(nrow, ncol)) * out.sigma_max;
    out.sigma_min_kept = std::numeric_limits<double>::infinity();
    for (int i = 0; i < singular_values.size(); ++i) {
        if (singular_values(i) > out.cutoff) {
            ++out.rank;
            out.sigma_min_kept = std::min(out.sigma_min_kept, singular_values(i));
        }
    }
    if (out.rank > 0 &&
        std::isfinite(out.sigma_min_kept) && out.sigma_min_kept > 0.0) {
        out.condition = out.sigma_max / out.sigma_min_kept;
    }
    return out;
}

second_order_svd_chart_result_t base_result(
    const Eigen::MatrixXd& local,
    int requested_dim,
    const Eigen::VectorXd& weights
) {
    second_order_svd_chart_result_t out;
    const int k = local.rows();
    const int n = local.cols();
    fill_quadratic_monomials(
        requested_dim,
        out.monomial_a,
        out.monomial_b,
        out.monomial_multiplier
    );
    out.selected_dim = requested_dim;
    out.effective_support = positive_weight_count(weights);
    out.quadratic_ncol = requested_dim * (requested_dim + 1) / 2;
    out.coordinates = nan_matrix(k, requested_dim);
    out.basis = nan_matrix(n, requested_dim);
    out.preliminary_basis = nan_matrix(n, requested_dim);
    out.rho_tangent = Eigen::MatrixXd(0, 0);
    out.corrected_residual = Eigen::MatrixXd(0, 0);
    out.curvature_coefficients = Eigen::MatrixXd(0, 0);
    out.first_singular_values = Eigen::VectorXd(0);
    out.second_singular_values = Eigen::VectorXd(0);
    out.design_condition = std::numeric_limits<double>::infinity();
    out.status = "ok";
    return out;
}

bool plain_pca_fallback_feasible(
    int k,
    int n,
    int requested_dim,
    const Eigen::VectorXd* original_weights,
    int effective_support
) {
    if (k < 1 || n < 1 || requested_dim < 1) {
        return false;
    }
    if (requested_dim > std::min(k, n)) {
        return false;
    }
    if (original_weights != nullptr && effective_support < 1) {
        return false;
    }
    return true;
}

second_order_svd_chart_result_t fallback_result(
    const Eigen::MatrixXd& local,
    const Eigen::RowVectorXd& anchor,
    int requested_dim,
    const std::string& center_mode,
    const Eigen::VectorXd* original_weights,
    const Eigen::VectorXd& sanitized_weights,
    bool rebase_to_anchor,
    bool orient_basis,
    const std::string& primary_failure_reason,
    const second_order_svd_chart_result_t* current = nullptr
) {
    second_order_svd_chart_result_t out = current == nullptr ?
        base_result(local, requested_dim, sanitized_weights) : *current;
    out.fallback_used = true;
    out.primary_failure_reason = primary_failure_reason;
    out.plain_pca_fallback_feasible = plain_pca_fallback_feasible(
        local.rows(),
        local.cols(),
        requested_dim,
        original_weights,
        out.effective_support
    );

    if (!out.plain_pca_fallback_feasible) {
        out.coordinates = nan_matrix(local.rows(), requested_dim);
        out.basis = nan_matrix(local.cols(), requested_dim);
        out.preliminary_basis = nan_matrix(local.cols(), requested_dim);
        out.fallback_reason = "plain_pca_fallback_not_feasible";
        out.status = "structured.failure";
        return out;
    }

    local_pca_chart_result_t chart = compute_local_pca_chart(
        local,
        anchor,
        requested_dim,
        "fixed",
        0.9,
        center_mode,
        original_weights,
        rebase_to_anchor,
        orient_basis
    );
    out.coordinates = chart.coordinates;
    out.basis = chart.basis;
    if (out.preliminary_basis.rows() != local.cols() ||
        out.preliminary_basis.cols() != requested_dim ||
        !matrix_is_finite(out.preliminary_basis)) {
        out.preliminary_basis = chart.basis;
    }
    if (out.first_singular_values.size() == 0) {
        out.first_singular_values = chart.singular_values;
    }
    if (current == nullptr) {
        out.first_rank = std::min(
            requested_dim,
            static_cast<int>(chart.singular_values.size())
        );
    }
    out.fallback_reason = primary_failure_reason;
    out.status = "pca.fallback";
    return out;
}

Eigen::MatrixXd build_quadratic_design(
    const Eigen::MatrixXd& rho,
    const std::vector<int>& monomial_a,
    const std::vector<int>& monomial_b,
    const std::vector<double>& monomial_multiplier
) {
    Eigen::MatrixXd A(rho.rows(), monomial_a.size());
    for (int j = 0; j < A.cols(); ++j) {
        const int a = monomial_a[j];
        const int b = monomial_b[j];
        const double mult = monomial_multiplier[j];
        A.col(j) = mult * rho.col(a).cwiseProduct(rho.col(b));
    }
    return A;
}

} // namespace

second_order_svd_chart_result_t compute_local_second_order_svd_chart(
    const Eigen::MatrixXd& local,
    const Eigen::RowVectorXd& anchor,
    int requested_dim,
    const std::string& center_mode,
    const Eigen::VectorXd* weights,
    double rank_tolerance,
    double rank_absolute_tolerance,
    double curvature_condition_max,
    double curvature_ridge,
    int min_curvature_support,
    bool rebase_to_anchor,
    bool orient_basis
) {
    if (local.rows() < 1 || local.cols() < 1) {
        Rcpp::stop("'X_support' must have positive dimensions.");
    }
    if (anchor.size() != local.cols()) {
        Rcpp::stop("'center' must have ncol(X_support) entries.");
    }
    if (!matrix_is_finite(local)) {
        Rcpp::stop("'X_support' must contain only finite values.");
    }
    if (!row_vector_is_finite(anchor)) {
        Rcpp::stop("'center' must contain only finite values.");
    }
    if (requested_dim < 1) {
        Rcpp::stop("'chart.dim' must be at least 1.");
    }
    if (requested_dim > local.cols()) {
        Rcpp::stop("'chart.dim' must be no larger than ncol(X_support).");
    }
    if (center_mode != "anchor" && center_mode != "mean") {
        Rcpp::stop("'center.mode' must be either 'anchor' or 'mean'.");
    }
    if (weights != nullptr && weights->size() != local.rows()) {
        Rcpp::stop("'weights' must have nrow(X_support) entries.");
    }
    if (rank_tolerance < 0.0 || !std::isfinite(rank_tolerance)) {
        Rcpp::stop("'rank.tolerance' must be a finite nonnegative scalar.");
    }
    if (rank_absolute_tolerance < 0.0 ||
        !std::isfinite(rank_absolute_tolerance)) {
        Rcpp::stop("'rank.absolute.tolerance' must be a finite nonnegative scalar.");
    }
    if (curvature_condition_max <= 0.0 ||
        std::isnan(curvature_condition_max)) {
        Rcpp::stop("'curvature.condition.max' must be positive.");
    }
    if (curvature_ridge != 0.0 || !std::isfinite(curvature_ridge)) {
        Rcpp::stop("'curvature.ridge' must be 0 for the H1 prototype.");
    }
    if (min_curvature_support < 0) {
        Rcpp::stop("'min.curvature.support' must be nonnegative.");
    }
    const Eigen::VectorXd sanitized_weights =
        sanitize_weights(weights, local.rows());

    if (requested_dim == local.cols()) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "chart_dim_not_less_than_ambient_dim"
        );
    }

    second_order_svd_chart_result_t out =
        base_result(local, requested_dim, sanitized_weights);

    if (out.effective_support < requested_dim + 1) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "too_few_effective_support",
            &out
        );
    }

    const Eigen::RowVectorXd svd_center =
        second_order_chart_center(local, anchor, center_mode);
    const Eigen::MatrixXd centered = local.rowwise() - svd_center;
    const Eigen::MatrixXd weighted_centered =
        apply_sqrt_weights(centered, sanitized_weights);

    Eigen::JacobiSVD<Eigen::MatrixXd> first_svd(
        weighted_centered,
        Eigen::ComputeThinV
    );
    out.first_singular_values = first_svd.singularValues();
    if (!vector_is_finite(out.first_singular_values) ||
        !matrix_is_finite(first_svd.matrixV())) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "first_svd_rank_deficient",
            &out
        );
    }
    const rank_result_t first_rank = relative_svd_rank(
        out.first_singular_values,
        weighted_centered.rows(),
        weighted_centered.cols(),
        rank_tolerance,
        rank_absolute_tolerance
    );
    out.first_rank = first_rank.rank;
    if (first_rank.rank < requested_dim ||
        first_svd.matrixV().cols() < requested_dim) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "first_svd_rank_deficient",
            &out
        );
    }

    const Eigen::MatrixXd preliminary_basis =
        first_svd.matrixV().leftCols(requested_dim);
    out.preliminary_basis = preliminary_basis;
    out.rho_tangent = centered * preliminary_basis;

    const Eigen::MatrixXd A = build_quadratic_design(
        out.rho_tangent,
        out.monomial_a,
        out.monomial_b,
        out.monomial_multiplier
    );
    const int required_curvature_support = std::max(
        std::max(requested_dim + 1, out.quadratic_ncol + 1),
        min_curvature_support
    );
    if (out.effective_support < required_curvature_support) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "curvature_under_determined",
            &out
        );
    }

    const Eigen::MatrixXd A_weighted =
        apply_sqrt_weights(A, sanitized_weights);
    Eigen::JacobiSVD<Eigen::MatrixXd> curvature_svd(
        A_weighted,
        Eigen::ComputeThinU | Eigen::ComputeThinV
    );
    const Eigen::VectorXd curvature_sing = curvature_svd.singularValues();
    if (!vector_is_finite(curvature_sing) ||
        !matrix_is_finite(curvature_svd.matrixU()) ||
        !matrix_is_finite(curvature_svd.matrixV())) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "curvature_solve_failure",
            &out
        );
    }

    const rank_result_t curvature_rank = relative_svd_rank(
        curvature_sing,
        A_weighted.rows(),
        A_weighted.cols(),
        rank_tolerance,
        rank_absolute_tolerance
    );
    out.design_rank = curvature_rank.rank;
    out.design_condition = curvature_rank.condition;
    if (curvature_rank.rank < out.quadratic_ncol) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "curvature_rank_deficient",
            &out
        );
    }
    if (!std::isfinite(curvature_rank.condition) ||
        curvature_rank.condition > curvature_condition_max) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "curvature_ill_conditioned",
            &out
        );
    }

    Eigen::VectorXd inv = Eigen::VectorXd::Zero(curvature_sing.size());
    for (int i = 0; i < curvature_sing.size(); ++i) {
        if (curvature_sing(i) > curvature_rank.cutoff) {
            inv(i) = 1.0 / curvature_sing(i);
        }
    }
    const Eigen::MatrixXd target_weighted =
        apply_sqrt_weights(2.0 * centered, sanitized_weights);
    out.curvature_coefficients =
        curvature_svd.matrixV() * inv.asDiagonal() *
        curvature_svd.matrixU().transpose() * target_weighted;
    if (!matrix_is_finite(out.curvature_coefficients)) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "curvature_solve_failure",
            &out
        );
    }
    out.fit_method = "svd";

    const Eigen::MatrixXd fitted = A * out.curvature_coefficients;
    const Eigen::MatrixXd fitted_weighted_residual =
        apply_sqrt_weights(fitted - 2.0 * centered, sanitized_weights);
    out.fit_residual_frobenius = fitted_weighted_residual.norm();
    const Eigen::MatrixXd curvature_fitted = 0.5 * fitted;
    out.curvature_fitted_frobenius = curvature_fitted.norm();
    out.corrected_residual = centered - curvature_fitted;
    out.corrected_residual_frobenius = out.corrected_residual.norm();
    if (!matrix_is_finite(out.corrected_residual) ||
        !std::isfinite(out.fit_residual_frobenius) ||
        !std::isfinite(out.curvature_fitted_frobenius) ||
        !std::isfinite(out.corrected_residual_frobenius)) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "curvature_solve_failure",
            &out
        );
    }

    const Eigen::MatrixXd weighted_corrected =
        apply_sqrt_weights(out.corrected_residual, sanitized_weights);
    Eigen::JacobiSVD<Eigen::MatrixXd> second_svd(
        weighted_corrected,
        Eigen::ComputeThinV
    );
    out.second_singular_values = second_svd.singularValues();
    if (!vector_is_finite(out.second_singular_values) ||
        !matrix_is_finite(second_svd.matrixV())) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "second_svd_failure",
            &out
        );
    }
    const rank_result_t second_rank = relative_svd_rank(
        out.second_singular_values,
        weighted_corrected.rows(),
        weighted_corrected.cols(),
        rank_tolerance,
        rank_absolute_tolerance
    );
    out.second_rank = second_rank.rank;
    if (second_rank.rank < requested_dim ||
        second_svd.matrixV().cols() < requested_dim) {
        return fallback_result(
            local,
            anchor,
            requested_dim,
            center_mode,
            weights,
            sanitized_weights,
            rebase_to_anchor,
            orient_basis,
            "second_svd_rank_deficient",
            &out
        );
    }

    out.basis = second_svd.matrixV().leftCols(requested_dim);
    if (orient_basis) {
        orient_second_order_basis_columns(out.basis);
    }
    const Eigen::RowVectorXd coord_origin =
        rebase_to_anchor ? anchor : svd_center;
    out.coordinates = (local.rowwise() - coord_origin) * out.basis;
    out.fallback_used = false;
    out.fallback_reason = "none";
    out.primary_failure_reason = "none";
    out.plain_pca_fallback_feasible = false;
    out.status = "ok";
    return out;
}

} // namespace geosmooth
