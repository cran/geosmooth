#include <geosmooth/local_second_order_svd_charts.h>

#include <Rcpp.h>

using Rcpp::DataFrame;
using Rcpp::List;
using Rcpp::Nullable;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

Eigen::MatrixXd numeric_matrix_to_eigen_second_order(const NumericMatrix& X) {
    Eigen::MatrixXd out(X.nrow(), X.ncol());
    for (int j = 0; j < X.ncol(); ++j) {
        for (int i = 0; i < X.nrow(); ++i) {
            out(i, j) = X(i, j);
        }
    }
    return out;
}

NumericMatrix eigen_matrix_to_numeric_second_order(const Eigen::MatrixXd& X) {
    NumericMatrix out(X.rows(), X.cols());
    for (int j = 0; j < X.cols(); ++j) {
        for (int i = 0; i < X.rows(); ++i) {
            out(i, j) = X(i, j);
        }
    }
    return out;
}

NumericVector eigen_vector_to_numeric_second_order(const Eigen::VectorXd& x) {
    NumericVector out(x.size());
    for (int i = 0; i < x.size(); ++i) out[i] = x(i);
    return out;
}

DataFrame monomials_to_frame(
    const geosmooth::second_order_svd_chart_result_t& chart
) {
    const int q = static_cast<int>(chart.monomial_a.size());
    Rcpp::IntegerVector component(q);
    Rcpp::IntegerVector a(q);
    Rcpp::IntegerVector b(q);
    Rcpp::NumericVector multiplier(q);
    for (int j = 0; j < q; ++j) {
        component[j] = j + 1;
        a[j] = chart.monomial_a[j] + 1;
        b[j] = chart.monomial_b[j] + 1;
        multiplier[j] = chart.monomial_multiplier[j];
    }
    return DataFrame::create(
        Rcpp::Named("component") = component,
        Rcpp::Named("a") = a,
        Rcpp::Named("b") = b,
        Rcpp::Named("multiplier") = multiplier
    );
}

DataFrame diagnostics_to_frame(
    const geosmooth::second_order_svd_chart_result_t& chart
) {
    return DataFrame::create(
        Rcpp::Named("effective.support") = chart.effective_support,
        Rcpp::Named("quadratic.ncol") = chart.quadratic_ncol,
        Rcpp::Named("design.rank") = chart.design_rank,
        Rcpp::Named("design.condition") = chart.design_condition,
        Rcpp::Named("fit.method") = chart.fit_method,
        Rcpp::Named("ridge.lambda") = chart.ridge_lambda,
        Rcpp::Named("fit.residual.frobenius") =
            chart.fit_residual_frobenius,
        Rcpp::Named("curvature.fitted.frobenius") =
            chart.curvature_fitted_frobenius,
        Rcpp::Named("corrected.residual.frobenius") =
            chart.corrected_residual_frobenius,
        Rcpp::Named("first.rank") = chart.first_rank,
        Rcpp::Named("second.rank") = chart.second_rank,
        Rcpp::Named("plain.pca.fallback.feasible") =
            chart.plain_pca_fallback_feasible,
        Rcpp::Named("primary.failure.reason") =
            chart.primary_failure_reason,
        Rcpp::Named("status") = chart.status
    );
}

} // namespace

//' Second-order local SVD chart
//'
//' Internal experimental C++ backend for Harlim-style second-order local SVD
//' chart construction.  This primitive is intentionally separate from the
//' shared plain local-PCA chart backend.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List rcpp_local_second_order_svd_chart(
        const NumericMatrix& X_support,
        const NumericVector& center,
        const int chart_dim,
        const std::string& center_mode = "anchor",
        Nullable<NumericVector> weights = R_NilValue,
        const double rank_tolerance = 1.4901161193847656e-8,
        const double rank_absolute_tolerance = 0.0,
        const double curvature_condition_max = 1e8,
        const double curvature_ridge = 0.0,
        const int min_curvature_support = 0,
        const bool rebase_to_anchor = true,
        const bool orient_basis = false) {
    if (center.size() != X_support.ncol()) {
        Rcpp::stop("'center' must have ncol(X_support) entries.");
    }
    Eigen::MatrixXd local =
        numeric_matrix_to_eigen_second_order(X_support);
    Eigen::RowVectorXd anchor(center.size());
    for (int j = 0; j < center.size(); ++j) anchor(j) = center[j];

    Eigen::VectorXd weight_vec;
    Eigen::VectorXd* weight_ptr = nullptr;
    if (weights.isNotNull()) {
        NumericVector w(weights);
        if (w.size() != X_support.nrow()) {
            Rcpp::stop("'weights' must have nrow(X_support) entries.");
        }
        weight_vec.resize(w.size());
        for (int i = 0; i < w.size(); ++i) weight_vec(i) = w[i];
        weight_ptr = &weight_vec;
    }

    geosmooth::second_order_svd_chart_result_t chart =
        geosmooth::compute_local_second_order_svd_chart(
            local,
            anchor,
            chart_dim,
            center_mode,
            weight_ptr,
            rank_tolerance,
            rank_absolute_tolerance,
            curvature_condition_max,
            curvature_ridge,
            min_curvature_support,
            rebase_to_anchor,
            orient_basis
        );

    return List::create(
        Rcpp::Named("coordinates") =
            eigen_matrix_to_numeric_second_order(chart.coordinates),
        Rcpp::Named("basis") =
            eigen_matrix_to_numeric_second_order(chart.basis),
        Rcpp::Named("preliminary.basis") =
            eigen_matrix_to_numeric_second_order(chart.preliminary_basis),
        Rcpp::Named("normal.basis") = R_NilValue,
        Rcpp::Named("rho.tangent") =
            eigen_matrix_to_numeric_second_order(chart.rho_tangent),
        Rcpp::Named("corrected.residual") =
            eigen_matrix_to_numeric_second_order(chart.corrected_residual),
        Rcpp::Named("first.singular.values") =
            eigen_vector_to_numeric_second_order(chart.first_singular_values),
        Rcpp::Named("second.singular.values") =
            eigen_vector_to_numeric_second_order(chart.second_singular_values),
        Rcpp::Named("chart.dim") = chart.selected_dim,
        Rcpp::Named("curvature.coefficients") =
            eigen_matrix_to_numeric_second_order(chart.curvature_coefficients),
        Rcpp::Named("curvature.monomials") =
            monomials_to_frame(chart),
        Rcpp::Named("curvature.diagnostics") =
            diagnostics_to_frame(chart),
        Rcpp::Named("fallback.used") = chart.fallback_used,
        Rcpp::Named("fallback.reason") = chart.fallback_reason,
        Rcpp::Named("primary.failure.reason") =
            chart.primary_failure_reason
    );
}
