#include <geosmooth/local_pca_charts.h>

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace geosmooth {

namespace {

Eigen::RowVectorXd chart_center(const Eigen::MatrixXd& local,
                                const Eigen::RowVectorXd& anchor,
                                const std::string& center_mode) {
    if (center_mode == "anchor") {
        return anchor;
    }
    if (center_mode == "mean") {
        return local.colwise().mean();
    }
    Rcpp::stop("'center.mode' must be either 'anchor' or 'mean'.");
}

void orient_basis_columns(Eigen::MatrixXd& basis) {
    for (int j = 0; j < basis.cols(); ++j) {
        Eigen::Index max_index = 0;
        basis.col(j).cwiseAbs().maxCoeff(&max_index);
        if (basis(max_index, j) < 0.0) {
            basis.col(j) *= -1.0;
        }
    }
}

int select_chart_dim_from_singular_values(const Eigen::VectorXd& sing,
                                          int requested_dim,
                                          const std::string& dim_rule,
                                          double eigen_tolerance,
                                          int local_cols) {
    const int max_dim =
        std::min(static_cast<int>(sing.size()), local_cols);
    if (max_dim < 1) {
        Rcpp::stop("Local PCA has no available singular vectors.");
    }

    double total_var = 0.0;
    for (int a = 0; a < sing.size(); ++a) {
        total_var += sing(a) * sing(a);
    }

    int m = requested_dim;
    if (dim_rule == "eigen.cumulative" && requested_dim <= 0) {
        m = 1;
        double cum = 0.0;
        for (int a = 0; a < max_dim; ++a) {
            cum += sing(a) * sing(a);
            if (total_var <= 0.0 || cum / total_var >= eigen_tolerance) {
                m = a + 1;
                break;
            }
        }
    } else if (dim_rule == "eigen.cumulative" && requested_dim > 0) {
        m = std::min(requested_dim, max_dim);
    }
    if (m < 1 || m > max_dim) {
        Rcpp::stop("The selected chart dimension must be between 1 and min(nrow(local), ncol(local)).");
    }
    return m;
}

double singular_variance(const Eigen::VectorXd& sing, int count) {
    double out = 0.0;
    const int n = std::min(count, static_cast<int>(sing.size()));
    for (int a = 0; a < n; ++a) {
        out += sing(a) * sing(a);
    }
    return out;
}

bool compute_chart_from_row_gram(const Eigen::MatrixXd& centered,
                                 int requested_dim,
                                 const std::string& dim_rule,
                                 double eigen_tolerance,
                                 Eigen::VectorXd& sing,
                                 Eigen::MatrixXd& basis,
                                 int& selected_dim) {
    if (centered.rows() >= centered.cols()) {
        return false;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(
        centered * centered.transpose()
    );
    if (eig.info() != Eigen::Success) {
        return false;
    }

    const Eigen::VectorXd evals = eig.eigenvalues();
    const Eigen::MatrixXd evecs = eig.eigenvectors();
    const int nval = static_cast<int>(evals.size());
    sing.resize(nval);
    for (int a = 0; a < nval; ++a) {
        const double lambda = std::max(0.0, evals(nval - 1 - a));
        sing(a) = std::sqrt(lambda);
    }

    selected_dim = select_chart_dim_from_singular_values(
        sing,
        requested_dim,
        dim_rule,
        eigen_tolerance,
        static_cast<int>(centered.cols())
    );

    const double max_sing = sing.size() > 0 ? sing.maxCoeff() : 0.0;
    const double min_usable =
        std::max(max_sing * 1e-10,
                 std::sqrt(std::numeric_limits<double>::epsilon()));
    if (max_sing <= 0.0 || sing(selected_dim - 1) <= min_usable) {
        return false;
    }

    basis.resize(centered.cols(), selected_dim);
    for (int a = 0; a < selected_dim; ++a) {
        const int idx = nval - 1 - a;
        Eigen::VectorXd v = centered.transpose() * evecs.col(idx);
        v /= sing(a);
        const double norm = v.norm();
        if (!std::isfinite(norm) || norm <= 0.0) {
            return false;
        }
        basis.col(a) = v / norm;
    }
    return true;
}

} // namespace

local_pca_chart_result_t compute_local_pca_chart(
    const Eigen::MatrixXd& local,
    const Eigen::RowVectorXd& anchor,
    int requested_dim,
    const std::string& dim_rule,
    double eigen_tolerance,
    const std::string& center_mode,
    const Eigen::VectorXd* weights,
    bool rebase_to_anchor,
    bool orient_basis
) {
    if (local.rows() < 1 || local.cols() < 1) {
        Rcpp::stop("'local' must have positive dimensions.");
    }
    if (anchor.size() != local.cols()) {
        Rcpp::stop("'anchor' must have ncol(local) entries.");
    }
    if (weights != nullptr && weights->size() != local.rows()) {
        Rcpp::stop("'weights' must have one entry per local row.");
    }
    if (dim_rule != "fixed" && dim_rule != "eigen.cumulative") {
        Rcpp::stop("'dim.rule' must be 'fixed' or 'eigen.cumulative'.");
    }

    const Eigen::RowVectorXd svd_center =
        chart_center(local, anchor, center_mode);
    Eigen::MatrixXd centered = local.rowwise() - svd_center;
    if (weights != nullptr) {
        for (int i = 0; i < centered.rows(); ++i) {
            const double wi = (*weights)(i);
            const double scale = (std::isfinite(wi) && wi > 0.0) ?
                std::sqrt(wi) : 0.0;
            centered.row(i) *= scale;
        }
    }

    Eigen::VectorXd sing;
    Eigen::MatrixXd basis;
    int m = 0;
    bool have_basis = compute_chart_from_row_gram(
        centered,
        requested_dim,
        dim_rule,
        eigen_tolerance,
        sing,
        basis,
        m
    );

    if (!have_basis) {
        Eigen::JacobiSVD<Eigen::MatrixXd> pca(
            centered, Eigen::ComputeThinU | Eigen::ComputeThinV
        );
        sing = pca.singularValues();
        m = select_chart_dim_from_singular_values(
            sing,
            requested_dim,
            dim_rule,
            eigen_tolerance,
            static_cast<int>(local.cols())
        );
        basis = pca.matrixV().leftCols(m);
    }

    if (orient_basis) {
        orient_basis_columns(basis);
    }
    const Eigen::RowVectorXd coord_origin =
        rebase_to_anchor ? anchor : svd_center;
    Eigen::MatrixXd coordinates = (local.rowwise() - coord_origin) * basis;

    const double total_var = singular_variance(
        sing,
        static_cast<int>(sing.size())
    );
    const double selected_var = singular_variance(sing, m);

    local_pca_chart_result_t out;
    out.coordinates = coordinates;
    out.basis = basis;
    out.singular_values = sing;
    out.selected_dim = m;
    out.total_variance = total_var;
    out.selected_variance_ratio = total_var > 0.0 ?
        selected_var / total_var : 1.0;
    return out;
}

} // namespace geosmooth
