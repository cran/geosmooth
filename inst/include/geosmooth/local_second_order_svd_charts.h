#ifndef GEOSMOOTH_LOCAL_SECOND_ORDER_SVD_CHARTS_H_
#define GEOSMOOTH_LOCAL_SECOND_ORDER_SVD_CHARTS_H_

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace geosmooth {

struct second_order_svd_chart_result_t {
    Eigen::MatrixXd coordinates;
    Eigen::MatrixXd basis;
    Eigen::MatrixXd preliminary_basis;
    Eigen::MatrixXd rho_tangent;
    Eigen::MatrixXd corrected_residual;
    Eigen::MatrixXd curvature_coefficients;
    Eigen::VectorXd first_singular_values;
    Eigen::VectorXd second_singular_values;
    std::vector<int> monomial_a;
    std::vector<int> monomial_b;
    std::vector<double> monomial_multiplier;
    int selected_dim = 0;
    int effective_support = 0;
    int quadratic_ncol = 0;
    int design_rank = 0;
    int first_rank = 0;
    int second_rank = 0;
    double design_condition = 0.0;
    double ridge_lambda = 0.0;
    double fit_residual_frobenius = 0.0;
    double curvature_fitted_frobenius = 0.0;
    double corrected_residual_frobenius = 0.0;
    std::string fit_method = "none";
    std::string status = "ok";
    bool fallback_used = false;
    bool plain_pca_fallback_feasible = false;
    std::string fallback_reason = "none";
    std::string primary_failure_reason = "none";
};

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
);

} // namespace geosmooth

#endif // GEOSMOOTH_LOCAL_SECOND_ORDER_SVD_CHARTS_H_
