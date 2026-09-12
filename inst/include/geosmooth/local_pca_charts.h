#ifndef GEOSMOOTH_LOCAL_PCA_CHARTS_H_
#define GEOSMOOTH_LOCAL_PCA_CHARTS_H_

#include <Eigen/Dense>

#include <string>

namespace geosmooth {

struct local_pca_chart_result_t {
    Eigen::MatrixXd coordinates;
    Eigen::MatrixXd basis;
    Eigen::VectorXd singular_values;
    int selected_dim = 0;
    double total_variance = 0.0;
    double selected_variance_ratio = 1.0;
};

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
);

} // namespace geosmooth

#endif // GEOSMOOTH_LOCAL_PCA_CHARTS_H_
