#include <Rcpp.h>

#include <geosmooth/local_pca_charts.h>

#include <R_ext/Lapack.h>

#include <ANN/ANN.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Rcpp::CharacterVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

using steady_clock = std::chrono::steady_clock;

double elapsed_seconds(const steady_clock::time_point& start,
                       const steady_clock::time_point& end) {
    return std::chrono::duration<double>(end - start).count();
}

enum class klp_kernel_t {
    gaussian,
    tricube,
    epanechnikov,
    triangular
};

klp_kernel_t parse_kernel(const std::string& kernel) {
    if (kernel == "gaussian") return klp_kernel_t::gaussian;
    if (kernel == "tricube") return klp_kernel_t::tricube;
    if (kernel == "epanechnikov") return klp_kernel_t::epanechnikov;
    if (kernel == "triangular") return klp_kernel_t::triangular;
    Rcpp::stop("Unsupported kernel: %s", kernel);
}

int design_ncol(const int degree, const int chart_dim) {
    if (degree == 0) return 1;
    if (degree == 1) return 1 + chart_dim;
    if (degree == 2) return 1 + chart_dim + chart_dim * (chart_dim + 1) / 2;
    Rcpp::stop("Unsupported local polynomial degree: %d", degree);
}

double kernel_weight(const double distance, const double bandwidth,
                     const klp_kernel_t kernel) {
    const double denom = bandwidth + std::sqrt(std::numeric_limits<double>::epsilon());
    const double u = distance / denom;
    if (!std::isfinite(u)) return 0.0;
    switch (kernel) {
    case klp_kernel_t::gaussian:
        return std::exp(-0.5 * u * u);
    case klp_kernel_t::tricube:
        if (u >= 1.0) return 0.0;
        return std::pow(1.0 - u * u * u, 3.0);
    case klp_kernel_t::epanechnikov:
        return std::max(0.0, 1.0 - u * u);
    case klp_kernel_t::triangular:
        return std::max(0.0, 1.0 - u);
    }
    return 0.0;
}

double weighted_mean(const std::vector<double>& y,
                     const std::vector<double>& weights,
                     const int n) {
    double sw = 0.0;
    double swy = 0.0;
    for (int i = 0; i < n; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (std::isfinite(wi) && wi > 0.0 && std::isfinite(yi)) {
            sw += wi;
            swy += wi * yi;
        }
    }
    if (sw <= 0.0 || !std::isfinite(sw)) {
        double sy = 0.0;
        int ny = 0;
        for (int i = 0; i < n; ++i) {
            const double yi = y[static_cast<size_t>(i)];
            if (std::isfinite(yi)) {
                sy += yi;
                ++ny;
            }
        }
        return ny > 0 ? sy / static_cast<double>(ny) : NA_REAL;
    }
    return swy / sw;
}

void fill_features(const NumericMatrix& X,
                   const std::vector<int>& original_index,
                   const std::vector<double>& center,
                   const int local_row,
                   const int degree,
                   std::vector<double>& features) {
    const int p = X.ncol();
    features[0] = 1.0;
    if (degree == 0) return;

        const int obs = original_index[static_cast<size_t>(local_row)];
        for (int j = 0; j < p; ++j) {
        features[static_cast<size_t>(1 + j)] = X(obs, j) - center[static_cast<size_t>(j)];
    }
    if (degree == 1) return;

    int col = 1 + p;
    for (int a = 0; a < p; ++a) {
        const double za = features[static_cast<size_t>(1 + a)];
        for (int b = a; b < p; ++b) {
            const double zb = features[static_cast<size_t>(1 + b)];
            features[static_cast<size_t>(col++)] = za * zb;
        }
    }
}

bool solve_spd(std::vector<double>& a, std::vector<double>& b, const int n) {
    int info = 0;
    int nrhs = 1;
    F77_CALL(dpotrf)("L", &n, a.data(), &n, &info FCONE);
    if (info != 0) return false;
    F77_CALL(dpotrs)("L", &n, &nrhs, a.data(), &n, b.data(), &n, &info FCONE);
    return info == 0 && std::isfinite(b[0]);
}

double solve_weighted_lm_qr(const std::vector<std::vector<double> >& features,
                            const std::vector<double>& y,
                            const std::vector<double>& weights,
                            const int n,
                            const int q) {
    int n_ok = 0;
    for (int i = 0; i < n; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (!std::isfinite(wi) || wi <= 0.0 || !std::isfinite(yi)) continue;
        bool finite_features = true;
        for (int j = 0; j < q; ++j) {
            if (!std::isfinite(features[static_cast<size_t>(i)][static_cast<size_t>(j)])) {
                finite_features = false;
                break;
            }
        }
        if (finite_features) ++n_ok;
    }
    if (n_ok < q) return NA_REAL;
    if (n_ok < static_cast<int>(std::ceil(1.25 * static_cast<double>(q)))) {
        return NA_REAL;
    }

    Eigen::MatrixXd Xw(n_ok, q);
    Eigen::VectorXd yw(n_ok);
    int row = 0;
    for (int i = 0; i < n; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (!std::isfinite(wi) || wi <= 0.0 || !std::isfinite(yi)) continue;
        bool finite_features = true;
        for (int j = 0; j < q; ++j) {
            if (!std::isfinite(features[static_cast<size_t>(i)][static_cast<size_t>(j)])) {
                finite_features = false;
                break;
            }
        }
        if (!finite_features) continue;
        const double sw = std::sqrt(wi);
        for (int j = 0; j < q; ++j) {
            Xw(row, j) = sw * features[static_cast<size_t>(i)][static_cast<size_t>(j)];
        }
        yw(row) = sw * yi;
        ++row;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Xw, Eigen::ComputeThinU);
    const Eigen::VectorXd singular_values = svd.singularValues();
    if (singular_values.size() < q) return NA_REAL;
    const double sv_max = singular_values.maxCoeff();
    const double sv_min = singular_values.minCoeff();
    const double tol = static_cast<double>(std::max(n_ok, q)) * sv_max *
        std::numeric_limits<double>::epsilon();
    int rank = 0;
    for (int j = 0; j < singular_values.size(); ++j) {
        if (singular_values(j) > tol) ++rank;
    }
    if (rank < q || sv_min <= 0.0 || !std::isfinite(sv_min)) {
        return NA_REAL;
    }
    const double condition = sv_max / sv_min;
    if (!std::isfinite(condition) || condition > 1e4) {
        return NA_REAL;
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Xw);
    qr.setThreshold(1e-7);
    if (qr.rank() < q) return NA_REAL;
    Eigen::VectorXd beta = qr.solve(yw);
    if (beta.size() < 1 || !std::isfinite(beta(0))) return NA_REAL;
    return beta(0);
}

double solve_weighted_lm_r_compatible(
        const std::vector<std::vector<double> >& features,
        const std::vector<double>& y,
        const std::vector<double>& weights,
        const int n,
        const int q) {
    int n_ok = 0;
    for (int i = 0; i < n; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (!std::isfinite(wi) || wi <= 0.0 || !std::isfinite(yi)) continue;
        bool finite_features = true;
        for (int j = 0; j < q; ++j) {
            if (!std::isfinite(features[static_cast<size_t>(i)][static_cast<size_t>(j)])) {
                finite_features = false;
                break;
            }
        }
        if (finite_features) ++n_ok;
    }
    if (n_ok < q) return NA_REAL;
    if (n_ok < static_cast<int>(std::ceil(1.25 * static_cast<double>(q)))) {
        return NA_REAL;
    }

    NumericMatrix design(n_ok, q);
    NumericVector response(n_ok);
    NumericVector fit_weights(n_ok);
    int row = 0;
    for (int i = 0; i < n; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (!std::isfinite(wi) || wi <= 0.0 || !std::isfinite(yi)) continue;
        bool finite_features = true;
        for (int j = 0; j < q; ++j) {
            if (!std::isfinite(features[static_cast<size_t>(i)][static_cast<size_t>(j)])) {
                finite_features = false;
                break;
            }
        }
        if (!finite_features) continue;
        for (int j = 0; j < q; ++j) {
            design(row, j) = features[static_cast<size_t>(i)][static_cast<size_t>(j)];
        }
        response[row] = yi;
        fit_weights[row] = wi;
        ++row;
    }

    Eigen::MatrixXd Xw(n_ok, q);
    for (int i = 0; i < n_ok; ++i) {
        const double sw = std::sqrt(fit_weights[i]);
        for (int j = 0; j < q; ++j) {
            Xw(i, j) = sw * design(i, j);
        }
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Xw, Eigen::ComputeThinU);
    const Eigen::VectorXd singular_values = svd.singularValues();
    if (singular_values.size() < q) return NA_REAL;
    const double sv_max = singular_values.maxCoeff();
    const double sv_min = singular_values.minCoeff();
    const double tol = static_cast<double>(std::max(n_ok, q)) * sv_max *
        std::numeric_limits<double>::epsilon();
    int rank = 0;
    for (int j = 0; j < singular_values.size(); ++j) {
        if (singular_values(j) > tol) ++rank;
    }
    if (rank < q || sv_min <= 0.0 || !std::isfinite(sv_min)) {
        return NA_REAL;
    }
    const double condition = sv_max / sv_min;
    if (!std::isfinite(condition) || condition > 1e4) {
        return NA_REAL;
    }

    try {
        Rcpp::Environment stats_env("package:stats");
        Rcpp::Function lm_wfit = stats_env["lm.wfit"];
        List fit = lm_wfit(design, response, fit_weights);
        NumericVector coefficients = fit["coefficients"];
        if (coefficients.size() < 1 || !std::isfinite(coefficients[0])) {
            return NA_REAL;
        }
        return coefficients[0];
    } catch (...) {
        return NA_REAL;
    }
}

double fit_intercept(const NumericMatrix& X,
                     const std::vector<int>& original_index,
                     const std::vector<double>& center,
                     const std::vector<double>& y,
                     const std::vector<double>& weights,
                     const int support_size,
                     const int degree) {
    const int p = X.ncol();
    const int q = design_ncol(degree, p);
    int n_ok = 0;
    for (int i = 0; i < support_size; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (std::isfinite(wi) && wi > 0.0 && std::isfinite(yi)) {
            ++n_ok;
        }
    }
    if (n_ok < q) {
        return weighted_mean(y, weights, support_size);
    }
    if (degree == 0) {
        return weighted_mean(y, weights, support_size);
    }

    std::vector<double> xtwx(static_cast<size_t>(q) * static_cast<size_t>(q), 0.0);
    std::vector<double> xtwy(static_cast<size_t>(q), 0.0);
    std::vector<double> f(static_cast<size_t>(q), 0.0);
    std::vector<std::vector<double> > qr_features(
        static_cast<size_t>(support_size),
        std::vector<double>(static_cast<size_t>(q), NA_REAL)
    );

    for (int i = 0; i < support_size; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (!std::isfinite(wi) || wi <= 0.0 || !std::isfinite(yi)) continue;
        fill_features(X, original_index, center, i, degree, f);
        qr_features[static_cast<size_t>(i)] = f;
        for (int a = 0; a < q; ++a) {
            xtwy[static_cast<size_t>(a)] += wi * f[static_cast<size_t>(a)] * yi;
            for (int b = 0; b <= a; ++b) {
                xtwx[static_cast<size_t>(a) + static_cast<size_t>(q) * static_cast<size_t>(b)] +=
                    wi * f[static_cast<size_t>(a)] * f[static_cast<size_t>(b)];
            }
        }
    }
    for (int a = 0; a < q; ++a) {
        for (int b = 0; b < a; ++b) {
            xtwx[static_cast<size_t>(b) + static_cast<size_t>(q) * static_cast<size_t>(a)] =
                xtwx[static_cast<size_t>(a) + static_cast<size_t>(q) * static_cast<size_t>(b)];
        }
    }
    const double r_intercept = solve_weighted_lm_r_compatible(
        qr_features, y, weights, support_size, q
    );
    if (std::isfinite(r_intercept)) return r_intercept;
    const double qr_intercept = solve_weighted_lm_qr(
        qr_features, y, weights, support_size, q
    );
    if (std::isfinite(qr_intercept)) return qr_intercept;
    if (solve_spd(xtwx, xtwy, q)) return xtwy[0];
    return weighted_mean(y, weights, support_size);
}

void fill_features_from_chart(const Eigen::MatrixXd& Z,
                              const int local_row,
                              const int degree,
                              const int chart_dim,
                              std::vector<double>& features) {
    features[0] = 1.0;
    if (degree == 0) return;

    for (int j = 0; j < chart_dim; ++j) {
        features[static_cast<size_t>(1 + j)] = Z(local_row, j);
    }
    if (degree == 1) return;

    int col = 1 + chart_dim;
    for (int a = 0; a < chart_dim; ++a) {
        const double za = Z(local_row, a);
        for (int b = a; b < chart_dim; ++b) {
            const double zb = Z(local_row, b);
            features[static_cast<size_t>(col++)] = za * zb;
        }
    }
}

double fit_intercept_from_chart(const Eigen::MatrixXd& Z,
                                const std::vector<double>& y,
                                const std::vector<double>& weights,
                                const int support_size,
                                const int degree,
                                const int chart_dim) {
    const int q = design_ncol(degree, chart_dim);
    int n_ok = 0;
    for (int i = 0; i < support_size; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (std::isfinite(wi) && wi > 0.0 && std::isfinite(yi)) {
            ++n_ok;
        }
    }
    if (n_ok < q || degree == 0) {
        return weighted_mean(y, weights, support_size);
    }

    std::vector<double> xtwx(static_cast<size_t>(q) * static_cast<size_t>(q), 0.0);
    std::vector<double> xtwy(static_cast<size_t>(q), 0.0);
    std::vector<double> f(static_cast<size_t>(q), 0.0);
    std::vector<std::vector<double> > qr_features(
        static_cast<size_t>(support_size),
        std::vector<double>(static_cast<size_t>(q), NA_REAL)
    );

    for (int i = 0; i < support_size; ++i) {
        const double wi = weights[static_cast<size_t>(i)];
        const double yi = y[static_cast<size_t>(i)];
        if (!std::isfinite(wi) || wi <= 0.0 || !std::isfinite(yi)) continue;
        fill_features_from_chart(Z, i, degree, chart_dim, f);
        bool finite_features = true;
        for (int a = 0; a < q; ++a) {
            if (!std::isfinite(f[static_cast<size_t>(a)])) {
                finite_features = false;
                break;
            }
        }
        if (!finite_features) continue;
        qr_features[static_cast<size_t>(i)] = f;
        for (int a = 0; a < q; ++a) {
            xtwy[static_cast<size_t>(a)] += wi * f[static_cast<size_t>(a)] * yi;
            for (int b = 0; b <= a; ++b) {
                xtwx[static_cast<size_t>(a) + static_cast<size_t>(q) * static_cast<size_t>(b)] +=
                    wi * f[static_cast<size_t>(a)] * f[static_cast<size_t>(b)];
            }
        }
    }
    for (int a = 0; a < q; ++a) {
        for (int b = 0; b < a; ++b) {
            xtwx[static_cast<size_t>(b) + static_cast<size_t>(q) * static_cast<size_t>(a)] =
                xtwx[static_cast<size_t>(a) + static_cast<size_t>(q) * static_cast<size_t>(b)];
        }
    }
    const double r_intercept = solve_weighted_lm_r_compatible(
        qr_features, y, weights, support_size, q
    );
    if (std::isfinite(r_intercept)) return r_intercept;
    return weighted_mean(y, weights, support_size);
}

geosmooth::local_pca_chart_result_t local_pca_chart_for_support(
        const NumericMatrix& X,
        const std::vector<int>& original_index,
        const std::vector<double>& center,
        const int support_size,
        const int chart_dim) {
    const int p = X.ncol();
    Eigen::MatrixXd local(support_size, p);
    Eigen::RowVectorXd anchor(p);
    for (int j = 0; j < p; ++j) {
        anchor(j) = center[static_cast<size_t>(j)];
    }
    for (int i = 0; i < support_size; ++i) {
        const int row = original_index[static_cast<size_t>(i)];
        for (int j = 0; j < p; ++j) {
            local(i, j) = X(row, j);
        }
    }
    return geosmooth::compute_local_pca_chart(
        local,
        anchor,
        chart_dim,
        "fixed",
        0.9,
        "anchor",
        nullptr,
        true,
        false
    );
}

double fit_intercept_local_pca(const NumericMatrix& X,
                               const std::vector<int>& original_index,
                               const std::vector<double>& center,
                               const std::vector<double>& y,
                               const std::vector<double>& weights,
                               const int support_size,
                               const int degree,
                               const int chart_dim) {
    if (chart_dim < 1 || chart_dim > X.ncol()) {
        Rcpp::stop("'chart_dim' must be between 1 and ncol(X).");
    }
    geosmooth::local_pca_chart_result_t chart =
        local_pca_chart_for_support(X, original_index, center,
                                    support_size, chart_dim);
    return fit_intercept_from_chart(chart.coordinates, y, weights,
                                    support_size, degree, chart_dim);
}

struct local_pca_chart_cache_entry_t {
    int support_size;
    int chart_dim;
    geosmooth::local_pca_chart_result_t chart;
};

int find_local_pca_chart_cache_entry(
        const std::vector<local_pca_chart_cache_entry_t>& cache,
        const int support_size,
        const int chart_dim) {
    for (int i = 0; i < static_cast<int>(cache.size()); ++i) {
        if (cache[static_cast<size_t>(i)].support_size == support_size &&
            cache[static_cast<size_t>(i)].chart_dim == chart_dim) {
            return i;
        }
    }
    return -1;
}

struct neighbor_candidate_t {
    int local_index;
    double squared_distance;
    int original_row;
};

double squared_distance_to_center(const NumericMatrix& X,
                                  const int row,
                                  const std::vector<double>& center) {
    double out = 0.0;
    for (int j = 0; j < X.ncol(); ++j) {
        const double diff = X(row, j) - center[static_cast<size_t>(j)];
        out += diff * diff;
    }
    return out;
}

void sort_neighbors_by_distance_row(std::vector<ANNidx>& nn_idx,
                                    std::vector<ANNdist>& nn_dist,
                                    const std::vector<int>& original_rows,
                                    const int n_neighbors) {
    std::vector<int> ord(static_cast<size_t>(n_neighbors));
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](const int a, const int b) {
        const double da = nn_dist[static_cast<size_t>(a)];
        const double db = nn_dist[static_cast<size_t>(b)];
        if (da < db) return true;
        if (db < da) return false;
        const int ia = nn_idx[static_cast<size_t>(a)];
        const int ib = nn_idx[static_cast<size_t>(b)];
        return original_rows[static_cast<size_t>(ia)] <
            original_rows[static_cast<size_t>(ib)];
    });

    std::vector<ANNidx> sorted_idx(static_cast<size_t>(n_neighbors));
    std::vector<ANNdist> sorted_dist(static_cast<size_t>(n_neighbors));
    for (int i = 0; i < n_neighbors; ++i) {
        const int src = ord[static_cast<size_t>(i)];
        sorted_idx[static_cast<size_t>(i)] = nn_idx[static_cast<size_t>(src)];
        sorted_dist[static_cast<size_t>(i)] = nn_dist[static_cast<size_t>(src)];
    }
    for (int i = 0; i < n_neighbors; ++i) {
        nn_idx[static_cast<size_t>(i)] = sorted_idx[static_cast<size_t>(i)];
        nn_dist[static_cast<size_t>(i)] = sorted_dist[static_cast<size_t>(i)];
    }
}

void recover_tie_complete_neighbors(const NumericMatrix& X,
                                    const std::vector<int>& original_rows,
                                    const std::vector<double>& center,
                                    std::vector<ANNidx>& nn_idx,
                                    std::vector<ANNdist>& nn_dist,
                                    const int n_neighbors) {
    sort_neighbors_by_distance_row(
        nn_idx,
        nn_dist,
        original_rows,
        n_neighbors
    );
    const double kth = nn_dist[static_cast<size_t>(n_neighbors - 1)];
    const double tolerance = 1e-12 * std::max(1.0, std::abs(kth));
    std::vector<neighbor_candidate_t> candidates;
    candidates.reserve(static_cast<size_t>(original_rows.size()));

    for (int local = 0; local < static_cast<int>(original_rows.size()); ++local) {
        const int original = original_rows[static_cast<size_t>(local)];
        const double dist = squared_distance_to_center(X, original, center);
        if (dist < kth || std::abs(dist - kth) <= tolerance) {
            candidates.push_back(neighbor_candidate_t{
                local,
                dist,
                original
            });
        }
    }
    if (static_cast<int>(candidates.size()) < n_neighbors) {
        for (int i = 0; i < n_neighbors; ++i) {
            const int local = nn_idx[static_cast<size_t>(i)];
            const int original = original_rows[static_cast<size_t>(local)];
            candidates.push_back(neighbor_candidate_t{
                local,
                nn_dist[static_cast<size_t>(i)],
                original
            });
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const neighbor_candidate_t& a,
                 const neighbor_candidate_t& b) {
                  if (a.squared_distance < b.squared_distance) return true;
                  if (b.squared_distance < a.squared_distance) return false;
                  return a.original_row < b.original_row;
              });
    candidates.erase(
        std::unique(candidates.begin(), candidates.end(),
                    [](const neighbor_candidate_t& a,
                       const neighbor_candidate_t& b) {
                        return a.local_index == b.local_index;
                    }),
        candidates.end()
    );
    if (static_cast<int>(candidates.size()) < n_neighbors) {
        Rcpp::stop("Internal error: tie-complete neighbor recovery returned fewer than k neighbors.");
    }
    for (int i = 0; i < n_neighbors; ++i) {
        nn_idx[static_cast<size_t>(i)] = candidates[static_cast<size_t>(i)].local_index;
        nn_dist[static_cast<size_t>(i)] = candidates[static_cast<size_t>(i)].squared_distance;
    }
}

class AnnTree {
public:
    AnnTree(const NumericMatrix& X, const std::vector<int>& rows)
        : n_(static_cast<int>(rows.size())), p_(X.ncol()),
          data_(nullptr), tree_(nullptr) {
        if (n_ <= 0 || p_ <= 0) {
            throw std::runtime_error("ANN tree needs positive dimensions");
        }
        data_ = annAllocPts(n_, p_);
        try {
            for (int i = 0; i < n_; ++i) {
                const int row = rows[static_cast<size_t>(i)];
                for (int j = 0; j < p_; ++j) {
                    data_[i][j] = X(row, j);
                }
            }
            tree_ = new ANNkd_tree(data_, n_, p_);
        } catch (...) {
            cleanup();
            throw;
        }
    }

    ~AnnTree() {
        cleanup();
    }

    void search(const std::vector<double>& center, const int k,
                std::vector<ANNidx>& nn_idx,
                std::vector<ANNdist>& nn_dist) const {
        ANNpoint query = annAllocPt(p_);
        try {
            for (int j = 0; j < p_; ++j) query[j] = center[static_cast<size_t>(j)];
            tree_->annkSearch(query, k, nn_idx.data(), nn_dist.data(), 0.0);
        } catch (...) {
            annDeallocPt(query);
            throw;
        }
        annDeallocPt(query);
    }

private:
    int n_;
    int p_;
    ANNpointArray data_;
    ANNkd_tree* tree_;

    void cleanup() {
        if (tree_ != nullptr) {
            delete tree_;
            tree_ = nullptr;
        }
        if (data_ != nullptr) {
            annDeallocPts(data_);
            data_ = nullptr;
        }
        annClose();
    }
};

NumericVector predict_coordinates_cpp(const NumericMatrix& X_train,
                                      const NumericVector& y_train,
                                      const NumericMatrix& X_eval,
                                      const IntegerVector& train_rows,
                                      const int support_size,
                                      const int degree,
                                      const klp_kernel_t kernel) {
    const int n_train = train_rows.size();
    const int n_eval = X_eval.nrow();
    const int k = std::min(support_size, n_train);
    if (k <= 0) Rcpp::stop("No training rows available.");

    std::vector<int> rows(static_cast<size_t>(n_train));
    for (int i = 0; i < n_train; ++i) {
        rows[static_cast<size_t>(i)] = train_rows[i] - 1;
    }
    AnnTree tree(X_train, rows);
    NumericVector out(n_eval);
    std::vector<ANNidx> nn_idx(static_cast<size_t>(k));
    std::vector<ANNdist> nn_dist(static_cast<size_t>(k));
    std::vector<double> local_y(static_cast<size_t>(k));
    std::vector<double> local_w(static_cast<size_t>(k));
    std::vector<int> original_index(static_cast<size_t>(k));
    std::vector<double> center(static_cast<size_t>(X_train.ncol()));

    for (int i = 0; i < n_eval; ++i) {
        for (int j = 0; j < X_train.ncol(); ++j) {
            center[static_cast<size_t>(j)] = X_eval(i, j);
        }
        tree.search(center, k, nn_idx, nn_dist);
        recover_tie_complete_neighbors(
            X_train,
            rows,
            center,
            nn_idx,
            nn_dist,
            k
        );
        double bandwidth = 0.0;
        for (int j = 0; j < k; ++j) {
            bandwidth = std::max(bandwidth, std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)])));
        }
        if (!std::isfinite(bandwidth) || bandwidth <= 0.0) bandwidth = 1.0;
        for (int j = 0; j < k; ++j) {
            const int local = nn_idx[static_cast<size_t>(j)];
            const int original = rows[static_cast<size_t>(local)];
            original_index[static_cast<size_t>(j)] = original;
            local_y[static_cast<size_t>(j)] = y_train[original];
            const double d = std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]));
            local_w[static_cast<size_t>(j)] = kernel_weight(d, bandwidth, kernel);
        }
        out[i] = fit_intercept(X_train, original_index, center, local_y,
                               local_w, k, degree);
    }
    return out;
}

} // namespace

//' Kernel local polynomial CV RMSE for ambient coordinates
//'
//' Internal C++ backend for `fit.lps()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector rcpp_kernel_local_polynomial_cv_coordinates(
        const NumericMatrix& X,
        const NumericVector& y,
        const IntegerVector& foldid,
        const IntegerVector& support_size,
        const IntegerVector& degree,
        const CharacterVector& kernel) {
    const int n = X.nrow();
    const int n_cand = support_size.size();
    if (y.size() != n || foldid.size() != n ||
        degree.size() != n_cand || kernel.size() != n_cand) {
        Rcpp::stop("Inconsistent input lengths.");
    }
    NumericVector sse(n_cand, 0.0);
    IntegerVector count(n_cand);

    std::vector<int> folds(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) folds[static_cast<size_t>(i)] = foldid[i];
    std::sort(folds.begin(), folds.end());
    folds.erase(std::unique(folds.begin(), folds.end()), folds.end());

    std::vector<klp_kernel_t> kernels(static_cast<size_t>(n_cand));
    int max_support = 0;
    for (int rr = 0; rr < n_cand; ++rr) {
        kernels[static_cast<size_t>(rr)] =
            parse_kernel(Rcpp::as<std::string>(kernel[rr]));
        max_support = std::max(max_support, support_size[rr]);
    }

    for (int ff = 0; ff < static_cast<int>(folds.size()); ++ff) {
        const int fold = folds[static_cast<size_t>(ff)];
        std::vector<int> train_zero;
        std::vector<int> test_zero;
        for (int i = 0; i < n; ++i) {
            if (foldid[i] == fold) {
                test_zero.push_back(i);
            } else {
                train_zero.push_back(i);
            }
        }
        if (train_zero.empty()) continue;
        const int fold_k = std::min(max_support, static_cast<int>(train_zero.size()));
        AnnTree tree(X, train_zero);
        std::vector<ANNidx> nn_idx(static_cast<size_t>(fold_k));
        std::vector<ANNdist> nn_dist(static_cast<size_t>(fold_k));
        std::vector<double> local_y(static_cast<size_t>(fold_k));
        std::vector<double> local_w(static_cast<size_t>(fold_k));
        std::vector<int> original_index(static_cast<size_t>(fold_k));
        std::vector<double> center(static_cast<size_t>(X.ncol()));

        for (const int target : test_zero) {
            for (int j = 0; j < X.ncol(); ++j) {
                center[static_cast<size_t>(j)] = X(target, j);
            }
            tree.search(center, fold_k, nn_idx, nn_dist);
            recover_tie_complete_neighbors(
                X,
                train_zero,
                center,
                nn_idx,
                nn_dist,
                fold_k
            );
            std::vector<double> distances(static_cast<size_t>(fold_k));
            for (int j = 0; j < fold_k; ++j) {
                const int local = nn_idx[static_cast<size_t>(j)];
                const int original = train_zero[static_cast<size_t>(local)];
                original_index[static_cast<size_t>(j)] = original;
                local_y[static_cast<size_t>(j)] = y[original];
                distances[static_cast<size_t>(j)] =
                    std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]));
            }
            std::vector<local_pca_chart_cache_entry_t> chart_cache;
            for (int rr = 0; rr < n_cand; ++rr) {
                const int k = std::min(support_size[rr], fold_k);
                double bandwidth = 0.0;
                for (int j = 0; j < k; ++j) {
                    bandwidth = std::max(bandwidth, distances[static_cast<size_t>(j)]);
                }
                if (!std::isfinite(bandwidth) || bandwidth <= 0.0) bandwidth = 1.0;
                for (int j = 0; j < k; ++j) {
                    local_w[static_cast<size_t>(j)] = kernel_weight(
                        distances[static_cast<size_t>(j)],
                        bandwidth,
                        kernels[static_cast<size_t>(rr)]
                    );
                }
                const double pred = fit_intercept(X, original_index, center,
                                                  local_y, local_w, k,
                                                  degree[rr]);
                const double err = pred - y[target];
                if (std::isfinite(err)) {
                    sse[rr] += err * err;
                    count[rr] += 1;
                }
            }
        }
    }

    NumericVector rmse(n_cand);
    for (int rr = 0; rr < n_cand; ++rr) {
        rmse[rr] = count[rr] > 0 ? std::sqrt(sse[rr] / count[rr]) : NA_REAL;
    }
    return rmse;
}

//' Kernel local polynomial CV RMSE for local-PCA coordinates
//'
//' Internal C++ prototype backend for `fit.lps()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector rcpp_kernel_local_polynomial_cv_local_pca(
        const NumericMatrix& X,
        const NumericVector& y,
        const IntegerVector& foldid,
        const IntegerVector& support_size,
        const IntegerVector& degree,
        const CharacterVector& kernel,
        const IntegerVector& chart_dim) {
    const int n = X.nrow();
    const int n_cand = support_size.size();
    if (y.size() != n || foldid.size() != n ||
        degree.size() != n_cand || kernel.size() != n_cand ||
        chart_dim.size() != n_cand) {
        Rcpp::stop("Inconsistent input lengths.");
    }
    NumericVector sse(n_cand, 0.0);
    IntegerVector count(n_cand);

    std::vector<int> folds(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) folds[static_cast<size_t>(i)] = foldid[i];
    std::sort(folds.begin(), folds.end());
    folds.erase(std::unique(folds.begin(), folds.end()), folds.end());

    std::vector<klp_kernel_t> kernels(static_cast<size_t>(n_cand));
    int max_support = 0;
    for (int rr = 0; rr < n_cand; ++rr) {
        kernels[static_cast<size_t>(rr)] =
            parse_kernel(Rcpp::as<std::string>(kernel[rr]));
        if (chart_dim[rr] < 1 || chart_dim[rr] > X.ncol()) {
            Rcpp::stop("'chart_dim' must be between 1 and ncol(X).");
        }
        max_support = std::max(max_support, support_size[rr]);
    }

    for (int ff = 0; ff < static_cast<int>(folds.size()); ++ff) {
        const int fold = folds[static_cast<size_t>(ff)];
        std::vector<int> train_zero;
        std::vector<int> test_zero;
        for (int i = 0; i < n; ++i) {
            if (foldid[i] == fold) {
                test_zero.push_back(i);
            } else {
                train_zero.push_back(i);
            }
        }
        if (train_zero.empty()) continue;
        const int fold_k = std::min(max_support, static_cast<int>(train_zero.size()));
        AnnTree tree(X, train_zero);
        std::vector<ANNidx> nn_idx(static_cast<size_t>(fold_k));
        std::vector<ANNdist> nn_dist(static_cast<size_t>(fold_k));
        std::vector<double> local_y(static_cast<size_t>(fold_k));
        std::vector<double> local_w(static_cast<size_t>(fold_k));
        std::vector<int> original_index(static_cast<size_t>(fold_k));
        std::vector<double> center(static_cast<size_t>(X.ncol()));

        for (const int target : test_zero) {
            for (int j = 0; j < X.ncol(); ++j) {
                center[static_cast<size_t>(j)] = X(target, j);
            }
            tree.search(center, fold_k, nn_idx, nn_dist);
            recover_tie_complete_neighbors(
                X,
                train_zero,
                center,
                nn_idx,
                nn_dist,
                fold_k
            );
            std::vector<double> distances(static_cast<size_t>(fold_k));
            for (int j = 0; j < fold_k; ++j) {
                const int local = nn_idx[static_cast<size_t>(j)];
                const int original = train_zero[static_cast<size_t>(local)];
                original_index[static_cast<size_t>(j)] = original;
                local_y[static_cast<size_t>(j)] = y[original];
                distances[static_cast<size_t>(j)] =
                    std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]));
            }
            std::vector<local_pca_chart_cache_entry_t> chart_cache;
            for (int rr = 0; rr < n_cand; ++rr) {
                const int k = std::min(support_size[rr], fold_k);
                double bandwidth = 0.0;
                for (int j = 0; j < k; ++j) {
                    bandwidth = std::max(bandwidth, distances[static_cast<size_t>(j)]);
                }
                if (!std::isfinite(bandwidth) || bandwidth <= 0.0) bandwidth = 1.0;
                for (int j = 0; j < k; ++j) {
                    local_w[static_cast<size_t>(j)] = kernel_weight(
                        distances[static_cast<size_t>(j)],
                        bandwidth,
                        kernels[static_cast<size_t>(rr)]
                    );
                }
                int cache_index = find_local_pca_chart_cache_entry(
                    chart_cache,
                    k,
                    chart_dim[rr]
                );
                if (cache_index < 0) {
                    chart_cache.push_back(local_pca_chart_cache_entry_t{
                        k,
                        chart_dim[rr],
                        local_pca_chart_for_support(
                            X,
                            original_index,
                            center,
                            k,
                            chart_dim[rr]
                        )
                    });
                    cache_index = static_cast<int>(chart_cache.size()) - 1;
                }
                const double pred = fit_intercept_from_chart(
                    chart_cache[static_cast<size_t>(cache_index)].chart.coordinates,
                    local_y,
                    local_w,
                    k,
                    degree[rr],
                    chart_dim[rr]
                );
                const double err = pred - y[target];
                if (std::isfinite(err)) {
                    sse[rr] += err * err;
                    count[rr] += 1;
                }
            }
        }
    }

    NumericVector rmse(n_cand);
    for (int rr = 0; rr < n_cand; ++rr) {
        rmse[rr] = count[rr] > 0 ? std::sqrt(sse[rr] / count[rr]) : NA_REAL;
    }
    return rmse;
}

//' Profile kernel local polynomial CV RMSE for local-PCA coordinates
//'
//' Internal diagnostic backend for `fit.lps()`. It mirrors
//' `rcpp_kernel_local_polynomial_cv_local_pca()` and returns timing/count
//' diagnostics for the native loop. It is intended for engineering reports,
//' not as a user-facing API.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List rcpp_kernel_local_polynomial_cv_local_pca_profile(
        const NumericMatrix& X,
        const NumericVector& y,
        const IntegerVector& foldid,
        const IntegerVector& support_size,
        const IntegerVector& degree,
        const CharacterVector& kernel,
        const IntegerVector& chart_dim) {
    const int n = X.nrow();
    const int n_cand = support_size.size();
    if (y.size() != n || foldid.size() != n ||
        degree.size() != n_cand || kernel.size() != n_cand ||
        chart_dim.size() != n_cand) {
        Rcpp::stop("Inconsistent input lengths.");
    }

    double t_fold_partition = 0.0;
    double t_tree_build = 0.0;
    double t_ann_search = 0.0;
    double t_tie_recovery = 0.0;
    double t_neighbor_extract = 0.0;
    double t_kernel_weights = 0.0;
    double t_chart_build = 0.0;
    double t_local_solve = 0.0;
    double t_accumulate = 0.0;
    double t_rmse = 0.0;

    int count_folds = 0;
    int count_trees = 0;
    int count_targets = 0;
    int count_ann_searches = 0;
    int count_tie_recoveries = 0;
    int count_candidate_evals = 0;
    int count_chart_builds = 0;
    int count_chart_cache_hits = 0;
    int count_local_solves = 0;

    NumericVector sse(n_cand, 0.0);
    IntegerVector count(n_cand);

    std::vector<int> folds(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) folds[static_cast<size_t>(i)] = foldid[i];
    std::sort(folds.begin(), folds.end());
    folds.erase(std::unique(folds.begin(), folds.end()), folds.end());

    std::vector<klp_kernel_t> kernels(static_cast<size_t>(n_cand));
    int max_support = 0;
    for (int rr = 0; rr < n_cand; ++rr) {
        kernels[static_cast<size_t>(rr)] =
            parse_kernel(Rcpp::as<std::string>(kernel[rr]));
        if (chart_dim[rr] < 1 || chart_dim[rr] > X.ncol()) {
            Rcpp::stop("'chart_dim' must be between 1 and ncol(X).");
        }
        max_support = std::max(max_support, support_size[rr]);
    }

    const steady_clock::time_point total_start = steady_clock::now();
    for (int ff = 0; ff < static_cast<int>(folds.size()); ++ff) {
        const int fold = folds[static_cast<size_t>(ff)];
        ++count_folds;

        steady_clock::time_point phase_start = steady_clock::now();
        std::vector<int> train_zero;
        std::vector<int> test_zero;
        for (int i = 0; i < n; ++i) {
            if (foldid[i] == fold) {
                test_zero.push_back(i);
            } else {
                train_zero.push_back(i);
            }
        }
        t_fold_partition += elapsed_seconds(phase_start, steady_clock::now());
        if (train_zero.empty()) continue;

        const int fold_k = std::min(max_support, static_cast<int>(train_zero.size()));
        phase_start = steady_clock::now();
        AnnTree tree(X, train_zero);
        t_tree_build += elapsed_seconds(phase_start, steady_clock::now());
        ++count_trees;

        std::vector<ANNidx> nn_idx(static_cast<size_t>(fold_k));
        std::vector<ANNdist> nn_dist(static_cast<size_t>(fold_k));
        std::vector<double> local_y(static_cast<size_t>(fold_k));
        std::vector<double> local_w(static_cast<size_t>(fold_k));
        std::vector<int> original_index(static_cast<size_t>(fold_k));
        std::vector<double> center(static_cast<size_t>(X.ncol()));

        for (const int target : test_zero) {
            ++count_targets;
            for (int j = 0; j < X.ncol(); ++j) {
                center[static_cast<size_t>(j)] = X(target, j);
            }

            phase_start = steady_clock::now();
            tree.search(center, fold_k, nn_idx, nn_dist);
            t_ann_search += elapsed_seconds(phase_start, steady_clock::now());
            ++count_ann_searches;

            phase_start = steady_clock::now();
            recover_tie_complete_neighbors(
                X,
                train_zero,
                center,
                nn_idx,
                nn_dist,
                fold_k
            );
            t_tie_recovery += elapsed_seconds(phase_start, steady_clock::now());
            ++count_tie_recoveries;

            phase_start = steady_clock::now();
            std::vector<double> distances(static_cast<size_t>(fold_k));
            for (int j = 0; j < fold_k; ++j) {
                const int local = nn_idx[static_cast<size_t>(j)];
                const int original = train_zero[static_cast<size_t>(local)];
                original_index[static_cast<size_t>(j)] = original;
                local_y[static_cast<size_t>(j)] = y[original];
                distances[static_cast<size_t>(j)] =
                    std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]));
            }
            t_neighbor_extract += elapsed_seconds(phase_start, steady_clock::now());

            std::vector<local_pca_chart_cache_entry_t> chart_cache;
            for (int rr = 0; rr < n_cand; ++rr) {
                ++count_candidate_evals;
                const int k = std::min(support_size[rr], fold_k);

                phase_start = steady_clock::now();
                double bandwidth = 0.0;
                for (int j = 0; j < k; ++j) {
                    bandwidth = std::max(bandwidth, distances[static_cast<size_t>(j)]);
                }
                if (!std::isfinite(bandwidth) || bandwidth <= 0.0) bandwidth = 1.0;
                for (int j = 0; j < k; ++j) {
                    local_w[static_cast<size_t>(j)] = kernel_weight(
                        distances[static_cast<size_t>(j)],
                        bandwidth,
                        kernels[static_cast<size_t>(rr)]
                    );
                }
                t_kernel_weights += elapsed_seconds(phase_start, steady_clock::now());

                int cache_index = find_local_pca_chart_cache_entry(
                    chart_cache,
                    k,
                    chart_dim[rr]
                );
                if (cache_index < 0) {
                    phase_start = steady_clock::now();
                    chart_cache.push_back(local_pca_chart_cache_entry_t{
                        k,
                        chart_dim[rr],
                        local_pca_chart_for_support(
                            X,
                            original_index,
                            center,
                            k,
                            chart_dim[rr]
                        )
                    });
                    t_chart_build += elapsed_seconds(phase_start, steady_clock::now());
                    ++count_chart_builds;
                    cache_index = static_cast<int>(chart_cache.size()) - 1;
                } else {
                    ++count_chart_cache_hits;
                }

                phase_start = steady_clock::now();
                const double pred = fit_intercept_from_chart(
                    chart_cache[static_cast<size_t>(cache_index)].chart.coordinates,
                    local_y,
                    local_w,
                    k,
                    degree[rr],
                    chart_dim[rr]
                );
                t_local_solve += elapsed_seconds(phase_start, steady_clock::now());
                ++count_local_solves;

                phase_start = steady_clock::now();
                const double err = pred - y[target];
                if (std::isfinite(err)) {
                    sse[rr] += err * err;
                    count[rr] += 1;
                }
                t_accumulate += elapsed_seconds(phase_start, steady_clock::now());
            }
        }
    }

    const steady_clock::time_point rmse_start = steady_clock::now();
    NumericVector rmse(n_cand);
    for (int rr = 0; rr < n_cand; ++rr) {
        rmse[rr] = count[rr] > 0 ? std::sqrt(sse[rr] / count[rr]) : NA_REAL;
    }
    t_rmse = elapsed_seconds(rmse_start, steady_clock::now());
    const double total_seconds = elapsed_seconds(total_start, steady_clock::now());

    CharacterVector phases = CharacterVector::create(
        "fold_partition",
        "tree_build",
        "ann_search",
        "tie_recovery",
        "neighbor_extract",
        "kernel_weights",
        "chart_build",
        "local_solve",
        "accumulate",
        "rmse"
    );
    NumericVector seconds = NumericVector::create(
        t_fold_partition,
        t_tree_build,
        t_ann_search,
        t_tie_recovery,
        t_neighbor_extract,
        t_kernel_weights,
        t_chart_build,
        t_local_solve,
        t_accumulate,
        t_rmse
    );
    NumericVector share(seconds.size());
    for (int i = 0; i < seconds.size(); ++i) {
        share[i] = total_seconds > 0.0 ? seconds[i] / total_seconds : NA_REAL;
    }
    List timing = List::create(
        Rcpp::Named("phase") = phases,
        Rcpp::Named("seconds") = seconds,
        Rcpp::Named("share.of.total") = share
    );
    List counts = List::create(
        Rcpp::Named("folds") = count_folds,
        Rcpp::Named("trees") = count_trees,
        Rcpp::Named("targets") = count_targets,
        Rcpp::Named("ann.searches") = count_ann_searches,
        Rcpp::Named("tie.recoveries") = count_tie_recoveries,
        Rcpp::Named("candidate.evals") = count_candidate_evals,
        Rcpp::Named("chart.builds") = count_chart_builds,
        Rcpp::Named("chart.cache.hits") = count_chart_cache_hits,
        Rcpp::Named("local.solves") = count_local_solves
    );

    return List::create(
        Rcpp::Named("rmse") = rmse,
        Rcpp::Named("timing") = timing,
        Rcpp::Named("counts") = counts,
        Rcpp::Named("total.seconds") = total_seconds
    );
}

//' Kernel local polynomial predictions for ambient coordinates
//'
//' Internal C++ backend for `fit.lps()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector rcpp_kernel_local_polynomial_predict_coordinates(
        const NumericMatrix& X_train,
        const NumericVector& y_train,
        const NumericMatrix& X_eval,
        const int support_size,
        const int degree,
        const std::string& kernel) {
    if (y_train.size() != X_train.nrow()) {
        Rcpp::stop("'y_train' must have length nrow(X_train).");
    }
    if (X_eval.ncol() != X_train.ncol()) {
        Rcpp::stop("'X_eval' must have ncol(X_train) columns.");
    }
    IntegerVector train_rows(X_train.nrow());
    for (int i = 0; i < X_train.nrow(); ++i) train_rows[i] = i + 1;
    return predict_coordinates_cpp(X_train, y_train, X_eval, train_rows,
                                   support_size, degree,
                                   parse_kernel(kernel));
}

//' Kernel local polynomial predictions for local-PCA coordinates
//'
//' Internal C++ prototype backend for `fit.lps()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector rcpp_kernel_local_polynomial_predict_local_pca(
        const NumericMatrix& X_train,
        const NumericVector& y_train,
        const NumericMatrix& X_eval,
        const int support_size,
        const int degree,
        const std::string& kernel,
        const int chart_dim) {
    if (y_train.size() != X_train.nrow()) {
        Rcpp::stop("'y_train' must have length nrow(X_train).");
    }
    if (X_eval.ncol() != X_train.ncol()) {
        Rcpp::stop("'X_eval' must have ncol(X_train) columns.");
    }
    const int n_train = X_train.nrow();
    const int n_eval = X_eval.nrow();
    const int k = std::min(support_size, n_train);
    if (k <= 0) Rcpp::stop("No training rows available.");
    if (chart_dim < 1 || chart_dim > X_train.ncol()) {
        Rcpp::stop("'chart_dim' must be between 1 and ncol(X_train).");
    }

    std::vector<int> rows(static_cast<size_t>(n_train));
    for (int i = 0; i < n_train; ++i) rows[static_cast<size_t>(i)] = i;
    AnnTree tree(X_train, rows);
    NumericVector out(n_eval);
    std::vector<ANNidx> nn_idx(static_cast<size_t>(k));
    std::vector<ANNdist> nn_dist(static_cast<size_t>(k));
    std::vector<double> local_y(static_cast<size_t>(k));
    std::vector<double> local_w(static_cast<size_t>(k));
    std::vector<int> original_index(static_cast<size_t>(k));
    std::vector<double> center(static_cast<size_t>(X_train.ncol()));
    const klp_kernel_t parsed_kernel = parse_kernel(kernel);

    for (int i = 0; i < n_eval; ++i) {
        for (int j = 0; j < X_train.ncol(); ++j) {
            center[static_cast<size_t>(j)] = X_eval(i, j);
        }
        tree.search(center, k, nn_idx, nn_dist);
        recover_tie_complete_neighbors(
            X_train,
            rows,
            center,
            nn_idx,
            nn_dist,
            k
        );
        double bandwidth = 0.0;
        for (int j = 0; j < k; ++j) {
            bandwidth = std::max(
                bandwidth,
                std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]))
            );
        }
        if (!std::isfinite(bandwidth) || bandwidth <= 0.0) bandwidth = 1.0;
        for (int j = 0; j < k; ++j) {
            const int local = nn_idx[static_cast<size_t>(j)];
            const int original = rows[static_cast<size_t>(local)];
            original_index[static_cast<size_t>(j)] = original;
            local_y[static_cast<size_t>(j)] = y_train[original];
            const double d = std::sqrt(static_cast<double>(nn_dist[static_cast<size_t>(j)]));
            local_w[static_cast<size_t>(j)] = kernel_weight(d, bandwidth, parsed_kernel);
        }
        out[i] = fit_intercept_local_pca(
            X_train,
            original_index,
            center,
            local_y,
            local_w,
            k,
            degree,
            chart_dim
        );
    }
    return out;
}

//' Inspect raw ANN and tie-complete neighbor order
//'
//' Internal diagnostic backend for tests. It exposes the raw ANN order and the
//' geosmooth tie-complete support order for one query point.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List rcpp_kernel_local_polynomial_neighbor_probe(
        const NumericMatrix& X,
        const NumericVector& center,
        const int k) {
    if (center.size() != X.ncol()) {
        Rcpp::stop("'center' must have length ncol(X).");
    }
    const int n = X.nrow();
    if (k < 1 || k > n) {
        Rcpp::stop("'k' must be between 1 and nrow(X).");
    }
    std::vector<int> rows(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) rows[static_cast<size_t>(i)] = i;
    std::vector<double> query(static_cast<size_t>(X.ncol()));
    for (int j = 0; j < X.ncol(); ++j) query[static_cast<size_t>(j)] = center[j];

    AnnTree tree(X, rows);
    std::vector<ANNidx> raw_idx(static_cast<size_t>(k));
    std::vector<ANNdist> raw_dist(static_cast<size_t>(k));
    tree.search(query, k, raw_idx, raw_dist);

    std::vector<ANNidx> tie_idx = raw_idx;
    std::vector<ANNdist> tie_dist = raw_dist;
    recover_tie_complete_neighbors(X, rows, query, tie_idx, tie_dist, k);

    std::vector<neighbor_candidate_t> reference;
    reference.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        reference.push_back(neighbor_candidate_t{
            i,
            squared_distance_to_center(X, i, query),
            i
        });
    }
    std::sort(reference.begin(), reference.end(),
              [](const neighbor_candidate_t& a,
                 const neighbor_candidate_t& b) {
                  if (a.squared_distance < b.squared_distance) return true;
                  if (b.squared_distance < a.squared_distance) return false;
                  return a.original_row < b.original_row;
              });

    IntegerVector raw_row(k);
    NumericVector raw_sqdist(k);
    IntegerVector tie_row(k);
    NumericVector tie_sqdist(k);
    IntegerVector reference_row(k);
    NumericVector reference_sqdist(k);
    for (int i = 0; i < k; ++i) {
        raw_row[i] = rows[static_cast<size_t>(raw_idx[static_cast<size_t>(i)])] + 1;
        raw_sqdist[i] = raw_dist[static_cast<size_t>(i)];
        tie_row[i] = rows[static_cast<size_t>(tie_idx[static_cast<size_t>(i)])] + 1;
        tie_sqdist[i] = tie_dist[static_cast<size_t>(i)];
        reference_row[i] = reference[static_cast<size_t>(i)].original_row + 1;
        reference_sqdist[i] = reference[static_cast<size_t>(i)].squared_distance;
    }

    return List::create(
        Rcpp::Named("raw.row") = raw_row,
        Rcpp::Named("raw.squared.distance") = raw_sqdist,
        Rcpp::Named("tie.complete.row") = tie_row,
        Rcpp::Named("tie.complete.squared.distance") = tie_sqdist,
        Rcpp::Named("reference.row") = reference_row,
        Rcpp::Named("reference.squared.distance") = reference_sqdist
    );
}
