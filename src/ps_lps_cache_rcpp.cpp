#include <geosmooth/local_pca_charts.h>

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

Eigen::MatrixXd numeric_matrix_to_eigen_rows(
        const NumericMatrix& X,
        const std::vector<int>& rows) {
    Eigen::MatrixXd out(rows.size(), X.ncol());
    for (int jj = 0; jj < X.ncol(); ++jj) {
        for (int ii = 0; ii < static_cast<int>(rows.size()); ++ii) {
            out(ii, jj) = X(rows[ii], jj);
        }
    }
    return out;
}

NumericMatrix eigen_matrix_to_numeric_local(const Eigen::MatrixXd& X) {
    NumericMatrix out(X.rows(), X.cols());
    for (int jj = 0; jj < X.cols(); ++jj) {
        for (int ii = 0; ii < X.rows(); ++ii) {
            out(ii, jj) = X(ii, jj);
        }
    }
    return out;
}

NumericVector ps_lps_kernel_weights(const std::vector<double>& distances,
                                    const std::string& kernel) {
    const int n = distances.size();
    NumericVector weights(n);
    if (n == 0) {
        return weights;
    }
    double h = 0.0;
    for (int ii = 0; ii < n; ++ii) {
        if (R_finite(distances[ii])) {
            h = std::max(h, distances[ii]);
        }
    }
    if (!R_finite(h) || h <= 0.0) {
        h = 1.0;
    }
    const double denom = h + std::sqrt(std::numeric_limits<double>::epsilon());
    for (int ii = 0; ii < n; ++ii) {
        const double u = distances[ii] / denom;
        double w = 0.0;
        if (kernel == "gaussian") {
            w = std::exp(-0.5 * u * u);
        } else if (kernel == "tricube") {
            w = u < 1.0 ? std::pow(1.0 - u * u * u, 3.0) : 0.0;
        } else if (kernel == "epanechnikov") {
            w = std::max(0.0, 1.0 - u * u);
        } else if (kernel == "triangular") {
            w = std::max(0.0, 1.0 - u);
        } else {
            stop("Unsupported PS-LPS kernel in native frame helper.");
        }
        weights[ii] = R_finite(w) ? w : 0.0;
    }
    bool any_positive = false;
    for (int ii = 0; ii < n; ++ii) {
        if (weights[ii] > 0.0) {
            any_positive = true;
            break;
        }
    }
    if (!any_positive) {
        for (int ii = 0; ii < n; ++ii) {
            weights[ii] = 1.0;
        }
    }
    return weights;
}

} // namespace

// [[Rcpp::export]]
List rcpp_ps_lps_local_pca_supports(
        const NumericMatrix& X,
        const int support_size,
        const IntegerVector& chart_dim_by_anchor,
        const std::string& kernel) {
    const int n = X.nrow();
    const int p = X.ncol();
    if (n < 1 || p < 1) {
        stop("'X' must have positive dimensions.");
    }
    if (chart_dim_by_anchor.size() != n) {
        stop("'chart_dim_by_anchor' must have one entry per row of X.");
    }
    const int k = std::min(std::max(1, support_size), n);
    List out(n);
    std::vector<std::pair<double, int> > order(n);
    for (int aa = 0; aa < n; ++aa) {
        for (int ii = 0; ii < n; ++ii) {
            double d2 = 0.0;
            for (int jj = 0; jj < p; ++jj) {
                const double delta = X(ii, jj) - X(aa, jj);
                d2 += delta * delta;
            }
            order[ii] = std::make_pair(d2, ii);
        }
        std::sort(order.begin(), order.end(),
                  [](const std::pair<double, int>& lhs,
                     const std::pair<double, int>& rhs) {
                      if (lhs.first < rhs.first) return true;
                      if (lhs.first > rhs.first) return false;
                      return lhs.second < rhs.second;
                  });
        IntegerVector index(k);
        NumericVector distances(k);
        std::vector<int> rows(k);
        std::vector<double> dvec(k);
        for (int rr = 0; rr < k; ++rr) {
            rows[rr] = order[rr].second;
            index[rr] = rows[rr] + 1;
            const double dist = std::sqrt(order[rr].first);
            distances[rr] = dist;
            dvec[rr] = dist;
        }
        NumericVector weights = ps_lps_kernel_weights(dvec, kernel);
        Eigen::MatrixXd local = numeric_matrix_to_eigen_rows(X, rows);
        Eigen::RowVectorXd anchor(p);
        for (int jj = 0; jj < p; ++jj) {
            anchor(jj) = X(aa, jj);
        }
        const int chart_dim = chart_dim_by_anchor[aa];
        geosmooth::local_pca_chart_result_t chart =
            geosmooth::compute_local_pca_chart(
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
        out[aa] = List::create(
            _["index"] = index,
            _["distances"] = distances,
            _["weights"] = weights,
            _["coordinates"] = eigen_matrix_to_numeric_local(chart.coordinates)
        );
    }
    return out;
}

// [[Rcpp::export]]
List rcpp_ps_lps_prepare_sync_rows(
        const List& frames,
        const int sync_neighbor_size,
        const std::string& overlap_weight) {
    const int n = frames.size();
    std::set<std::pair<int, int> > pairs;
    for (int ii = 0; ii < n; ++ii) {
        const List fr = frames[ii];
        const IntegerVector index = fr["index"];
        int added = 0;
        for (int rr = 0; rr < index.size(); ++rr) {
            const int jj = index[rr] - 1;
            if (jj == ii) {
                continue;
            }
            pairs.insert(std::make_pair(std::min(ii, jj), std::max(ii, jj)));
            ++added;
            if (added >= sync_neighbor_size) {
                break;
            }
        }
    }
    List out(pairs.size());
    int out_pos = 0;
    for (const auto& pr : pairs) {
        const int ii = pr.first;
        const int jj = pr.second;
        const List fi = frames[ii];
        const List fj = frames[jj];
        const IntegerVector index_i = fi["index"];
        const IntegerVector index_j = fj["index"];
        const NumericVector weights_i = fi["weights"];
        const NumericVector weights_j = fj["weights"];
        std::map<int, int> row_j_by_point;
        for (int rr = 0; rr < index_j.size(); ++rr) {
            row_j_by_point[index_j[rr]] = rr + 1;
        }
        std::vector<int> point;
        std::vector<int> row_i;
        std::vector<int> row_j;
        std::vector<double> prod;
        for (int rr = 0; rr < index_i.size(); ++rr) {
            const int pp = index_i[rr];
            auto it = row_j_by_point.find(pp);
            if (it == row_j_by_point.end()) {
                continue;
            }
            const int rj = it->second;
            const double value = weights_i[rr] * weights_j[rj - 1];
            point.push_back(pp);
            row_i.push_back(rr + 1);
            row_j.push_back(rj);
            prod.push_back(value);
        }
        if (point.empty()) {
            continue;
        }
        double prod_sum = 0.0;
        for (double value : prod) {
            prod_sum += value;
        }
        std::vector<int> keep_point;
        std::vector<int> keep_row_i;
        std::vector<int> keep_row_j;
        std::vector<double> keep_omega;
        for (int rr = 0; rr < static_cast<int>(point.size()); ++rr) {
            double omega = prod[rr];
            if (overlap_weight == "normalized.product") {
                omega = point.size() * prod[rr] /
                    (prod_sum + std::sqrt(std::numeric_limits<double>::epsilon()));
            }
            if (R_finite(omega) && omega > 0.0) {
                keep_point.push_back(point[rr]);
                keep_row_i.push_back(row_i[rr]);
                keep_row_j.push_back(row_j[rr]);
                keep_omega.push_back(omega);
            }
        }
        if (keep_point.empty()) {
            continue;
        }
        out[out_pos] = List::create(
            _["i"] = ii + 1,
            _["j"] = jj + 1,
            _["point"] = wrap(keep_point),
            _["row.i"] = wrap(keep_row_i),
            _["row.j"] = wrap(keep_row_j),
            _["omega"] = wrap(keep_omega)
        );
        ++out_pos;
    }
    if (out_pos < out.size()) {
        List trimmed(out_pos);
        for (int ii = 0; ii < out_pos; ++ii) {
            trimmed[ii] = out[ii];
        }
        return trimmed;
    }
    return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_ps_lps_rhs_matrix(const List& frames,
                                     const NumericMatrix& y_mat) {
    const int ncoef = as<int>(frames.attr("ncoef"));
    const int ncol_y = y_mat.ncol();
    NumericMatrix rhs(ncoef, ncol_y);
    for (int ii = 0; ii < frames.size(); ++ii) {
        const List fr = frames[ii];
        const IntegerVector index = fr["index"];
        const NumericVector weights = fr["weights"];
        const NumericMatrix design = fr["design"];
        const int q = as<int>(fr["q"]);
        const int offset = as<int>(fr["offset"]);
        for (int rr = 0; rr < index.size(); ++rr) {
            const int point = index[rr] - 1;
            const double w = weights[rr];
            if (!R_finite(w)) {
                continue;
            }
            for (int qq = 0; qq < q; ++qq) {
                const double x = design(rr, qq) * w;
                for (int cc = 0; cc < ncol_y; ++cc) {
                    rhs(offset + qq, cc) += x * y_mat(point, cc);
                }
            }
        }
    }
    return rhs;
}

// [[Rcpp::export]]
NumericMatrix rcpp_ps_lps_fitted_matrix(const List& frames,
                                        const NumericMatrix& beta) {
    const int n = frames.size();
    const int ncol_beta = beta.ncol();
    NumericMatrix fitted(n, ncol_beta);
    for (int ii = 0; ii < n; ++ii) {
        const List fr = frames[ii];
        const NumericMatrix anchor_design = fr["anchor.design"];
        const int q = as<int>(fr["q"]);
        const int offset = as<int>(fr["offset"]);
        for (int cc = 0; cc < ncol_beta; ++cc) {
            double value = 0.0;
            for (int qq = 0; qq < q; ++qq) {
                value += anchor_design(0, qq) * beta(offset + qq, cc);
            }
            fitted(ii, cc) = value;
        }
    }
    return fitted;
}

// [[Rcpp::export]]
NumericMatrix rcpp_ps_lps_independent_fitted_matrix(
        const List& frames,
        const NumericMatrix& y_mat,
        const double ridge_multiplier) {
    const int n = frames.size();
    const int ncol_y = y_mat.ncol();
    NumericMatrix fitted(n, ncol_y);
    for (int ii = 0; ii < n; ++ii) {
        const List fr = frames[ii];
        bool active = true;
        if (fr.containsElementNamed("active")) {
            active = as<bool>(fr["active"]);
        }
        if (!active) {
            for (int cc = 0; cc < ncol_y; ++cc) {
                fitted(ii, cc) = 0.0;
            }
            continue;
        }
        const IntegerVector index = fr["index"];
        const NumericVector weights = fr["weights"];
        const NumericMatrix design = fr["design"];
        const NumericMatrix anchor_design = fr["anchor.design"];
        const int q = as<int>(fr["q"]);
        Eigen::MatrixXd cross = Eigen::MatrixXd::Zero(q, q);
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(q, ncol_y);
        for (int rr = 0; rr < index.size(); ++rr) {
            const int point = index[rr] - 1;
            const double w = weights[rr];
            if (!R_finite(w) || w <= 0.0) {
                continue;
            }
            for (int a = 0; a < q; ++a) {
                const double xa = design(rr, a);
                const double wxa = w * xa;
                for (int b = 0; b < q; ++b) {
                    cross(a, b) += wxa * design(rr, b);
                }
                for (int cc = 0; cc < ncol_y; ++cc) {
                    rhs(a, cc) += wxa * y_mat(point, cc);
                }
            }
        }
        double scale = 0.0;
        for (int a = 0; a < q; ++a) {
            scale = std::max(scale, cross(a, a));
        }
        if (!R_finite(scale) || scale <= 0.0) {
            scale = 1.0;
        }
        const double ridge = ridge_multiplier * scale;
        Eigen::MatrixXd normal = cross;
        for (int a = 0; a < q; ++a) {
            normal(a, a) += ridge;
        }
        Eigen::MatrixXd beta(q, ncol_y);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(normal);
        if (ldlt.info() == Eigen::Success) {
            beta = ldlt.solve(rhs);
        } else {
            beta = normal.colPivHouseholderQr().solve(rhs);
        }
        for (int cc = 0; cc < ncol_y; ++cc) {
            double value = 0.0;
            for (int a = 0; a < q; ++a) {
                value += anchor_design(0, a) * beta(a, cc);
            }
            fitted(ii, cc) = value;
        }
    }
    return fitted;
}

// [[Rcpp::export]]
List rcpp_ps_lps_assemble_cached_system(
    const List& cache,
    const NumericVector& y,
    const NumericVector& response_weights,
    const double lambda_sync) {

    const List data_blocks = cache["data.blocks"];
    const List sync_blocks = cache["sync.blocks"];
    const int ncoef = as<int>(cache["ncoef"]);

    int data_nrows = 0;
    int data_nnz = 0;
    for (int ii = 0; ii < data_blocks.size(); ++ii) {
        const List db = data_blocks[ii];
        const IntegerVector index = db["index"];
        const NumericVector weights = db["weights"];
        const int q = as<int>(db["q"]);
        int n_ok = 0;
        for (int aa = 0; aa < index.size(); ++aa) {
            const int point = index[aa] - 1;
            const double rw = response_weights[point];
            const double ww = rw * weights[aa];
            if (R_finite(rw) && R_finite(ww) && ww > 0.0) {
                ++n_ok;
            }
        }
        data_nrows += n_ok;
        data_nnz += n_ok * q;
    }

    int sync_nrows = 0;
    int sync_nnz = 0;
    if (lambda_sync > 0.0 && sync_blocks.size() > 0) {
        for (int ss = 0; ss < sync_blocks.size(); ++ss) {
            const List sb = sync_blocks[ss];
            const IntegerVector point = sb["point"];
            const int q_i = as<int>(sb["q.i"]);
            const int q_j = as<int>(sb["q.j"]);
            const int n_sb = point.size();
            sync_nrows += n_sb;
            sync_nnz += n_sb * (q_i + q_j);
        }
    }

    const int rr = data_nrows + sync_nrows;
    const int nnz = data_nnz + sync_nnz;
    if (rr <= 0 || nnz <= 0) {
        stop("PS-LPS cached native system has no rows.");
    }

    IntegerVector rows(nnz);
    IntegerVector cols(nnz);
    NumericVector vals(nnz);
    NumericVector rhs(rr);

    int row_pos = 0;
    int nz_pos = 0;
    for (int ii = 0; ii < data_blocks.size(); ++ii) {
        const List db = data_blocks[ii];
        const IntegerVector index = db["index"];
        const NumericVector weights = db["weights"];
        const NumericMatrix design = db["design"];
        const IntegerVector db_cols = db["cols"];
        const int q = as<int>(db["q"]);
        for (int aa = 0; aa < index.size(); ++aa) {
            const int point = index[aa] - 1;
            const double rw = response_weights[point];
            const double ww = rw * weights[aa];
            if (!(R_finite(rw) && R_finite(ww) && ww > 0.0)) {
                continue;
            }
            const double sw = std::sqrt(ww);
            ++row_pos;
            for (int qq = 0; qq < q; ++qq) {
                rows[nz_pos] = row_pos;
                cols[nz_pos] = db_cols[qq];
                vals[nz_pos] = sw * design(aa, qq);
                ++nz_pos;
            }
            rhs[row_pos - 1] = sw * y[point];
        }
    }

    if (lambda_sync > 0.0 && sync_blocks.size() > 0) {
        const double lambda_sqrt = std::sqrt(lambda_sync);
        for (int ss = 0; ss < sync_blocks.size(); ++ss) {
            const List sb = sync_blocks[ss];
            const IntegerVector point = sb["point"];
            const NumericVector omega_sqrt = sb["omega.sqrt"];
            const IntegerVector cols_i = sb["cols.i"];
            const IntegerVector cols_j = sb["cols.j"];
            const NumericMatrix values_i = sb["values.i"];
            const NumericMatrix values_j = sb["values.j"];
            const int q_i = as<int>(sb["q.i"]);
            const int q_j = as<int>(sb["q.j"]);
            const int n_sb = point.size();
            for (int aa = 0; aa < n_sb; ++aa) {
                const double scale = lambda_sqrt * omega_sqrt[aa];
                ++row_pos;
                for (int qq = 0; qq < q_i; ++qq) {
                    rows[nz_pos] = row_pos;
                    cols[nz_pos] = cols_i[qq];
                    vals[nz_pos] = scale * values_i(aa, qq);
                    ++nz_pos;
                }
                for (int qq = 0; qq < q_j; ++qq) {
                    rows[nz_pos] = row_pos;
                    cols[nz_pos] = cols_j[qq];
                    vals[nz_pos] = -scale * values_j(aa, qq);
                    ++nz_pos;
                }
            }
        }
    }

    if (row_pos != rr || nz_pos != nnz) {
        stop("Internal PS-LPS cached native sparse assembly count mismatch.");
    }

    return List::create(
        _["rows"] = rows,
        _["cols"] = cols,
        _["vals"] = vals,
        _["rhs"] = rhs,
        _["nrow"] = rr,
        _["ncol"] = ncoef,
        _["nnz"] = nnz,
        _["data_nrow"] = data_nrows,
        _["sync_nrow"] = sync_nrows
    );
}
