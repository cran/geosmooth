#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

std::vector<std::vector<int>> parse_adj_list(const List& adj_list) {
    const int n = adj_list.size();
    std::vector<std::vector<int>> adj(static_cast<size_t>(n));

    for (int i = 0; i < n; ++i) {
        IntegerVector neighbors(adj_list[i]);
        adj[static_cast<size_t>(i)].reserve(static_cast<size_t>(neighbors.size()));
        for (int k = 0; k < neighbors.size(); ++k) {
            const int v = neighbors[k];
            if (v == NA_INTEGER || v < 0 || v >= n) {
                Rcpp::stop("adj.list contains an out-of-range 0-based index.");
            }
            adj[static_cast<size_t>(i)].push_back(v);
        }
    }

    return adj;
}

std::vector<std::vector<double>> parse_weight_list(const List& weight_list, int n) {
    if (weight_list.size() != n) {
        Rcpp::stop("weight.list must have the same length as adj.list.");
    }

    std::vector<std::vector<double>> weights(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        NumericVector w(weight_list[i]);
        weights[static_cast<size_t>(i)].reserve(static_cast<size_t>(w.size()));
        for (int k = 0; k < w.size(); ++k) {
            const double wk = w[k];
            if (!std::isfinite(wk) || wk <= 0.0) {
                Rcpp::stop("weight.list entries must be finite positive values.");
            }
            weights[static_cast<size_t>(i)].push_back(wk);
        }
    }

    return weights;
}

std::unordered_set<int> parse_region_vertices(const IntegerVector& region_vertices, int n) {
    std::unordered_set<int> region;
    region.reserve(static_cast<size_t>(region_vertices.size()));

    for (int i = 0; i < region_vertices.size(); ++i) {
        const int v_r = region_vertices[i];
        if (v_r == NA_INTEGER || v_r < 1 || v_r > n) {
            Rcpp::stop("region.vertices contains an out-of-range 1-based index.");
        }
        region.insert(v_r - 1);
    }

    return region;
}

IntegerMatrix extrema_matrix(const std::vector<std::vector<int>>& adj,
                             const std::vector<double>& values) {
    std::vector<int> vertices;
    std::vector<int> is_maximum;

    const int n = static_cast<int>(adj.size());
    for (int i = 0; i < n; ++i) {
        if (adj[static_cast<size_t>(i)].empty()) {
            continue;
        }
        bool is_min = true;
        bool is_max = true;
        bool strict_min = false;
        bool strict_max = false;
        const double yi = values[static_cast<size_t>(i)];
        for (const int j : adj[static_cast<size_t>(i)]) {
            const double yj = values[static_cast<size_t>(j)];
            if (yi > yj) {
                is_min = false;
                strict_max = true;
            } else if (yi < yj) {
                is_max = false;
                strict_min = true;
            }
        }
        if (is_min && strict_min) {
            vertices.push_back(i + 1);
            is_maximum.push_back(0);
        }
        if (is_max && strict_max) {
            vertices.push_back(i + 1);
            is_maximum.push_back(1);
        }
    }

    IntegerMatrix out(static_cast<int>(vertices.size()), 2);
    for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
        out(i, 0) = vertices[static_cast<size_t>(i)];
        out(i, 1) = is_maximum[static_cast<size_t>(i)];
    }
    Rcpp::CharacterVector names = Rcpp::CharacterVector::create("evertex", "is_max");
    out.attr("dimnames") = Rcpp::List::create(R_NilValue, names);
    return out;
}

double extrema_difference(const IntegerMatrix& a, const IntegerMatrix& b) {
    const int na = a.nrow();
    const int nb = b.nrow();
    if (na == 0 && nb == 0) {
        return 0.0;
    }
    if (na == 0 || nb == 0) {
        return 1.0;
    }

    std::set<std::pair<int, int>> aa;
    std::set<std::pair<int, int>> bb;
    for (int i = 0; i < na; ++i) {
        aa.insert(std::make_pair(a(i, 0), a(i, 1)));
    }
    for (int i = 0; i < nb; ++i) {
        bb.insert(std::make_pair(b(i, 0), b(i, 1)));
    }

    int common = 0;
    for (const auto& item : aa) {
        if (bb.find(item) != bb.end()) {
            ++common;
        }
    }
    return 1.0 - (2.0 * static_cast<double>(common)) /
        static_cast<double>(aa.size() + bb.size());
}

bool last_window_is_stable(const std::vector<double>& diffs,
                           int stability_window,
                           double stability_threshold) {
    if (stability_window <= 0 ||
        static_cast<int>(diffs.size()) < stability_window) {
        return false;
    }
    const int start = static_cast<int>(diffs.size()) - stability_window;
    for (int i = start; i < static_cast<int>(diffs.size()); ++i) {
        if (!std::isfinite(diffs[static_cast<size_t>(i)]) ||
            diffs[static_cast<size_t>(i)] > stability_threshold) {
            return false;
        }
    }
    return true;
}

} // namespace

//' Native graph harmonic smoothing backend
//'
//' Internal C++ backend for `perform.harmonic.smoothing()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List rcpp_perform_harmonic_smoothing(
        const List& adj_list,
        const List& weight_list,
        const NumericVector& values,
        const IntegerVector& region_vertices,
        const int max_iterations,
        const double tolerance) {
    const int n = adj_list.size();
    if (values.size() != n) {
        Rcpp::stop("values must have the same length as adj.list.");
    }
    if (max_iterations <= 0) {
        Rcpp::stop("max.iterations must be a positive integer.");
    }
    if (!std::isfinite(tolerance) || tolerance <= 0.0) {
        Rcpp::stop("tolerance must be a finite positive scalar.");
    }

    std::vector<std::vector<int>> adj = parse_adj_list(adj_list);
    std::vector<std::vector<double>> weights = parse_weight_list(weight_list, n);
    for (int i = 0; i < n; ++i) {
        if (adj[static_cast<size_t>(i)].size() != weights[static_cast<size_t>(i)].size()) {
            Rcpp::stop("adj.list and weight.list entries must have matching lengths.");
        }
    }

    std::vector<double> f(values.begin(), values.end());
    std::unordered_set<int> region = parse_region_vertices(region_vertices, n);

    std::unordered_set<int> boundary;
    boundary.reserve(region.size());
    for (const int v : region) {
        const std::vector<int>& nbrs = adj[static_cast<size_t>(v)];
        for (const int u : nbrs) {
            if (region.count(u) == 0) {
                boundary.insert(v);
                break;
            }
        }
    }

    std::unordered_set<int> interior;
    interior.reserve(region.size());
    for (const int v : region) {
        if (boundary.count(v) == 0) {
            interior.insert(v);
        }
    }

    std::vector<double> max_change;
    std::vector<double> max_residual;
    max_change.reserve(static_cast<size_t>(max_iterations));
    max_residual.reserve(static_cast<size_t>(max_iterations));

    bool converged = false;
    int num_iterations = 0;

    if (!boundary.empty() && !interior.empty()) {
        const double eps = 1e-10;
        std::vector<double> next_f = f;

        auto weighted_average = [&](int v, const std::vector<double>& current, double& wsum) {
            double sum = 0.0;
            wsum = 0.0;
            const std::vector<int>& nbrs = adj[static_cast<size_t>(v)];
            const std::vector<double>& ws = weights[static_cast<size_t>(v)];
            for (size_t k = 0; k < nbrs.size(); ++k) {
                const int u = nbrs[k];
                if (region.count(u) == 0) {
                    continue;
                }
                const double conductance = 1.0 / (ws[k] + eps);
                sum += conductance * current[static_cast<size_t>(u)];
                wsum += conductance;
            }
            return wsum > 0.0 ? sum / wsum : current[static_cast<size_t>(v)];
        };

        for (int iter = 0; iter < max_iterations; ++iter) {
            double iter_change = 0.0;
            for (const int v : interior) {
                double wsum = 0.0;
                const double proposed = weighted_average(v, f, wsum);
                const double change = std::abs(proposed - f[static_cast<size_t>(v)]);
                iter_change = std::max(iter_change, change);
                next_f[static_cast<size_t>(v)] = proposed;
            }

            for (const int v : interior) {
                f[static_cast<size_t>(v)] = next_f[static_cast<size_t>(v)];
            }

            double iter_residual = 0.0;
            for (const int v : interior) {
                double wsum = 0.0;
                const double avg = weighted_average(v, f, wsum);
                if (wsum > 0.0) {
                    iter_residual = std::max(
                        iter_residual,
                        std::abs(f[static_cast<size_t>(v)] - avg)
                    );
                }
            }

            max_change.push_back(iter_change);
            max_residual.push_back(iter_residual);
            num_iterations = iter + 1;

            if (iter_change < tolerance && iter_residual < tolerance) {
                converged = true;
                break;
            }
        }
    }

    NumericVector harmonic_predictions(n);
    for (int i = 0; i < n; ++i) {
        harmonic_predictions[i] = f[static_cast<size_t>(i)];
    }

    return List::create(
        Named("harmonic_predictions") = harmonic_predictions,
        Named("converged") = converged,
        Named("num_region") = static_cast<int>(region.size()),
        Named("num_boundary") = static_cast<int>(boundary.size()),
        Named("num_interior") = static_cast<int>(interior.size()),
        Named("num_iterations") = num_iterations,
        Named("max_change") = NumericVector(max_change.begin(), max_change.end()),
        Named("max_residual") = NumericVector(max_residual.begin(), max_residual.end())
    );
}

//' Native graph harmonic smoothing backend with topology tracking
//'
//' Internal C++ backend for `harmonic.smoother()`.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List rcpp_harmonic_smoother(
        const List& adj_list,
        const List& weight_list,
        const NumericVector& values,
        const IntegerVector& region_vertices,
        const int max_iterations,
        const double tolerance,
        const int record_frequency,
        const int stability_window,
        const double stability_threshold) {
    const int n = adj_list.size();
    if (values.size() != n) {
        Rcpp::stop("values must have the same length as adj.list.");
    }
    if (max_iterations <= 0) {
        Rcpp::stop("max.iterations must be a positive integer.");
    }
    if (!std::isfinite(tolerance) || tolerance <= 0.0) {
        Rcpp::stop("tolerance must be a finite positive scalar.");
    }
    if (record_frequency <= 0) {
        Rcpp::stop("record.frequency must be a positive integer.");
    }
    if (stability_window <= 0) {
        Rcpp::stop("stability.window must be a positive integer.");
    }
    if (!std::isfinite(stability_threshold) ||
        stability_threshold < 0.0 ||
        stability_threshold > 1.0) {
        Rcpp::stop("stability.threshold must be in [0, 1].");
    }

    std::vector<std::vector<int>> adj = parse_adj_list(adj_list);
    std::vector<std::vector<double>> weights = parse_weight_list(weight_list, n);
    for (int i = 0; i < n; ++i) {
        if (adj[static_cast<size_t>(i)].size() !=
            weights[static_cast<size_t>(i)].size()) {
            Rcpp::stop("adj.list and weight.list entries must have matching lengths.");
        }
    }

    std::vector<double> f(values.begin(), values.end());
    std::unordered_set<int> region = parse_region_vertices(region_vertices, n);

    std::unordered_set<int> boundary;
    boundary.reserve(region.size());
    for (const int v : region) {
        for (const int u : adj[static_cast<size_t>(v)]) {
            if (region.count(u) == 0) {
                boundary.insert(v);
                break;
            }
        }
    }

    std::vector<int> interior;
    interior.reserve(region.size());
    for (const int v : region) {
        if (boundary.count(v) == 0) {
            interior.push_back(v);
        }
    }

    std::vector<std::vector<double>> prediction_history;
    std::vector<IntegerMatrix> basin_history;
    std::vector<double> topology_differences;
    prediction_history.push_back(f);
    basin_history.push_back(extrema_matrix(adj, f));

    bool converged = interior.empty();
    int stable_iteration = max_iterations;
    bool stability_seen = false;
    const double eps = 1e-10;
    std::vector<double> next_f = f;

    for (int iter = 0; iter < max_iterations && !converged; ++iter) {
        double iter_change = 0.0;
        for (const int v : interior) {
            double sum = 0.0;
            double wsum = 0.0;
            const std::vector<int>& nbrs = adj[static_cast<size_t>(v)];
            const std::vector<double>& ws = weights[static_cast<size_t>(v)];
            for (size_t k = 0; k < nbrs.size(); ++k) {
                const double conductance = 1.0 / (ws[k] + eps);
                sum += conductance * f[static_cast<size_t>(nbrs[k])];
                wsum += conductance;
            }
            const double proposed = wsum > 0.0 ?
                sum / wsum :
                f[static_cast<size_t>(v)];
            iter_change = std::max(
                iter_change,
                std::abs(proposed - f[static_cast<size_t>(v)])
            );
            next_f[static_cast<size_t>(v)] = proposed;
        }

        for (const int v : interior) {
            f[static_cast<size_t>(v)] = next_f[static_cast<size_t>(v)];
        }

        converged = iter_change <= tolerance;
        if (((iter + 1) % record_frequency == 0) || converged) {
            prediction_history.push_back(f);
            IntegerMatrix current_basins = extrema_matrix(adj, f);
            topology_differences.push_back(
                extrema_difference(basin_history.back(), current_basins)
            );
            basin_history.push_back(current_basins);
            if (!stability_seen &&
                last_window_is_stable(
                    topology_differences,
                    stability_window,
                    stability_threshold
                )) {
                stable_iteration = iter + 1;
                stability_seen = true;
            }
        }
    }

    if (!stability_seen && converged) {
        stable_iteration = static_cast<int>(
            std::min<size_t>(
                static_cast<size_t>(max_iterations),
                topology_differences.size() * static_cast<size_t>(record_frequency)
            )
        );
    }

    NumericVector harmonic_predictions(n);
    for (int i = 0; i < n; ++i) {
        harmonic_predictions[i] = f[static_cast<size_t>(i)];
    }

    const int n_records = static_cast<int>(prediction_history.size());
    NumericMatrix i_harmonic_predictions(n, n_records);
    for (int col = 0; col < n_records; ++col) {
        for (int row = 0; row < n; ++row) {
            i_harmonic_predictions(row, col) =
                prediction_history[static_cast<size_t>(col)][static_cast<size_t>(row)];
        }
    }

    List i_basins(static_cast<int>(basin_history.size()));
    for (int i = 0; i < static_cast<int>(basin_history.size()); ++i) {
        i_basins[i] = basin_history[static_cast<size_t>(i)];
    }

    NumericVector diffs(
        topology_differences.begin(),
        topology_differences.end()
    );

    return List::create(
        Named("harmonic_predictions") = harmonic_predictions,
        Named("i_harmonic_predictions") = i_harmonic_predictions,
        Named("i_basins") = i_basins,
        Named("stable_iteration") = stable_iteration,
        Named("topology_differences") = diffs,
        Named("basin_cx_differences") = diffs,
        Named("converged") = converged,
        Named("num_region") = static_cast<int>(region.size()),
        Named("num_boundary") = static_cast<int>(boundary.size()),
        Named("num_interior") = static_cast<int>(interior.size())
    );
}
