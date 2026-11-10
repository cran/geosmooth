#include <Rcpp.h>
#include <ANN/ANN.h>
#include <dgraphs/linf_simplex_distance.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
using namespace Rcpp;

namespace {
struct Neighbor { double distance; int row; };
struct TrainingGroups {
    std::vector<std::vector<int> > rows;
    std::vector<int> discordant;
    int people = 0, cases = 0;
    TrainingGroups(const IntegerVector& train, const IntegerVector& group,
                   const NumericVector& y, int ng) : rows(ng + 1) {
        for (int j : train) rows[group[j]].push_back(j);
        for (int g = 1; g <= ng; ++g) if (!rows[g].empty()) {
            ++people;
            bool same = true;
            for (int j : rows[g]) if (y[j] != y[rows[g][0]]) same = false;
            if (same) cases += y[rows[g][0]] == 1;
            else discordant.push_back(g);
        }
    }
};
// ANN uses shared search scratch state: this backend is deliberately serial.
// Keep the shared trivial leaf alive until the last CV tree is destroyed.
int live_trees = 0;
struct NeighborAnnIndex {
    ANNpointArray points = nullptr;
    ANNkd_tree* tree = nullptr;
    std::vector<int> rows;
    TrainingGroups groups;
    int dimension;
    NeighborAnnIndex(const NumericMatrix& X, const IntegerVector& train,
                     const IntegerVector& priority, const IntegerVector& group,
                     const NumericVector& y, int ng)
        : rows(train.begin(), train.end()), groups(train, group, y, ng),
          dimension(X.ncol()) {
        std::sort(rows.begin(), rows.end(), [&](int a, int b) {
            return priority[a] < priority[b];
        });
        points = annAllocPts(rows.size(), dimension);
        try {
            for (size_t i = 0; i < rows.size(); ++i)
                for (int c = 0; c < dimension; ++c) points[i][c] = X(rows[i], c);
            tree = new ANNkd_tree(points, rows.size(), dimension);
            ++live_trees;
        } catch (...) {
            annDeallocPts(points);
            throw;
        }
    }
    ~NeighborAnnIndex() {
        delete tree;
        annDeallocPts(points);
        if (--live_trees == 0) annClose();
    }
};
dgraphs_simplex::Row surface_row(const NumericMatrix& X, int i) {
    dgraphs_simplex::Row row;
    double maximum = 0;
    for (int c = 0; c < X.ncol(); ++c) {
        const double value = X(i,c);
        if (!R_finite(value) || value < 0) stop("Invalid intrinsic surface coordinates.");
        maximum = std::max(maximum, value);
    }
    if (maximum != 1) stop("Native intrinsic coordinates must be max-normalized.");
    for (int c = 0; c < X.ncol(); ++c) if (X(i,c) > 0) {
        row.id.push_back(c);
        row.value.push_back(X(i,c));
    }
    return row;
}
struct NeighborIntrinsicIndex {
    std::vector<dgraphs_simplex::Row> points;
    std::vector<int> training;
    int dimension;
    NeighborIntrinsicIndex(const NumericMatrix& X, const IntegerVector& train)
        : points(X.nrow()), training(train.begin(), train.end()), dimension(X.ncol()) {
        for (int j : train) {
            if (j < 0 || j >= X.nrow()) stop("Invalid intrinsic training row.");
            if (j % 128 == 0) checkUserInterrupt();
            points[j] = surface_row(X, j);
        }
    }
};
void validate_training(const NumericVector& y, const IntegerVector& train,
                       const IntegerVector& group, const IntegerVector& priority) {
    int n = y.size();
    if (group.size() != n || priority.size() != n) stop("Invalid native training dimensions.");
    for (int j = 0; j < n; ++j)
        if (group[j] < 1 || group[j] > n || priority[j] < 0 || priority[j] >= n)
            stop("Invalid native group or priority.");
    for (int j : train)
        if (j < 0 || j >= n || !R_finite(y[j]) || (y[j] != 0 && y[j] != 1))
            stop("Native training rows must have binary observed labels.");
}
}

// [[Rcpp::export]]
SEXP rcpp_neighbor_ann_index(const NumericMatrix& X, const IntegerVector& train,
                             const IntegerVector& priority, const IntegerVector& group,
                             const NumericVector& y) {
    validate_training(y, train, group, priority);
    if (X.nrow() != y.size() || X.ncol() < 1 || train.size() < 1)
        stop("Invalid ANN training coordinates.");
    int ng = *std::max_element(group.begin(), group.end());
    XPtr<NeighborAnnIndex> ptr(new NeighborAnnIndex(X, train, priority, group, y, ng), true);
    return ptr;
}

// [[Rcpp::export]]
SEXP rcpp_neighbor_intrinsic_index(const NumericMatrix& X, const IntegerVector& train) {
    if (X.nrow() < 1 || X.ncol() < 1 || train.size() < 1)
        stop("Invalid intrinsic training dimensions.");
    XPtr<NeighborIntrinsicIndex> ptr(new NeighborIntrinsicIndex(X, train), true);
    return ptr;
}

// [[Rcpp::export]]
void rcpp_neighbor_release_index(SEXP neighbor_index, bool intrinsic) {
    if (Rf_isNull(neighbor_index)) return;
    if (intrinsic) XPtr<NeighborIntrinsicIndex>(neighbor_index).release();
    else XPtr<NeighborAnnIndex>(neighbor_index).release();
}

// One query neighborhood at a time; no all-pairs coordinate distance matrix.
// [[Rcpp::export]]
List rcpp_neighbor_probability(const NumericMatrix& X,
                               const NumericMatrix& query,
                               const NumericMatrix& distance,
                               const NumericVector& y,
                               const IntegerVector& train,
                               const IntegerVector& group,
                               const IntegerVector& query_group,
                               const IntegerVector& priority,
                               const IntegerVector& sizes,
                               int method, int metric, bool precomputed,
                               bool diagnostics, SEXP neighbor_index, double eps) {
    const int n = y.size(), nq = query_group.size(), nk = sizes.size();
    validate_training(y, train, group, priority);
    if (method < 0 || method > 2 || metric < 0 || metric > 2 || nk < 1 ||
        !R_finite(eps) || eps < 0 || !R_finite((1 + eps) * (1 + eps)))
        stop("Invalid native neighbor arguments.");
    if (precomputed) {
        if (distance.nrow() != nq || distance.ncol() != n)
            stop("Invalid query-by-training distance dimensions.");
    } else if (X.nrow() != n || query.nrow() != nq ||
               X.ncol() != query.ncol() || X.ncol() < 1)
        stop("Invalid coordinate dimensions.");
    for (int j = 0; j < nk; ++j)
        if (sizes[j] < (method == 2 ? 2 : 1) || (j > 0 && sizes[j] <= sizes[j-1]))
            stop("Native sizes must be increasing positive integers.");
    int ng = *std::max_element(group.begin(), group.end());
    for (int g : query_group) if (g < 0 || g > ng) stop("Invalid native query group.");
    NeighborAnnIndex* index = nullptr;
    NeighborIntrinsicIndex* surface = nullptr;
    if (!Rf_isNull(neighbor_index) && metric == 2 && !precomputed) {
        XPtr<NeighborIntrinsicIndex> ptr(neighbor_index);
        surface = ptr.get();
        if (!surface || surface->dimension != X.ncol() ||
            surface->points.size() != static_cast<size_t>(n) ||
            !std::equal(train.begin(), train.end(), surface->training.begin(), surface->training.end()))
            stop("Invalid intrinsic index.");
    } else if (!Rf_isNull(neighbor_index)) {
        XPtr<NeighborAnnIndex> ptr(neighbor_index);
        index = ptr.get();
        if (!index || precomputed || metric != 0 || index->dimension != X.ncol() ||
            index->rows.size() != static_cast<size_t>(train.size()))
            stop("Invalid ANN index.");
    }
    if (metric == 2 && !precomputed && train.size() && !surface)
        stop("Intrinsic coordinate search needs a prepared dgraphs index.");
    if (eps > 0 && !index && train.size()) stop("Positive eps requires ANN.");
    std::unique_ptr<TrainingGroups> brute_groups;
    if (!index) brute_groups.reset(new TrainingGroups(train, group, y, ng));
    const TrainingGroups& groups = index ? index->groups : *brute_groups;

    NumericMatrix probability(nq, nk);
    std::fill(probability.begin(), probability.end(), NA_REAL);
    IntegerMatrix neighbors(diagnostics ? nq : 0, diagnostics ? nk : 0);
    IntegerMatrix cases(diagnostics ? nq : 0, diagnostics ? nk : 0);
    NumericMatrix radius(diagnostics ? nq : 0, diagnostics ? nk : 0);
    IntegerMatrix ties(diagnostics ? nq : 0, diagnostics ? nk : 0);
    LogicalMatrix split(diagnostics ? nq : 0, diagnostics ? nk : 0);
    if (diagnostics) {
        std::fill(neighbors.begin(), neighbors.end(), NA_INTEGER);
        std::fill(cases.begin(), cases.end(), NA_INTEGER);
        std::fill(radius.begin(), radius.end(), NA_REAL);
        std::fill(ties.begin(), ties.end(), NA_INTEGER);
        std::fill(split.begin(), split.end(), NA_LOGICAL);
    }
    IntegerVector available(nq), positive(nq), searched(nq);
    std::vector<Neighbor> ordered;
    ordered.reserve(train.size());
    std::vector<int> seen(ng + 1);
    std::vector<double> center(X.ncol());
    const int target = sizes[nk - 1];
    for (int q = 0; q < nq; ++q) {
        checkUserInterrupt();
        const int excluded = query_group[q];
        dgraphs_simplex::Row surface_query;
        if (surface) surface_query = surface_row(query, q);
        auto dist = [&](int j) {
            double d = 0;
            if (precomputed) {
                d = distance(q, j);
                if (ISNAN(d) || d < 0) stop("Invalid native distance.");
            } else if (surface) {
                d = dgraphs_simplex::solve(surface_query, surface->points[j]).distance;
                if (!R_finite(d)) stop("Nonfinite intrinsic distance from dgraphs.");
            } else {
                for (int c = 0; c < X.ncol(); ++c) {
                    const double delta = query(q,c) - X(j,c);
                    if (metric == 0) d += delta * delta;
                    else d = std::max(d, std::abs(delta));
                }
                if (metric == 0) d = std::sqrt(d);
                if (!R_finite(d)) stop("Coordinate distance overflow.");
            }
            return d;
        };
        auto less = [&](const Neighbor& a, const Neighbor& b) {
            if (a.distance != b.distance) return a.distance < b.distance;
            return priority[a.row] < priority[b.row];
        };
        auto append = [&](int j) {
            if (excluded > 0 && group[j] == excluded) return;
            double d = dist(j);
            if (R_finite(d)) ordered.push_back({d,j});
        };
        auto sort_people = [&]() {
            std::sort(ordered.begin(), ordered.end(), less);
            std::fill(seen.begin(), seen.end(), 0);
            size_t keep = 0;
            for (size_t j = 0; j < ordered.size(); ++j) {
                const int g = group[ordered[j].row];
                if (!seen[g]) { seen[g] = 1; ordered[keep++] = ordered[j]; }
            }
            ordered.resize(keep);
        };
        auto stopping_rank = [&]() {
            int count = 0;
            for (size_t j = 0; j < ordered.size(); ++j)
                if ((method == 0 || y[ordered[j].row] == 1) && ++count == target)
                    return static_cast<int>(j + 1);
            return 0;
        };
        ordered.clear();
        if (!index) {
            int evaluated = 0;
            for (int j : train) {
                if ((++evaluated % 128) == 0) checkUserInterrupt();
                append(j);
            }
            searched[q] = ordered.size();
            sort_people();
        } else {
            for (int c = 0; c < X.ncol(); ++c) center[c] = query(q,c);
            const int nt = index->rows.size();
            int request = std::min(nt, std::max(64, target < nt ? target + 1 : nt));
            while (true) {
                checkUserInterrupt();
                ordered.clear();
                if (request == nt) {
                    for (int j : index->rows) append(j);
                    searched[q] = ordered.size();
                    sort_people();
                    break;
                }
                std::vector<ANNidx> ids(request);
                std::vector<ANNdist> ds(request);
                index->tree->annkSearch(center.data(), request, ids.data(), ds.data(), eps);
                double farthest = 0;
                std::fill(seen.begin(), seen.end(), 0);
                int candidates = 0;
                for (int i : ids) {
                    if (i < 0 || i >= nt) stop("ANN returned an invalid neighbor.");
                    int j = index->rows[i];
                    farthest = std::max(farthest, dist(j));
                    if (eps > 0) {
                        // Preserve closest-answering-visit semantics even with ANN omissions.
                        int g = group[j];
                        if (seen[g]) continue;
                        seen[g] = 1;
                        for (int v : groups.rows[g]) { append(v); ++candidates; }
                    } else { append(j); ++candidates; }
                }
                searched[q] = candidates;
                sort_people();
                int stop_at = stopping_rank();
                if (stop_at && eps == 0) {
                    const double cutoff = ordered[stop_at - 1].distance;
                    // Include diagnostic near-ties and guard floating point pruning/sqrt.
                    const double pad = (diagnostics ? 1e-12 : 0) +
                        64 * std::numeric_limits<double>::epsilon() * X.ncol() *
                        std::max(cutoff, std::numeric_limits<double>::min());
                    double bound = std::nextafter(cutoff + pad,
                                                  std::numeric_limits<double>::infinity());
                    if (farthest <= bound) {
                        double squared = std::nextafter(bound * bound,
                                                        std::numeric_limits<double>::infinity());
                        int count = index->tree->annkFRSearch(center.data(), squared, 0,
                                                             nullptr, nullptr, 0.0);
                        std::vector<ANNidx> shell(count);
                        std::vector<ANNdist> shell_dist(count);
                        index->tree->annkFRSearch(center.data(), squared, count,
                                                 shell.data(), shell_dist.data(), 0.0);
                        ordered.clear();
                        for (int i : shell) {
                            if (i < 0 || i >= nt) stop("ANN returned an invalid radius neighbor.");
                            append(index->rows[i]);
                        }
                        searched[q] = count;
                        sort_people();
                        stop_at = stopping_rank();
                    }
                }
                if (stop_at) break;
                request = request > nt / 2 ? nt : request * 2;
            }
        }
        std::vector<int> cumulative(ordered.size()), stops;
        int count = 0;
        for (size_t j = 0; j < ordered.size(); ++j) {
            if (y[ordered[j].row] == 1) { ++count; stops.push_back(j+1); }
            cumulative[j] = count;
        }
        if (!index) {
            available[q] = ordered.size();
            positive[q] = count;
        } else {
            available[q] = groups.people;
            positive[q] = groups.cases;
            for (int g : groups.discordant) {
                Neighbor best{std::numeric_limits<double>::infinity(), groups.rows[g][0]};
                for (int j : groups.rows[g]) {
                    Neighbor candidate{dist(j), j};
                    if (less(candidate, best)) best = candidate;
                }
                positive[q] += y[best.row] == 1;
            }
            if (excluded > 0 && !groups.rows[excluded].empty()) {
                --available[q];
                Neighbor best{std::numeric_limits<double>::infinity(), groups.rows[excluded][0]};
                for (int j : groups.rows[excluded]) {
                    Neighbor candidate{dist(j), j};
                    if (less(candidate, best)) best = candidate;
                }
                positive[q] -= y[best.row] == 1;
            }
        }
        for (int z = 0; z < nk; ++z) {
            const int k = sizes[z];
            if ((method == 0 && static_cast<size_t>(k) > ordered.size()) ||
                (method != 0 && k > count)) continue;
            const int t = method == 0 ? k : stops[k-1];
            probability(q,z) = method == 2 ? double(k-1)/(t-1) : double(cumulative[t-1])/t;
            if (!diagnostics) continue;
            neighbors(q,z) = t;
            cases(q,z) = cumulative[t-1];
            const double rad = ordered[t-1].distance;
            radius(q,z) = rad;
            int shell = 0; bool outside = false;
            for (size_t j = 0; j < ordered.size(); ++j)
                if (std::abs(ordered[j].distance - rad) <= 1e-12) {
                    ++shell; if (int(j) >= t) outside = true;
                }
            ties(q,z) = shell;
            split(q,z) = outside;
        }
    }
    return List::create(_["probability"] = probability,
                        _["neighbors"] = neighbors, _["cases"] = cases,
                        _["radius"] = radius, _["boundary.ties"] = ties,
                        _["split.tie"] = split, _["available.people"] = available,
                        _["available.cases"] = positive, _["search.candidates"] = searched);
}
