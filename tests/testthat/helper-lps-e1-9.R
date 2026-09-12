# E1.9 shared fixtures: deterministic pinning configurations for the
# bandwidth-multiplier gates.
#
# These constructors are the single source of truth for the b = 1
# backward-compatibility GATE's data and fit.lps argument lists. They are
# sourced by BOTH dev/methods/lps/ci/e1_9_pin_reference_fits.R (which generated the
# pre-change reference values in helper-lps-e1-9-reference.R at the commit
# recorded there) and tests/testthat/test-lps-bandwidth-multiplier.R (the
# GATE itself). Editing any constructor invalidates the pinned references;
# regenerate them from the pre-change commit if that ever becomes necessary.

# Plan G1 truth, D = 2, p = 2 (lps_experimental_plan_2026-06-09.tex, sec:dgp).
e19.g1.truth <- function(U) {
    0.5 + 1.0 * U[, 1L] - 0.7 * U[, 2L] + 0.4 * U[, 1L]^2 +
        0.3 * U[, 1L] * U[, 2L] - 0.6 * U[, 2L]^2
}

# Fixed ambient G1 dataset (D = 2) with one realized noise draw.
e19.pin.data.ambient <- function() {
    n <- 60L
    set.seed(4101)
    X <- matrix(stats::runif(n * 2L, -1, 1), ncol = 2L)
    truth <- e19.g1.truth(X)
    set.seed(4102)
    y <- truth + 0.1 * stats::rnorm(n)
    list(X = X, y = y, truth = truth)
}

# Fixed flat embedded dataset (intrinsic d = 2 in ambient D = 4) for the
# local-PCA code path.
e19.pin.data.embedded <- function() {
    n <- 50L
    set.seed(4103)
    U <- matrix(stats::runif(n * 2L, -1, 1), ncol = 2L)
    set.seed(4104)
    Q <- qr.Q(qr(matrix(stats::rnorm(4L * 2L), nrow = 4L)))[, 1:2,
                                                            drop = FALSE]
    X <- U %*% t(Q)
    truth <- e19.g1.truth(U)
    set.seed(4105)
    y <- truth + 0.05 * stats::rnorm(n)
    list(X = X, y = y, truth = truth)
}

# Configuration A: multi-candidate ambient selection (2 supports x 3 degrees
# x 4 kernels = 24 candidates) under the guarded default solve settings.
e19.fit.A <- function(...) {
    d <- e19.pin.data.ambient()
    fit.lps(
        X = d$X,
        y = d$y,
        foldid = rep(1:5, length.out = nrow(d$X)),
        support.grid = c(12L, 20L),
        degree.grid = 0:2,
        kernel.grid = c("gaussian", "tricube", "epanechnikov", "triangular"),
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10, 1e-8),
        ridge.condition.max = 1e12,
        unstable.action = "na",
        ...
    )
}

# Configuration B: strict single-candidate ambient fit (Tier-0 style:
# unguarded exact solve, compact kernel).
e19.fit.B <- function(...) {
    d <- e19.pin.data.ambient()
    fit.lps(
        X = d$X,
        y = d$y,
        foldid = rep(1:2, length.out = nrow(d$X)),
        support.grid = 18L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na",
        ...
    )
}

# Configuration C: single-candidate local-PCA chart fit on the flat
# embedded dataset.
e19.fit.C <- function(...) {
    d <- e19.pin.data.embedded()
    fit.lps(
        X = d$X,
        y = d$y,
        foldid = rep(1:5, length.out = nrow(d$X)),
        support.grid = 16L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 2L,
        local.chart.method = "pca",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na",
        ...
    )
}

# Deterministic fixed distance vector for the E1.9(a) characterization GATE:
# a 2-D K-NN-like profile (K = 20, max distance 1, no RNG).
e19.characterization.distances <- function() {
    sqrt(seq_len(20L) / 20L)
}
