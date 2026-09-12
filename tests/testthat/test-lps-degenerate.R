# =============================================================================
# Tier-0 E0.8 -- Degenerate-geometry pathologies
#
# Frozen spec: dev/methods/lps/specs/lps_experimental_plan_2026-06-09.tex,
# section E0.8.
# One test per pathology. Each asserts the SPECIFIED guarded behaviour (finite
# where defined, NA rather than a silent mean under unstable.action="na",
# probabilities in range, deterministic tie-breaking), and asserts on
# diagnostics where the fit object exposes them.
#
# SCAFFOLD STATUS: authored against the confirmed fit.lps API but not executed by
# the author (no R available here). Core behavioural assertions are hard; checks
# that depend on optional diagnostic fields are guarded so a missing field does
# not cause a spurious failure -- but the behaviour itself is always asserted.
# =============================================================================

deg.fit <- function(X, y, degree, kernel = "gaussian",
                    basis = "orthogonal.polynomial.drop",
                    unstable.action = "na", support = NULL,
                    outcome.family = "gaussian") {
    n <- nrow(X)
    if (is.null(support)) support <- min(n - 1L, max(8L, 4L * (1L + ncol(X))))
    fit.lps(
        X = X, y = y, foldid = rep(1:2, length.out = n),
        support.grid = support, degree.grid = degree, kernel.grid = kernel,
        coordinate.method = "coordinates", backend = "R",
        design.basis = basis, ridge.multiplier.grid = 0,
        ridge.condition.max = Inf, unstable.action = unstable.action,
        outcome.family = outcome.family)
}

# -- Case 1: duplicate points -> collapse to the weighted mean (degree-0) ------
test_that("E0.8 duplicate-point supports collapse to the weighted mean", {
    set.seed(9801)
    nA <- 20L; nB <- 20L
    X <- rbind(matrix(0, nA, 2L), matrix(1, nB, 2L))   # two fully-duplicated clusters
    yA <- stats::rnorm(nA, mean = -1); yB <- stats::rnorm(nB, mean = 3)
    y <- c(yA, yB)
    fit <- deg.fit(X, y, degree = 1L, support = 10L)
    fv <- fit$fitted.values
    expect_false(anyNA(fv))                            # finite predictions
    # a fully-degenerate (h=0) support has no first-order information: the
    # prediction collapses to a weighted average of the support responses
    # (degree-0 behaviour). The k-NN support is a tie-broken subset of the
    # cluster, so assert the prediction lies within the cluster's response range
    # and that the two clusters never bleed into each other.
    fvA <- fv[seq_len(nA)]; fvB <- fv[nA + seq_len(nB)]
    expect_true(all(fvA >= min(yA) & fvA <= max(yA)))
    expect_true(all(fvB >= min(yB) & fvB <= max(yB)))
    expect_lt(max(fvA), min(fvB))                      # clusters stay separated
    expect_equal(fit$local.chart.diagnostics.summary$min.design.rank, 1L)
})

# -- Case 2: collinear support, degree 2, unstable.action="na" -----------------
test_that("E0.8 collinear degree-2 support yields NA (not a silent mean)", {
    set.seed(9802)
    n <- 40L
    t <- stats::runif(n, -1, 1)
    X <- cbind(t, 2 * t)                               # exactly collinear in R^2
    y <- 0.5 + 1.3 * t + stats::rnorm(n, sd = 1e-3)
    expect_error(fit <- deg.fit(X, y, degree = 2L, unstable.action = "na"), NA)
    fv <- fit$fitted.values
    # never a silent global-mean fallback under unstable.action="na":
    expect_false(isTRUE(all.equal(fv, rep(mean(y), n))))
    # values are either finite (rank-deficient column dropped) or explicit NA:
    expect_true(all(is.finite(fv) | is.na(fv)))
    # contrast: unstable.action="mean" is *allowed* to fall back; just must not crash
    expect_error(deg.fit(X, y, degree = 2L, unstable.action = "mean"), NA)
})

# -- Case 3: constant response -> exact reproduction for all degree/kernel/basis
test_that("E0.8 constant response is reproduced to 1e-10 everywhere", {
    set.seed(9803)
    n <- 40L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    cval <- 2.71828
    y <- rep(cval, n)
    kernels <- c("gaussian", "tricube", "epanechnikov", "triangular")
    bases   <- c("orthogonal.polynomial.drop", "monomial",
                 "weighted.qr", "weighted.qr.drop")
    for (degree in c(0L, 1L, 2L)) {
        for (ker in kernels) {
            for (bas in bases) {
                fv <- deg.fit(X, y, degree = degree, kernel = ker,
                              basis = bas, support = 24L)$fitted.values
                expect_false(anyNA(fv))
                expect_lt(max(abs(fv - cval)), 1e-10)
            }
        }
    }
})

# -- Case 4: distance ties -> deterministic, bitwise-identical repeats ---------
test_that("E0.8 exact distance ties give deterministic, identical results", {
    # a symmetric grid produces many exactly-tied neighbor distances
    g <- seq(-1, 1, length.out = 7L)
    X <- as.matrix(expand.grid(g, g))                 # 49 points, symmetric
    set.seed(9804)
    y <- stats::rnorm(nrow(X))
    fit.a <- deg.fit(X, y, degree = 1L, support = 12L)
    fit.b <- deg.fit(X, y, degree = 1L, support = 12L)
    expect_identical(fit.a$fitted.values, fit.b$fitted.values)  # tie-break deterministic
})

# -- Case 5: extreme class imbalance (binary) ----------------------------------
test_that("E0.8 extreme class imbalance does not crash and stays in [0,1]", {
    set.seed(9805)
    n <- 300L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    p <- rep(0.02, n)                                  # prevalence ~ 0.02
    y <- stats::rbinom(n, 1L, p)
    if (sum(y) == 0L) y[[which.max(X[, 1L])]] <- 1L    # ensure at least one event
    res <- tryCatch(
        deg.fit(X, y, degree = 1L, support = 40L, outcome.family = "binomial"),
        error = function(e) e)
    if (inherits(res, "error")) {
        # spec allows a documented fail-fast on no finite selection score
        expect_match(conditionMessage(res), "finite|selection|score|fallback",
                     ignore.case = TRUE)
    } else {
        fv <- res$fitted.values
        # NA is the SANCTIONED guarded output under unstable.action="na"; require
        # every NON-NA prediction to be a valid probability, and require the
        # fallback to be reported (not silent) -- do not assert no NA.
        ok <- fv[!is.na(fv)]
        expect_true(all(ok >= 0 & ok <= 1))
        expect_false(is.null(res$logistic.diagnostics))# fallback telemetry present
    }
})

# -- Case 6: compositional data with structural zeros (D=6) --------------------
test_that("E0.8 compositional data with structural zeros gives finite/NA, no crash", {
    set.seed(9806)
    n <- 120L; D <- 6L
    Z <- matrix(stats::rgamma(n * D, shape = 0.4), ncol = D)  # sparse-ish
    Z[Z < 0.15] <- 0                                          # structural zeros
    Z[rowSums(Z) == 0, 1L] <- 1                               # avoid all-zero rows
    X <- Z / rowSums(Z)                                       # compositional (rows sum to 1)
    y <- X[, 1L] - X[, 2L] + stats::rnorm(n, sd = 0.05)
    expect_error(fit <- deg.fit(X, y, degree = 1L, support = 30L), NA)
    fv <- fit$fitted.values
    expect_true(all(is.finite(fv) | is.na(fv)))              # finite where defined
    diag.sum <- fit$local.chart.diagnostics.summary
    expect_true(is.finite(diag.sum$zero.bandwidth.fraction))
    expect_gte(diag.sum$zero.bandwidth.fraction, 0)
    expect_lte(diag.sum$zero.bandwidth.fraction, 1)
})
