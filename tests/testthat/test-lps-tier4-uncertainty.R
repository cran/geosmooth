# E4.1 Part A — deterministic unit GATE for the pointwise-variance/band
# machinery (Tier 4, uncertainty). No DGP-library dependency: fixtures are
# seeded deterministic point clouds.
#
# The implemented route (lps.smoother.matrix / lps.pointwise.band) extracts S
# analytically from the local WLS algebra. The reference route below extracts
# S independently, column by column, through the public API as fit(e_j) — the
# E0.2 protocol (linearity makes finite differencing exact). The GATE asserts
# the two routes agree to the program's algebraic tolerance 1e-10, per point,
# and that df = tr S, sigma.hat^2 = RSS / (n - tr S), and the band endpoints
# match their defining formulas computed from the independent S.

e41.tol <- 1e-10

e41.fixed.lps.fit <- function(X, y, coordinate.method = "coordinates",
                              chart.dim = NULL, support.size = 18L,
                              X.eval = NULL) {
    fit.lps(
        X = X,
        y = y,
        foldid = rep(1:2, length.out = nrow(X)),
        support.grid = support.size,
        degree.grid = 1L,
        kernel.grid = "tricube",
        X.eval = X.eval,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        local.chart.method = "pca",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
}

# Independent S extraction: column j is the fitted vector of the j-th unit
# response through the full public pipeline (E0.2 protocol).
e41.probe.smoother.matrix <- function(X, coordinate.method = "coordinates",
                                      chart.dim = NULL,
                                      support.size = 18L) {
    n <- nrow(X)
    S <- matrix(NA_real_, nrow = n, ncol = n)
    for (j in seq_len(n)) {
        e.j <- numeric(n)
        e.j[[j]] <- 1
        fit <- e41.fixed.lps.fit(
            X = X,
            y = e.j,
            coordinate.method = coordinate.method,
            chart.dim = chart.dim,
            support.size = support.size
        )
        S[, j] <- fit$fitted.values
    }
    S
}

# Seeded quadratic-surface fixture in D = 3 with intrinsic dimension 2 (a
# deterministic unit fixture for the chart.dim = 2 local-PCA configuration;
# not a DGP-library generator and consumed by no STUDY).
e41.surface.fixture <- function(n, seed) {
    set.seed(seed)
    u <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    cbind(
        u[, 1L],
        u[, 2L],
        0.25 * u[, 1L]^2 + 0.15 * u[, 2L]^2 - 0.1 * u[, 1L] * u[, 2L]
    )
}

e41.check.band.against.probe <- function(fit, S.probe, sigma0) {
    n <- nrow(S.probe)
    y <- fit$y
    raw <- fit$fitted.values.raw
    expect_false(anyNA(fit$fitted.values))
    expect_equal(nrow(fit$cv.table), 1L)

    ## The probe S satisfies the linear-smoother identity on this fit.
    expect_lt(max(abs(as.numeric(S.probe %*% y) - raw)), e41.tol)

    ## Implemented analytic S agrees with the independent probe S entrywise.
    S.impl <- lps.smoother.matrix(fit)
    expect_lt(max(abs(S.impl - S.probe)), e41.tol)

    ## Known-sigma variance: sigma0^2 * sum_j S_ij^2 from the probe S.
    band <- lps.pointwise.band(fit, sigma = sigma0)
    var.ref <- sigma0^2 * rowSums(S.probe^2)
    expect_lt(max(abs(band$variance - var.ref)), e41.tol)
    expect_lt(max(abs(band$se - sqrt(var.ref))), e41.tol)
    expect_identical(band$sigma.source, "known")
    expect_identical(band$sigma, sigma0)
    expect_identical(band$fitted, raw)

    ## df = tr S from the probe S.
    expect_lt(abs(band$df - sum(diag(S.probe))), e41.tol)

    ## Band endpoints: fitted +/- z_{0.975} * se.
    z <- stats::qnorm(0.975)
    expect_lt(max(abs(band$lower - (raw - z * band$se))), e41.tol)
    expect_lt(max(abs(band$upper - (raw + z * band$se))), e41.tol)
    expect_identical(band$z, z)
    expect_identical(band$level, 0.95)

    ## Plug-in variant: sigma.hat^2 = RSS / (n - tr S), all from the probe S.
    band.plugin <- lps.pointwise.band(fit)
    df.ref <- sum(diag(S.probe))
    rss.ref <- sum((y - raw)^2)
    expect_lt(abs(band.plugin$rss - rss.ref), e41.tol)
    expect_lt(abs(band.plugin$sigma.hat^2 - rss.ref / (n - df.ref)), e41.tol)
    expect_lt(
        max(abs(band.plugin$variance -
                    band.plugin$sigma.hat^2 * rowSums(S.probe^2))),
        e41.tol
    )
    expect_identical(band.plugin$sigma.source, "plug.in")
    expect_identical(band.plugin$sigma, NA_real_)

    ## sigma.hat is reported (identically) in known mode too.
    expect_lt(abs(band$sigma.hat - band.plugin$sigma.hat), e41.tol)

    invisible(NULL)
}

test_that("E4.1 variance/band machinery matches an independently extracted S (ambient coordinates)", {
    set.seed(9401)
    n <- 42L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    y <- stats::rnorm(n)

    fit <- e41.fixed.lps.fit(X, y)
    S.probe <- e41.probe.smoother.matrix(X)
    e41.check.band.against.probe(fit, S.probe, sigma0 = 0.37)
})

test_that("E4.1 variance/band machinery matches an independently extracted S (local PCA chart, chart.dim = 2)", {
    n <- 34L
    X <- e41.surface.fixture(n, seed = 9402)
    set.seed(9403)
    y <- stats::rnorm(n)

    fit <- e41.fixed.lps.fit(
        X, y,
        coordinate.method = "local.pca",
        chart.dim = 2L
    )
    S.probe <- e41.probe.smoother.matrix(
        X,
        coordinate.method = "local.pca",
        chart.dim = 2L
    )
    e41.check.band.against.probe(fit, S.probe, sigma0 = 0.37)
})

test_that("E4.1 smoother extraction supports rectangular X.eval while the band requires the square training fit", {
    set.seed(9404)
    n <- 36L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    y <- stats::rnorm(n)
    X.eval <- X[seq_len(11L), , drop = FALSE]

    fit <- e41.fixed.lps.fit(X, y, X.eval = X.eval)
    S <- lps.smoother.matrix(fit)
    expect_identical(dim(S), c(11L, n))
    expect_lt(max(abs(as.numeric(S %*% y) - fit$fitted.values.raw)), e41.tol)
    expect_error(lps.pointwise.band(fit), "X.eval identical to X")
})

test_that("E4.1 extraction reproduces rank-dropped local fits (support below design size)", {
    set.seed(9409)
    n <- 24L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    y <- stats::rnorm(n)

    ## support.size = 2 < ncol(design) = 3 at degree 1 in D = 2: the
    ## orthogonal.polynomial.drop basis drops rank-deficient directions and
    ## still fits. The extraction must reproduce that path exactly; its
    ## internal self-guard (max |S y - fitted| <= 1e-10) enforces it.
    fit.drop <- e41.fixed.lps.fit(X, y, support.size = 2L)
    expect_false(anyNA(fit.drop$fitted.values))
    S.drop <- lps.smoother.matrix(fit.drop)
    expect_lt(max(abs(as.numeric(S.drop %*% y) - fit.drop$fitted.values.raw)),
              e41.tol)
})

test_that("E4.1 machinery rejects configurations outside the linear-smoother premise", {
    set.seed(9405)
    n <- 40L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    y <- stats::rnorm(n)

    ## Non-singleton grids: data-driven selection breaks linearity.
    fit.grid <- fit.lps(
        X, y,
        foldid = rep(1:2, length.out = n),
        support.grid = c(12L, 18L),
        degree.grid = 1L,
        kernel.grid = "tricube",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
    expect_error(lps.smoother.matrix(fit.grid), "singleton")
    expect_error(lps.pointwise.band(fit.grid), "singleton")

    ## Auto chart dimension: a separate map, excluded by the spec.
    X3 <- e41.surface.fixture(n, seed = 9406)
    set.seed(9407)
    y3 <- stats::rnorm(n)
    fit.auto <- e41.fixed.lps.fit(
        X3, y3,
        coordinate.method = "local.pca",
        chart.dim = "auto"
    )
    expect_error(lps.smoother.matrix(fit.auto), "chart")

    ## local.pca with the NULL default chart dimension: not explicit.
    fit.null.dim <- e41.fixed.lps.fit(
        X3, y3,
        coordinate.method = "local.pca",
        chart.dim = NULL
    )
    expect_error(lps.smoother.matrix(fit.null.dim), "chart")

    ## Non-gaussian outcome family.
    y.bin <- as.numeric(stats::runif(n) < 0.5)
    fit.bern <- fit.lps(
        X, y.bin,
        foldid = rep(1:2, length.out = n),
        support.grid = 18L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na",
        outcome.family = "bernoulli"
    )
    expect_error(lps.smoother.matrix(fit.bern), "gaussian")

    ## Unsupported design basis.
    fit.monomial <- fit.lps(
        X, y,
        foldid = rep(1:2, length.out = n),
        support.grid = 18L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        backend = "R",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
    expect_error(lps.smoother.matrix(fit.monomial),
                 "orthogonal.polynomial.drop")

    ## Band argument validation on a valid fit.
    fit.ok <- e41.fixed.lps.fit(X, y)
    expect_error(lps.pointwise.band(fit.ok, sigma = -1), "positive")
    expect_error(lps.pointwise.band(fit.ok, level = 1.2), "level")
})
