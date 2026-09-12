# E2.13 — ridge-penalty structure alignment (Tier 2 GATEs).
#
# Contract S C / E2.13 + the S G4 resolution
# (dev/methods/lps/audit_contracts/lps_g4_ridge_resolution_2026-06-12.md):
# the aligned ridge
# is OPT-IN via ridge.shrinkage.target = "local.mean" (weighted-centering
# reparametrization, gaussian/least-squares solve only); the default
# "zero" preserves the historical penalty structure bit-for-bit.
#
# GATEs in this file:
#   (a) aligned-mode GATE  -- with "local.mean", design.basis =
#       "orthogonal.polynomial.drop", ridge.condition.max = Inf, singleton
#       grids, rho in {0, 1e-8, 1e-2, 1, 1e2}, paired across rho on the
#       same data: |f(rho=1e2) - ybar_w| is small relative to
#       |f(rho=1e2) - 0| (shrinks to the LOCAL WEIGHTED MEAN, not 0), and
#       |f(rho=1e-8) - f(rho=0)| < 1e-6 (tiny-ridge invariance).
#   (b) legacy regression test -- documents the shrink-to-ZERO behavior of
#       the default "zero" arm at rho = 1e2, so the contrast is intentional
#       and visible.
#   (c) S A2 backward-compat GATE -- default-configuration fits (default
#       ridge grid c(0, 1e-10, 1e-8), default basis, no
#       ridge.shrinkage.target supplied) reproduce the committed
#       pre-change reference values EXACTLY (17-significant-digit round
#       trip; reference generated with R/lps.R in its pre-E2.13 state at
#       HEAD c621e2f and committed in c796408, by
#       dev/methods/lps/ci/e2_13_pin_reference_fits.R into
#       dev/methods/lps/runs/tier2/e2_13_reference_fits.csv). This is the pin
#       that protects
#       Tier-0 and E1.9.
#
# DGP (plan S E2.13): G1, D = 2, n = 150, mild noise (sd 0.1), deterministic
# seed. Realized quantities at implementation (this machine, vecLib BLAS):
# tiny-ridge invariance 6.95e-09 (threshold 1e-6, ~144x margin); at
# rho = 1e2 the median |f - ybar_w| / |f - 0| ratio is 0.0018 (threshold
# 0.05, ~28x) and the max over anchors with |ybar_w| > 0.2 is 0.0114
# (threshold 0.1, ~8.7x); the legacy arm's mean |f(1e2)| is 0.0060 against
# mean |ybar_w| = 0.50 (threshold ratio 0.1, ~8x).

e213.g1 <- function(n = 150L, seed = 1301L, noise.sd = 0.1) {
    set.seed(seed)
    X <- matrix(stats::runif(n * 2L, -1, 1), ncol = 2L)
    truth <- 0.5 + 1.0 * X[, 1L] - 0.7 * X[, 2L] + 0.4 * X[, 1L]^2 +
        0.3 * X[, 1L] * X[, 2L] - 0.6 * X[, 2L]^2
    y <- truth + stats::rnorm(n, sd = noise.sd)
    list(X = X, y = y, truth = truth)
}

# Local weighted mean per anchor, computed through the actual package
# internals (the same support ordering and kernel weights the fit uses).
e213.local.weighted.means <- function(X, y, support.size, kernel) {
    vapply(seq_len(nrow(X)), function(i) {
        ordered <- geosmooth:::.klp.local.order(X, X[i, ], support.size)
        w <- geosmooth:::.klp.kernel.weights(ordered$distances, kernel)
        sum(w * y[ordered$index]) / sum(w)
    }, numeric(1L))
}

e213.fit <- function(g, rho, target, degree = 1L) {
    fit.lps(
        X = g$X,
        y = g$y,
        foldid = rep(1:5, length.out = length(g$y)),
        support.grid = 30L,
        degree.grid = degree,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = rho,
        ridge.condition.max = Inf,
        unstable.action = "na",
        ridge.shrinkage.target = target
    )$fitted.values
}

test_that("E2.13 aligned ridge shrinks to the local weighted mean, not zero", {
    g <- e213.g1()
    ybar.w <- e213.local.weighted.means(g$X, g$y, 30L, "gaussian")

    f0 <- e213.fit(g, 0, "local.mean")
    f.tiny <- e213.fit(g, 1e-8, "local.mean")
    f.mid1 <- e213.fit(g, 1e-2, "local.mean")
    f.mid2 <- e213.fit(g, 1, "local.mean")
    f.big <- e213.fit(g, 1e2, "local.mean")

    expect_true(all(is.finite(c(f0, f.tiny, f.mid1, f.mid2, f.big))))

    # Tiny-ridge prediction-invariance (contract threshold 1e-6):
    expect_lt(max(abs(f.tiny - f0)), 1e-6)

    # At rho = 0 the penalty vanishes, so the two shrinkage targets must
    # coincide exactly (the aligned branch defers to the legacy solve):
    expect_identical(f0, e213.fit(g, 0, "zero"))

    # Large ridge shrinks toward the local weighted mean, not toward 0:
    # per-anchor relative comparison |f - ybar_w| / |f - 0|.
    ratio <- abs(f.big - ybar.w) / abs(f.big)
    expect_lt(stats::median(ratio), 0.05)
    expect_lt(max(ratio[abs(ybar.w) > 0.2]), 0.1)
    expect_gt(sum(abs(ybar.w) > 0.2), 100L)   # the restriction is not vacuous
    # ... and the fitted level matches the weighted-mean level, not zero:
    expect_lt(abs(mean(abs(f.big)) - mean(abs(ybar.w))), 0.05)
    expect_gt(mean(abs(f.big)), 0.4)

    # The shrinkage path is monotone toward the mean: larger rho gets
    # closer to ybar_w.
    expect_lt(mean(abs(f.big - ybar.w)), mean(abs(f.mid2 - ybar.w)))
    expect_lt(mean(abs(f.mid2 - ybar.w)), mean(abs(f.mid1 - ybar.w)))

    # Degree-0 exactness: with only the constant in the design, the aligned
    # ridge prediction IS the local weighted mean for any rho > 0.
    f.deg0 <- e213.fit(g, 1, "local.mean", degree = 0L)
    expect_lt(max(abs(f.deg0 - ybar.w)), 1e-10)
})

test_that("E2.13 legacy default arm shrinks to zero (documented pre-fix behavior)", {
    g <- e213.g1()
    ybar.w <- e213.local.weighted.means(g$X, g$y, 30L, "gaussian")

    z.big <- e213.fit(g, 1e2, "zero")
    f.big <- e213.fit(g, 1e2, "local.mean")

    # The legacy penalty acts on the constant direction too, so large
    # ridge collapses the prediction toward 0 -- the statistically wrong
    # shrinkage target the plan documents:
    expect_lt(mean(abs(z.big)), 0.1 * mean(abs(ybar.w)))
    expect_lt(max(abs(z.big)), 0.1)

    # The two arms differ materially at large rho (the change is real and
    # visible, not cosmetic):
    expect_gt(max(abs(f.big - z.big)), 0.5)
})

test_that("E2.13 S A2: the default path reproduces pre-change fits", {
    # Reference values were generated by
    # dev/methods/lps/ci/e2_13_pin_reference_fits.R
    # with R/lps.R in its pre-E2.13 state (HEAD c621e2f at generation; the
    # reference was committed immediately after in c796408). The stored values
    # retain 17 significant digits; comparisons allow only floating-point
    # roundoff across R and BLAS implementations.
    ref.path <- "../../dev/methods/lps/runs/tier2/e2_13_reference_fits.csv"
    skip_if_not(file.exists(ref.path),
                "source-tree reference fits are not included in the package tarball")
    ref <- utils::read.csv(
        ref.path,
        colClasses = c(value17 = "character"),
        stringsAsFactors = FALSE)
    expect_identical(unique(ref$generated.at.commit),
                     "c621e2ff311a776ca09130cd41d189a05e5d9a1a")

    g <- e213.g1()   # same construction as the pinning script
    fit.gaussian <- fit.lps(
        X = g$X, y = g$y, foldid = rep(1:5, length.out = 150L),
        support.grid = c(20L, 30L), degree.grid = c(0L, 1L),
        kernel.grid = "gaussian", coordinate.method = "coordinates",
        backend = "R"
        # defaults: orthogonal.polynomial.drop, ridge grid c(0,1e-10,1e-8),
        # ridge.condition.max 1e12, unstable.action "na",
        # ridge.shrinkage.target "zero" (not supplied)
    )
    expect_equal(
        fit.gaussian$fitted.values,
        as.numeric(ref$value17[ref$family == "gaussian" &
                                   ref$quantity == "fitted"]),
        tolerance = 1e-14
    )
    expect_equal(
        fit.gaussian$cv.table$cv.rmse.observed,
        as.numeric(ref$value17[ref$family == "gaussian" &
                                   ref$quantity == "cv.rmse.observed"]),
        tolerance = 1e-14
    )

    set.seed(1302L)
    p <- pmin(pmax(stats::plogis(2 * g$X[, 1L]), 0.05), 0.95)
    y.bin <- stats::rbinom(150L, 1L, p)
    fit.bernoulli <- fit.lps(
        X = g$X, y = y.bin, foldid = rep(1:5, length.out = 150L),
        support.grid = c(20L, 30L), degree.grid = c(0L, 1L),
        kernel.grid = "gaussian", coordinate.method = "coordinates",
        backend = "R", outcome.family = "bernoulli"
    )
    expect_equal(
        fit.bernoulli$fitted.values,
        as.numeric(ref$value17[ref$family == "bernoulli" &
                                   ref$quantity == "fitted"]),
        tolerance = 1e-14
    )
    expect_equal(
        fit.bernoulli$cv.table$cv.brier.observed,
        as.numeric(ref$value17[ref$family == "bernoulli" &
                                   ref$quantity == "cv.brier.observed"]),
        tolerance = 1e-14
    )
})

test_that("E2.13 scope: binomial mode warns that the alignment does not apply", {
    set.seed(1303L)
    X <- matrix(stats::runif(80L * 2L, -1, 1), ncol = 2L)
    y <- stats::rbinom(80L, 1L, pmin(pmax(stats::plogis(2 * X[, 1L]),
                                          0.05), 0.95))
    expect_warning(
        fit.lps(
            X = X, y = y, foldid = rep(1:5, length.out = 80L),
            support.grid = 30L, degree.grid = 1L, kernel.grid = "gaussian",
            coordinate.method = "coordinates", backend = "R",
            outcome.family = "binomial",
            ridge.shrinkage.target = "local.mean"
        ),
        "has no effect"
    )
})
