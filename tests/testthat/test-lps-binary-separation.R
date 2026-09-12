# E2.14 — local logistic robustness under separation (Tier 2 GATE).
#
# Contract (dev/methods/lps/audit_contracts/lps_tiers1to4 via
# dev/methods/lps/specs, frozen spec
# S E2.14): with IRLS step-halving the deviance trajectory is non-increasing
# to within 1e-8 per step; the final fitted probability lies strictly in
# (0, 1) and the coefficients are finite; non-convergence within the
# iteration cap triggers the documented, telemetered fallback; EXACT
# separation hits the fallback (no unbounded loop, no NaN). Assertions are
# on the deviance TRAJECTORY, not only the endpoint, so oscillation is
# caught.
#
# DGP (frozen spec): a single constructed support, z in R, y = 1{z > 0},
# one flipped label for the near-separable case, kernel weights from
# "gaussian"; the solver is called directly (degree 1,
# outcome.family = "binomial", backend = "R"), singleton ridge rho = 0 with
# ridge.condition.max = Inf so exactly one IRLS attempt is traced.
#
# Fixture geometry (deterministic, no RNG): ten equispaced negatives near
# the evaluation center plus one far positive. The flipped label sits at
# the EDGE of the negative cluster (index 2). This placement matters: the
# pre-fix plain-Newton iteration on this support overshoots (measured
# deviance increase of about 4.9e2 at step 2, eight increases > 1e-8 in 50
# iterations, no convergence; replication script:
# dev/methods/lps/ci/e2_14_prefix_newton_overshoot.R), while a mid-cluster flip
# (e.g. index 5) yields a monotone plain-Newton trajectory and would make
# the step-halving mutation undetectable. Step-halving recovers a genuine
# converged fit on this support (one halving, convergence at iteration 9).

e214.fixture <- function(flip = TRUE) {
    z <- c(seq(-0.20, -0.04, by = 0.02), 6)
    y <- as.numeric(z > 0)
    if (flip) {
        y[[2L]] <- 1
    }
    weights <- geosmooth:::.klp.kernel.weights(abs(z), "gaussian")
    list(
        design = cbind(1, z),
        y = y,
        weights = weights,
        prediction.row = matrix(c(1, 0), nrow = 1L)
    )
}

e214.solve <- function(fx) {
    geosmooth:::.klp.solve.local.logistic(
        design = fx$design,
        y = fx$y,
        weights = fx$weights,
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-8,
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        prediction.row = fx$prediction.row
    )
}

test_that("E2.14 near-separable support: deviance trajectory is monotone, fit is bounded", {
    fx <- e214.fixture(flip = TRUE)
    # Two genuinely distinct classes, not linearly separable (the flipped
    # positive lies inside the negative cluster), separable once unflipped.
    expect_setequal(unique(fx$y), c(0, 1))
    solved <- e214.solve(fx)

    expect_true(isTRUE(solved$ok))
    expect_identical(solved$status, "ok")
    expect_true(isTRUE(solved$converged))

    trace <- solved$deviance.trace
    expect_true(is.numeric(trace) && length(trace) >= 2L)
    expect_true(all(is.finite(trace)))
    # The GATE assertion: per-step deviance non-increasing within 1e-8,
    # asserted on the whole trajectory.
    expect_true(all(diff(trace) <= 1e-8))

    expect_true(all(is.finite(solved$coefficients)))
    expect_true(is.finite(solved$prediction[[1L]]))
    expect_gt(solved$prediction[[1L]], 0)
    expect_lt(solved$prediction[[1L]], 1)

    # Non-vacuousness: the step-halving actually engaged on this support.
    # Disabling it (the contract's named mutation) makes the plain-Newton
    # trajectory non-monotone here, reddening the trajectory assertion above.
    expect_gte(solved$step.halvings, 1L)
    expect_identical(solved$iterations, length(trace) - 1L)

    # Determinism: the solve is a fixed-point computation with no RNG.
    solved.again <- e214.solve(fx)
    expect_identical(solved$deviance.trace, solved.again$deviance.trace)
    expect_identical(solved$coefficients, solved.again$coefficients)
})

test_that("E2.14 exact separation: bounded non-convergence, documented fallback, no NaN", {
    fx <- e214.fixture(flip = FALSE)
    expect_setequal(unique(fx$y), c(0, 1))
    solved <- e214.solve(fx)

    # The unpenalized logistic MLE does not exist under exact separation;
    # the solve must stop at the iteration cap (bounded loop), with a
    # finite, monotone trajectory and no NaN anywhere.
    expect_false(isTRUE(solved$ok))
    expect_identical(solved$status, "not_converged")
    expect_false(isTRUE(solved$converged))
    trace <- solved$deviance.trace
    expect_true(all(is.finite(trace)))
    expect_true(all(diff(trace) <= 1e-8))
    expect_identical(solved$iterations, 50L)
    expect_identical(length(trace), 51L)

    # The documented fallback fires at the fitting layer and is recorded in
    # the telemetry: weighted event rate under unstable.action = "mean".
    telemetry.mean <- geosmooth:::.klp.logistic.telemetry.new("binomial")
    fitted.mean <- geosmooth:::.klp.fit.logistic.prob.design(
        design = fx$design,
        y = fx$y,
        weights = fx$weights,
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-8,
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean",
        logistic.telemetry = telemetry.mean
    )
    expect_true(is.finite(fitted.mean))
    expect_gte(fitted.mean, 0)
    expect_lte(fitted.mean, 1)
    expect_equal(
        fitted.mean,
        stats::weighted.mean(fx$y, fx$weights),
        tolerance = 1e-12
    )
    summary.mean <- geosmooth:::.klp.logistic.telemetry.summary(telemetry.mean)
    expect_identical(summary.mean$attempted, 1L)
    expect_identical(summary.mean$converged, 0L)
    expect_identical(summary.mean$fallback.path.count, 1L)
    expect_identical(summary.mean$event.rate.fallback.count, 1L)
    expect_identical(summary.mean$na.failure.count, 0L)
    expect_identical(
        as.integer(summary.mean$status.counts[["not_converged"]]), 1L
    )

    # ... and to NA under unstable.action = "na" (the sanctioned guarded
    # output), likewise telemetered.
    telemetry.na <- geosmooth:::.klp.logistic.telemetry.new("binomial")
    fitted.na <- geosmooth:::.klp.fit.logistic.prob.design(
        design = fx$design,
        y = fx$y,
        weights = fx$weights,
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-8,
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na",
        logistic.telemetry = telemetry.na
    )
    expect_true(is.na(fitted.na))
    summary.na <- geosmooth:::.klp.logistic.telemetry.summary(telemetry.na)
    expect_identical(summary.na$fallback.path.count, 1L)
    expect_identical(summary.na$event.rate.fallback.count, 0L)
    expect_identical(summary.na$na.failure.count, 1L)
})

test_that("E2.14 well-behaved supports keep plain-IRLS numerics (zero halvings)", {
    # Regression pin for the always-on change: when every full Newton step
    # already satisfies the deviance slack, the step-halving check accepts
    # the identical Newton candidate, so the iterate sequence -- and hence
    # the returned solve -- is numerically identical to the pre-change
    # solver. Verified here against an independent plain-IRLS replication
    # (no halving) of the documented update equations on a benign support.
    z <- matrix(seq(-1, 1, length.out = 30L), ncol = 1L)
    y <- as.numeric(z[, 1L] > 0)
    y[[18L]] <- 0   # one flipped label well inside the positive run
    weights <- geosmooth:::.klp.kernel.weights(abs(z[, 1L]), "gaussian")
    design <- cbind(1, z)
    prediction.row <- matrix(c(1, 0), nrow = 1L)

    solved <- geosmooth:::.klp.solve.local.logistic(
        design = design,
        y = y,
        weights = weights,
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-8,
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        prediction.row = prediction.row
    )
    expect_true(isTRUE(solved$ok))
    expect_true(isTRUE(solved$converged))
    expect_identical(solved$step.halvings, 0L)

    # Plain-IRLS replication: identical update equations, no step control.
    transformed <- geosmooth:::.klp.orthogonal.polynomial.transform(
        design = design,
        weights = weights,
        prediction.rows = prediction.row,
        design.drop.tol = 1e-8
    )
    X <- transformed$design
    ybar <- stats::weighted.mean(y, weights)
    ybar <- min(1 - 1e-6, max(1e-6, ybar))
    beta <- c(stats::qlogis(ybar), rep(0, ncol(X) - 1L))
    converged <- FALSE
    for (iter in seq_len(50L)) {
        eta <- as.numeric(X %*% beta)
        mu <- stats::plogis(pmax(-35, pmin(35, eta)))
        variance <- pmax(mu * (1 - mu), 1e-8)
        ww <- weights * variance
        zz <- eta + (y - mu) / variance
        xw <- X * sqrt(ww)
        beta.new <- as.numeric(solve(crossprod(xw), crossprod(xw, zz * sqrt(ww))))
        step.converged <- max(abs(beta.new - beta)) < 1e-7 * (1 + max(abs(beta)))
        beta <- beta.new
        if (step.converged) {
            converged <- TRUE
            break
        }
    }
    expect_true(converged)
    expect_identical(solved$coefficients, beta)
    expect_identical(
        solved$prediction[[1L]],
        stats::plogis(as.numeric(transformed$prediction.rows %*% beta))
    )
})
