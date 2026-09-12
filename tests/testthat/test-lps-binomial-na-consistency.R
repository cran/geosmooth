# E2.15 — binomial selection NA-consistency (Tier 2 GATE).
#
# Amendment (dev/methods/lps/audit_contracts/
# lps_e2_15_binomial_na_consistency_amendment_2026-06-12.md, promoting
# spec-memo item 11): gaussian and bernoulli
# selection score any candidate with a non-finite CV prediction as Inf
# (the .klp.rmse convention), but binomial selection previously delegated
# to .klp.logloss, which DROPS non-finite pairs -- so a candidate
# predicting NA on most points was scored only where it happened to
# succeed and could win. Post-fix, the binomial selection column
# cv.logloss.observed is Inf for any candidate with a non-finite CV
# prediction; the .klp.logloss function itself is unchanged and still
# backs the logloss.clipped probability DIAGNOSTIC (reporting over the
# observed pairs is intended there).
#
# Fixture (fully deterministic, no RNG): n = 120 points on a 1-D grid in
# [-1, 1], labels alternating in ten stripes of twelve points, plus ONE
# flipped label at index 66 (the center of stripe 5). Under
# unstable.action = "na":
#   - candidate A (support 8): windows straddling a stripe boundary
#     without the flip are EXACTLY SEPARABLE -> the logistic solve cannot
#     converge -> NA. Realized: 55.8% of A's out-of-fold predictions are
#     NA -- but its retained points are dominated by confident-correct
#     stripe interiors, so under the old drop-NA rule A scores 0.3016;
#   - candidate B (support 110): every window contains the flipped label,
#     so no window is separable; B predicts everywhere (0% NA) and scores
#     0.6818.
# Under the old rule A WINS by margin 0.38 (a real flip, not a tie) --
# rewarding the candidate for failing on its hard points. Post-fix A is
# unselectable (Inf) and B is selected.

e215.fixture <- function() {
    n <- 120L
    x <- seq(-1, 1, length.out = n)
    stripe <- pmin(floor((x + 1) / 0.2), 9)
    y <- as.integer(stripe %% 2 == 1)
    y[66L] <- 1L - y[66L]
    list(X = matrix(x, ncol = 1L), y = y,
         foldid = rep(1:5, length.out = n))
}

e215.fit <- function(g) {
    fit.lps(
        X = g$X,
        y = g$y,
        foldid = g$foldid,
        support.grid = c(8L, 110L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na",
        outcome.family = "binomial",
        keep.cv.predictions = TRUE
    )
}

test_that("E2.15 an NA-heavy binomial candidate is unselectable", {
    g <- e215.fixture()
    fit <- e215.fit(g)
    ct <- fit$cv.table
    pred <- fit$cv.predictions
    a <- which(ct$support.size == 8L)
    b <- which(ct$support.size == 110L)

    # Fixture validity: A fails to predict on a MAJORITY of eval points;
    # B predicts everywhere.
    na.fraction <- colMeans(!is.finite(pred))
    expect_gt(na.fraction[[a]], 0.5)
    expect_identical(na.fraction[[b]], 0)

    # Pre-fix demonstration (documented from the scored predictions, the
    # way E2.12 documents its raw-vs-clipped flip): under the old drop-NA
    # rule -- .klp.logloss, unchanged, scoring only the finite pairs --
    # the NA-heavy candidate WINS, and not by a tie.
    old.scores <- vapply(
        seq_len(ncol(pred)),
        function(j) geosmooth:::.klp.logloss(g$y, pred[, j]),
        numeric(1L)
    )
    expect_identical(which.min(old.scores), a)
    expect_gt(abs(old.scores[[b]] - old.scores[[a]]), 0.1)

    # (1) Post-fix: the NA-heavy candidate's selection score is Inf and it
    # is not selected.
    expect_identical(ct$cv.logloss.observed[[a]], Inf)
    expect_false(fit$selected$support.size[[1L]] == 8L)

    # (2) The complete-prediction candidate is selected instead, on a
    # finite score that equals the (unchanged) observed-pairs log loss.
    expect_identical(fit$selected$support.size[[1L]], 110L)
    expect_true(is.finite(ct$cv.logloss.observed[[b]]))
    expect_equal(ct$cv.logloss.observed[[b]],
                 geosmooth:::.klp.logloss(g$y, pred[, b]),
                 tolerance = 1e-12)

    # (3) Cross-family consistency: the Inf-on-any-non-finite rule now
    # holds across both binomial selection-facing columns, matching the
    # gaussian/bernoulli convention.
    expect_identical(ct$cv.brier.observed[[a]], Inf)
    expect_true(is.finite(ct$cv.brier.observed[[b]]))
})

test_that("E2.15 healthy-data binomial selection is unchanged (regression pin)", {
    # On all-finite predictions the Inf-on-any-non-finite guard is inert:
    # the selection score is computed by the same unchanged .klp.logloss
    # call as before the fix, so healthy-data selection is bit-for-bit
    # identical. Pinned here on a no-NA fixture.
    set.seed(7001L)
    X <- matrix(stats::runif(400L * 2L, -1, 1), ncol = 2L)
    eta <- 6 * tanh(15 * X[, 1L])
    p <- pmin(pmax(stats::plogis(eta), 0.05), 0.95)
    y <- stats::rbinom(400L, 1L, p)
    fit <- fit.lps(
        X = X, y = y, foldid = rep(1:5, length.out = 400L),
        support.grid = 60L, degree.grid = c(0L, 1L),
        kernel.grid = "gaussian", coordinate.method = "coordinates",
        backend = "R", design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0, ridge.condition.max = Inf,
        unstable.action = "mean", outcome.family = "binomial",
        keep.cv.predictions = TRUE
    )
    pred <- fit$cv.predictions
    expect_true(all(is.finite(pred)))
    # The selection column equals the unchanged observed-pairs log loss on
    # every candidate (no Inf substitution fired), and the selected
    # candidate is the finite argmin.
    expect_equal(
        fit$cv.table$cv.logloss.observed,
        vapply(seq_len(ncol(pred)),
               function(j) geosmooth:::.klp.logloss(y, pred[, j]),
               numeric(1L)),
        tolerance = 1e-12
    )
    expect_true(all(is.finite(fit$cv.table$cv.logloss.observed)))
    expect_identical(
        which.min(fit$cv.table$cv.logloss.observed),
        which(fit$cv.table$support.size ==
                  fit$selected$support.size[[1L]] &
                  fit$cv.table$degree == fit$selected$degree[[1L]])
    )
})
