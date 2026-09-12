# E1.9 GATEs -- decouple bandwidth from support size.
#
# Contract: dev/methods/lps/audit_contracts/lps_tiers1to4 (sB / E1.9);
# frozen spec dev/methods/lps/specs/lps_experimental_plan_2026-06-09.tex
# (E1.9).
# Fixtures: helper-lps-e1-9.R (configurations) and
# helper-lps-e1-9-reference.R (pre-change pinned values; provenance in its
# header).
#
# Sub-item (a): ESS/K + last-weight characterization GATE. Pins the current
#   kernel-weight semantics (bandwidth tied to the K-th NN distance) on the
#   ACTUAL internal routine geosmooth:::.klp.kernel.weights -- not a
#   re-implementation -- per resolved contract question G1.
# Sub-item (b): b = 1 backward-compatibility GATE. Fits with the default
#   bandwidth.multiplier.grid = 1 must equal the pre-change fits within
#   tau_alg = 1e-10 (they are bit-identical at authorship; the tolerance is
#   the contract's frozen threshold).
# Sub-item (c) (benefit STUDY over G3a/G3d) is NOT here: it is a
#   dev/methods/lps/ci study deferred until Amendment 1 binds the consolidated DGP
#   generators.

test_that("E1.9a characterization: ESS/K and last-weight pin current kernel-weight behavior", {
    distances <- e19.characterization.distances()
    K <- length(distances)
    expect_identical(K, 20L)

    kish.ess <- function(w) sum(w)^2 / sum(w^2)

    weights <- lapply(
        c(gaussian = "gaussian", tricube = "tricube",
          epanechnikov = "epanechnikov", triangular = "triangular"),
        function(kernel) .klp.kernel.weights(distances, kernel)
    )

    # All kernels return one finite nonnegative weight per distance, and the
    # weights are nonincreasing in distance (sanity of the fixture).
    for (w in weights) {
        expect_length(w, K)
        expect_true(all(is.finite(w)))
        expect_true(all(w >= 0))
        expect_true(all(diff(w) <= 0))
    }

    # Frozen thresholds (contract sB / E1.9):
    # gaussian is near-flat at the K-NN bandwidth: ESS/K > 0.9.
    expect_gt(kish.ess(weights$gaussian) / K, 0.9)
    # tricube concentrates: ESS/K < 0.85.
    expect_lt(kish.ess(weights$tricube) / K, 0.85)
    # compact kernels effectively zero the K-th neighbor:
    # w_(K) / max_j w_j < 1e-6.
    for (kernel in c("tricube", "epanechnikov", "triangular")) {
        w <- weights[[kernel]]
        expect_lt(w[[K]] / max(w), 1e-6)
    }
})

test_that("E1.9b exactness: default fits equal pre-change pinned fits within 1e-10", {
    tau.alg <- 1e-10

    configs <- list(
        A = e19.fit.A(),
        B = e19.fit.B(),
        C = e19.fit.C()
    )
    for (name in names(configs)) {
        fit <- configs[[name]]
        ref <- e19.reference[[name]]

        expect_false(anyNA(fit$fitted.values))
        expect_lt(max(abs(fit$fitted.values - ref$fitted.values)), tau.alg)
        expect_lt(
            max(abs(fit$cv.table$cv.rmse.observed - ref$cv.rmse.observed)),
            tau.alg
        )
        # Candidate-table alignment: adding the bandwidth axis must not
        # reorder the historical (support, degree, kernel) grid.
        expect_identical(as.integer(fit$cv.table$support.size),
                         ref$cand.support.size)
        expect_identical(as.integer(fit$cv.table$degree), ref$cand.degree)
        expect_identical(as.character(fit$cv.table$kernel), ref$cand.kernel)
        # Selection is unchanged.
        expect_identical(as.integer(fit$selected$support.size[[1L]]),
                         ref$selected.support.size)
        expect_identical(as.integer(fit$selected$degree[[1L]]),
                         ref$selected.degree)
        expect_identical(as.character(fit$selected$kernel[[1L]]),
                         ref$selected.kernel)
        # The new diagnostic fields exist and carry the default multiplier.
        expect_identical(unique(fit$cv.table$bandwidth.multiplier), 1)
        expect_identical(fit$selected$bandwidth.multiplier[[1L]], 1)
        expect_identical(fit$bandwidth.multiplier.grid, 1)
    }
})

test_that("E1.9b exactness: explicit bandwidth.multiplier.grid = 1 is bit-identical to the default", {
    for (fitter in list(e19.fit.A, e19.fit.B, e19.fit.C)) {
        fit.default <- fitter()
        fit.b1 <- fitter(bandwidth.multiplier.grid = 1)
        expect_identical(fit.default$fitted.values, fit.b1$fitted.values)
        expect_identical(fit.default$fitted.values.raw,
                         fit.b1$fitted.values.raw)
        expect_identical(fit.default$cv.table, fit.b1$cv.table)
        expect_identical(fit.default$selected, fit.b1$selected)
    }
})

test_that("E1.9 multiplier semantics: large b approaches the unweighted local mean", {
    # As b -> Inf, u -> 0 and every kernel weight -> its maximum (flat
    # weighting), so a degree-0 fit at an interior point converges to the
    # plain mean of the support responses -- a closed-form target that does
    # not consult the weight routine.
    d <- e19.pin.data.ambient()
    K <- 12L
    center.idx <- 1L
    dist <- sqrt(colSums((t(d$X) - d$X[center.idx, ])^2))
    support <- order(dist, seq_along(dist))[seq_len(K)]
    flat.target <- mean(d$y[support])

    fit <- fit.lps(
        X = d$X,
        y = d$y,
        foldid = rep(1:2, length.out = nrow(d$X)),
        support.grid = K,
        degree.grid = 0L,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na",
        bandwidth.multiplier.grid = 1e6
    )
    expect_lt(abs(fit$fitted.values[[center.idx]] - flat.target), 1e-9)

    # And b != 1 genuinely changes the fits relative to the pinned b = 1
    # reference (the multiplier is live, not cosmetic).
    fit.b2 <- e19.fit.B(bandwidth.multiplier.grid = 2)
    expect_gt(
        max(abs(fit.b2$fitted.values - e19.reference$B$fitted.values)),
        1e-6
    )
})

test_that("E1.9 selection grid: multiplier candidates expand the CV table and are selectable", {
    grid <- c(0.5, 1, 2, 4)
    fit <- e19.fit.B(bandwidth.multiplier.grid = grid)

    expect_identical(nrow(fit$cv.table), length(grid))
    expect_identical(fit$cv.table$bandwidth.multiplier, grid)
    expect_true(all(is.finite(fit$cv.table$cv.rmse.observed)))
    expect_identical(fit$bandwidth.multiplier.grid, grid)
    # The selected multiplier achieves the minimum CV score.
    expect_identical(
        fit$selected$cv.rmse.observed[[1L]],
        min(fit$cv.table$cv.rmse.observed)
    )
    expect_true(fit$selected$bandwidth.multiplier[[1L]] %in% grid)
    # The final refit honors the selected multiplier: refitting with the
    # selected multiplier as a singleton grid reproduces fitted.values.
    fit.singleton <- e19.fit.B(
        bandwidth.multiplier.grid = fit$selected$bandwidth.multiplier[[1L]]
    )
    expect_identical(fit$fitted.values, fit.singleton$fitted.values)
})

test_that("E1.9 backend coupling: b != 1 requires the R backend", {
    d <- e19.pin.data.ambient()
    base.args <- list(
        X = d$X,
        y = d$y,
        foldid = rep(1:2, length.out = nrow(d$X)),
        support.grid = 10L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )

    fit.auto.default <- do.call(fit.lps, c(base.args, list(backend = "auto")))
    expect_identical(fit.auto.default$backend.used, "cpp")

    fit.auto.b2 <- do.call(
        fit.lps,
        c(base.args, list(backend = "auto", bandwidth.multiplier.grid = 2))
    )
    expect_identical(fit.auto.b2$backend.used, "R")

    expect_error(
        do.call(
            fit.lps,
            c(base.args, list(backend = "cpp", bandwidth.multiplier.grid = 2))
        ),
        "bandwidth multipliers"
    )
})

test_that("E1.9 grid validation mirrors the ridge-grid cleaner", {
    expect_error(
        e19.fit.B(bandwidth.multiplier.grid = numeric(0)),
        "bandwidth.multiplier.grid"
    )
    expect_error(
        e19.fit.B(bandwidth.multiplier.grid = c(-1, NA)),
        "bandwidth.multiplier.grid"
    )
    # Duplicates and unsorted input are cleaned deterministically.
    fit <- e19.fit.B(bandwidth.multiplier.grid = c(2, 1, 2))
    expect_identical(fit$bandwidth.multiplier.grid, c(1, 2))
    expect_identical(fit$cv.table$bandwidth.multiplier, c(1, 2))
})

test_that("E1.9 predict.lps applies the selected multiplier and tolerates pre-change objects", {
    fit <- e19.fit.B(bandwidth.multiplier.grid = c(0.5, 1, 2, 4))
    # In-sample predict reproduces fitted.values under the selected b.
    expect_identical(predict(fit), fit$fitted.values)

    # An object lacking the bandwidth.multiplier field (fitted before this
    # change) predicts with b = 1.
    fit.b1 <- e19.fit.B()
    legacy <- fit.b1
    legacy$selected$bandwidth.multiplier <- NULL
    legacy$bandwidth.multiplier.grid <- NULL
    expect_identical(predict(legacy), predict(fit.b1))
})
