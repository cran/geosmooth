# E2.12 — binary selection-metric consistency & log-loss clipping (Tier 2
# GATEs).
#
# Contract (S C / E2.12; frozen spec S E2.12):
#   (a) Bernoulli-mode selection scores the DEPLOYED (clipped) metric; the
#       pre-fix discrepancy (selection on raw, unclipped predictions while
#       deployment clips to [0,1]) is demonstrated in this same file as a
#       documented motivating case. Post-fix equality within comparison
#       tolerance 1e-6.
#   (b) The log-loss probability truncation is pinned at 1e-6 (from 1e-15);
#       the 1e-15 score's clip-sensitive component is demonstrably one
#       deliberately confident-wrong held-out point (dominance definition
#       per spec-memo item 11b: unique binder; > 99.9% of the cross-clip
#       score change; largest single contributor by > 5x on the selected
#       candidate). Cross-clip selection stability over {1e-6, 1e-3} is a
#       STUDY (dev/methods/lps/ci/e2_12_crossclip_stability_study.R), NOT gated
#       here.
#
# DGP: constructed G6 cases, n = 400, deterministic seeds, no replication
# (contract smoke/full sizing). G6 with named overrides (plan S sec:dgp
# allows overriding named parameters): log-odds eta(x) = 6 * tanh(15 * x1)
# on x in [-1,1]^2 -- a sharp cliff at x1 = 0 forcing some raw bernoulli
# CV predictions outside [0,1] -- with p = expit(eta) clipped to
# [0.05, 0.95] by the G6 construction (alpha = 0; the surface is symmetric,
# so target prevalence 0.5).

e212.g6 <- function(n = 400L, seed, sharpness = 15, amplitude = 6) {
    set.seed(seed)
    X <- matrix(stats::runif(n * 2L, -1, 1), ncol = 2L)
    eta <- amplitude * tanh(sharpness * X[, 1L])
    p <- pmin(pmax(stats::plogis(eta), 0.05), 0.95)
    y <- stats::rbinom(n, 1L, p)
    list(X = X, p = p, y = y, foldid = rep(1:5, length.out = n))
}

test_that("E2.12a bernoulli selection scores the deployed clipped metric", {
    g <- e212.g6(seed = 2121L)
    fit <- fit.lps(
        X = g$X,
        y = g$y,
        foldid = g$foldid,
        support.grid = c(8L, 60L),
        degree.grid = c(0L, 2L),
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean",
        outcome.family = "bernoulli",
        keep.cv.predictions = TRUE
    )
    ct <- fit$cv.table
    pred <- fit$cv.predictions
    expect_identical(dim(pred), c(length(g$y), nrow(ct)))
    expect_true(all(is.finite(pred)))

    # Fixture validity: the clip does real work -- raw out-of-fold
    # predictions genuinely leave [0,1] on the sharp surface.
    expect_gt(sum(pred < 0 | pred > 1), 0L)

    # The selection score IS the deployed clipped metric, recomputed here
    # from the actual CV predictions without package helpers: the observed
    # Brier score of the [0,1]-clipped predictions.
    deployed <- vapply(
        seq_len(ncol(pred)),
        function(j) mean((g$y - pmin(1, pmax(0, pred[, j])))^2),
        numeric(1L)
    )
    expect_equal(ct$cv.brier.observed, deployed, tolerance = 1e-12)
    # Contract acceptance tolerance for the equality (1e-6) on the selected
    # candidate's score:
    selected.row <- which(
        ct$support.size == fit$selected$support.size[[1L]] &
            ct$degree == fit$selected$degree[[1L]] &
            ct$kernel == fit$selected$kernel[[1L]]
    )
    expect_lt(
        abs(fit$selected$cv.brier.observed[[1L]] - deployed[[selected.row]]),
        1e-6
    )

    # Documented motivating case (the pre-fix discrepancy): the raw metric
    # (the OLD selection rule scored cv.rmse.observed, i.e. unclipped
    # predictions) and the deployed clipped metric rank the candidates
    # DIFFERENTLY on this fixture.
    raw.idx <- which.min(ct$cv.rmse.observed)
    clip.idx <- which.min(ct$cv.brier.observed)
    expect_true(raw.idx != clip.idx)
    # The old rule would have selected (support 8, degree 0); the deployed
    # metric selects (support 60, degree 2):
    expect_identical(ct$support.size[[raw.idx]], 8L)
    expect_identical(ct$degree[[raw.idx]], 0L)
    expect_identical(ct$support.size[[clip.idx]], 60L)
    expect_identical(ct$degree[[clip.idx]], 2L)
    # The post-fix fit deploys the clipped winner:
    expect_identical(fit$selected$support.size[[1L]], 60L)
    expect_identical(fit$selected$degree[[1L]], 2L)
    expect_identical(selected.row, clip.idx)

    # The flip is real, not a tie (plan safeguard: verify by hand): the two
    # candidates differ by > 100x the 1e-6 comparison tolerance on BOTH
    # metrics.
    expect_gt(
        abs(ct$cv.rmse.observed[[raw.idx]]^2 - ct$cv.rmse.observed[[clip.idx]]^2),
        1e-4
    )
    expect_gt(
        abs(ct$cv.brier.observed[[raw.idx]] - ct$cv.brier.observed[[clip.idx]]),
        1e-4
    )

    # And the discrepancy is visible inside the returned table itself:
    # cv.brier.observed is NOT the raw rmse^2 wherever clipping engaged.
    expect_gt(
        max(abs(ct$cv.brier.observed - ct$cv.rmse.observed^2)),
        1e-4
    )
})

test_that("E2.12 keep.cv.predictions default reproduces the prior fit object exactly", {
    # S A2 regression pin for the new opt-in argument: at the default
    # (FALSE) the returned object carries no cv.predictions element and
    # every other component is identical to the keep.cv.predictions = TRUE
    # fit (the argument only appends the matrix; `call` necessarily
    # differs and is excluded).
    g <- e212.g6(n = 60L, seed = 99L)
    fit.args <- list(
        X = g$X,
        y = g$y,
        foldid = rep(1:5, length.out = 60L),
        support.grid = c(6L, 8L),
        degree.grid = c(0L, 1L),
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean",
        outcome.family = "bernoulli"
    )
    fit.default <- do.call(fit.lps, fit.args)
    fit.keep <- do.call(fit.lps, c(fit.args, list(keep.cv.predictions = TRUE)))
    expect_false("cv.predictions" %in% names(fit.default))
    expect_true("cv.predictions" %in% names(fit.keep))
    expect_identical(
        setdiff(names(fit.keep), names(fit.default)),
        "cv.predictions"
    )
    shared <- setdiff(names(fit.default), "call")
    expect_identical(fit.default[shared], fit.keep[shared])
    # And the stored matrix is the per-candidate out-of-fold predictions:
    expect_identical(dim(fit.keep$cv.predictions),
                     c(60L, nrow(fit.keep$cv.table)))
})

test_that("E2.12 bernoulli always uses the R CV path (no raw-metric backend corner)", {
    # Audit-required fix (dev/methods/lps/audits/tier2_audit_2026-06-11.md):
    # a legal
    # bernoulli configuration (monomial basis, singleton ridge 0,
    # ridge.condition.max = Inf, coordinates) used to reach the C++ CV
    # kernels, which return only the aggregate raw RMSE, silently reverting
    # selection to the pre-fix raw metric. Bernoulli now resolves
    # backend = "auto" to the R backend (exactly like "binomial") and
    # rejects an explicit C++ backend, making the deployed-metric
    # invariant unconditional over legal bernoulli configurations.
    g <- e212.g6(n = 60L, seed = 99L)
    cpp.eligible <- list(
        X = g$X,
        y = g$y,
        foldid = rep(1:5, length.out = 60L),
        support.grid = c(6L, 8L),
        degree.grid = c(0L, 1L),
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )
    fit.auto <- do.call(fit.lps, c(cpp.eligible, list(
        backend = "auto",
        outcome.family = "bernoulli",
        keep.cv.predictions = TRUE
    )))
    expect_identical(fit.auto$backend.used, "R")
    expect_false(is.null(fit.auto$cv.predictions))
    # The selection column is the deployed clipped metric, computed from
    # actual per-point CV predictions (never decorated from raw RMSE):
    expect_equal(
        fit.auto$cv.table$cv.brier.observed,
        vapply(
            seq_len(ncol(fit.auto$cv.predictions)),
            function(j) {
                mean((g$y - pmin(1, pmax(0, fit.auto$cv.predictions[, j])))^2)
            },
            numeric(1L)
        ),
        tolerance = 1e-12
    )
    # An explicit C++ backend with bernoulli is an error, mirroring the
    # binomial constraint:
    expect_error(
        do.call(fit.lps, c(cpp.eligible, list(
            backend = "cpp",
            outcome.family = "bernoulli"
        ))),
        "currently uses the R backend"
    )
    # The forcing is family-scoped: the same configuration in gaussian
    # mode still resolves to the C++ fast path.
    fit.gauss <- do.call(fit.lps, c(cpp.eligible, list(
        backend = "auto",
        outcome.family = "gaussian"
    )))
    expect_identical(fit.gauss$backend.used, "cpp")
})

test_that("E2.12b log-loss clip is pinned at 1e-6 and the 1e-15 instability is one point", {
    # -- The pin ------------------------------------------------------------
    expect_identical(
        formals(geosmooth:::.klp.clip.probability)$eps,
        1e-6
    )
    expect_identical(
        geosmooth:::.klp.clip.probability(c(0, 1)),
        c(1e-6, 1 - 1e-6)
    )
    # Deployed log loss of two maximally confident-wrong predictions is
    # -log(1e-6): the clip, not the raw probabilities, bounds the score.
    # Tolerance 1e-9 because the y = 0 term is log1p(-(1 - 1e-6)), which
    # differs from log(1e-6) by ~3e-11 in double precision.
    expect_equal(
        geosmooth:::.klp.logloss(c(1, 0), c(0, 1)),
        -log(1e-6),
        tolerance = 1e-9
    )

    # -- The motivating case: G6 with ONE deliberately confident-wrong
    # -- held-out point ------------------------------------------------------
    g <- e212.g6(seed = 7001L)
    fit.bin <- function(y) fit.lps(
        X = g$X,
        y = y,
        foldid = g$foldid,
        support.grid = 60L,
        degree.grid = c(0L, 1L),
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean",
        outcome.family = "binomial",
        keep.cv.predictions = TRUE
    )
    # The deliberate flip: observation 7 sits deep in the p = 0.05 plateau
    # (x1 = -0.702); its label is drawn 0 and every candidate's out-of-fold
    # prediction for it is ~0 (all-zero 60-point neighborhood). Setting its
    # label to 1 makes it confidently wrong when held out.
    flip <- 7L
    expect_identical(g$y[[flip]], 0L)
    expect_lt(g$X[flip, 1L], -0.5)
    y.flipped <- g$y
    y.flipped[[flip]] <- 1
    fit <- fit.bin(y.flipped)
    ct <- fit$cv.table
    pred <- fit$cv.predictions
    expect_true(all(is.finite(pred)))

    # Deployed-metric equality for the binomial selection score: the
    # cv.logloss.observed column equals the clip-1e-6 log loss recomputed
    # from the actual CV predictions without package helpers.
    manual.logloss <- function(p, y, eps) {
        p <- pmin(1 - eps, pmax(eps, p))
        -mean(y * log(p) + (1 - y) * log1p(-p))
    }
    expect_equal(
        ct$cv.logloss.observed,
        vapply(seq_len(ncol(pred)),
               function(j) manual.logloss(pred[, j], y.flipped, 1e-6),
               numeric(1L)),
        tolerance = 1e-12
    )

    # (i) The flipped point is the UNIQUE point whose wrong-side clip binds
    # at 1e-6, for every candidate; its realized out-of-fold prediction is
    # exactly the all-zero-neighborhood event-rate fallback, 0.
    for (j in seq_len(ncol(pred))) {
        binders <- which(
            (y.flipped == 1 & pred[, j] < 1e-6) |
                (y.flipped == 0 & pred[, j] > 1 - 1e-6)
        )
        expect_identical(binders, flip)
        expect_lt(pred[flip, j], 1e-12)
    }

    # (ii) The single point accounts for > 99.9% of the score change
    # between clips 1e-15 and 1e-6, for every candidate (spec-memo item 11b
    # dominance definition).
    for (j in seq_len(ncol(pred))) {
        total.swing <- manual.logloss(pred[, j], y.flipped, 1e-15) -
            manual.logloss(pred[, j], y.flipped, 1e-6)
        point.swing <- (-log(pmax(pred[flip, j], 1e-15)) +
                            log(pmax(pred[flip, j], 1e-6))) / length(y.flipped)
        expect_gt(total.swing, 0)
        expect_gt(point.swing / total.swing, 0.999)
    }

    # (iii) At clip 1e-15 the point is the largest single contributor to
    # the selected candidate's score by a factor > 5 (realized: 34.54 vs
    # 4.03).
    selected.row <- which(
        ct$support.size == fit$selected$support.size[[1L]] &
            ct$degree == fit$selected$degree[[1L]] &
            ct$kernel == fit$selected$kernel[[1L]]
    )
    contrib15 <- function(p, y) {
        p <- pmin(1 - 1e-15, pmax(1e-15, p))
        -(y * log(p) + (1 - y) * log1p(-p))
    }
    contribs <- contrib15(pred[, selected.row], y.flipped)
    expect_equal(contribs[[flip]], -log(1e-15), tolerance = 1e-9)
    expect_gt(contribs[[flip]] / max(contribs[-flip]), 5)
})
