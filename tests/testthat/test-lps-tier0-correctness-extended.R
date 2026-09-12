# =============================================================================
# Tier-0 statistical-correctness battery, continued: E0.3a, E0.4, E0.5, E0.6, E0.7
#
# Companion to test-lps-tier0-correctness.R (E0.1, E0.2). These gates implement
# the FROZEN spec in dev/methods/lps/specs/lps_experimental_plan_2026-06-09.tex
# (sections E0.3a..E0.7). Acceptance thresholds below are copied verbatim from
# that spec; do not relax them without a spec change (post-hoc threshold changes
# invalidate the gate).
#
# SCAFFOLD STATUS: authored against the confirmed fit.lps / predict.lps API but
# NOT yet executed by the author (no R in the authoring environment). E0.3a/E0.4/
# E0.7 are exact deterministic checks. E0.5/E0.6 are accuracy studies; they run
# at SMOKE size by default and at full frozen-spec size when LPS_TIER0_FULL is
# set. The slope CI is fit over PER-REPLICATE points (not the spec's averaged
# points) purely so the CI is stable at smoke size; this is stricter, not looser.
# If a smoke CI is still wide on first run, raise R rather than loosening -0.1.
#
# Self-contained on purpose: testthat does not share helpers across test files,
# so each gate calls fit.lps / predict.lps directly.
# =============================================================================

# -- E0.3a --------------------------------------------------------------------
test_that("E0.3a response-removal LOO shortcut matches the code path, or defers", {
    # Spec E0.3a: compare the code-reported per-point LOO/GCV residual to the
    # analytic shortcut (y_i - yhat_i)/(1 - S_ii) built from an INDEPENDENT S.
    # If fit.lps exposes no separate GCV/LOO path, the test has no independent
    # target and is DEFERRED (folded into E0.2). It must NOT be reconstructed
    # from S and checked against itself -- that is tautological and prohibited.
    set.seed(9301)
    n <- 36L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    y <- stats::rnorm(n)
    fixed.fit <- function(yy) fit.lps(
        X = X, y = yy, foldid = rep(1:2, length.out = n),
        support.grid = 18L, degree.grid = 1L, kernel.grid = "tricube",
        coordinate.method = "coordinates", backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0, ridge.condition.max = Inf,
        unstable.action = "na")
    fit <- fixed.fit(y)

    loo.field <- intersect(
        c("loo.residuals", "loo", "gcv.residuals", "press.residuals",
          "cv.residuals", "leave.one.out", "loocv.residuals"),
        names(fit))
    if (!length(loo.field)) {
        skip(paste0(
            "E0.3a deferred: fit.lps exposes no independent per-point GCV/LOO ",
            "residual path (checked object names). Per frozen spec E0.3a this ",
            "folds into E0.2; reconstructing the shortcut from S and checking ",
            "it against itself would be tautological."))
    }

    # Future-proof branch: an independent LOO/GCV path appeared -> verify it.
    S <- matrix(NA_real_, n, n)
    for (j in seq_len(n)) {
        ej <- numeric(n); ej[[j]] <- 1
        S[, j] <- fixed.fit(ej)$fitted.values
    }
    sii <- diag(S)
    keep <- sii < 1 - 1e-6                       # exclude near-interpolation rows
    shortcut <- (y - fit$fitted.values) / (1 - sii)
    reported <- as.numeric(fit[[loo.field[[1L]]]])
    rel.err <- max(abs(reported[keep] - shortcut[keep]) /
                       (abs(shortcut[keep]) + 1e-12))
    expect_lt(rel.err, 1e-8)
})

# -- E0.4 ---------------------------------------------------------------------
test_that("E0.4 local-linear is boundary-bias-free where local-constant is not", {
    # Spec E0.4: degree-1 is first-order boundary-bias-free; degree-0 carries
    # O(h) boundary bias. Noiseless f(x)=e^x on [0,1], n=300, K=30, gaussian.
    set.seed(9401)
    n <- 300L
    grid <- seq(0, 1, length.out = n)
    jitter <- stats::runif(n, -1, 1) * (0.2 / n)    # fixed small jitter (seed set)
    x <- pmin(pmax(grid + jitter, 0), 1)
    X <- matrix(x, ncol = 1L)
    f <- exp(x)                                     # nonzero slope & curvature at edges
    y <- f                                          # noiseless: isolates bias
    K <- 30L
    fit.deg <- function(p) fit.lps(
        X = X, y = y, foldid = rep(1:2, length.out = n),
        support.grid = K, degree.grid = p, kernel.grid = "gaussian",
        coordinate.method = "coordinates", backend = "R",
        design.basis = "monomial", ridge.multiplier.grid = 0,
        ridge.condition.max = Inf, unstable.action = "na")
    e0 <- abs(fit.deg(0L)$fitted.values - f)
    e1 <- abs(fit.deg(1L)$fitted.values - f)

    # h = realized K-NN bandwidth at the edge (K-th smallest distance from x_min)
    h <- sort(abs(x - min(x)))[K]
    B <- x < h | x > 1 - h
    I <- !B
    eB0 <- mean(e0[B]); eB1 <- mean(e1[B])
    eI0 <- mean(e0[I]); eI1 <- mean(e1[I])

    # Frozen spec acceptance (E0.4) -- verbatim thresholds:
    expect_lt(eB1 / eB0, 0.5)        # degree-1 halves (or better) the boundary error
    expect_lt(eB1 / eI1, 3)          # degree-1 boundary error ~ interior error
    expect_gt(eB0 / eI0, 5)          # degree-0 boundary error blows up vs interior
})

# -- E0.5 ---------------------------------------------------------------------
test_that("E0.5 Truth-RMSE is consistent (slope < -0.1) on a curved 2-manifold", {
    # Spec E0.5: with k(n) on the optimal growth schedule, RMSE decreases at the
    # local-linear manifold rate. Required smoke criterion: the 95% CI for the
    # log-RMSE-vs-log-n slope lies entirely below -0.1.
    full <- nzchar(Sys.getenv("LPS_TIER0_FULL"))
    ns <- if (full) c(200L, 400L, 800L, 1600L, 3200L) else c(200L, 400L, 800L, 1600L)
    R  <- if (full) 30L else 12L
    d <- 2L
    c.k <- 20 / (200^(4 / (d + 4)))                 # pin k(200)=20; k(n)=c*n^{2/3}

    gen <- function(n, seed) {
        set.seed(seed)
        u <- matrix(stats::runif(n * 2L, -1, 1), ncol = 2L)   # intrinsic coords
        amb <- cbind(u[, 1L], u[, 2L], 0.4 * (u[, 1L]^2 - u[, 2L]^2))  # curved saddle
        f <- sin(pi * u[, 1L]) * cos(0.5 * pi * u[, 2L])      # smooth truth
        list(X = amb, f = f, y = f + stats::rnorm(n, sd = 0.1))
    }
    one.rmse <- function(n, r) {
        k <- max(8L, as.integer(round(c.k * n^(4 / (d + 4)))))
        g <- gen(n, 50000L + n + r)
        fit <- fit.lps(
            X = g$X, y = g$y, foldid = rep(1:2, length.out = n),
            support.grid = k, degree.grid = 1L, kernel.grid = "gaussian",
            coordinate.method = "local.pca", chart.dim = 2L,
            local.chart.method = "pca", backend = "R",
            design.basis = "orthogonal.polynomial.drop",
            ridge.multiplier.grid = 0, ridge.condition.max = Inf,
            unstable.action = "na")
        sqrt(mean((fit$fitted.values - g$f)^2, na.rm = TRUE))
    }
    # per-replicate points -> stable slope CI even at smoke size
    pts <- do.call(rbind, lapply(ns, function(n) {
        rs <- vapply(seq_len(R), function(r) one.rmse(n, r), numeric(1L))
        data.frame(logn = log(n), logr = log(rs))
    }))
    lmf <- stats::lm(logr ~ logn, data = pts)
    b  <- stats::coef(lmf)[["logn"]]
    se <- summary(lmf)$coefficients["logn", 2L]
    expect_lt(b + 1.96 * se, -0.1)   # frozen spec smoke criterion
})

# -- E0.6 ---------------------------------------------------------------------
test_that("E0.6 binary modes recover and calibrate probabilities", {
    # Spec E0.6: both binary modes give probability estimates that are (i)
    # consistent (RMSE_p decreasing in n) and (ii) calibrated at large n.
    #
    # Tier-0 amendment (2026-06-12, orchestrator-adjudicated; raised by the
    # Tier-2 implementer in dev/methods/lps/audit_contracts/tiers1to4/
    # e2_15_e06_interaction_raise_2026-06-12.md): the BINOMIAL arm runs
    # unstable.action = "mean" (deployed event-rate fallback predictions at
    # every point) instead of "na". Rationale: (i) under E2.15 a candidate
    # with any non-finite CV prediction is unselectable, and at small
    # supports / low prevalence every candidate can have NA predictions
    # (logistic exact separation), erroring the fit -- 7 of 72 smoke cells
    # pre-amendment; (ii) the accepted "na"-era statistics were computed
    # over selections E2.15 classifies as invalid (candidates scored only
    # on the points they happened to predict, failing on up to ~13%);
    # (iii) scoring the DEPLOYED predictions at every point is the same
    # train/deploy-consistency principle as E2.12, with the fallback
    # fraction telemetered per E2.14. The bernoulli arm keeps "na" (its
    # least-squares solve produces no NA on these fixtures; realized
    # na.fraction 0 in the accepted bundles).
    full <- nzchar(Sys.getenv("LPS_TIER0_FULL"))
    ns <- if (full) c(500L, 1000L, 2000L, 4000L) else c(500L, 1000L, 2000L)
    R  <- if (full) 40L else 8L
    prevalence.grid <- c(0.1, 0.3, 0.5)
    d <- 2L
    # Binary probability recovery needs a slightly wider local averaging
    # scale than the continuous E0.5 regression gate.  The rate is still the
    # frozen nonparametric schedule k(n) proportional to n^{4/(d+4)}; the
    # constant is set so the n=500 reference support is 30 rather than 20,
    # avoiding the full-size calibration overconfidence seen with too-small
    # supports.
    c.k <- 30 / (500^(4 / (d + 4)))
    support.grid.for <- function(n) {
        k <- max(12L, as.integer(round(c.k * n^(4 / (d + 4)))))
        sort(unique(pmax(8L, as.integer(round(k * c(0.75, 1, 1.25))))))
    }

    gen <- function(n, seed, prev) {
        set.seed(seed)
        X <- matrix(stats::runif(n * 2L, -1, 1), ncol = 2L)
        lin <- 1.2 * X[, 1L] - 0.8 * X[, 2L] + 0.6 * X[, 1L] * X[, 2L]
        a0 <- stats::uniroot(function(a) mean(stats::plogis(a + lin)) - prev,
                             c(-12, 12))$root
        p <- pmin(pmax(stats::plogis(a0 + lin), 0.05), 0.95)   # spec: p in [0.05,0.95]
        list(X = X, p = p, y = stats::rbinom(n, 1L, p))
    }
    fit.bin <- function(g, fam) fit.lps(
        X = g$X, y = g$y, foldid = rep(1:5, length.out = length(g$y)),
        support.grid = support.grid.for(length(g$y)), degree.grid = c(0L, 1L),
        kernel.grid = "tricube", coordinate.method = "coordinates",
        backend = "R", design.basis = "orthogonal.polynomial.drop",
        outcome.family = fam, ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = if (identical(fam, "binomial")) "mean" else "na")

    for (fam in c("bernoulli", "binomial")) {
        for (prev in prevalence.grid) {
            pts <- do.call(rbind, lapply(ns, function(n) {
                rows <- lapply(seq_len(R), function(r) {
                    g <- gen(n, 60000L + 1000L * round(10 * prev) + n + r,
                             prev = prev)
                    fit <- fit.bin(g, fam)
                    ph <- fit$fitted.values
                    finite <- is.finite(ph)
                    expect_true(all(ph[finite] >= 0 & ph[finite] <= 1))
                    fallback.frac <- if (identical(fam, "binomial")) {
                        fit$logistic.diagnostics$final$fallback.path.fraction
                    } else {
                        NA_real_
                    }
                    data.frame(
                        family = fam,
                        prevalence = prev,
                        n = n,
                        logn = log(n),
                        rmse = sqrt(mean((ph[finite] - g$p[finite])^2)),
                        na.fraction = mean(!finite),
                        fallback.fraction = fallback.frac
                    )
                })
                do.call(rbind, rows)
            }))
            expect_lt(max(pts$na.fraction, na.rm = TRUE), 0.25)
            pts <- pts[is.finite(pts$rmse) & pts$na.fraction < 0.25, , drop = FALSE]
            lmf <- stats::lm(log(rmse) ~ logn, data = pts)
            ci.hi <- stats::coef(lmf)[["logn"]] +
                1.96 * summary(lmf)$coefficients["logn", 2L]
            median.fallback <- stats::median(pts$fallback.fraction, na.rm = TRUE)
            cat(sprintf(
                "E0.6 family=%s prevalence=%.1f support=%s slope=%.4f ci_hi=%.4f max_na=%.4f median_fallback=%.4f\n",
                fam, prev, paste(vapply(ns, function(n) {
                    paste(support.grid.for(n), collapse = "/")
                }, character(1L)), collapse = ","),
                stats::coef(lmf)[["logn"]], ci.hi,
                max(pts$na.fraction, na.rm = TRUE),
                median.fallback
            ))
            expect_lt(ci.hi, -0.1)    # consistency slope CI below -0.1
            # E0.6 fallback-fraction bound (BINOMIAL arm only). The binomial arm
            # runs unstable.action = "mean" (Tier-0 amendment above), so every
            # non-converged logistic solve deploys an event-rate fallback
            # prediction; the per-fit fallback fraction is telemetered in
            # fit$logistic.diagnostics$final$fallback.path.fraction. This bound
            # promotes the already-printed median_fallback telemetry into a gate:
            # it asserts the deployed fallback stays a small minority of final
            # fits. Realized per-prevalence medians are ~0.0155 (smoke) /
            # ~0.0020 (full-size) at prevalence 0.1 and 0 elsewhere; the 0.3 bound
            # clears the real code by ~20x while reddening a degenerate
            # all-fallback regime (median 1.0). The bernoulli arm records NA
            # fallback (least-squares solve, no logistic path), so the assertion
            # is binomial-only -- a bernoulli median over all-NA is NA, which is
            # not a meaningful bound.
            if (identical(fam, "binomial")) {
                expect_lt(median.fallback, 0.3)
            }
        }
    }

    # Calibration on a held-out half (spec: slope/intercept of y ~ logit(phat)).
    nbig <- ns[length(ns)]
    g <- gen(nbig, 69999L, prev = 0.3)
    half <- (seq_len(nbig) %% 2L) == 0L
    fit <- fit.bin(list(X = g$X[!half, , drop = FALSE], y = g$y[!half]), "bernoulli")
    ph <- as.numeric(predict(fit, newdata = g$X[half, , drop = FALSE]))
    finite <- is.finite(ph)
    expect_gt(mean(finite), 0.9)
    ph <- ph[finite]
    y.cal <- g$y[half][finite]
    ph <- pmin(pmax(ph, 1e-4), 1 - 1e-4)
    cal <- stats::glm(y.cal ~ stats::qlogis(ph), family = stats::binomial())
    slope <- stats::coef(cal)[[2L]]
    intercept <- stats::coef(cal)[[1L]]
    if (full) {
        # frozen spec acceptance is stated at n=4000:
        expect_gt(slope, 0.8);  expect_lt(slope, 1.25)
        expect_gt(intercept, -0.25); expect_lt(intercept, 0.25)
    } else {
        # smoke: calibration is meaningful and positively sloped (looser band).
        expect_gt(slope, 0.3)
        expect_lt(abs(intercept), 0.6)
    }
})

# -- E0.7 ---------------------------------------------------------------------
test_that("E0.7 cross-validation has no held-out response leakage", {
    # Spec E0.7: perturbing a held-out y_i leaves i's own out-of-fold prediction
    # exactly unchanged; as a positive control it DOES change a j in another fold
    # whose training support contains i (proving the test is not vacuous).
    set.seed(9701)
    n <- 100L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    y <- as.numeric(cbind(1, X) %*% c(0.7, 0.4, -0.3)) + stats::rnorm(n, sd = 0.1)
    foldid <- rep(1:5, length.out = n)

    oof.pred <- function(yvec, idx) {
        f <- foldid[idx]
        tr <- which(foldid != f)
        fit <- fit.lps(
            X = X[tr, , drop = FALSE], y = yvec[tr],
            foldid = rep(1:2, length.out = length(tr)),
            support.grid = 24L, degree.grid = 1L, kernel.grid = "tricube",
            coordinate.method = "coordinates", backend = "R",
            design.basis = "orthogonal.polynomial.drop",
            ridge.multiplier.grid = 0, ridge.condition.max = Inf,
            unstable.action = "na")
        as.numeric(predict(fit, newdata = X[idx, , drop = FALSE]))
    }

    i <- 7L
    y.pert <- y; y.pert[[i]] <- y.pert[[i]] + 10   # Delta = 10
    base.i <- oof.pred(y, i)
    pert.i <- oof.pred(y.pert, i)
    expect_lt(abs(pert.i - base.i), 1e-12)         # no leakage into i's own OOF pred

    # choose j as i's NEAREST different-fold neighbour, so i is in j's k-NN
    # support and perturbing y_i must move j's OOF prediction (non-vacuous control)
    diff.fold <- which(foldid != foldid[[i]])
    d2 <- rowSums(sweep(X[diff.fold, , drop = FALSE], 2L, X[i, ], "-")^2)
    j <- diff.fold[which.min(d2)]
    base.j <- oof.pred(y, j)                        # j's training set includes i
    pert.j <- oof.pred(y.pert, j)
    expect_gt(abs(pert.j - base.j), 0)             # positive control: not vacuous
})
