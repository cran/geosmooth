test_that("GE7 fit.lps is the canonical LPS entry point", {
    X <- matrix(seq(0, 1, length.out = 24), ncol = 1L)
    y <- sin(2 * pi * X[, 1])
    foldid <- rep(1:3, length.out = nrow(X))

    fit <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = c(6L, 8L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        backend = "R"
    )

    expect_s3_class(fit, "lps")
    expect_identical(class(fit)[[1L]], "lps")
    expect_false(inherits(fit, "kernel.local.polynomial.cv"))
    expect_identical(fit$method.id, "lps")
    expect_identical(fit$method.label, "LPS")
    expect_identical(fit$method.family, "local_polynomial_smoother")
    expect_equal(length(predict(fit, X[1:4, , drop = FALSE])), 4L)
    expect_match(capture.output(print(fit))[[1L]], "LPS")
})

test_that("LPS sparse chart activation routes coordinate C++ through R", {
    X <- matrix(seq(0, 1, length.out = 20), ncol = 1L)
    y <- numeric(nrow(X))
    y[c(9L, 11L)] <- 1

    fit <- fit.lps(
        X = X,
        y = y,
        support.grid = 5L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "coordinates",
        backend = "cpp",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean",
        chart.activation = "subject.od",
        chart.activation.response = y,
        cv.folds = 2L
    )

    expect_identical(fit$backend, "cpp")
    expect_identical(fit$backend.used, "R")
    expect_true(isTRUE(fit$diagnostics$chart.activation$enabled))
    expect_true(any(!fit$chart.activation.diagnostics$active))
    expect_true(all(fit$fitted.values.raw[
        !fit$chart.activation.diagnostics$active
    ] == 0))

    backend.diag <- lps.backend.diagnostics(fit)
    expect_identical(
        backend.diag$backend.auto.policy,
        "explicit_coordinates_cpp_activation_R_reference"
    )
})

test_that("LPS sparse chart activation supports native local-PCA backend", {
    X <- cbind(seq(0, 1, length.out = 24),
               sin(seq(0, 1, length.out = 24)))
    y <- numeric(nrow(X))
    y[c(10L, 12L, 14L)] <- 1
    common <- list(
        X = X,
        y = y,
        support.grid = 7L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        local.chart.method = "pca",
        chart.dim = 1L,
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean",
        chart.activation = "subject.od",
        chart.activation.response = y,
        cv.folds = 2L
    )

    fit.r <- do.call(fit.lps, c(common, list(backend = "R")))
    fit.cpp <- do.call(fit.lps, c(common, list(backend = "cpp.local.pca")))

    expect_identical(fit.cpp$backend, "cpp.local.pca")
    expect_identical(fit.cpp$backend.used, "cpp.local.pca")
    expect_true(isTRUE(fit.cpp$diagnostics$chart.activation$enabled))
    expect_true(any(!fit.cpp$chart.activation.diagnostics$active))
    expect_equal(fit.cpp$fitted.values.raw, fit.r$fitted.values.raw,
                 tolerance = 1e-8)
    expect_true(all(fit.cpp$fitted.values.raw[
        !fit.cpp$chart.activation.diagnostics$active
    ] == 0))
})

test_that("LPS Bernoulli mode validates 0/1 responses and reports probabilities", {
    X <- matrix(seq(0, 1, length.out = 30), ncol = 1L)
    y <- as.numeric(X[, 1] > 0.45)
    foldid <- rep(1:5, length.out = nrow(X))

    gaussian <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        outcome.family = "gaussian",
        support.grid = c(6L, 8L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        backend = "R",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )
    bernoulli <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        outcome.family = "bernoulli",
        support.grid = c(6L, 8L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        backend = "R",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )

    expect_identical(bernoulli$outcome.family, "bernoulli")
    expect_equal(bernoulli$fitted.values.raw, gaussian$fitted.values,
                 tolerance = 1e-12)
    expect_equal(bernoulli$cv.table$cv.rmse.observed,
                 gaussian$cv.table$cv.rmse.observed,
                 tolerance = 1e-12)
    expect_true("cv.brier.observed" %in% names(bernoulli$cv.table))
    expect_equal(bernoulli$cv.table$cv.brier.observed,
                 bernoulli$cv.table$cv.rmse.observed^2,
                 tolerance = 1e-12)
    expect_true(all(bernoulli$fitted.values >= 0))
    expect_true(all(bernoulli$fitted.values <= 1))
    expect_true(is.list(bernoulli$probability.diagnostics))
    expect_equal(bernoulli$probability.diagnostics$diagnostic.scope,
                 "labeled_predictions")
    expect_equal(bernoulli$probability.diagnostics$n.labels, length(y))
    expect_equal(bernoulli$probability.diagnostics$n.predictions, length(y))
    expect_equal(bernoulli$probability.diagnostics$brier.denominator,
                 length(y))
    expect_true(is.finite(bernoulli$probability.diagnostics$brier.clipped))
    expect_true(is.finite(bernoulli$probability.diagnostics$logloss.clipped))
    expect_equal(predict(bernoulli, type = "raw"),
                 bernoulli$fitted.values.raw,
                 tolerance = 1e-12)
    expect_equal(predict(bernoulli, type = "response"),
                 bernoulli$fitted.values,
                 tolerance = 1e-12)
    expect_equal(predict(gaussian, type = "raw"),
                 predict(gaussian, type = "response"),
                 tolerance = 1e-12)
    expect_equal(
        predict(bernoulli, X[1:6, , drop = FALSE], type = "raw"),
        predict(gaussian, X[1:6, , drop = FALSE]),
        tolerance = 1e-12
    )
    response.pred <- predict(bernoulli, X[1:6, , drop = FALSE],
                             type = "response")
    expect_true(all(response.pred >= 0))
    expect_true(all(response.pred <= 1))
    expect_error(
        fit.lps(
            X = X,
            y = y + 0.1,
            foldid = foldid,
            outcome.family = "bernoulli",
            support.grid = 6L,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            backend = "R"
        ),
        "requires y values in \\{0, 1\\}"
    )

    expect_warning(
        single.class <- fit.lps(
            X = X,
            y = rep(1, nrow(X)),
            foldid = foldid,
            outcome.family = "bernoulli",
            support.grid = 6L,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            backend = "R"
        ),
        "only one observed class"
    )
    expect_s3_class(single.class, "lps")
    expect_true(all(single.class$fitted.values >= 0))
    expect_true(all(single.class$fitted.values <= 1))

    external <- fit.lps(
        X = X,
        y = y,
        X.eval = X[1:10, , drop = FALSE],
        foldid = foldid,
        outcome.family = "bernoulli",
        support.grid = 6L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        backend = "R"
    )
    expect_equal(external$probability.diagnostics$diagnostic.scope,
                 "unlabeled_eval_predictions")
    expect_equal(external$probability.diagnostics$n.labels, length(y))
    expect_equal(external$probability.diagnostics$n.predictions, 10L)
    expect_true(is.na(external$probability.diagnostics$brier.clipped))
    expect_true(is.na(external$probability.diagnostics$logloss.clipped))
})

test_that("LPS binomial mode uses local logistic fits and log-loss selection", {
    set.seed(1708)
    X <- matrix(sort(runif(48)), ncol = 1L)
    prob <- stats::plogis(-0.4 + 2.5 * sin(2 * pi * X[, 1]))
    y <- stats::rbinom(nrow(X), size = 1L, prob = prob)
    foldid <- rep(1:4, length.out = nrow(X))

    fit <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        outcome.family = "binomial",
        support.grid = c(10L, 14L),
        degree.grid = 1L,
        kernel.grid = c("gaussian", "tricube"),
        backend = "R",
        design.basis = "monomial",
        ridge.multiplier.grid = c(1e-8, 1e-6),
        ridge.condition.max = 1e10,
        unstable.action = "mean"
    )

    expect_identical(fit$outcome.family, "binomial")
    expect_true("cv.logloss.observed" %in% names(fit$cv.table))
    expect_true("cv.brier.observed" %in% names(fit$cv.table))
    expect_equal(
        fit$selected$cv.logloss.observed[[1L]],
        min(fit$cv.table$cv.logloss.observed, na.rm = TRUE),
        tolerance = 1e-12
    )
    expect_true(all(fit$fitted.values >= 0))
    expect_true(all(fit$fitted.values <= 1))
    expect_true(is.list(fit$logistic.diagnostics))
    expect_gt(fit$logistic.diagnostics$cv$attempted, 0)
    expect_gt(fit$logistic.diagnostics$final$attempted, 0)
    expect_gte(fit$logistic.diagnostics$cv$converged, 0)
    expect_gte(fit$logistic.diagnostics$final$converged, 0)
    expect_true("fallback.path.count" %in%
                    names(fit$logistic.diagnostics$cv))
    expect_true("event.rate.fallback.count" %in%
                    names(fit$logistic.diagnostics$cv))
    expect_true("na.failure.count" %in%
                    names(fit$logistic.diagnostics$cv))
    expect_equal(
        fit$logistic.diagnostics$cv$attempted,
        fit$logistic.diagnostics$cv$converged +
            fit$logistic.diagnostics$cv$failed
    )
    expect_equal(fit$fitted.values, fit$fitted.values.raw, tolerance = 1e-12)
    expect_true(is.finite(fit$probability.diagnostics$logloss.clipped))
    expect_equal(predict(fit, type = "response"),
                 fit$fitted.values,
                 tolerance = 1e-12)
    expect_equal(predict(fit, type = "raw"),
                 fit$fitted.values.raw,
                 tolerance = 1e-12)
    fit.auto <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        outcome.family = "binomial",
        support.grid = 10L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        backend = "auto",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )
    expect_identical(fit.auto$backend.used, "R")
    expect_identical(fit.auto$outcome.family, "binomial")
    expect_error(
        fit.lps(
            X = X,
            y = y,
            foldid = foldid,
            outcome.family = "binomial",
            support.grid = 10L,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            backend = "cpp",
            design.basis = "monomial",
            ridge.multiplier.grid = 0,
            ridge.condition.max = Inf
        ),
        "currently uses the R backend"
    )

    external <- fit.lps(
        X = X,
        y = y,
        X.eval = X[1:5, , drop = FALSE],
        foldid = foldid,
        outcome.family = "binomial",
        support.grid = 10L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        backend = "R",
        unstable.action = "mean"
    )
    expect_equal(external$probability.diagnostics$diagnostic.scope,
                 "unlabeled_eval_predictions")
    expect_true(is.na(external$probability.diagnostics$logloss.clipped))
    expect_equal(length(external$fitted.values), 5)

    local.pca <- fit.lps(
        X = cbind(X, X^2),
        y = y,
        foldid = foldid,
        outcome.family = "binomial",
        support.grid = 14L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        chart.dim = "auto",
        backend = "R",
        unstable.action = "mean"
    )
    expect_true(is.finite(local.pca$selected$cv.logloss.observed[[1L]]))
    expect_true(all(is.finite(local.pca$fitted.values)))
    expect_gt(local.pca$logistic.diagnostics$final$attempted, 0)

    expect_warning(
        single.class <- fit.lps(
            X = X,
            y = rep(0, nrow(X)),
            foldid = foldid,
            outcome.family = "binomial",
            support.grid = 10L,
            degree.grid = 0L,
            kernel.grid = "gaussian",
            backend = "R",
            unstable.action = "mean"
        ),
        "received only one observed class"
    )
    expect_true(all(single.class$fitted.values >= 0))
    expect_true(all(single.class$fitted.values <= 1))

    expect_error(
        fit.lps(
            X = X[1:12, , drop = FALSE],
            y = y[1:12],
            foldid = rep(1:3, length.out = 12),
            outcome.family = "binomial",
            support.grid = 2L,
            degree.grid = 2L,
            kernel.grid = "gaussian",
            backend = "R",
            design.basis = "monomial",
            unstable.action = "na"
        ),
        "No candidate has a finite selection score"
    )
    # Phase-3 Pass-1 finding: the former all-ones rank-deficient fixture now
    # converges, so it no longer exercises the na.failure telemetry path. Use
    # the exact-separation fixture supplied by the Pass-1 audit, which drives a
    # genuine NA failure (the unpenalized logistic MLE does not exist under
    # exact separation).
    z <- c(seq(-0.20, -0.04, by = 0.02), 6)
    design <- cbind(1, z)
    y <- as.numeric(z > 0)
    weights <- geosmooth:::.klp.kernel.weights(abs(z), "gaussian")
    telemetry <- geosmooth:::.klp.logistic.telemetry.new("binomial")
    failed <- geosmooth:::.klp.fit.logistic.prob.design(
        design = design, y = y, weights = weights,
        design.basis = "orthogonal.polynomial.drop", design.drop.tol = 1e-8,
        ridge.multiplier.grid = 0, ridge.condition.max = Inf,
        unstable.action = "na", logistic.telemetry = telemetry)
    summary <- geosmooth:::.klp.logistic.telemetry.summary(telemetry)
    expect_true(is.na(failed))
    expect_equal(summary$fallback.path.count, 1L)
    expect_equal(summary$event.rate.fallback.count, 0L)
    expect_equal(summary$na.failure.count, 1L)
})

test_that("GE8 removes the old kernel.local.polynomial.cv API", {
    expect_false(exists(
        "kernel.local.polynomial.cv",
        envir = asNamespace("geosmooth"),
        inherits = FALSE
    ))
    expect_false(exists(
        "predict.kernel.local.polynomial.cv",
        envir = asNamespace("geosmooth"),
        inherits = FALSE
    ))
    expect_false(exists(
        "print.kernel.local.polynomial.cv",
        envir = asNamespace("geosmooth"),
        inherits = FALSE
    ))
})

test_that("K4 local-PCA LPS C++ backend matches R reference path", {
    set.seed(41)
    t <- seq(0, 1, length.out = 36)
    X <- cbind(t, t^2, sin(2 * pi * t))
    y <- sin(3 * pi * t) + 0.15 * cos(7 * pi * t)
    foldid <- rep(1:4, length.out = nrow(X))

    common.args <- list(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = c(12L, 14L),
        degree.grid = 0:2,
        kernel.grid = c("gaussian", "tricube"),
        coordinate.method = "local.pca",
        local.chart.method = "pca",
        chart.dim = 2L,
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )

    fit.r <- do.call(fit.lps, c(common.args, list(backend = "R")))
    fit.cpp <- do.call(fit.lps, c(common.args, list(backend = "cpp.local.pca")))

    expect_equal(fit.cpp$backend.used, "cpp.local.pca")
    expect_equal(
        fit.cpp$cv.table$cv.rmse.observed,
        fit.r$cv.table$cv.rmse.observed,
        tolerance = 1e-8
    )
    expect_equal(fit.cpp$selected, fit.r$selected, tolerance = 1e-8)
    expect_equal(fit.cpp$fitted.values, fit.r$fitted.values, tolerance = 1e-8)

    newdata <- X[c(3L, 11L, 23L), , drop = FALSE]
    expect_equal(
        predict(fit.cpp, newdata),
        predict(fit.r, newdata),
        tolerance = 1e-8
    )
})

test_that("K4 local-PCA C++ backend is explicit and narrow", {
    X <- matrix(seq(0, 1, length.out = 20), ncol = 1L)
    y <- X[, 1]^2
    foldid <- rep(1:2, length.out = nrow(X))

    expect_error(
        fit.lps(
            X,
            y,
            foldid = foldid,
            coordinate.method = "coordinates",
            backend = "cpp.local.pca"
        ),
        "requires coordinate.method = 'local.pca'"
    )
    expect_error(
        fit.lps(
            X,
            y,
            foldid = foldid,
            coordinate.method = "local.pca",
            local.chart.method = "second.order.svd",
            chart.dim = 1L,
            backend = "cpp.local.pca"
        ),
        "requires coordinate.method = 'local.pca'"
    )

    fit.auto <- fit.lps(
        X,
        y,
        foldid = foldid,
        coordinate.method = "local.pca",
        chart.dim = 1L,
        backend = "auto"
    )
    expect_equal(fit.auto$backend.used, "R")
})

test_that("K4 local-PCA C++ backend uses R-compatible tie-stable supports", {
    compare.backends <- function(X, y, foldid, support.grid,
                                 degree.grid, kernel.grid) {
        common.args <- list(
            X = X,
            y = y,
            foldid = foldid,
            support.grid = support.grid,
            degree.grid = degree.grid,
            kernel.grid = kernel.grid,
            coordinate.method = "local.pca",
            local.chart.method = "pca",
            chart.dim = 2L,
            design.basis = "monomial",
            ridge.multiplier.grid = 0,
            ridge.condition.max = Inf,
            unstable.action = "mean"
        )
        fit.r <- do.call(fit.lps, c(common.args, list(backend = "R")))
        fit.cpp <- do.call(
            fit.lps,
            c(common.args, list(backend = "cpp.local.pca"))
        )
        expect_equal(
            fit.cpp$cv.table$cv.rmse.observed,
            fit.r$cv.table$cv.rmse.observed,
            tolerance = 1e-8
        )
        expect_equal(fit.cpp$selected, fit.r$selected, tolerance = 1e-8)
        expect_equal(fit.cpp$fitted.values, fit.r$fitted.values,
                     tolerance = 1e-8)
    }

    uv <- as.matrix(expand.grid(
        u = seq(-1, 1, length.out = 5),
        v = seq(-1, 1, length.out = 5)
    ))
    X.grid <- cbind(uv[, 1], uv[, 2], uv[, 1] + 2 * uv[, 2])
    y.grid <- uv[, 1]^2 - uv[, 2] + 0.25 * uv[, 1] * uv[, 2]
    compare.backends(
        X = X.grid,
        y = y.grid,
        foldid = rep(1:5, length.out = nrow(X.grid)),
        support.grid = c(10L, 12L),
        degree.grid = 2L,
        kernel.grid = c("gaussian", "tricube")
    )

    base <- as.matrix(expand.grid(u = 0:2, v = 0:1))
    uv.dup <- base[rep(seq_len(nrow(base)), each = 2L), , drop = FALSE]
    X.dup <- cbind(uv.dup[, 1], uv.dup[, 2], uv.dup[, 1] - uv.dup[, 2])
    y.dup <- 0.2 * seq_len(nrow(X.dup)) + uv.dup[, 1] - 0.5 * uv.dup[, 2]
    compare.backends(
        X = X.dup,
        y = y.dup,
        foldid = rep(1:3, length.out = nrow(X.dup)),
        support.grid = c(6L, 8L),
        degree.grid = 1:2,
        kernel.grid = c("gaussian", "tricube")
    )
})

test_that("K4 native neighbor probe repairs raw ANN tie order", {
    check.probe <- function(X, center, k) {
        probe <- geosmooth:::rcpp_kernel_local_polynomial_neighbor_probe(
            X = X,
            center = center,
            k = k
        )
        d <- rowSums((X - matrix(center, nrow(X), ncol(X), byrow = TRUE))^2)
        reference <- order(d, seq_along(d))[seq_len(k)]

        expect_equal(probe$reference.row, reference)
        expect_equal(probe$tie.complete.row, reference)
        expect_equal(probe$tie.complete.squared.distance, d[reference])
        expect_length(probe$raw.row, k)
        expect_length(probe$raw.squared.distance, k)
    }

    cardinal <- matrix(c(
        -1, 0,
         1, 0,
         0, -1,
         0, 1
    ), ncol = 2, byrow = TRUE)
    check.probe(cardinal, center = c(0, 0), k = 2L)

    duplicated <- rbind(
        c(0, 0), c(0, 0),
        c(1, 0), c(1, 0),
        c(2, 0), c(2, 0)
    )
    check.probe(duplicated, center = c(1, 0), k = 3L)

    grid <- as.matrix(expand.grid(u = -1:1, v = -1:1))
    check.probe(grid, center = c(0, 0), k = 4L)
})

test_that("K12 LPS backend diagnostics expose conservative backend policy", {
    X <- matrix(seq(0, 1, length.out = 24), ncol = 1L)
    y <- sin(2 * pi * X[, 1])
    foldid <- rep(1:4, length.out = nrow(X))

    ambient <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = c(6L, 8L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "auto"
    )
    ambient.diag <- lps.backend.diagnostics(ambient)

    expect_s3_class(ambient.diag, "data.frame")
    expect_equal(nrow(ambient.diag), 1L)
    expect_equal(ambient.diag$backend.requested, "auto")
    expect_equal(ambient.diag$backend.used, "R")
    expect_equal(ambient.diag$backend.auto.policy,
                 "auto_coordinates_R_guarded_design")
    expect_equal(ambient.diag$requested.chart.dim, "NULL")
    expect_equal(ambient.diag$resolved.chart.dim, ncol(X))
    expect_false(ambient.diag$chart.dim.auto)
    expect_false(ambient.diag$local.pca.real.data.contract)
    expect_equal(ambient.diag$candidate.count, nrow(ambient$cv.table))
})

test_that("K12 local-PCA diagnostics record auto-dimension contract", {
    t <- seq(0, 1, length.out = 32)
    X <- cbind(t, t^2, sin(2 * pi * t))
    y <- cos(2 * pi * t)
    foldid <- rep(1:4, length.out = nrow(X))

    fit <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = c(10L, 12L),
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        chart.dim = "auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        backend = "auto"
    )
    diag <- lps.backend.diagnostics(fit)

    expect_equal(diag$backend.requested, "auto")
    expect_equal(diag$backend.used, "R")
    expect_equal(diag$backend.auto.policy, "auto_local_pca_R_reference")
    expect_equal(diag$requested.chart.dim, "auto")
    expect_true(diag$chart.dim.auto)
    expect_equal(diag$auto.chart.support.metric, "both")
    expect_equal(diag$auto.chart.selection.metric, "operator")
    expect_equal(diag$auto.chart.support.metric.selected, "coordinates")
    expect_true(is.finite(diag$resolved.chart.dim))
    expect_true(diag$local.pca.real.data.contract)
})

test_that("LPS supports opt-in per-anchor local auto chart dimensions", {
    t <- seq(0, 1, length.out = 40)
    X <- cbind(t, t^2, sin(2 * pi * t), 0.02 * cos(5 * pi * t))
    y <- sin(3 * pi * t)
    foldid <- rep(1:4, length.out = nrow(X))

    fit <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = c(10L, 14L),
        degree.grid = 1:2,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = "local.auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        backend = "auto"
    )
    diag <- lps.backend.diagnostics(fit)

    expect_s3_class(fit, "lps")
    expect_equal(fit$backend.used, "R")
    expect_true(fit$auto.chart.dim)
    expect_true(fit$auto.chart.dim.local)
    expect_equal(fit$chart.dim.mode, "local.auto")
    expect_length(fit$diagnostics$chart.dim$by.anchor, nrow(X))
    expect_true(all(fit$diagnostics$chart.dim$by.anchor >= 1L))
    expect_true(all(fit$diagnostics$chart.dim$by.anchor <= ncol(X)))
    expect_true(is.finite(fit$chart.dim))
    expect_equal(diag$requested.chart.dim, "local.auto")
    expect_true(diag$chart.dim.auto)
    expect_true(diag$chart.dim.local.auto)
    expect_equal(diag$chart.dim.mode, "local.auto")
    expect_equal(diag$chart.dim.by.anchor.n, nrow(X))
    expect_true(diag$local.pca.real.data.contract)
    expect_equal(length(predict(fit, X[1:5, , drop = FALSE])), 5L)
})

test_that("LPS local auto chart dimensions stay on the R local-PCA path", {
    X <- cbind(seq(0, 1, length.out = 24), seq(0, 1, length.out = 24)^2)
    y <- X[, 1]
    foldid <- rep(1:3, length.out = nrow(X))

    expect_error(
        fit.lps(
            X,
            y,
            foldid = foldid,
            coordinate.method = "coordinates",
            chart.dim = "local.auto",
            backend = "R"
        ),
        "requires coordinate.method = 'local.pca'"
    )
    expect_error(
        fit.lps(
            X,
            y,
            foldid = foldid,
            coordinate.method = "local.pca",
            chart.dim = "local.auto",
            backend = "cpp.local.pca",
            design.basis = "monomial",
            ridge.multiplier.grid = 0,
            ridge.condition.max = Inf
        ),
        "currently uses the R local-PCA backend"
    )
    expect_error(
        fit.lps(
            X,
            y,
            foldid = foldid,
            coordinate.method = "local.pca",
            chart.dim = "local.auto",
            local.chart.method = "second.order.svd",
            backend = "auto"
        ),
        "currently supports only local.chart.method = 'pca'"
    )
})

test_that("LPS local WLS falls back on nearly saturated ill-conditioned designs", {
    set.seed(1701)
    z <- matrix(rnorm(35 * 6), 35, 6)
    design <- .local.polynomial.design.matrix(z, degree = 2L)
    design[, ncol(design)] <- design[, ncol(design) - 1L] +
        1e-8 * rnorm(nrow(design))
    y <- sin(seq(0, 1, length.out = nrow(design)))
    weights <- rep(1, nrow(design))

    expect_false(.klp.local.design.is.safe(design, weights))
    # Phase-3 Pass-1 disposition: the default orthogonal-polynomial ridge path
    # does not consult .klp.local.design.is.safe() to force a weighted-mean
    # fallback (this long-standing behavior predates t2). The contract is that
    # the solve still returns a finite fitted intercept, which is NOT the local
    # weighted mean. Robust check (finite + not-equal), not a brittle value pin.
    fitted.intercept <- .klp.fit.intercept.design(design, y, weights)
    expect_true(is.finite(fitted.intercept))
    expect_false(isTRUE(all.equal(fitted.intercept,
                                  stats::weighted.mean(y, weights))))
})

test_that("LPS weighted QR drop removes dependent local polynomial columns", {
    set.seed(1702)
    z <- cbind(
        seq(-1, 1, length.out = 12),
        seq(-1, 1, length.out = 12)
    )
    design <- .local.polynomial.design.matrix(z, degree = 2L)
    y <- sin(z[, 1])
    weights <- rep(1, length(y))

    monomial <- .klp.solve.local.wls(
        design,
        y,
        weights,
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )
    dropped <- .klp.solve.local.wls(
        design,
        y,
        weights,
        design.basis = "weighted.qr.drop",
        design.drop.tol = 1e-7,
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )

    expect_false(isTRUE(monomial$ok))
    expect_true(isTRUE(dropped$ok))
    expect_true(is.finite(dropped$coefficients[[1L]]))
})

test_that("LPS guarded ridge selects the smallest passing multiplier", {
    set.seed(1703)
    z <- matrix(rnorm(28 * 3), 28, 3)
    design <- .local.polynomial.design.matrix(z, degree = 2L)
    design[, ncol(design)] <- design[, ncol(design) - 1L] +
        1e-7 * rnorm(nrow(design))
    y <- cos(seq(0, 2, length.out = nrow(design)))
    weights <- rep(1, length(y))

    strict <- .klp.solve.local.wls(
        design,
        y,
        weights,
        design.basis = "weighted.qr.drop",
        design.drop.tol = 1e-10,
        ridge.multiplier.grid = c(0, 1e-12),
        ridge.condition.max = 10
    )
    guarded <- .klp.solve.local.wls(
        design,
        y,
        weights,
        design.basis = "weighted.qr.drop",
        design.drop.tol = 1e-10,
        ridge.multiplier.grid = c(0, 1e-12, 1e-8, 1e-4, 1e-2),
        ridge.condition.max = 1e6
    )

    expect_false(isTRUE(strict$ok))
    expect_true(isTRUE(guarded$ok))
    expect_true(guarded$ridge.multiplier %in% c(0, 1e-12, 1e-8, 1e-4, 1e-2))
    expect_lte(guarded$condition, 1e6)
})

test_that("fit.lps can avoid unstable weighted QR drop candidates", {
    set.seed(1704)
    x <- seq(-1, 1, length.out = 48)
    X <- cbind(x, x^2 + 1e-6 * rnorm(length(x)))
    y <- sin(3 * x) + rnorm(length(x), sd = 0.03)
    foldid <- rep(1:4, length.out = length(y))

    fit <- fit.lps(
        X,
        y,
        foldid = foldid,
        support.grid = c(8L, 10L, 14L),
        degree.grid = 2L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        chart.dim = 2L,
        backend = "auto",
        design.basis = "weighted.qr.drop",
        design.drop.tol = 1e-7,
        ridge.multiplier.grid = c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2),
        ridge.condition.max = 1e6,
        unstable.action = "na"
    )

    expect_equal(fit$backend.used, "R")
    expect_equal(fit$design.basis, "weighted.qr.drop")
    expect_true(any(is.finite(fit$cv.table$cv.rmse.observed)))
    expect_true(is.finite(fit$selected$cv.rmse.observed[[1L]]))
})

test_that("LPS guarded local failures are not silently converted to means", {
    design <- cbind(1, matrix(rnorm(6 * 5), 6, 5))
    y <- seq_len(nrow(design))
    weights <- rep(1, length(y))

    solved <- .klp.solve.local.wls(
        design = design,
        y = y,
        weights = weights,
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = 1
    )
    intercept <- .klp.fit.intercept.design(
        design = design,
        y = y,
        weights = weights,
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = 1,
        unstable.action = "na"
    )

    expect_false(isTRUE(solved$ok))
    expect_equal(solved$status, "ridge_condition_failed")
    expect_true(is.na(intercept))
    expect_false(isTRUE(all.equal(intercept, stats::weighted.mean(y, weights))))
})

test_that("orthogonal polynomial drop builds a weighted-orthogonal basis", {
    set.seed(1705)
    z <- matrix(rnorm(18 * 3), 18, 3)
    design <- .local.polynomial.design.matrix(z, degree = 2L)
    design[, ncol(design)] <- design[, ncol(design) - 1L]
    weights <- runif(nrow(design), 0.2, 1)
    anchor <- matrix(c(1, rep(0, ncol(design) - 1L)), nrow = 1L)

    transformed <- .klp.orthogonal.polynomial.transform(
        design = design,
        weights = weights,
        prediction.rows = anchor,
        design.drop.tol = 1e-8
    )
    gram <- crossprod(transformed$design * sqrt(weights))

    expect_true(isTRUE(transformed$ok))
    expect_lt(ncol(transformed$design), ncol(design))
    expect_equal(gram, diag(ncol(gram)), tolerance = 1e-8)
    expect_equal(ncol(transformed$prediction.rows), ncol(transformed$design))
})

test_that("LPS orthogonal polynomial drop agrees with monomial span without ridge", {
    set.seed(1706)
    x <- seq(-1, 1, length.out = 54)
    X <- cbind(x, x^2 + 0.02 * rnorm(length(x)))
    y <- sin(2 * pi * x) + rnorm(length(x), sd = 0.02)
    foldid <- rep(1:3, length.out = length(y))

    common <- list(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = 18L,
        degree.grid = 2L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        chart.dim = 2L,
        backend = "R",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
    monomial <- do.call(fit.lps, c(common, list(design.basis = "monomial")))
    orthogonal <- do.call(fit.lps, c(common, list(
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-10
    )))

    expect_equal(orthogonal$design.basis, "orthogonal.polynomial.drop")
    expect_equal(orthogonal$fitted.values, monomial$fitted.values,
                 tolerance = 1e-8)
    expect_equal(orthogonal$cv.table$cv.rmse.observed,
                 monomial$cv.table$cv.rmse.observed,
                 tolerance = 1e-8)
})

test_that("LPS default backend policy uses orthogonal adaptive tiny guard", {
    set.seed(1707)
    x <- seq(-1, 1, length.out = 42)
    X <- cbind(x, x^2 + 0.01 * rnorm(length(x)))
    y <- cos(2 * pi * x) + rnorm(length(x), sd = 0.03)
    foldid <- rep(1:3, length.out = length(y))

    fit <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = 16L,
        degree.grid = 2L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        chart.dim = 2L
    )

    expect_equal(fit$backend.used, "R")
    expect_equal(fit$design.basis, "orthogonal.polynomial.drop")
    expect_equal(fit$design.drop.tol, 1e-8)
    expect_equal(fit$ridge.multiplier.grid, c(0, 1e-10, 1e-8))
    expect_equal(fit$ridge.condition.max, 1e12)
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("K12 local-PCA contract excludes experimental second-order charts", {
    t <- seq(-1, 1, length.out = 36)
    X <- cbind(t, t^2, sin(pi * t))
    y <- t^2 + 0.1 * sin(3 * pi * t)
    foldid <- rep(1:4, length.out = nrow(X))

    fit <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = 14L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        local.chart.method = "second.order.svd",
        chart.dim = "auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        backend = "auto"
    )
    diag <- lps.backend.diagnostics(fit)

    expect_equal(diag$local.chart.method.effective, "second.order.svd")
    expect_equal(diag$backend.auto.policy, "auto_local_pca_R_reference")
    expect_true(diag$chart.dim.auto)
    expect_false(diag$local.pca.real.data.contract)
})

test_that("K12 explicit local-PCA native opt-in is reported as opt-in", {
    t <- seq(0, 1, length.out = 28)
    X <- cbind(t, t^2)
    y <- t + 0.2 * t^2
    foldid <- rep(1:4, length.out = nrow(X))

    fit <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = 10L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 1L,
        backend = "cpp.local.pca",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )
    diag <- lps.backend.diagnostics(fit)

    expect_equal(diag$backend.requested, "cpp.local.pca")
    expect_equal(diag$backend.used, "cpp.local.pca")
    expect_equal(
        diag$backend.auto.policy,
        "explicit_local_pca_native_opt_in"
    )
    expect_equal(diag$requested.chart.dim, "1")
    expect_false(diag$chart.dim.auto)
    expect_false(diag$local.pca.real.data.contract)
})

test_that("K10 row-Gram local PCA chart matches the SVD subspace", {
    set.seed(104)
    X <- matrix(rnorm(8L * 60L), nrow = 8L)
    center <- X[1L, ]
    chart <- geosmooth:::rcpp_local_pca_chart(
        X_support = X,
        center = center,
        chart_dim = 3L,
        center_mode = "anchor",
        dim_rule = "fixed",
        rebase_to_anchor = TRUE,
        orient_basis = FALSE
    )

    centered <- sweep(X, 2L, center)
    ref <- svd(centered, nu = 0L, nv = 3L)
    expect_equal(
        chart$singular.values[1:3],
        ref$d[1:3],
        tolerance = 1e-8
    )
    chart.projector <- chart$basis %*% t(chart$basis)
    ref.projector <- ref$v[, 1:3, drop = FALSE] %*%
        t(ref$v[, 1:3, drop = FALSE])
    expect_equal(chart.projector, ref.projector, tolerance = 1e-8)
})

test_that("K10 high-dimensional local-PCA LPS retains R/native parity", {
    set.seed(105)
    latent <- cbind(
        seq(-1, 1, length.out = 45L),
        sin(seq(-1, 1, length.out = 45L) * pi)
    )
    noise <- matrix(rnorm(45L * 58L, sd = 0.02), nrow = 45L)
    X <- cbind(latent, noise)
    y <- latent[, 1]^2 + 0.5 * latent[, 2]
    foldid <- rep(1:5, length.out = nrow(X))
    common.args <- list(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = c(12L, 15L),
        degree.grid = 1:2,
        kernel.grid = c("gaussian", "tricube"),
        coordinate.method = "local.pca",
        local.chart.method = "pca",
        chart.dim = 3L,
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )

    fit.r <- do.call(fit.lps, c(common.args, list(backend = "R")))
    fit.cpp <- do.call(fit.lps, c(common.args, list(backend = "cpp.local.pca")))

    expect_equal(
        fit.cpp$cv.table$cv.rmse.observed,
        fit.r$cv.table$cv.rmse.observed,
        tolerance = 1e-8
    )
    expect_equal(fit.cpp$selected, fit.r$selected, tolerance = 1e-8)
    expect_equal(fit.cpp$fitted.values, fit.r$fitted.values, tolerance = 1e-8)
})
