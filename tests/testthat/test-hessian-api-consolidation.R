test_that("quadratic selection rejects ambiguous inputs before building geometry", {
    fit <- function(...) fit.ssrhe.hessian.regression(...)
    expect_error(fit(lambda1 = 1, lambda1.grid = 1), "Penalty grids require")
    expect_error(fit(), "lambda1 is required")
    expect_error(fit(lambda.selection = "cv"), "lambda1.grid is required")
    expect_error(fit(lambda.selection = "gcv", lambda1.grid = 1, lambda1 = 1),
                 "not scalar")
    expect_error(fit(lambda.selection = "cv", lambda1.grid = 1, lambda2 = 0),
                 "not scalar")
    expect_error(fit(cv.control = list(cv.folds = 2)), "cv.control requires")
    expect_error(fit(lambda.selection = "cv", gcv.control = list(trace.seed = 1)),
                 "gcv.control requires")
    expect_error(fit(lambda.selection = "gcv", cv.control = list(cv.folds = 2)),
                 "cv.control requires")
    expect_error(fit(cv.control = list(fold.id = 1:2)), "Unknown cv.control.*fold.id")
    expect_error(fit(gcv.control = list(trace.sead = 1)), "Unknown gcv.control.*trace.sead")
    expect_error(fit(cv.control = list(2)), "named list")
    expect_error(fit(cv.control = list(loss = "mse", loss = "mae")), "unique")
    expect_error(fit(gcv.control = 1), "named list")
    expect_error(fit(cv.control = data.frame(cv.folds = 2)), "named list")
})

test_that("unified CV scores agree with independent fixed-operator fold solves", {
    X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
    y <- sin(2 * pi * X[, 1])
    y[12] <- NA_real_
    foldid <- c(rep(1:2, length.out = 11), 0L)
    weights <- seq(1, 2, length.out = 12)
    fit <- fit.ssrhe.hessian.regression(
        X, y, k = 6L, tangent.dim = 1L,
        lambda.selection = "cv", lambda1.grid = c(0.01, 0.1),
        cv.control = list(foldid = foldid), weights = weights, ridge = 1e-8
    )
    B <- as.matrix(fit$operator$B)
    target <- replace(y, is.na(y), 0)
    observed.weights <- replace(weights, is.na(y), 0)
    expected <- vapply(fit$cv.table$lambda1, function(lambda) {
        mean(vapply(1:2, function(fold) {
            held <- which(foldid == fold)
            train.weights <- replace(observed.weights, held, 0)
            prediction <- solve(diag(train.weights) + lambda * B + diag(1e-8, 12),
                                train.weights * target)
            sum(weights[held] * (prediction[held] - y[held])^2) / sum(weights[held])
        }, numeric(1)))
    }, numeric(1))
    expect_equal(fit$cv.table$cv.mean, expected, tolerance = 1e-8)
    expect_identical(fit$lambda.selection, "cv")
    expect_equal(fit$selection$selected.index, which.min(expected))
    expect_identical(fit$fold.id, foldid)
    expect_identical(attr(fit, "call")[[1]], quote(fit.ssrhe.hessian.regression))
})

test_that("selection enables a requested stabilizer and preserves fixed refitting", {
    X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
    y <- sin(2 * pi * X[, 1])
    for (mode in c("cv", "gcv")) {
        fit <- fit.ssrhe.hessian.regression(
            X, y, k = 6L, tangent.dim = 1L, lambda.selection = mode,
            lambda1.grid = c(0.01, 0.1), lambda2.grid = c(0, 0.02), ridge = 1e-8
        )
        expect_true(fit$operator$parameters$stabilizer)
        expect_identical(fit$lambda.selection, mode)
        expect_s3_class(fit, paste0("ssrhe.hessian.", mode, ".fit"))
        expected <- fit.ssrhe.hessian.regression(
            X, y, k = 6L, tangent.dim = 1L, lambda1 = fit$lambda$lambda1,
            lambda2 = fit$lambda$lambda2, ridge = 1e-8
        )
        expect_equal(fit$fitted.values, expected$fitted.values, tolerance = 1e-9)
        expect_equal(refit(fit)$fitted.values,
                     fit$fitted.values, tolerance = 1e-9)
        expect_error(fit.ssrhe.hessian.regression(
            X, y, k = 6L, tangent.dim = 1L, lambda.selection = mode,
            lambda1.grid = 0.1, lambda2.grid = 0.02, stabilizer = FALSE
        ), "require stabilizer")
    }
    expect_error(fit.ssrhe.hessian.regression(
        X, y, k = 6L, tangent.dim = 1L, lambda.selection = "gcv",
        lambda1.grid = 0.1, weights = c(0, rep(1, 11))
    ), "strictly positive weights")
    expect_error(fit.ssrhe.hessian.regression(
        X, cbind(y, y), k = 6L, tangent.dim = 1L, lambda.selection = "gcv",
        lambda1.grid = 0.1
    ), "one response vector")
})
