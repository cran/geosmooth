test_that("concise summaries describe actual fits and preserve their objects", {
    X <- matrix(seq(0, 1, length.out = 18), ncol = 1)
    y <- X[, 1]^2
    fits <- list(
        fit.chart.kernel(X, y, support.size = 7L),
        fit.local.likelihood(X, rep(0, 18), support.size = 7L, likelihood.family = "bernoulli"),
        fit.ps.lps(X, y, foldid = rep(1:3, 6), support.size = 7L,
                   degree = 1L, chart.dim = 1L, lambda.sync.grid = 0,
                   lambda.ridge = 0),
        normalize.density(c(-.1, .5, .6))
    )
    if (requireNamespace("genlasso", quietly = TRUE)) fits <- c(fits, list(
        fit.lpl.tf(X, y, degree = 1L, support.type = "knn", support.size = 7L,
                    lambda = .1, lambda.selection = "fixed"),
        fit.slpl.tf(X, y, degree = 1L, support.type = "knn", support.size = 7L,
                     lambda1 = .1, lambda2 = .1)))
    for (fit in fits) {
        before <- serialize(fit, NULL)
        printed <- capture.output(result <- withVisible(print(fit)))
        expect_false(result$visible)
        expect_identical(result$value, fit)
        expect_lt(length(printed), 20L)
        sm <- summary(fit)
        expect_s3_class(sm, "summary.geosmooth_fit")
        values <- if (inherits(fit, "density_fit")) fit$rho else fit$fitted.values
        expect_equal(sm$n.evaluation, NROW(values))
        expect_equal(sm$nonfinite, sum(!is.finite(values)))
        expect_identical(serialize(fit, NULL), before)
    }
    expect_equal(summary(fits[[2L]])$diagnostics$fallback.count,
                 fits[[2L]]$diagnostics$fallback.count)
    expect_equal(summary(fits[[4L]])$mass$total, 1)
    expect_identical(summary(fits[[4L]])$accounting, fits[[4L]]$accounting)
    failed <- fits[[1L]]
    failed$fitted.values[1] <- NA_real_
    failed$solver <- list(converged = FALSE)
    expect_output(print(failed), "NOT CONVERGED")
    expect_equal(summary(failed)$nonfinite, 1)
    larger <- failed
    larger$fitted.values <- rep(failed$fitted.values, 100)
    expect_equal(length(capture.output(print(failed))), length(capture.output(print(larger))))
})

test_that("standard accessors respect evaluation coordinates and response scale", {
    X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
    y <- X[, 1]^2
    fit <- fit.lps(X, y, foldid = rep(1:4, 5), support.grid = 8L,
                    degree.grid = 1L, kernel.grid = "gaussian")
    expect_equal(fitted(fit), predict(fit))
    expect_equal(residuals(fit), y - predict(fit))
    new <- fit.lps(X, y, X.eval = X[20:1, , drop = FALSE],
                    foldid = rep(1:4, 5), support.grid = 8L,
                    degree.grid = 1L, kernel.grid = "gaussian")
    expect_equal(fitted(new), predict(new))
    expect_error(residuals(new), "evaluation coordinates identical")
    expect_error(fitted(fit, newdata = X), "Unused arguments")
    expect_error(residuals(fit, type = "deviance"), "Unused arguments")
    averaged <- fit.malps(X, y, degree = 1L, support.type = "knn", support.size = 8L)
    updated <- refit(averaged, 1 - y)
    expect_equal(fitted(updated), predict(updated))
    expect_equal(residuals(updated), 1 - y - fitted(updated))
})

test_that("Hessian accessors preserve multiple responses and missing-label masks", {
    X <- matrix(seq(0, 1, length.out = 18), ncol = 1)
    Y <- cbind(linear = X[, 1], quadratic = X[, 1]^2)
    Y[c(3, 8), 1] <- NA_real_
    fit <- fit.ssrhe.hessian.regression(X, Y, k = 6L, tangent.dim = 1L,
                                        lambda1 = .01)
    expect_identical(dim(fitted(fit)), c(18L, 2L))
    expect_identical(colnames(fitted(fit)), colnames(Y))
    expect_identical(is.na(residuals(fit)), is.na(Y))
    expect_equal(residuals(fit), Y - fitted(fit))
    updated <- refit(fit, 2 * Y)
    expect_equal(residuals(updated), 2 * Y - fitted(updated))
    expect_equal(fitted(updated), updated$fitted.values)
})

test_that("tracked harmonic values have an accessor without changing the basic list", {
    adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
    weights <- list(1, c(1,1), c(1,1), 1)
    tracked <- harmonic.smoother(adj, weights, c(0,2,-1,1), region.vertices = 1:3)
    expect_identical(fitted(tracked), tracked$harmonic_predictions)
    basic <- perform.harmonic.smoothing(adj, weights, c(0,2,-1,1), region.vertices = 1:3)
    expect_identical(class(basic), "list")
})
