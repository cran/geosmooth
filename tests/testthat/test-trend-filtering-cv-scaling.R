test_that("weighted CV agrees with independently normalized mean-loss fits", {
    x <- seq(-1, 2, length.out = 24)
    w <- rep(c(1, 2, 5, 3), 6)
    z <- matrix(cos(4*x), ncol = 1)
    folds <- rep(1:3, 8)
    lambda <- c(Inf, .02)
    for (family in c("gaussian", "binomial")) {
        fitter <- if (family == "gaussian") fit.gaussian.trend.filtering else
            fit.logistic.trend.filtering
        y <- if (family == "gaussian") sin(3*x) + .2*z[, 1] else
            plogis(sin(3*x) + .2*z[, 1])
        for (backend in c("dense", "banded")) {
            fit <- fitter(x, y, lambda, weights = w, covariates = z,
                nuisance.ridge = .3, curve.ridge = .2, backend = backend,
                lambda.selection = "cv", foldid = folds)
            expect_true(all(fit$cv.table$eligible))
            for (k in 1:3) {
                train <- folds != k
                # Normalize the data loss itself instead of scaling raw penalties.
                ref <- fitter(x[train], y[train], lambda / sum(w),
                    weights = w[train] / sum(w[train]), covariates = z[train, , drop = FALSE],
                    nuisance.ridge = .3 / sum(w), curve.ridge = .2 / sum(w),
                    backend = backend)
                eta <- predict(ref, x[!train], covariates = z[!train, , drop = FALSE],
                               type = "link")
                loss <- if (family == "gaussian") (y[!train] - eta)^2 else
                    pmax(eta, 0) - y[!train]*eta + log1p(exp(-abs(eta)))
                expect_equal(unname(fit$cv.fold.loss.sums[k, ]),
                    unname(colSums(loss*w[!train])), tolerance = 2e-7)
                expect_equal(unname(fit$cv.fold.lambda[k, ]),
                    lambda * sum(w[train])/sum(w))
            }
            expect_identical(fit$cv.penalty.scale, "training_mass")
        }
    }
})

test_that("legacy CV, zero penalties, and fixed paths retain their semantics", {
    x <- seq(0, 1, length.out = 24); y <- sin(5*x)
    folds <- rep(1:3, 8); lambda <- c(Inf, .01, 0)
    scaled <- fit.gaussian.trend.filtering(x, y, lambda, lambda.selection = "cv",
        foldid = folds, backend = "banded")
    legacy <- fit.gaussian.trend.filtering(x, y, lambda, lambda.selection = "cv",
        foldid = folds, backend = "banded", cv.penalty.scale = "none")
    for (k in 1:3) {
        train <- folds != k
        ref <- fit.gaussian.trend.filtering(x[train], y[train], lambda, backend = "banded")
        expect_equal(unname(legacy$cv.fold.loss.sums[k, ]),
            unname(colSums((y[!train] - predict(ref, x[!train]))^2)), tolerance = 1e-9)
    }
    expect_equal(scaled$fits, legacy$fits)
    expect_true(all(is.infinite(scaled$cv.fold.lambda[, 1])))
    expect_equal(scaled$cv.fold.lambda[, 3], rep(0, 3))
    expect_equal(unname(legacy$cv.fold.penalty.factors), rep(1, 3))
    expect_error(fit.gaussian.trend.filtering(x, y, .1, cv.penalty.scale = "typo"),
                 "arg")
})

test_that("graph and local-polynomial CV satisfy the balanced-fold scale identity", {
    skip_if_not_installed("genlasso")
    n <- 18L; x <- seq(0, 1, length.out = n); y <- sin(5*x)
    folds <- rep(1:3, 6); lambda <- c(.03, .1, .3)
    adj <- lapply(seq_len(n), function(i) as.integer(intersect(c(i-1L, i+1L), seq_len(n))))
    g <- fit.graph.trend.filtering(adj, y = y, order = 1L,
        lambda.grid = lambda, foldid = folds)
    old <- fit.graph.trend.filtering(adj, y = y, order = 1L,
        lambda.grid = lambda*2/3, foldid = folds, cv.penalty.scale = "none")
    expect_equal(g$cv$fold.errors, old$cv$fold.errors, ignore_attr = TRUE,
                 tolerance = 1e-10)
    expect_equal(unname(g$cv$fold.penalty.factors), rep(2/3, 3))
    op <- lpl.tf.operator(matrix(x), degree = 1L, support.type = "knn", support.size = 7L)
    a <- fit.lpl.tf(y = y, operator = op, lambda.grid = lambda, foldid = folds)
    b <- fit.lpl.tf(y = y, operator = op, lambda.grid = lambda*2/3,
                    foldid = folds, cv.penalty.scale = "none")
    expect_equal(a$cv$fold.errors, b$cv$fold.errors, ignore_attr = TRUE,
                 tolerance = 1e-10)
    expect_equal(unname(a$cv$fold.lambda), outer(rep(2/3, 3), lambda))
})

test_that("LPL-TF never averages over a missing fold", {
    skip_if_not_installed("genlasso")
    n <- 12L; x <- seq(0, 1, length.out = n)
    op <- lpl.tf.operator(matrix(x), degree = 1L, support.type = "knn", support.size = 5L)
    original <- geosmooth:::.lpl.tf.genlasso.path
    calls <- 0L
    testthat::local_mocked_bindings(.lpl.tf.genlasso.path = function(...) {
        calls <<- calls + 1L
        if (calls == 2L) stop("deliberate fold failure")
        original(...)
    }, .package = "geosmooth")
    expect_error(fit.lpl.tf(y = sin(x), operator = op, lambda.grid = c(.01, .1),
        foldid = rep(1:3, 4)), "All LPL-TF CV fits failed")
    displayed <- geosmooth:::.lpl.tf.decorate.cv(
        list(fold.errors = matrix(c(1, NA, 4, 9, 16, 25), nrow = 3)), "rmse")
    expect_true(is.infinite(displayed$mean.error[1]))
    expect_true(is.infinite(displayed$se[1]))
    expect_equal(displayed$mean.error[2], 4)
})

test_that("nested refits carry the convention and scale both training levels", {
    x <- seq(0, 1, length.out = 24); y <- sin(5*x)
    object <- fit.gaussian.trend.filtering(x, y, c(Inf, .02), curve.ridge = .1,
                                         backend = "banded")
    rows <- seq_len(18); folds <- rep(1:3, 6)
    refit <- geosmooth:::.univariate.tf.assessment.refit(object, y, rows, folds,
                                                        summarize = FALSE, gaussian = TRUE)
    expect_equal(refit$fit$lambda, object$lambda*18/24)
    expect_equal(refit$fit$curve.ridge, .1*18/24)
    expect_equal(unname(refit$fit$cv.fold.lambda), outer(rep(12/24, 3), object$lambda))
    expect_identical(refit$record$cv.penalty.scale, "training_mass")
    object$cv.penalty.scale <- "none"
    legacy <- geosmooth:::.univariate.tf.assessment.refit(object, y, rows, folds,
                                                         summarize = FALSE, gaussian = TRUE)
    expect_equal(legacy$fit$lambda, object$lambda)
    expect_equal(legacy$fit$curve.ridge, .1)
    object$cv.penalty.scale <- NULL
    historical <- geosmooth:::.univariate.tf.assessment.refit(object, y, rows, folds,
                                                             summarize = FALSE, gaussian = TRUE)
    expect_equal(historical$fit$cv.table, legacy$fit$cv.table)
    expect_identical(historical$record$cv.penalty.scale, "none")
})
