test_that("shared model operations dispatch by class and reject unsupported inputs", {
    expect_named(formals(refit), c("object", "y", "..."))
    expect_named(formals(smoother.matrix), c("object", "..."))
    for (class in c("malps", "lpl_tf", "slpl_tf", "metric.graph.lowpass.fit",
                    "ssrhe.hessian.fit", "ssrhe.hessian.l1.fit")) {
        expect_true(is.function(getS3method("refit", class)))
    }
    for (class in c("lps", "malps")) {
        expect_true(is.function(getS3method("smoother.matrix", class)))
    }
    expect_error(refit(list(), 1:3), "No refit.*method for class")
    expect_error(smoother.matrix(list()), "No smoother.matrix.*method for class")
})

test_that("MALPS refits and matrices agree for new responses and retain guards", {
    X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
    fit <- fit.malps(X, sin(2 * pi * X[, 1]), degree = 1L,
                     support.type = "knn", support.size = 8L)
    y <- cos(2 * pi * X[, 1])
    updated <- refit(object = fit, y = y)
    S <- smoother.matrix(fit)
    expect_s3_class(updated, "malps")
    expect_equal(as.vector(S %*% y), predict(updated), tolerance = 1e-10)
    expect_equal(predict(refit(fit)), predict(fit), tolerance = 1e-10)
    expect_equal(predict(refit(updated, fit$y)), predict(fit), tolerance = 1e-10)
    expect_equal(malps.gcv(fit), malps.gcv(fit, smoother.matrix = S))
    expect_error(refit(fit, y, reuse.selection = FALSE), "reuse.selection")
    expect_error(refit(fit, y, lambda = 0.1), "Unused arguments")
    expect_error(smoother.matrix(fit, max.n = 10L), "max.n")
    expect_error(smoother.matrix(fit, check.tol = 1e-8), "Unused arguments")

    robust <- fit.malps(X, fit$y, degree = 1L, support.type = "knn",
                        support.size = 8L, robust.iterations = 2L)
    expect_error(smoother.matrix(robust), "fixed-weight")
    frozen <- smoother.matrix(robust, allow.robust = TRUE)
    expect_equal(as.vector(frozen %*% robust$y), predict(robust), tolerance = 1e-10)
})

test_that("lifting refit methods reuse operators and allow new fixed penalties", {
    skip_if_not_installed("genlasso")
    X <- matrix(seq(0, 1, length.out = 18), ncol = 1)
    y <- sin(2 * pi * X[, 1])
    fit <- fit.lpl.tf(X, X[, 1]^2, degree = 1L, support.type = "knn",
                      support.size = 7L, lambda = 0.1, lambda.selection = "fixed")
    updated <- refit(fit, y, lambda = 0.2)
    expected <- fit.lpl.tf(X, y, operator = fit$operator,
                           lambda = 0.2, lambda.selection = "fixed")
    expect_s3_class(updated, "lpl_tf")
    expect_identical(updated$operator, fit$operator)
    expect_equal(predict(updated), predict(expected), tolerance = 1e-10)
    expect_equal(refit(fit, y)$lambda, fit$lambda)
    expect_error(refit(fit, y, reuse.lambda = FALSE), "Provide 'lambda'")
    expect_error(refit(fit, y, lambda2 = 1), "Unused arguments")

    sync <- fit.slpl.tf(X, X[, 1]^2, degree = 1L, support.type = "knn",
                        support.size = 7L, lambda1 = 0.1, lambda2 = 0.1)
    updated <- refit(sync, y, lambda1 = 0.2, lambda2 = 0.05)
    expected <- fit.slpl.tf(X, y, operator = sync$operator,
                            lambda1 = 0.2, lambda2 = 0.05)
    expect_s3_class(updated, "slpl_tf")
    expect_identical(updated$operator, sync$operator)
    expect_equal(predict(updated), predict(expected), tolerance = 1e-10)
    expect_equal(refit(sync, y)$lambda2, sync$lambda2)
    expect_error(refit(sync, y, lambda1 = 0.2, reuse.lambda = FALSE), "lambda2")
    expect_error(refit(sync, y, misspelled = 0.2), "Unused arguments")
})

test_that("Hessian refits dispatch through CV/GCV classes and preserve missing labels", {
    X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
    for (selection in c("cv", "gcv")) {
        fit <- fit.ssrhe.hessian.regression(
            X, X[, 1]^2, k = 6L, tangent.dim = 1L,
            lambda.selection = selection, lambda1.grid = c(0.01, 0.1),
            cv.control = if (selection == "cv") list(foldid = rep(1:2, 6)) else list(),
            gcv.control = if (selection == "gcv") list(trace.method = "exact") else list()
        )
        expect_equal(refit(fit)$fitted.values, fit$fitted.values, tolerance = 1e-10)
        Y <- cbind(cubic = X[, 1]^3, linear = X[, 1])
        Y[2, 1] <- NA_real_
        updated <- refit(object = fit, y = Y)
        expected <- fit.ssrhe.hessian.regression(
            X, Y, k = 6L, tangent.dim = 1L, lambda1 = fit$lambda$lambda1,
            lambda2 = fit$lambda$lambda2, ridge = fit$lambda$ridge
        )
        expect_equal(updated$fitted.values, expected$fitted.values, tolerance = 1e-10)
        expect_identical(updated$operator, fit$operator)
        expect_error(refit(updated, Y), "No refit.*method for class")
        expect_error(refit(fit, y = Y, y.new = Y), "Unused arguments")
    }
})

test_that("Hessian L1 refits preserve one-response restriction and are reusable", {
    X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
    fit <- fit.ssrhe.hessian.l1.regression(
        X, X[, 1]^2, k = 6L, tangent.dim = 1L,
        lambda.grid = 0.05, lambda.selection = "fixed", solver = "admm"
    )
    updated <- refit(fit, X[, 1]^3, solver = "admm")
    expect_s3_class(updated, "ssrhe.hessian.l1.refit")
    expect_identical(updated$operator, fit$operator)
    expect_equal(refit(updated, solver = "admm")$fitted.values,
                 updated$fitted.values, tolerance = 1e-10)
    expect_error(refit(fit, cbind(X, X), solver = "admm"), "one response vector")
    expect_error(refit(fit, X[, 1], check.tol = 1e-8), "Unused arguments")
})

test_that("metric low-pass refit retains matrix response support and rejects typos", {
    adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
    lengths <- lapply(adj, function(x) rep(1, length(x)))
    fit <- fit.metric.graph.lowpass(adj, lengths, 1:4, n.eigenpairs = 4L,
                                    eta.grid = c(0.1, 1), eigen.solver = "dense")
    Y <- cbind(a = 4:1, b = c(1, 3, 2, 4))
    updated <- refit(object = fit, y = Y, block.size = 1L)
    expect_identical(dimnames(updated$fitted.values), dimnames(Y))
    expect_equal(updated$fitted.values[, 1], refit(fit, Y[, 1])$fitted.values)
    expect_error(refit(fit, Y, per.colum.gcv = TRUE), "Unused arguments")
    expect_error(refit(updated, Y), "No refit.*method for class")
})

test_that("LPS smoother matrix maps to evaluation points and rejects unknown controls", {
    X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
    fit <- fit.lps(
        X, X[, 1]^2, X.eval = matrix(c(0.15, 0.45, 0.75), ncol = 1),
        foldid = rep(1:2, 10), support.grid = 8L, degree.grid = 1L,
        kernel.grid = "tricube", backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0, ridge.condition.max = Inf, unstable.action = "na"
    )
    S <- smoother.matrix(fit)
    expect_identical(dim(S), c(3L, 20L))
    expect_equal(as.vector(S %*% fit$y), fit$fitted.values.raw, tolerance = 1e-10)
    expect_error(smoother.matrix(fit, max.n = 100L), "Unused arguments")
    expect_error(refit(fit, fit$y), "No refit.*method for class")
})
