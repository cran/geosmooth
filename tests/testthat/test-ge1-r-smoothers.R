test_that("GE1 LPS coordinate R backend fits and predicts", {
    set.seed(11)
    X <- matrix(seq(0, 1, length.out = 24), ncol = 1)
    y <- sin(2 * pi * X[, 1]) + stats::rnorm(nrow(X), sd = 0.02)

    fit <- fit.lps(
        X, y,
        support.grid = c(6L, 8L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        cv.folds = 3L,
        backend = "R"
    )

    expect_s3_class(fit, "lps")
    expect_identical(fit$local.chart.method, "pca")
    expect_identical(fit$local.chart.method.effective, "none")
    expect_identical(fit$local.chart.diagnostics.summary$local.chart.method,
                     "none")
    expect_equal(length(fit$fitted.values), nrow(X))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_equal(length(predict(fit, X[1:3, , drop = FALSE])), 3L)
})

test_that("GE2 LPS coordinate C++ backend fits and predicts", {
    set.seed(12)
    X <- matrix(seq(0, 1, length.out = 24), ncol = 1)
    y <- cos(2 * pi * X[, 1]) + stats::rnorm(nrow(X), sd = 0.02)

    fit <- fit.lps(
        X, y,
        support.grid = c(6L, 8L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        cv.folds = 3L,
        backend = "cpp",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )

    expect_s3_class(fit, "lps")
    expect_identical(fit$backend.used, "cpp")
    expect_equal(length(fit$fitted.values), nrow(X))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_equal(length(predict(fit, X[1:3, , drop = FALSE])), 3L)
})

test_that("GE2 local PCA chart backend supports LPS and LPL-TF", {
    X <- cbind(seq(0, 1, length.out = 18), seq(0, 1, length.out = 18)^2)
    y <- sin(2 * pi * X[, 1])

    fit <- fit.lps(
        X, y,
        support.grid = 8L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        cv.folds = 3L,
        coordinate.method = "local.pca",
        chart.dim = 1L,
        backend = "R"
    )
    lpl <- lpl.tf.operator(
        X,
        degree = 1L,
        support.type = "knn",
        support.size = 8L,
        kernel = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 1L
    )

    expect_s3_class(fit, "lps")
    expect_identical(fit$coordinate.method, "local.pca")
    expect_identical(fit$local.chart.method, "pca")
    expect_identical(fit$local.chart.method.effective, "pca")
    expect_identical(fit$local.chart.diagnostics.summary$local.chart.method,
                     "pca")
    expect_identical(fit$local.chart.diagnostics.summary$fallback.count, 0L)
    expect_s3_class(lpl, "lpl_tf_operator")
    expect_equal(ncol(lpl$A), nrow(X))
})

test_that("H5 LPS default local PCA matches explicit PCA chart method", {
    X <- cbind(seq(0, 1, length.out = 20),
               sin(seq(0, 1, length.out = 20)))
    y <- X[, 1]^2 + 0.1 * X[, 2]
    foldid <- rep(1:5, length.out = nrow(X))

    default <- fit.lps(
        X, y,
        foldid = foldid,
        support.grid = c(8L, 10L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 1L,
        backend = "R"
    )
    explicit <- fit.lps(
        X, y,
        foldid = foldid,
        support.grid = c(8L, 10L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 1L,
        local.chart.method = "pca",
        backend = "R"
    )

    expect_identical(default$local.chart.method, "pca")
    expect_identical(default$local.chart.method.effective, "pca")
    expect_equal(default$fitted.values, explicit$fitted.values,
                 tolerance = 1e-12)
    expect_equal(default$cv.table$cv.rmse.observed,
                 explicit$cv.table$cv.rmse.observed,
                 tolerance = 1e-12)
})

test_that("H5 ambient-coordinate LPS reports no effective chart method", {
    X <- matrix(seq(0, 1, length.out = 24), ncol = 1)
    y <- X[, 1]^2
    foldid <- rep(1:4, length.out = nrow(X))

    default <- fit.lps(
        X, y,
        foldid = foldid,
        support.grid = c(6L, 8L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R"
    )
    explicit <- fit.lps(
        X, y,
        foldid = foldid,
        support.grid = c(6L, 8L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        local.chart.method = "pca",
        backend = "R"
    )

    expect_identical(default$local.chart.method, "pca")
    expect_identical(default$local.chart.method.effective, "none")
    expect_identical(
        default$local.chart.diagnostics.summary$local.chart.method,
        "none"
    )
    expect_equal(default$fitted.values, explicit$fitted.values,
                 tolerance = 1e-12)
    expect_equal(default$cv.table$cv.rmse.observed,
                 explicit$cv.table$cv.rmse.observed,
                 tolerance = 1e-12)
})

test_that("H4 LPS second-order chart method is opt-in and diagnosable", {
    grid <- expand.grid(u = seq(-0.5, 0.5, length.out = 5L),
                        v = seq(-0.5, 0.5, length.out = 5L))
    X <- cbind(grid$u, grid$v, grid$u^2 + 0.5 * grid$v^2)
    y <- sin(grid$u) + 0.5 * grid$v

    fit <- fit.lps(
        X, y,
        support.grid = 12L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        cv.folds = 5L,
        coordinate.method = "local.pca",
        chart.dim = 2L,
        local.chart.method = "second.order.svd",
        backend = "R"
    )

    expect_s3_class(fit, "lps")
    expect_identical(fit$coordinate.method, "local.pca")
    expect_identical(fit$local.chart.method, "second.order.svd")
    expect_identical(fit$local.chart.method.effective, "second.order.svd")
    expect_equal(length(fit$fitted.values), nrow(X))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_equal(nrow(fit$local.chart.diagnostics), nrow(X))
    expect_identical(
        fit$local.chart.diagnostics.summary$local.chart.method,
        "second.order.svd"
    )
    expect_true(is.numeric(fit$local.chart.diagnostics.summary$fallback.rate))
    expect_true(is.logical(
        fit$local.chart.diagnostics.summary$any.pca.fallback.used
    ))
    pred <- predict(fit, X[1:3, , drop = FALSE])
    expect_true(is.numeric(pred))
    expect_false(is.list(pred))
    expect_equal(length(pred), 3L)
})

test_that("H4 second-order LPS chart rejects ambient coordinates", {
    X <- cbind(seq(0, 1, length.out = 12))
    y <- X[, 1]

    expect_error(
        fit.lps(
            X, y,
            support.grid = 5L,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            coordinate.method = "coordinates",
            local.chart.method = "second.order.svd",
            backend = "R"
        ),
        "requires coordinate.method = 'local.pca'"
    )
})

test_that("H4 second-order LPS matches PCA closely on a flat plane", {
    grid <- expand.grid(u = seq(-1, 1, length.out = 5L),
                        v = seq(-1, 1, length.out = 5L))
    X <- cbind(grid$u, grid$v, 0)
    y <- grid$u - 2 * grid$v
    foldid <- rep(1:5, length.out = nrow(X))

    pca <- fit.lps(
        X, y,
        foldid = foldid,
        support.grid = 12L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 2L,
        local.chart.method = "pca",
        backend = "R"
    )
    second.order <- fit.lps(
        X, y,
        foldid = foldid,
        support.grid = 12L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 2L,
        local.chart.method = "second.order.svd",
        backend = "R"
    )

    expect_equal(second.order$fitted.values, pca$fitted.values,
                 tolerance = 1e-8)
    expect_equal(
        .klp.rmse(second.order$fitted.values, pca$fitted.values),
        0,
        tolerance = 1e-8
    )
})

test_that("GE1 MALPS coordinate fit works", {
    X <- cbind(seq(0, 1, length.out = 20))
    y <- X[, 1]^2

    fit <- fit.malps(
        X, y,
        degree = 1L,
        support.type = "knn",
        support.size = 7L,
        kernel = "tricube",
        support.selection = "fixed",
        coordinate.method = "coordinates"
    )

    expect_s3_class(fit, "malps")
    expect_equal(length(fit$fitted.values), nrow(X))
    expect_true(all(is.finite(fit$fitted.values)))
})

test_that("GE1 LPL-TF and SLPLiFT coordinate operators build", {
    X <- cbind(seq(0, 1, length.out = 18))

    lpl <- lpl.tf.operator(
        X,
        degree = 1L,
        support.type = "knn",
        support.size = 7L,
        kernel = "gaussian",
        coordinate.method = "coordinates"
    )
    slpl <- slpl.tf.operator(
        X,
        degree = 1L,
        support.type = "knn",
        support.size = 7L,
        kernel = "gaussian",
        coordinate.method = "coordinates"
    )

    expect_s3_class(lpl, "lpl_tf_operator")
    expect_s3_class(slpl, "slpl_tf_operator")
    expect_equal(ncol(lpl$A), nrow(X))
    expect_equal(ncol(slpl$A_LPL), nrow(X))
})

test_that("GE2 unsupported native combinations fail informatively", {
    X <- cbind(seq(0, 1, length.out = 12))
    y <- X[, 1]

    expect_error(
        fit.lps(
            X, y,
            support.grid = 5L,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            coordinate.method = "local.pca",
            chart.dim = 1L,
            backend = "cpp"
        ),
        "currently supports only"
    )
})
