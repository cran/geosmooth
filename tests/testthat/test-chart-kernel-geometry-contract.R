test_that("chart-kernel PCA feasibility uses rank rather than polynomial degree", {
    set.seed(260726)
    X <- matrix(rnorm(80), nrow = 10L, ncol = 8L)

    info <- .chart.kernel.resolve.chart.dim(
        X = X,
        support.size = 4L,
        coordinate.method = "local.pca",
        chart.dim = 8L
    )
    expect_identical(info$chart.dim, 3L)

    plan <- .coupled.kd.chart.candidate.spec(
        X = X,
        support.grid = 4L,
        kernel.grid = "tricube",
        bandwidth.multiplier.grid = 0.1,
        chart.dim.grid = 8L,
        coordinate.method = "local.pca",
        auto.chart.support.metric = "coordinates",
        auto.chart.selection.metric = "coordinates"
    )
    expect_identical(plan$candidates$chart.dim.clipped, 3L)
    expect_identical(plan$candidates$geometry.rank.cap, 3L)
    expect_identical(plan$candidates$feasibility.contract, "chart_kernel")
    expect_false(any(c("degree", "design.ncol", "design.margin") %in%
                     names(plan$candidates)))
})

test_that("finite low-bandwidth Nadaraya-Watson fits are not design-pruned", {
    X <- cbind(seq_len(12L), seq_len(12L)^2)
    y <- seq_len(12L) / sum(seq_len(12L))

    fit <- fit.chart.kernel(
        X = X,
        y = y,
        support.size = 6L,
        kernel = "tricube",
        bandwidth.multiplier = 1e-8,
        coordinate.method = "local.pca",
        chart.dim = 2L,
        denominator.floor = 1e-12
    )

    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(fit$diagnostics$per.eval$effective.support == 1L))
    expect_false(any(fit$diagnostics$per.eval$used.denominator.floor))
    expect_true(all(fit$diagnostics$per.eval$raw.denominator > 0))
})

test_that("local likelihood retains polynomial degree candidate controls", {
    set.seed(260727)
    X <- cbind(
        x = seq(-1, 1, length.out = 24L),
        z = rnorm(24L)
    )
    y <- as.integer(X[, "x"] > 0)

    fit <- fit.local.likelihood(
        X = X,
        y = y,
        likelihood.family = "bernoulli",
        support.grid = c(9L, 13L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim.grid = 1L,
        lambda.ridge.grid = 1e-6,
        cv.folds = 3L,
        cv.seed = 260727L
    )

    expect_setequal(unique(fit$cv.table$degree), 0:1)
    expect_true(all(c("design.ncol", "design.margin") %in%
                    names(fit$cv.table)))
})
