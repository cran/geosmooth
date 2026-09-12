make.odcv1.path.graph <- function(n) {
    adj <- vector("list", n)
    wt <- vector("list", n)
    for (i in seq_len(n - 1L)) {
        j <- i + 1L
        adj[[i]] <- c(adj[[i]], j)
        wt[[i]] <- c(wt[[i]], 1)
        adj[[j]] <- c(adj[[j]], i)
        wt[[j]] <- c(wt[[j]], 1)
    }
    list(
        adj.list = lapply(adj, as.integer),
        weight.list = lapply(wt, as.double)
    )
}

make.odcv1.curved.X <- function(n = 36L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t^2, sin(pi * t))
}

test_that("OD-CV1 chart-kernel supports global auto local-PCA dimension", {
    X <- make.odcv1.curved.X()
    y <- dnorm(X[, 1], mean = 0.1, sd = 0.35)
    y <- y / sum(y)

    fit <- fit.chart.kernel(
        X = X,
        y = y,
        support.size = 11L,
        kernel = "tricube",
        coordinate.method = "local.pca",
        chart.dim = "auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator"
    )

    expect_s3_class(fit, "chart_kernel")
    expect_true(all(is.finite(fit$fitted.values)))
    expect_identical(fit$selected$requested.chart.dim, "auto")
    expect_true(isTRUE(fit$selected$auto.chart.dim))
    expect_false(isTRUE(fit$selected$auto.chart.dim.local))
    expect_identical(fit$selected$chart.dim.mode, "global.auto")
    expect_true(is.finite(fit$selected$chart.dim))
    expect_true(fit$selected$chart.dim >= 1L)
    expect_true(is.list(fit$diagnostics$chart.dim$auto.diagnostics))
    expect_true(is.data.frame(fit$diagnostics$per.eval))
    expect_true(all(fit$diagnostics$chart.dim$by.anchor ==
                    fit$selected$chart.dim))
})

test_that("OD-CV1 local likelihood supports local.auto local-PCA dimension", {
    X <- make.odcv1.curved.X()
    y <- as.numeric(X[, 1] > 0)

    fit <- fit.local.likelihood(
        X = X,
        y = y,
        likelihood.family = "bernoulli",
        support.size = 13L,
        degree = 1L,
        kernel = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = "local.auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        lambda.ridge = 1e-6,
        optimizer = "newton"
    )

    expect_s3_class(fit, "local_likelihood")
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(fit$fitted.values >= 0 & fit$fitted.values <= 1))
    expect_identical(fit$selected$requested.chart.dim, "local.auto")
    expect_true(isTRUE(fit$selected$auto.chart.dim))
    expect_true(isTRUE(fit$selected$auto.chart.dim.local))
    expect_identical(fit$selected$chart.dim.mode, "local.auto")
    expect_true(is.data.frame(fit$diagnostics$per.eval))
    expect_true(all(is.finite(fit$diagnostics$chart.dim$by.anchor)))
    expect_true(all(fit$diagnostics$chart.dim$by.anchor >= 1L))
})

test_that("OD-CV1 chart-kernel visit CV preserves auto chart-dimension policy", {
    n <- 34L
    X <- make.odcv1.curved.X(n)
    subject.index <- c(2L, 5L, 7L, 9L, 11L, 16L, 20L, 23L, 27L, 31L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = make.odcv1.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        kernel.grid = c("gaussian", "tricube"),
        coordinate.method = "local.pca",
        chart.dim = "auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator"
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), 4L)
    expect_identical(fit$diagnostics$selection$requested.chart.dim, "auto")
    expect_true(isTRUE(fit$diagnostics$selection$auto.chart.dim))
    expect_identical(fit$diagnostics$selection$chart.dim.mode,
                     "global.auto")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})

test_that("OD-CV1 Bernoulli local likelihood visit CV preserves local.auto policy", {
    n <- 36L
    X <- make.odcv1.curved.X(n)
    subject.index <- c(3L, 6L, 8L, 13L, 13L, 17L, 22L, 26L, 30L, 34L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_bernoulli",
        graph = make.odcv1.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        lambda.ridge.grid = c(0, 1e-6),
        coordinate.method = "local.pca",
        chart.dim = "local.auto",
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        optimizer = "newton"
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), 8L)
    expect_identical(fit$diagnostics$selection$requested.chart.dim,
                     "local.auto")
    expect_true(isTRUE(fit$diagnostics$selection$auto.chart.dim))
    expect_true(isTRUE(fit$diagnostics$selection$auto.chart.dim.local))
    expect_identical(fit$diagnostics$selection$chart.dim.mode,
                     "local.auto")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})
