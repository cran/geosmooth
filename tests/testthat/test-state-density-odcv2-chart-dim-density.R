make.odcv2.path.graph <- function(n) {
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

make.odcv2.curved.X <- function(n = 38L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t^2, sin(pi * t))
}

test_that("OD-CV2 chart-kernel visit CV searches chart-dimension policies", {
    n <- 38L
    X <- make.odcv2.curved.X(n)
    subject.index <- c(3L, 5L, 8L, 11L, 14L, 17L, 21L, 24L, 28L, 33L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = make.odcv2.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = 11L,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = c("1", "auto", "local.auto"),
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator"
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), 3L)
    expect_setequal(fit$visit.cv.table$chart.dim,
                    c("1", "auto", "local.auto"))
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_true("chart.dim.rank" %in% names(fit$visit.cv.table))
    expect_true(fit$diagnostics$od.visit.cv.selection$chart.dim %in%
                c("1", "auto", "local.auto"))
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})

test_that("OD-CV2 local-likelihood density supports visit CV and chart-dim grid", {
    n <- 40L
    X <- make.odcv2.curved.X(n)
    subject.index <- c(2L, 4L, 7L, 10L, 12L, 15L, 18L, 22L, 26L, 31L, 35L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_density",
        graph = make.odcv2.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(11L, 13L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        lambda.ridge.grid = c(1e-8, 1e-6),
        coordinate.method = "local.pca",
        chart.dim.grid = c("1", "auto"),
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator"
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$method.id, "local_likelihood_density")
    expect_identical(fit$theta$od.cv, "visit")
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), 16L)
    expect_true(all(c("support.size", "degree", "kernel",
                      "bandwidth.multiplier", "lambda.ridge",
                      "chart.dim", "chart.dim.rank") %in%
                    names(fit$visit.cv.table)))
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "auto"))
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_identical(fit$diagnostics$likelihood.family, "density")
    expect_true(fit$diagnostics$od.visit.cv.selection$chart.dim %in%
                c("1", "auto"))
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})

test_that("OD-CV2 Bernoulli local likelihood can compare fixed and local-auto charts", {
    n <- 36L
    X <- make.odcv2.curved.X(n)
    subject.index <- c(4L, 6L, 9L, 12L, 16L, 19L, 23L, 27L, 30L, 34L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_bernoulli",
        graph = make.odcv2.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = 11L,
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        lambda.ridge.grid = 1e-6,
        coordinate.method = "local.pca",
        chart.dim.grid = c("1", "local.auto"),
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        optimizer = "newton"
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$method.id, "local_likelihood_bernoulli")
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), 4L)
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "local.auto"))
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_true(fit$diagnostics$od.visit.cv.selection$chart.dim %in%
                c("1", "local.auto"))
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})
