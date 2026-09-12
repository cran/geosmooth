make.csd3.curved.X <- function(n = 28L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t^2, sin(2 * t), cos(2 * t))
}

make.csd3.path.graph <- function(n) {
    adj <- vector("list", n)
    wt <- vector("list", n)
    for (ii in seq_len(n - 1L)) {
        jj <- ii + 1L
        adj[[ii]] <- c(adj[[ii]], jj)
        wt[[ii]] <- c(wt[[ii]], 1)
        adj[[jj]] <- c(adj[[jj]], ii)
        wt[[jj]] <- c(wt[[jj]], 1)
    }
    list(
        adj.list = lapply(adj, as.integer),
        weight.list = lapply(wt, as.double)
    )
}

with.csd3.local.pca.call.counter <- function(expr) {
    expr <- substitute(expr)
    old.calls <- getOption("geosmooth.test.csd3.local.pca.calls")
    options(geosmooth.test.csd3.local.pca.calls = 0L)
    on.exit(options(geosmooth.test.csd3.local.pca.calls = old.calls),
            add = TRUE)
    trace(
        "rcpp_ps_lps_local_pca_supports",
        tracer = quote(options(
            geosmooth.test.csd3.local.pca.calls =
                getOption("geosmooth.test.csd3.local.pca.calls", 0L) + 1L
        )),
        where = asNamespace("geosmooth"),
        print = FALSE
    )
    on.exit(
        suppressWarnings(untrace(
            "rcpp_ps_lps_local_pca_supports",
            where = asNamespace("geosmooth")
        )),
        add = TRUE
    )
    value <- eval(expr, parent.frame())
    list(
        value = value,
        calls = getOption("geosmooth.test.csd3.local.pca.calls",
                          NA_integer_)
    )
}

expect.csd3.sparse.chart.cv <- function(fit, table.name = "cv.table") {
    tab <- fit[[table.name]]
    expect_true(is.data.frame(tab))
    expect_equal(nrow(tab), 9L)
    expect_equal(sort(unique(tab$support.size)), c(7L, 10L, 13L))
    expect_equal(sort(unique(as.character(tab$chart.dim))), c("1", "2", "4"))
    expect_true(is.list(fit$diagnostics$coupled.kd.selection))
    expect_identical(
        fit$diagnostics$coupled.kd.selection$selection.strategy,
        "sparse_kd"
    )
    expect_true(isTRUE(
        fit$diagnostics$coupled.kd.selection$coupled.chart.dim.search
    ))
}

test_that("CSD3 direct chart-kernel and local-likelihood fits use sparse k-d grids", {
    X <- make.csd3.curved.X(30L)
    t <- seq(-1, 1, length.out = nrow(X))
    y <- sin(pi * t)
    density.y <- abs(y) + 0.05
    binary.y <- as.integer(y > 0)

    chart.fit <- fit.chart.kernel(
        X = X,
        y = y,
        support.grid = 7:13,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = 1:4,
        selection.strategy = "sparse_kd",
        cv.folds = 3L,
        cv.seed = 31L
    )
    expect_s3_class(chart.fit, "chart_kernel")
    expect.csd3.sparse.chart.cv(chart.fit)
    expect_true(is.finite(chart.fit$selected$cv.rmse.observed))

    density.fit <- fit.local.likelihood(
        X = X,
        y = density.y,
        likelihood.family = "density",
        support.grid = 7:13,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        lambda.ridge.grid = 1e-8,
        coordinate.method = "local.pca",
        chart.dim.grid = 1:4,
        selection.strategy = "sparse_kd",
        cv.folds = 3L,
        cv.seed = 32L
    )
    expect_s3_class(density.fit, "local_likelihood")
    expect.csd3.sparse.chart.cv(density.fit)
    expect_true(is.finite(density.fit$selected$cv.rmse.observed))

    bernoulli.fit <- fit.local.likelihood(
        X = X,
        y = binary.y,
        likelihood.family = "bernoulli",
        support.grid = 7:13,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        lambda.ridge.grid = 1e-8,
        coordinate.method = "local.pca",
        chart.dim.grid = 1:4,
        selection.strategy = "sparse_kd",
        cv.folds = 3L,
        cv.seed = 33L
    )
    expect_s3_class(bernoulli.fit, "local_likelihood")
    expect.csd3.sparse.chart.cv(bernoulli.fit)
    expect_true(is.finite(bernoulli.fit$selected$cv.brier.observed))
})

test_that("CSD3 OD visit CV applies sparse k-d grids to chart OD methods", {
    X <- make.csd3.curved.X(28L)
    subject.index <- c(2L, 4L, 6L, 8L, 11L, 14L, 18L, 20L, 22L, 26L)
    foldid <- rep(1:5, length.out = length(subject.index))

    for (method in c("chart_kernel", "local_likelihood_density",
                     "local_likelihood_bernoulli")) {
        fit <- fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = method,
            graph = make.csd3.path.graph(nrow(X)),
            od.cv = "visit",
            visit.foldid = foldid,
            support.grid = 7:13,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            bandwidth.multiplier.grid = 1,
            lambda.ridge.grid = 1e-8,
            coordinate.method = "local.pca",
            chart.dim.grid = 1:4,
            selection.strategy = "sparse_kd"
        )

        expect_s3_class(fit, "density_fit")
        expect_identical(fit$theta$od.cv, "visit")
        expect.csd3.sparse.chart.cv(fit, table.name = "visit.cv.table")
        expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
        expect_true(is.data.frame(fit$coupled.kd.candidate.plan))
    }
})

test_that("CSD3 OD chart methods reuse max-dimension PCA coordinates by support", {
    X <- make.csd3.curved.X(28L)
    subject.index <- c(2L, 4L, 6L, 8L, 11L, 14L, 18L, 20L, 22L, 26L)
    foldid <- rep(1:5, length.out = length(subject.index))

    counted <- with.csd3.local.pca.call.counter({
        fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = "chart_kernel",
            graph = make.csd3.path.graph(nrow(X)),
            od.cv = "visit",
            visit.foldid = foldid,
            support.grid = 7:13,
            kernel.grid = "gaussian",
            bandwidth.multiplier.grid = 1,
            coordinate.method = "local.pca",
            chart.dim.grid = 1:4,
            selection.strategy = "sparse_kd"
        )
    })

    expect_equal(counted$calls, 3L)
    expect_equal(nrow(counted$value$visit.cv.table), 9L)
})
