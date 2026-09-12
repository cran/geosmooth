make.csd2.curved.X <- function(n = 28L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t^2, sin(2 * t), cos(2 * t))
}

make.csd2.path.graph <- function(n) {
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

test_that("CSD2 fit.lps evaluates candidate-specific sparse k-d grids", {
    X <- make.csd2.curved.X(30L)
    y <- sin(seq(-1, 1, length.out = nrow(X)) * pi)

    fit <- fit.lps(
        X = X,
        y = y,
        support.grid = 7:13,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = 1:4,
        selection.strategy = "sparse_kd",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = Inf,
        cv.folds = 3L,
        cv.seed = 11L
    )

    expect_s3_class(fit, "lps")
    expect_identical(fit$selection.strategy, "sparse_kd")
    expect_equal(sort(unique(fit$cv.table$support.size)), c(7L, 10L, 13L))
    expect_equal(sort(unique(fit$cv.table$chart.dim)), c(1L, 2L, 4L))
    expect_equal(nrow(fit$cv.table), 9L)
    expect_true(is.list(fit$diagnostics$coupled.kd.selection))
    expect_true(isTRUE(
        fit$diagnostics$coupled.kd.selection$coupled.chart.dim.search
    ))
    expect_equal(
        fit$diagnostics$coupled.kd.selection$evaluated.candidates,
        nrow(fit$cv.table)
    )
    expect_true(is.data.frame(fit$coupled.kd.candidate.plan))
})

test_that("CSD2 full grid remains available as an LPS reference", {
    X <- make.csd2.curved.X(24L)
    y <- X[, 1] + X[, 2]

    fit <- fit.lps(
        X = X,
        y = y,
        support.grid = c(7L, 9L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = 1:3,
        selection.strategy = "grid",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = Inf,
        cv.folds = 3L,
        cv.seed = 12L
    )

    expect_equal(nrow(fit$cv.table), 6L)
    expect_equal(sort(unique(fit$cv.table$support.size)), c(7L, 9L))
    expect_equal(sort(unique(fit$cv.table$chart.dim)), 1:3)
    expect_identical(
        fit$diagnostics$coupled.kd.selection$selection.strategy,
        "grid"
    )
})

test_that("CSD2 OD visit CV uses sparse k-d grids for LPS count and binary", {
    X <- make.csd2.curved.X(28L)
    subject.index <- c(2L, 4L, 6L, 8L, 11L, 14L, 18L, 20L, 22L, 26L)
    foldid <- rep(1:5, length.out = length(subject.index))

    for (method in c("lps_count", "lps_logistic_binary")) {
        fit <- fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = method,
            graph = make.csd2.path.graph(nrow(X)),
            od.cv = "visit",
            visit.foldid = foldid,
            support.grid = 7:13,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            bandwidth.multiplier.grid = 1,
            coordinate.method = "local.pca",
            chart.dim.grid = 1:4,
            selection.strategy = "sparse_kd",
            backend = "R",
            design.basis = "orthogonal.polynomial.drop",
            ridge.multiplier.grid = c(0, 1e-10),
            ridge.condition.max = Inf
        )

        expect_s3_class(fit, "density_fit")
        expect_identical(fit$theta$od.cv, "visit")
        expect_equal(nrow(fit$visit.cv.table), 9L)
        expect_equal(sort(unique(fit$visit.cv.table$support.size)),
                     c(7L, 10L, 13L))
        expect_equal(sort(unique(fit$visit.cv.table$chart.dim)),
                     c("1", "2", "4"))
        expect_true(is.list(fit$diagnostics$coupled.kd.selection))
        expect_identical(
            fit$diagnostics$coupled.kd.selection$selection.strategy,
            "sparse_kd"
        )
        expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    }
})
