make.csd4.curved.X <- function(n = 28L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t^2, sin(2 * t), cos(3 * t))
}

make.csd4.path.graph <- function(n) {
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

with.csd4.local.pca.call.counter <- function(expr) {
    expr <- substitute(expr)
    old.calls <- getOption("geosmooth.test.csd4.local.pca.calls")
    options(geosmooth.test.csd4.local.pca.calls = 0L)
    on.exit(options(geosmooth.test.csd4.local.pca.calls = old.calls),
            add = TRUE)
    trace(
        "rcpp_ps_lps_local_pca_supports",
        tracer = quote(options(
            geosmooth.test.csd4.local.pca.calls =
                getOption("geosmooth.test.csd4.local.pca.calls", 0L) + 1L
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
        calls = getOption("geosmooth.test.csd4.local.pca.calls",
                          NA_integer_)
    )
}

csd4.ps.lps.args <- function() {
    list(
        support.grid = 7:13,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        chart.dim.grid = 1:4,
        selection.strategy = "sparse_kd",
        local.candidate.search = "full",
        lambda.sync.grid = c(0, 0.1),
        lambda.sync.search = "grid",
        lambda.ridge = 1e-8,
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = 1e10,
        sync.neighbor.size = 3L
    )
}

test_that("CSD4 fit.ps.lps evaluates sparse coupled k-d candidates", {
    X <- make.csd4.curved.X(28L)
    y <- sin(pi * seq(-1, 1, length.out = nrow(X)))
    args <- csd4.ps.lps.args()

    fit <- do.call(
        fit.ps.lps,
        c(list(X = X, y = y, foldid = rep(1:4, length.out = nrow(X))),
          args)
    )

    expect_s3_class(fit, "ps_lps")
    expect_identical(fit$selection.strategy, "sparse_kd")
    expect_identical(
        fit$selection.contract,
        "sparse_kd_coupled_support_chart_dim_with_lambda_sync"
    )
    expect_equal(nrow(fit$local.candidate.table), 9L)
    expect_equal(sort(unique(fit$local.candidate.table$support.size)),
                 c(7L, 10L, 13L))
    expect_equal(sort(unique(as.integer(fit$local.candidate.table$chart.dim))),
                 c(1L, 2L, 4L))
    expect_equal(nrow(fit$lambda.cv.table), 18L)
    expect_true(isTRUE(
        fit$diagnostics$coupled.kd.selection$coupled.chart.dim.search
    ))
    expect_equal(
        fit$diagnostics$coupled.kd.selection$evaluated.candidates,
        nrow(fit$local.candidate.table)
    )
    expect_true(is.data.frame(fit$coupled.kd.candidate.plan))

    default.args <- args
    default.args$local.candidate.search <- NULL
    default.fit <- do.call(
        fit.ps.lps,
        c(list(X = X, y = y, foldid = rep(1:4, length.out = nrow(X))),
          default.args)
    )
    expect_identical(default.fit$local.candidate.search, "full")
    expect_equal(nrow(default.fit$local.candidate.table), 9L)
})

test_that("CSD4 direct PS-LPS reuses max-dimension PCA supports", {
    X <- make.csd4.curved.X(28L)
    y <- cos(pi * seq(-1, 1, length.out = nrow(X)))
    args <- csd4.ps.lps.args()

    counted <- with.csd4.local.pca.call.counter({
        do.call(
            fit.ps.lps,
            c(list(X = X, y = y, foldid = rep(1:4, length.out = nrow(X))),
              args)
        )
    })

    expect_equal(counted$calls, 3L)
    expect_equal(nrow(counted$value$local.candidate.table), 9L)
})

test_that("CSD4 OD ps_lps_count uses sparse coupled k-d visit CV", {
    X <- make.csd4.curved.X(28L)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 25L)
    args <- csd4.ps.lps.args()

    fit <- do.call(
        fit.subject.od,
        c(list(
            X = X,
            subject.index = subject.index,
            method = "ps_lps_count",
            graph = make.csd4.path.graph(nrow(X)),
            od.cv = "visit",
            visit.foldid = rep(1:4, length.out = length(subject.index))
        ), args)
    )

    expect_s3_class(fit, "density_fit")
    expect_equal(nrow(fit$visit.cv.table), 18L)
    expect_equal(sort(unique(fit$visit.cv.table$support.size)),
                 c(7L, 10L, 13L))
    expect_equal(sort(unique(fit$visit.cv.table$chart.dim)),
                 c("1", "2", "4"))
    expect_equal(sort(unique(fit$visit.cv.table$lambda.sync)), c(0, 0.1))
    expect_true(isTRUE(
        fit$diagnostics$coupled.kd.selection$coupled.chart.dim.search
    ))
    expect_identical(
        fit$diagnostics$coupled.kd.selection$selection.strategy,
        "sparse_kd"
    )
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})

test_that("CSD4 OD PS-LPS sparse visit CV reuses PCA supports by support/kernel", {
    X <- make.csd4.curved.X(28L)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 25L)
    args <- csd4.ps.lps.args()

    counted <- with.csd4.local.pca.call.counter({
        do.call(
            fit.subject.od,
            c(list(
                X = X,
                subject.index = subject.index,
                method = "ps_lps_count",
                graph = make.csd4.path.graph(nrow(X)),
                od.cv = "visit",
                visit.foldid = rep(1:4, length.out = length(subject.index))
            ), args)
        )
    })

    expect_equal(counted$calls, 3L)
    expect_equal(nrow(counted$value$visit.cv.table), 18L)
})
