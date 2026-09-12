make.csd1.curved.rank.deficient.X <- function(n = 24L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t, t^2, sin(2 * t))
}

make.csd1.path.graph <- function(n) {
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

with.csd1.local.pca.call.counter <- function(expr) {
    expr <- substitute(expr)
    old.calls <- getOption("geosmooth.test.csd1.local.pca.calls")
    options(geosmooth.test.csd1.local.pca.calls = 0L)
    on.exit(options(geosmooth.test.csd1.local.pca.calls = old.calls),
            add = TRUE)
    trace(
        "rcpp_ps_lps_local_pca_supports",
        tracer = quote(options(
            geosmooth.test.csd1.local.pca.calls =
                getOption("geosmooth.test.csd1.local.pca.calls", 0L) + 1L
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
        calls = getOption("geosmooth.test.csd1.local.pca.calls",
                          NA_integer_)
    )
}

test_that("CSD1 reuse plans group candidates by the requested cache contract", {
    candidates <- data.frame(
        candidate.id = seq_len(8L),
        support.size = c(9L, 9L, 9L, 9L, 11L, 11L, 11L, 11L),
        degree = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
        kernel = c("gaussian", "gaussian", "tricube", "tricube",
                   "gaussian", "gaussian", "tricube", "tricube"),
        chart.dim = c("1", "3", "2", "auto", "1", "2", "4", "local.auto"),
        feasible = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
        stringsAsFactors = FALSE
    )

    weighted <- .coupled.kd.reuse.plan(candidates, reuse.type = "weighted")
    expect_equal(nrow(weighted), 3L)
    expect_equal(
        weighted$max.chart.dim[
            weighted$support.size == 9L & weighted$kernel == "gaussian"
        ],
        3L
    )
    expect_equal(
        weighted$n.candidates[
            weighted$support.size == 9L & weighted$kernel == "tricube"
        ],
        1L
    )
    expect_true(all(c("reuse.key", "reuse.chart.dim.max",
                      "candidate.ids") %in% names(weighted)))

    chart <- .coupled.kd.reuse.plan(candidates, reuse.type = "chart")
    expect_equal(nrow(chart), 2L)
    expect_equal(chart$max.chart.dim[chart$support.size == 9L], 3L)
    expect_equal(chart$max.chart.dim[chart$support.size == 11L], 2L)
    expect_false("kernel" %in% names(chart))
})

test_that("CSD1 support cache builds one local PCA object per reuse group", {
    X <- make.csd1.curved.rank.deficient.X(22L)
    candidates <- data.frame(
        candidate.id = 1:4,
        support.size = c(9L, 9L, 9L, 11L),
        kernel = "gaussian",
        chart.dim = c("1", "2", "3", "2"),
        stringsAsFactors = FALSE
    )

    counted <- with.csd1.local.pca.call.counter({
        cache <- .coupled.kd.local.pca.support.cache(
            X = X,
            candidates = candidates,
            reuse.type = "weighted"
        )
        s1 <- cache$get(support.size = 9L, chart.dim = 1L,
                        kernel = "gaussian")
        s2 <- cache$get(support.size = 9L, chart.dim = 3L,
                        kernel = "gaussian")
        s3 <- cache$get(support.size = 11L, chart.dim = 2L,
                        kernel = "gaussian")
        list(cache = cache, s1 = s1, s2 = s2, s3 = s3)
    })

    expect_equal(counted$calls, 2L)
    expect_identical(counted$value$s1, counted$value$s2)
    expect_equal(ncol(counted$value$s1[[1L]]$coordinates), 3L)
    expect_equal(ncol(counted$value$s3[[1L]]$coordinates), 2L)
})

test_that("CSD1 chart-method cache ignores kernel and reuses chart coordinates", {
    X <- make.csd1.curved.rank.deficient.X(22L)
    candidates <- data.frame(
        candidate.id = 1:4,
        support.size = 9L,
        kernel = c("gaussian", "tricube", "gaussian", "tricube"),
        chart.dim = c("1", "1", "3", "3"),
        stringsAsFactors = FALSE
    )

    counted <- with.csd1.local.pca.call.counter({
        cache <- .coupled.kd.local.pca.support.cache(
            X = X,
            candidates = candidates,
            reuse.type = "chart"
        )
        s1 <- cache$get(support.size = 9L, chart.dim = 1L,
                        kernel = "gaussian")
        s2 <- cache$get(support.size = 9L, chart.dim = 3L,
                        kernel = "tricube")
        list(s1 = s1, s2 = s2)
    })

    expect_equal(counted$calls, 1L)
    expect_identical(counted$value$s1, counted$value$s2)
    expect_equal(ncol(counted$value$s1[[1L]]$coordinates), 3L)
})

test_that("CSD1 cached and uncached LPS visit predictions match", {
    X <- make.csd1.curved.rank.deficient.X(24L)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 23L)
    foldid <- rep(1:4, length.out = length(subject.index))
    candidates <- data.frame(
        candidate.id = 1:4,
        support.size = 9L,
        degree = c(1L, 1L, 2L, 2L),
        kernel = c("gaussian", "tricube", "gaussian", "tricube"),
        chart.dim = c("1", "2", "1", "2"),
        stringsAsFactors = FALSE
    )
    cache <- .coupled.kd.local.pca.support.cache(
        X = X,
        candidates = candidates,
        reuse.type = "weighted"
    )

    for (ii in seq_len(nrow(candidates))) {
        cand <- candidates[ii, , drop = FALSE]
        dots <- list(
            support.grid = cand$support.size,
            degree.grid = cand$degree,
            kernel.grid = cand$kernel,
            bandwidth.multiplier.grid = 1,
            coordinate.method = "local.pca",
            chart.dim = as.integer(cand$chart.dim),
            backend = "R",
            design.basis = "orthogonal.polynomial.drop",
            ridge.multiplier.grid = c(0, 1e-10),
            ridge.condition.max = Inf
        )
        local.pca.supports <- cache$get(
            support.size = cand$support.size,
            chart.dim = as.integer(cand$chart.dim),
            kernel = cand$kernel
        )
        cached <- .state.density.lps.fixed.visit.predictions(
            X = X,
            subject.index = subject.index,
            foldid = foldid,
            dots = dots,
            od.control = list(),
            local.pca.supports = local.pca.supports
        )
        uncached <- .state.density.lps.fixed.visit.predictions(
            X = X,
            subject.index = subject.index,
            foldid = foldid,
            dots = dots,
            od.control = list(),
            local.pca.supports = NULL
        )
        expect_equal(cached, uncached, tolerance = 1e-8)
    }
})
