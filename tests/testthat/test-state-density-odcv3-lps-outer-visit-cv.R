make.odcv3.path.graph <- function(n) {
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

make.odcv3.curved.X <- function(n = 34L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t^2)
}

expect.odcv3.visit.fit <- function(fit, n.visits, n.candidates) {
    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), n.candidates)
    expect_equal(dim(fit$visit.cv.predicted.mass),
                 c(n.visits, n.candidates))
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_equal(
        fit$diagnostics$od.visit.cv.selection$visit.cv.neg.log.rho,
        min(fit$visit.cv.table$visit.cv.neg.log.rho),
        tolerance = 1e-12
    )
}

with.odcv3.local.pca.call.counter <- function(expr) {
    expr <- substitute(expr)
    old.calls <- getOption("geosmooth.test.ps.lps.local.pca.calls")
    options(geosmooth.test.ps.lps.local.pca.calls = 0L)
    on.exit(options(geosmooth.test.ps.lps.local.pca.calls = old.calls),
            add = TRUE)
    trace(
        "rcpp_ps_lps_local_pca_supports",
        tracer = quote(options(
            geosmooth.test.ps.lps.local.pca.calls =
                getOption("geosmooth.test.ps.lps.local.pca.calls", 0L) + 1L
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
        calls = getOption("geosmooth.test.ps.lps.local.pca.calls",
                          NA_integer_)
    )
}

test_that("OD-CV3 LPS count visit CV searches outer scalar candidates", {
    n <- 34L
    X <- make.odcv3.curved.X(n)
    subject.index <- c(3L, 5L, 8L, 8L, 11L, 15L, 19L, 24L, 29L, 32L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "lps_count",
        graph = make.odcv3.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = c(1, 1.2),
        coordinate.method = "local.pca",
        chart.dim.grid = c("1", "auto"),
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )

    expect.odcv3.visit.fit(fit, length(subject.index), 16L)
    expect_identical(fit$method.id, "lps_count")
    expect_true(all(c("support.size", "degree", "kernel",
                      "bandwidth.multiplier", "chart.dim",
                      "chart.dim.rank") %in% names(fit$visit.cv.table)))
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "auto"))
    expect_true(fit$diagnostics$od.visit.cv.selection$chart.dim %in%
                c("1", "auto"))
    expect_true(fit$diagnostics$selection$support.size %in% c(9L, 11L))
    expect_true(fit$diagnostics$selection$degree %in% 0:1)
    expect_true(fit$diagnostics$selection$bandwidth.multiplier %in% c(1, 1.2))
})

test_that("OD-CV3 LPS Bernoulli visit CV preserves binary OD workflow", {
    n <- 32L
    X <- make.odcv3.curved.X(n)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 16L, 20L, 24L, 28L, 31L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "lps_logistic_binary",
        graph = make.odcv3.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        degree.grid = 0:1,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = c("1", "local.auto"),
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )

    expect.odcv3.visit.fit(fit, length(subject.index), 8L)
    expect_identical(fit$method.id, "lps_logistic_binary")
    expect_identical(fit$diagnostics$outcome.family, "bernoulli")
    expect_identical(
        fit$diagnostics$binary.workflow$probability.link,
        "identity_lps_least_squares_clipped"
    )
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "local.auto"))
})

test_that("OD-CV3 PS-LPS count visit CV searches lambda.sync candidates", {
    n <- 26L
    X <- make.odcv3.curved.X(n)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 25L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "ps_lps_count",
        graph = make.odcv3.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:4, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        chart.dim.grid = c("1", "auto"),
        auto.chart.support.metric = "both",
        auto.chart.selection.metric = "operator",
        lambda.sync.grid = c(0, 0.1),
        lambda.ridge = 1e-8,
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = 1e10,
        sync.neighbor.size = 3L
    )

    expect.odcv3.visit.fit(fit, length(subject.index), 8L)
    expect_identical(fit$method.id, "ps_lps_count")
    expect_true(all(c("support.size", "degree", "kernel", "lambda.sync",
                      "chart.dim", "chart.dim.rank") %in%
                    names(fit$visit.cv.table)))
    expect_setequal(fit$visit.cv.table$lambda.sync, c(0, 0.1))
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "auto"))
    expect_true(fit$diagnostics$od.visit.cv.selection$lambda.sync %in%
                c(0, 0.1))
    expect_true(fit$diagnostics$selection$lambda.sync %in% c(0, 0.1))
    expect_identical(fit$diagnostics$source.method, "fit.ps.lps")
})

test_that("OD-CV3 PS-LPS numeric chart-dimension grids reuse max PCA supports", {
    n <- 28L
    t <- seq(-1, 1, length.out = n)
    X <- cbind(t, t^2, sin(2 * t), cos(3 * t))
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 25L)

    counted <- with.odcv3.local.pca.call.counter(fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "ps_lps_count",
        graph = make.odcv3.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:4, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        chart.dim.grid = c(1L, 2L, 4L),
        lambda.sync.grid = 0,
        lambda.ridge = 1e-8,
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = 1e10,
        sync.neighbor.size = 3L,
        chart.activation = "none"
    ))

    fit <- counted$value
    expect.odcv3.visit.fit(fit, length(subject.index), 6L)
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "2", "4"))
    expect_equal(counted$calls, 2L)
})

test_that("OD-CV3 LPS numeric chart-dimension grids reuse max PCA supports", {
    n <- 28L
    t <- seq(-1, 1, length.out = n)
    X <- cbind(t, t^2, sin(2 * t), cos(3 * t))
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 25L)

    counted <- with.odcv3.local.pca.call.counter(fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "lps_count",
        graph = make.odcv3.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:4, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        degree.grid = 1L,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = c(1L, 2L, 4L),
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = 1e10
    ))

    fit <- counted$value
    expect.odcv3.visit.fit(fit, length(subject.index), 6L)
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "2", "4"))
    expect_equal(counted$calls, 2L)
})

test_that("OD-CV3 chart-kernel numeric chart-dimension grids reuse PCA coordinates", {
    n <- 28L
    t <- seq(-1, 1, length.out = n)
    X <- cbind(t, t^2, sin(2 * t), cos(3 * t))
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 25L)

    counted <- with.odcv3.local.pca.call.counter(fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = make.odcv3.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:4, length.out = length(subject.index)),
        support.grid = c(9L, 11L),
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = c(1L, 2L, 4L)
    ))

    fit <- counted$value
    expect.odcv3.visit.fit(fit, length(subject.index), 6L)
    expect_setequal(fit$visit.cv.table$chart.dim, c("1", "2", "4"))
    expect_equal(counted$calls, 2L)
})

test_that("OD-CV3 local-likelihood numeric chart-dimension grids reuse PCA coordinates", {
    n <- 28L
    t <- seq(-1, 1, length.out = n)
    X <- cbind(t, t^2, sin(2 * t), cos(3 * t))
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 25L)

    for (method in c("local_likelihood_density",
                     "local_likelihood_bernoulli")) {
        counted <- with.odcv3.local.pca.call.counter(fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = method,
            graph = make.odcv3.path.graph(n),
            od.cv = "visit",
            visit.foldid = rep(1:4, length.out = length(subject.index)),
            support.grid = c(9L, 11L),
            degree.grid = 1L,
            kernel.grid = "gaussian",
            bandwidth.multiplier.grid = 1,
            coordinate.method = "local.pca",
            chart.dim.grid = c(1L, 2L, 4L),
            lambda.ridge.grid = 1e-8,
            fallback = "degree0"
        ))

        fit <- counted$value
        expect.odcv3.visit.fit(fit, length(subject.index), 6L)
        expect_setequal(fit$visit.cv.table$chart.dim, c("1", "2", "4"))
        expect_equal(counted$calls, 2L)
    }
})

test_that("OD-CV3 PS-LPS fixed-candidate fast path matches fold loop", {
    n <- 24L
    X <- make.odcv3.curved.X(n)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 23L)
    foldid <- rep(1:4, length.out = length(subject.index))
    graph <- make.odcv3.path.graph(n)

    for (lambda.sync in c(0, 0.1)) {
        dots <- list(
            support.size = 9L,
            degree = 1L,
            kernel = "gaussian",
            chart.dim = 1L,
            auto.chart.support.metric = "coordinates",
            auto.chart.selection.metric = "coordinates",
            lambda.sync.grid = lambda.sync,
            lambda.sync.search = "grid",
            lambda.sync.selection = "fixed",
            lambda.ridge = 1e-8,
            design.basis = "orthogonal.polynomial.drop",
            ridge.multiplier.grid = c(0, 1e-10),
            ridge.condition.max = 1e10,
            sync.neighbor.size = 3L
        )
        cache <- .state.density.ps.lps.geometry.cache(X, dots)
        fast <- .state.density.ps.lps.fixed.visit.predictions(
            X = X,
            subject.index = subject.index,
            foldid = foldid,
            geometry.cache = cache,
            dots = dots,
            od.control = list(),
            graph = graph
        )
        slow <- rep(NA_real_, length(subject.index))
        for (fold in sort(unique(foldid))) {
            test.pos <- which(foldid == fold)
            train.pos <- which(foldid != fold)
            fit <- do.call(
                fit.subject.od,
                c(
                    list(
                        X = X,
                        subject.index = subject.index[train.pos],
                        method = "ps_lps_count",
                        graph = graph,
                        od.cv = "none",
                        return.details = FALSE
                    ),
                    dots
                )
            )
            slow[test.pos] <- fit$rho[subject.index[test.pos]]
        }
        expect_equal(fast, slow, tolerance = 1e-8)
    }
})

test_that("OD-CV3 LPS count fixed-candidate fast path matches fold loop", {
    n <- 24L
    X <- make.odcv3.curved.X(n)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 21L, 23L)
    foldid <- rep(1:4, length.out = length(subject.index))
    graph <- make.odcv3.path.graph(n)

    for (bandwidth.multiplier in c(1, 1.2)) {
        dots <- list(
            support.grid = 9L,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            bandwidth.multiplier.grid = bandwidth.multiplier,
            coordinate.method = "local.pca",
            chart.dim = 1L,
            auto.chart.support.metric = "coordinates",
            auto.chart.selection.metric = "coordinates",
            backend = "R",
            design.basis = "orthogonal.polynomial.drop",
            ridge.multiplier.grid = c(0, 1e-10),
            ridge.condition.max = Inf
        )
        fast <- .state.density.lps.fixed.visit.predictions(
            X = X,
            subject.index = subject.index,
            foldid = foldid,
            dots = dots,
            od.control = list()
        )
        slow <- rep(NA_real_, length(subject.index))
        for (fold in sort(unique(foldid))) {
            test.pos <- which(foldid == fold)
            train.pos <- which(foldid != fold)
            fit <- do.call(
                fit.subject.od,
                c(
                    list(
                        X = X,
                        subject.index = subject.index[train.pos],
                        method = "lps_count",
                        graph = graph,
                        od.cv = "none",
                        return.details = FALSE
                    ),
                    dots
                )
            )
            slow[test.pos] <- fit$rho[subject.index[test.pos]]
        }
        expect_equal(fast, slow, tolerance = 1e-8)
    }
})

test_that("OD-CV3 LPS Bernoulli fixed-candidate fast path matches fold loop", {
    n <- 24L
    X <- make.odcv3.curved.X(n)
    subject.index <- c(2L, 4L, 4L, 7L, 10L, 13L, 17L, 21L, 23L)
    foldid <- rep(1:3, length.out = length(subject.index))
    graph <- make.odcv3.path.graph(n)

    for (bandwidth.multiplier in c(1, 1.2)) {
        dots <- list(
            support.grid = 9L,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            bandwidth.multiplier.grid = bandwidth.multiplier,
            coordinate.method = "local.pca",
            chart.dim = 1L,
            auto.chart.support.metric = "coordinates",
            auto.chart.selection.metric = "coordinates",
            backend = "R",
            design.basis = "orthogonal.polynomial.drop",
            ridge.multiplier.grid = c(0, 1e-10),
            ridge.condition.max = Inf
        )
        fast <- .state.density.lps.fixed.visit.predictions(
            X = X,
            subject.index = subject.index,
            foldid = foldid,
            dots = dots,
            od.control = list(),
            outcome.family = "bernoulli"
        )
        slow <- rep(NA_real_, length(subject.index))
        for (fold in sort(unique(foldid))) {
            test.pos <- which(foldid == fold)
            train.pos <- which(foldid != fold)
            fit <- do.call(
                fit.subject.od,
                c(
                    list(
                        X = X,
                        subject.index = subject.index[train.pos],
                        method = "lps_logistic_binary",
                        graph = graph,
                        od.cv = "none",
                        return.details = FALSE
                    ),
                    dots
                )
            )
            slow[test.pos] <- fit$rho[subject.index[test.pos]]
        }
        expect_equal(fast, slow, tolerance = 1e-8)
    }
})

test_that("OD-CV3 chart-kernel fixed-candidate fast path matches fold loop", {
    n <- 24L
    t <- seq(-1, 1, length.out = n)
    X <- cbind(t, t^2, sin(t))
    subject.index <- c(2L, 4L, 4L, 7L, 10L, 13L, 17L, 21L, 23L)
    foldid <- rep(1:3, length.out = length(subject.index))
    graph <- make.odcv3.path.graph(n)
    dots <- list(
        support.size = 9L,
        kernel = "gaussian",
        bandwidth.multiplier = 1,
        coordinate.method = "local.pca",
        chart.dim = 2L,
        auto.chart.support.metric = "coordinates",
        auto.chart.selection.metric = "coordinates"
    )

    fast <- .state.density.chart.kernel.fixed.visit.predictions(
        X = X,
        subject.index = subject.index,
        foldid = foldid,
        dots = dots,
        od.control = list()
    )
    slow <- rep(NA_real_, length(subject.index))
    for (fold in sort(unique(foldid))) {
        test.pos <- which(foldid == fold)
        train.pos <- which(foldid != fold)
        fit <- do.call(
            fit.subject.od,
            c(
                list(
                    X = X,
                    subject.index = subject.index[train.pos],
                    method = "chart_kernel",
                    graph = graph,
                    od.cv = "none",
                    return.details = FALSE
                ),
                dots
            )
        )
        slow[test.pos] <- fit$rho[subject.index[test.pos]]
    }
    expect_equal(fast, slow, tolerance = 1e-8)
})

test_that("OD-CV3 local-likelihood fixed-candidate fast paths match fold loop", {
    n <- 24L
    t <- seq(-1, 1, length.out = n)
    X <- cbind(t, t^2, sin(t))
    subject.index <- c(2L, 4L, 4L, 7L, 10L, 13L, 17L, 21L, 23L)
    foldid <- rep(1:3, length.out = length(subject.index))
    graph <- make.odcv3.path.graph(n)
    dots <- list(
        support.size = 9L,
        degree = 1L,
        kernel = "gaussian",
        bandwidth.multiplier = 1,
        coordinate.method = "local.pca",
        chart.dim = 2L,
        auto.chart.support.metric = "coordinates",
        auto.chart.selection.metric = "coordinates",
        lambda.ridge = 1e-8,
        fallback = "degree0"
    )

    for (method in c("local_likelihood_density",
                     "local_likelihood_bernoulli")) {
        fast <- if (identical(method, "local_likelihood_density")) {
            .state.density.local.likelihood.density.fixed.visit.predictions(
                X = X,
                subject.index = subject.index,
                foldid = foldid,
                dots = dots,
                od.control = list()
            )
        } else {
            .state.density.local.likelihood.bernoulli.fixed.visit.predictions(
                X = X,
                subject.index = subject.index,
                foldid = foldid,
                dots = dots,
                od.control = list()
            )
        }
        slow <- rep(NA_real_, length(subject.index))
        for (fold in sort(unique(foldid))) {
            test.pos <- which(foldid == fold)
            train.pos <- which(foldid != fold)
            fit <- do.call(
                fit.subject.od,
                c(
                    list(
                        X = X,
                        subject.index = subject.index[train.pos],
                        method = method,
                        graph = graph,
                        od.cv = "none",
                        return.details = FALSE
                    ),
                    dots
                )
            )
            slow[test.pos] <- fit$rho[subject.index[test.pos]]
        }
        expect_equal(fast, slow, tolerance = 1e-8)
    }
})

test_that("OD-CV3 PS-LPS visit CV accepts explicit fixed lambda selection", {
    n <- 22L
    X <- make.odcv3.curved.X(n)
    subject.index <- c(2L, 4L, 7L, 10L, 13L, 17L, 20L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "ps_lps_count",
        od.cv = "visit",
        visit.foldid = rep(1:3, length.out = length(subject.index)),
        support.grid = 9L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        chart.dim = 1L,
        lambda.sync.grid = 0.1,
        lambda.sync.selection = "fixed",
        lambda.ridge = 1e-8,
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 1e-10,
        ridge.condition.max = Inf,
        sync.neighbor.size = 3L
    )

    expect.odcv3.visit.fit(fit, length(subject.index), 1L)
    expect_equal(fit$visit.cv.table$lambda.sync, 0.1)
    expect_identical(fit$diagnostics$selection$lambda.sync, 0.1)
})

test_that("OD-CV3 PS-LPS visit CV requires an explicit chart dimension policy", {
    X <- make.odcv3.curved.X(18L)
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(2L, 4L, 7L, 11L),
            method = "ps_lps_count",
            od.cv = "visit",
            visit.cv.folds = 2L
        ),
        "requires 'chart.dim'"
    )
})
