make.odcv0.path.graph <- function(n) {
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

expect.od.visit.cv.contract <- function(fit, n.visits, n.candidates) {
    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_true(is.list(fit$diagnostics$od.visit.cv))
    expect_identical(fit$diagnostics$od.visit.cv$score.column,
                     "visit.cv.neg.log.rho")
    expect_true(is.list(fit$diagnostics$od.visit.cv.selection))
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), n.candidates)
    required <- c(
        "candidate.id",
        "visit.cv.neg.log.rho",
        "visit.cv.mean.heldout.rho",
        "visit.cv.nonfinite.count",
        "visit.cv.zero.count",
        "visit.cv.status",
        "visit.cv.error.message"
    )
    expect_true(all(required %in% names(fit$visit.cv.table)))
    expect_equal(dim(fit$visit.cv.predicted.mass),
                 c(n.visits, n.candidates))
    expect_equal(length(fit$visit.foldid), n.visits)
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_equal(
        fit$diagnostics$od.visit.cv.selection$visit.cv.neg.log.rho,
        min(fit$visit.cv.table$visit.cv.neg.log.rho),
        tolerance = 1e-12
    )
}

test_that("OD-CV0 chart-kernel visit CV exposes the frozen telemetry contract", {
    n <- 30L
    X <- cbind(seq(0, 1, length.out = n),
               sin(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(3L, 5L, 7L, 8L, 8L, 13L, 18L, 22L, 26L, 29L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = make.odcv0.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(7L, 9L),
        kernel.grid = c("gaussian", "tricube"),
        bandwidth.multiplier.grid = c(1, 1.2),
        coordinate.method = "local.pca",
        chart.dim = 1L
    )

    expect.od.visit.cv.contract(
        fit = fit,
        n.visits = length(subject.index),
        n.candidates = 8L
    )
    expect_true(all(c("support.size", "kernel", "bandwidth.multiplier") %in%
                    names(fit$visit.cv.table)))
})

test_that("OD-CV0 Bernoulli local likelihood visit CV exposes the frozen telemetry contract", {
    n <- 31L
    X <- cbind(seq(0, 1, length.out = n),
               cos(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(2L, 6L, 9L, 11L, 11L, 14L, 19L, 23L, 27L, 30L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_bernoulli",
        graph = make.odcv0.path.graph(n),
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(7L, 9L),
        degree.grid = 0:1,
        kernel.grid = c("gaussian", "tricube"),
        bandwidth.multiplier.grid = c(1, 1.1),
        lambda.ridge.grid = c(0, 1e-6),
        coordinate.method = "local.pca",
        chart.dim = 1L,
        optimizer = "newton"
    )

    expect.od.visit.cv.contract(
        fit = fit,
        n.visits = length(subject.index),
        n.candidates = 32L
    )
    expect_true(all(c("support.size", "degree", "kernel",
                      "bandwidth.multiplier", "lambda.ridge") %in%
                    names(fit$visit.cv.table)))
})

test_that("OD-CV0 compact visit CV keeps selected OD metadata without heavy tables", {
    n <- 24L
    X <- matrix(seq(0, 1, length.out = n), ncol = 1L)
    subject.index <- c(3L, 5L, 7L, 9L, 14L, 20L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        od.cv = "visit",
        visit.foldid = rep(1:3, length.out = length(subject.index)),
        support.grid = c(5L, 7L),
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        return.details = FALSE
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_true(is.list(fit$diagnostics$od.visit.cv))
    expect_true(is.list(fit$diagnostics$od.visit.cv.selection))
    expect_null(fit$visit.cv.table)
    expect_null(fit$visit.cv.predicted.mass)
    expect_null(fit$visit.foldid)
})

test_that("OD-CV0 failed candidates do not report finite mean held-out mass", {
    n <- 18L
    X <- matrix(seq(0, 1, length.out = n), ncol = 1L)
    subject.index <- c(2L, 4L, 7L, 11L, 14L, 17L)
    candidates <- data.frame(
        candidate.id = 1:2,
        support.size = c(5L, 5L),
        kernel = c("gaussian", "not_a_kernel"),
        bandwidth.multiplier = c(1, 1),
        stringsAsFactors = FALSE
    )

    cv <- .state.density.visit.cv.table(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = NULL,
        graph.control = list(),
        od.control = list(),
        foldid = rep(1:3, length.out = length(subject.index)),
        candidates = candidates,
        base.dots = list(coordinate.method = "coordinates"),
        epsilon = 1e-15
    )

    expect_identical(cv$cv.table$visit.cv.status,
                     c("ok", "failed"))
    expect_true(is.finite(cv$cv.table$visit.cv.mean.heldout.rho[[1L]]))
    expect_true(is.na(cv$cv.table$visit.cv.mean.heldout.rho[[2L]]))
    expect_identical(cv$cv.table$visit.cv.nonfinite.count[[2L]],
                     length(subject.index))
    expect_identical(cv$cv.table$visit.cv.neg.log.rho[[2L]], Inf)
})

test_that("OD-CV0 rejects ambiguous or unsupported visit-CV requests", {
    X <- matrix(seq(0, 1, length.out = 10), ncol = 1L)
    subject.index <- c(2L, 4L, 6L, 8L)

    expect_error(
        fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = "empirical",
            od.cv = "visit"
        ),
        "currently implemented"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = "chart_kernel",
            od.cv = "visit",
            visit.foldid = c(1L, 2L)
        ),
        "length\\(subject.index\\)"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = subject.index,
            method = "chart_kernel",
            od.cv = "visit",
            foldid = rep(1:2, length.out = nrow(X))
        ),
        "visit.foldid"
    )
})
