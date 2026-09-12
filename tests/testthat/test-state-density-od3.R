make.od3.path.graph <- function(n) {
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

test_that("OD3 fit.chart.kernel returns a finite fitted field", {
    X <- matrix(seq(0, 1, length.out = 21), ncol = 1L)
    y <- c(rep(0, 5), 0.2, 0.5, 1, 0.4, rep(0, 12))

    fit <- fit.chart.kernel(
        X = X,
        y = y,
        support.size = 7L,
        kernel = "gaussian",
        coordinate.method = "coordinates"
    )

    expect_s3_class(fit, "chart_kernel")
    expect_identical(fit$method.id, "chart_kernel")
    expect_length(fit$fitted.values, nrow(X))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(fit$fitted.values >= 0))
    expect_identical(fit$selected$support.size, 7L)
    expect_identical(fit$selected$kernel, "gaussian")
    expect_identical(fit$diagnostics$denominator.floor.count, 0L)
    expect_true(is.data.frame(fit$diagnostics$per.eval))
})

test_that("OD3 chart-kernel fit normalizes through density accounting", {
    X <- matrix(seq(0, 1, length.out = 25), ncol = 1L)
    subject.index <- c(3L, 5L, 8L, 8L, 14L, 21L, 24L)
    weights <- tabulate(subject.index, nbins = nrow(X))
    empirical <- weights / sum(weights)
    graph <- make.od3.path.graph(nrow(X))

    fit <- fit.chart.kernel(
        X = X,
        y = empirical,
        support.size = 9L,
        kernel = "tricube",
        coordinate.method = "coordinates"
    )
    fit$empirical.rho <- empirical
    density <- normalize.density(
        fit,
        X = X,
        adj.list = graph$adj.list,
        keep.source.fit = FALSE
    )

    expect_s3_class(density, "density_fit")
    expect_identical(density$status, "ok")
    expect_identical(density$method.id, "normalized_chart_kernel")
    expect_equal(sum(density$rho), 1, tolerance = 1e-12)
    expect_true(all(density$rho >= -1e-12))
    expect_equal(density$empirical.rho, empirical, tolerance = 1e-12)
    expect_true(is.finite(density$smoothness$n.local.maxima))
    expect_identical(density$diagnostics$source.class, "chart_kernel")
})

test_that("OD3 subject wrapper dispatches chart_kernel workflow", {
    n <- 28L
    X <- cbind(seq(0, 1, length.out = n),
               sin(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(2L, 4L, 7L, 13L, 13L, 18L, 24L, 27L)
    graph <- make.od3.path.graph(n)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = graph,
        support.size = 9L,
        kernel = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 1L
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$status, "ok")
    expect_identical(fit$method.id, "chart_kernel")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_equal(fit$empirical.rho, tabulate(subject.index, nbins = n) /
                     length(subject.index), tolerance = 1e-12)
    expect_identical(fit$diagnostics$source.method, "fit.chart.kernel")
    expect_identical(fit$diagnostics$response.summary$type,
                     "normalized_count_mass")
    expect_identical(fit$diagnostics$selection$coordinate.method, "local.pca")
    expect_true(is.finite(fit$smoothness$n.local.maxima))
    expect_equal(fit$subject$max.multiplicity, 2)
})

test_that("OD3 chart-kernel CV selects support/kernel/bandwidth candidates", {
    n <- 30L
    X <- matrix(seq(0, 1, length.out = n), ncol = 1L)
    y <- dnorm(X[, 1], mean = 0.42, sd = 0.12)
    y <- y / sum(y)
    foldid <- rep(1:3, length.out = n)

    fit <- fit.chart.kernel(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = c(7L, 11L),
        kernel.grid = c("gaussian", "tricube"),
        bandwidth.multiplier.grid = c(0.8, 1.2),
        coordinate.method = "coordinates"
    )

    expect_s3_class(fit, "chart_kernel")
    expect_identical(fit$foldid, as.integer(foldid))
    expect_true(is.data.frame(fit$cv.table))
    expect_equal(nrow(fit$cv.table), 8L)
    expect_true(all(is.finite(fit$cv.table$cv.rmse.observed)))
    expect_true(fit$selected$support.size %in% c(7L, 11L))
    expect_true(fit$selected$kernel %in% c("gaussian", "tricube"))
    expect_true(fit$selected$bandwidth.multiplier %in% c(0.8, 1.2))
    expect_equal(
        fit$selected$cv.rmse.observed,
        min(fit$cv.table$cv.rmse.observed),
        tolerance = 1e-12
    )
    expect_equal(dim(fit$cv.predictions), c(n, 8L))
})

test_that("OD3 subject chart-kernel workflow preserves CV selection diagnostics", {
    n <- 32L
    X <- cbind(seq(0, 1, length.out = n),
               cos(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(3L, 5L, 9L, 12L, 12L, 17L, 22L, 30L)
    graph <- make.od3.path.graph(n)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = graph,
        foldid = rep(1:4, length.out = n),
        support.grid = c(7L, 9L),
        kernel.grid = c("gaussian", "tricube"),
        bandwidth.multiplier.grid = c(1, 1.3),
        coordinate.method = "local.pca",
        chart.dim = 1L
    )

    expect_s3_class(fit, "density_fit")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_true(is.data.frame(fit$diagnostics$source.fit$cv.table))
    expect_equal(nrow(fit$diagnostics$source.fit$cv.table), 8L)
    expect_true(is.finite(fit$diagnostics$selection$cv.rmse.observed))
})

test_that("OD3 subject chart-kernel supports OD-level visit CV", {
    n <- 34L
    X <- cbind(seq(0, 1, length.out = n),
               sin(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(3L, 5L, 7L, 9L, 9L, 15L, 19L, 23L, 29L, 32L)
    graph <- make.od3.path.graph(n)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "chart_kernel",
        graph = graph,
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(7L, 9L),
        kernel.grid = c("gaussian", "tricube"),
        bandwidth.multiplier.grid = c(0.9, 1.2),
        coordinate.method = "local.pca",
        chart.dim = 1L
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), 8L)
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_equal(dim(fit$visit.cv.predicted.mass),
                 c(length(subject.index), 8L))
    expect_true(fit$diagnostics$selection$support.size %in% c(7L, 9L))
    expect_true(fit$diagnostics$selection$kernel %in%
                    c("gaussian", "tricube"))
    expect_identical(fit$diagnostics$od.visit.cv$score,
                     "negative_log_heldout_mass")
    expect_equal(
        fit$diagnostics$od.visit.cv.selection$visit.cv.neg.log.rho,
        min(fit$visit.cv.table$visit.cv.neg.log.rho),
        tolerance = 1e-12
    )
})

test_that("OD3 chart-kernel validates density-relevant inputs", {
    X <- cbind(seq(0, 1, length.out = 6), 0)

    expect_error(
        fit.chart.kernel(X = X, y = c(1, 2, 3)),
        "length nrow"
    )
    expect_error(
        fit.chart.kernel(X = X, y = rep(1, 6), quadrature.weights = rep(0, 6)),
        "positive finite"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(1L, 3L),
            method = "chart_kernel",
            y = rep(0, nrow(X))
        ),
        "reserved argument"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(1L, 3L, 5L),
            method = "chart_kernel",
            od.cv = "visit",
            foldid = rep(1:2, length.out = nrow(X))
        ),
        "visit.foldid"
    )
})

test_that("OD3 local-likelihood density returns finite nonnegative values", {
    X <- matrix(seq(0, 1, length.out = 31), ncol = 1L)
    y <- dnorm(X[, 1], mean = 0.45, sd = 0.14)
    y <- y / sum(y)

    fit <- fit.local.likelihood(
        X = X,
        y = y,
        likelihood.family = "density",
        support.size = 11L,
        degree = 1L,
        kernel = "gaussian",
        coordinate.method = "coordinates",
        optimizer = "optim"
    )

    expect_s3_class(fit, "local_likelihood")
    expect_identical(fit$method.id, "local_likelihood")
    expect_identical(fit$likelihood.family, "density")
    expect_length(fit$fitted.values, nrow(X))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(fit$fitted.values >= 0))
    expect_true(is.data.frame(fit$diagnostics$per.eval))
    expect_equal(fit$diagnostics$fallback.fraction, 0, tolerance = 1e-12)
})

test_that("OD3 local-likelihood density normalizes through density accounting", {
    X <- matrix(seq(0, 1, length.out = 27), ncol = 1L)
    subject.index <- c(4L, 8L, 8L, 12L, 17L, 21L, 24L)
    empirical <- tabulate(subject.index, nbins = nrow(X)) / length(subject.index)
    graph <- make.od3.path.graph(nrow(X))

    fit <- fit.local.likelihood(
        X = X,
        y = empirical,
        likelihood.family = "density",
        support.size = 9L,
        degree = 0L,
        kernel = "tricube",
        coordinate.method = "coordinates"
    )
    fit$empirical.rho <- empirical
    density <- normalize.density(
        fit,
        X = X,
        adj.list = graph$adj.list,
        keep.source.fit = FALSE
    )

    expect_s3_class(density, "density_fit")
    expect_identical(density$status, "ok")
    expect_identical(density$method.id, "normalized_local_likelihood")
    expect_equal(sum(density$rho), 1, tolerance = 1e-12)
    expect_true(all(is.finite(density$rho)))
    expect_true(all(density$rho >= -1e-12))
    expect_equal(density$empirical.rho, empirical, tolerance = 1e-12)
    expect_true(is.finite(density$smoothness$n.local.maxima))
    expect_identical(density$diagnostics$source.class, "local_likelihood")
})

test_that("OD3 local-likelihood supports local PCA charts in dimension greater than one", {
    theta <- seq(0, 2 * pi, length.out = 36)
    X <- cbind(cos(theta), sin(theta), 0.15 * cos(2 * theta))
    y <- exp(-8 * (sin(theta) - 0.2)^2)
    y <- y / sum(y)

    fit <- fit.local.likelihood(
        X = X,
        y = y,
        likelihood.family = "density",
        support.size = 9L,
        degree = 1L,
        kernel = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 2L,
        optimizer = "optim"
    )

    expect_s3_class(fit, "local_likelihood")
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(fit$fitted.values >= 0))
    expect_true(all(fit$diagnostics$chart.dim$by.anchor == 2L))
})

test_that("OD3 local-likelihood surfaces sparse-mass fallbacks", {
    X <- matrix(seq(0, 1, length.out = 30), ncol = 1L)
    y <- numeric(nrow(X))
    y[seq_len(3L)] <- c(0.3, 0.4, 0.3)

    fit <- fit.local.likelihood(
        X = X,
        y = y,
        likelihood.family = "density",
        support.size = 5L,
        degree = 1L,
        kernel = "tricube",
        coordinate.method = "coordinates",
        min.local.mass = sqrt(.Machine$double.eps)
    )

    expect_true(any(fit$diagnostics$per.eval$status == "zero_mass_fallback"))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(fit$fitted.values >= 0))
    expect_gt(fit$diagnostics$fallback.count, 0)
})

test_that("OD3 local-likelihood Bernoulli returns fitted probabilities", {
    X <- matrix(seq(0, 1, length.out = 31), ncol = 1L)
    eta <- -2 + 5 * X[, 1]
    p <- 1 / (1 + exp(-eta))
    y <- as.numeric(p > 0.5)

    fit <- fit.local.likelihood(
        X = X,
        y = y,
        likelihood.family = "bernoulli",
        support.size = 11L,
        degree = 1L,
        kernel = "gaussian",
        coordinate.method = "coordinates",
        optimizer = "optim"
    )

    expect_s3_class(fit, "local_likelihood")
    expect_identical(fit$likelihood.family, "bernoulli")
    expect_length(fit$fitted.values, nrow(X))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true(all(fit$fitted.values >= 0))
    expect_true(all(fit$fitted.values <= 1))
    expect_true(is.data.frame(fit$diagnostics$per.eval))
})

test_that("OD3 local-likelihood Bernoulli CV selects local candidates", {
    n <- 36L
    X <- matrix(seq(0, 1, length.out = n), ncol = 1L)
    p <- plogis(-3 + 7 * X[, 1])
    y <- as.numeric(p > 0.45)
    foldid <- rep(1:4, length.out = n)

    fit <- fit.local.likelihood(
        X = X,
        y = y,
        likelihood.family = "bernoulli",
        foldid = foldid,
        support.grid = c(9L, 13L),
        degree.grid = 0:1,
        kernel.grid = c("gaussian", "tricube"),
        bandwidth.multiplier.grid = c(0.9, 1.1),
        lambda.ridge.grid = c(0, 1e-6),
        coordinate.method = "coordinates",
        optimizer = "newton"
    )

    expect_s3_class(fit, "local_likelihood")
    expect_identical(fit$likelihood.family, "bernoulli")
    expect_identical(fit$foldid, as.integer(foldid))
    expect_true(is.data.frame(fit$cv.table))
    expect_equal(nrow(fit$cv.table), 32L)
    expect_true(all(is.finite(fit$cv.table$cv.brier.observed)))
    expect_true(all(is.finite(fit$cv.table$cv.logloss.observed)))
    expect_true(fit$selected$support.size %in% c(9L, 13L))
    expect_true(fit$selected$degree %in% 0:1)
    expect_true(fit$selected$kernel %in% c("gaussian", "tricube"))
    expect_true(fit$selected$lambda.ridge %in% c(0, 1e-6))
    expect_equal(
        fit$selected$cv.brier.observed,
        min(fit$cv.table$cv.brier.observed),
        tolerance = 1e-12
    )
    expect_equal(dim(fit$cv.predictions), c(n, 32L))
    expect_true(all(fit$fitted.values >= 0 & fit$fitted.values <= 1))
})

test_that("OD3 local-likelihood density and Bernoulli OD wrappers share fixture", {
    n <- 29L
    X <- cbind(seq(0, 1, length.out = n),
               cos(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(4L, 6L, 9L, 9L, 14L, 18L, 24L, 27L)
    empirical <- tabulate(subject.index, nbins = n) / length(subject.index)
    graph <- make.od3.path.graph(n)

    density.fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_density",
        graph = graph,
        support.size = 9L,
        degree = 0L,
        kernel = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 1L
    )
    bernoulli.fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_bernoulli",
        graph = graph,
        support.size = 9L,
        degree = 0L,
        kernel = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 1L
    )

    expect_s3_class(density.fit, "density_fit")
    expect_s3_class(bernoulli.fit, "density_fit")
    expect_identical(density.fit$method.id, "local_likelihood_density")
    expect_identical(bernoulli.fit$method.id, "local_likelihood_bernoulli")
    expect_equal(sum(density.fit$rho), 1, tolerance = 1e-12)
    expect_equal(sum(bernoulli.fit$rho), 1, tolerance = 1e-12)
    expect_equal(density.fit$empirical.rho, empirical, tolerance = 1e-12)
    expect_equal(bernoulli.fit$empirical.rho, empirical, tolerance = 1e-12)
    expect_identical(density.fit$diagnostics$source.method,
                     "fit.local.likelihood")
    expect_identical(bernoulli.fit$diagnostics$source.method,
                     "fit.local.likelihood")
    expect_identical(density.fit$diagnostics$likelihood.family, "density")
    expect_identical(bernoulli.fit$diagnostics$likelihood.family, "bernoulli")
    expect_true(is.finite(density.fit$smoothness$n.local.maxima))
    expect_true(is.finite(bernoulli.fit$smoothness$n.local.maxima))
})

test_that("OD3 local-likelihood density fixed visit-CV path matches fold loop", {
    n <- 30L
    X <- cbind(seq(0, 1, length.out = n),
               cos(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(3L, 5L, 5L, 8L, 11L, 14L, 18L, 22L, 25L, 28L)
    foldid <- rep(1:4, length.out = length(subject.index))
    graph <- make.od3.path.graph(n)

    for (degree in 0:1) {
        dots <- list(
            support.size = 9L,
            degree = degree,
            kernel = "gaussian",
            bandwidth.multiplier = 1.1,
            lambda.ridge = 1e-8,
            coordinate.method = "local.pca",
            chart.dim = 1L,
            optimizer = "newton"
        )
        fast <- .state.density.local.likelihood.density.fixed.visit.predictions(
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
                        method = "local_likelihood_density",
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

test_that("OD3 subject local-likelihood Bernoulli workflow preserves CV diagnostics", {
    n <- 31L
    X <- cbind(seq(0, 1, length.out = n),
               sin(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(2L, 7L, 9L, 10L, 15L, 19L, 24L, 29L)
    graph <- make.od3.path.graph(n)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_bernoulli",
        graph = graph,
        foldid = rep(1:3, length.out = n),
        support.grid = c(7L, 9L),
        degree.grid = 0:1,
        kernel.grid = c("gaussian", "tricube"),
        lambda.ridge.grid = c(0, 1e-6),
        coordinate.method = "local.pca",
        chart.dim = 1L,
        optimizer = "newton"
    )

    expect_s3_class(fit, "density_fit")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_true(is.data.frame(fit$diagnostics$source.fit$cv.table))
    expect_equal(nrow(fit$diagnostics$source.fit$cv.table), 16L)
    expect_true(is.finite(fit$diagnostics$selection$cv.brier.observed))
    expect_identical(fit$diagnostics$likelihood.family, "bernoulli")
})

test_that("OD3 subject local-likelihood Bernoulli supports OD-level visit CV", {
    n <- 33L
    X <- cbind(seq(0, 1, length.out = n),
               cos(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(4L, 5L, 8L, 12L, 12L, 16L, 20L, 25L, 28L, 31L)
    graph <- make.od3.path.graph(n)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "local_likelihood_bernoulli",
        graph = graph,
        od.cv = "visit",
        visit.foldid = rep(1:5, length.out = length(subject.index)),
        support.grid = c(7L, 9L),
        degree.grid = 0:1,
        kernel.grid = c("gaussian", "tricube"),
        bandwidth.multiplier.grid = c(1, 1.2),
        lambda.ridge.grid = c(0, 1e-6),
        coordinate.method = "local.pca",
        chart.dim = 1L,
        optimizer = "newton"
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$theta$od.cv, "visit")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_true(is.data.frame(fit$visit.cv.table))
    expect_equal(nrow(fit$visit.cv.table), 32L)
    expect_true(all(is.finite(fit$visit.cv.table$visit.cv.neg.log.rho)))
    expect_true(fit$diagnostics$selection$support.size %in% c(7L, 9L))
    expect_true(fit$diagnostics$selection$degree %in% 0:1)
    expect_true(fit$diagnostics$selection$lambda.ridge %in% c(0, 1e-6))
    expect_identical(fit$diagnostics$od.visit.cv$score.column,
                     "visit.cv.neg.log.rho")
    expect_equal(
        fit$diagnostics$od.visit.cv.selection$visit.cv.neg.log.rho,
        min(fit$visit.cv.table$visit.cv.neg.log.rho),
        tolerance = 1e-12
    )
})

test_that("OD3 local-likelihood Bernoulli fixed visit-CV path matches fold loop", {
    n <- 28L
    X <- cbind(seq(0, 1, length.out = n),
               sin(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(3L, 5L, 5L, 8L, 11L, 15L, 18L, 22L, 25L, 27L)
    foldid <- rep(1:4, length.out = length(subject.index))
    graph <- make.od3.path.graph(n)

    for (bandwidth.multiplier in c(1, 1.2)) {
        dots <- list(
            support.size = 9L,
            degree = 1L,
            kernel = "gaussian",
            bandwidth.multiplier = bandwidth.multiplier,
            lambda.ridge = 1e-6,
            coordinate.method = "local.pca",
            chart.dim = 1L,
            optimizer = "newton"
        )
        fast <- .state.density.local.likelihood.bernoulli.fixed.visit.predictions(
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
                        method = "local_likelihood_bernoulli",
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

test_that("OD3 local-likelihood validates Bernoulli responses", {
    X <- matrix(seq(0, 1, length.out = 10), ncol = 1L)

    expect_error(
        fit.local.likelihood(
            X = X,
            y = rep(2, nrow(X)),
            likelihood.family = "bernoulli"
        ),
        "values in \\[0, 1\\]"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(1L, 3L),
            method = "local_likelihood_bernoulli",
            likelihood.family = "density"
        ),
        "reserved argument"
    )
    density.cv <- fit.local.likelihood(
        X = X,
        y = rep(1, nrow(X)),
        likelihood.family = "density",
        support.grid = c(5L, 7L),
        cv.folds = 2L,
        cv.seed = 13L
    )
    expect_s3_class(density.cv, "local_likelihood")
    expect_equal(nrow(density.cv$cv.table), 2L)
    expect_true(all(is.finite(density.cv$cv.table$cv.rmse.observed)))

    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(1L, 3L, 5L),
            method = "empirical",
            od.cv = "visit"
        ),
        "currently implemented"
    )
})
