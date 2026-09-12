test_that("OD0 empirical density normalizes repeated subject visits", {
    X <- matrix(seq(0, 1, length.out = 5), ncol = 1L)
    fit <- fit.subject.od(
        X = X,
        subject.index = c(1L, 2L, 2L, 5L),
        method = "empirical"
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$method.id, "empirical")
    expect_identical(fit$status, "ok")
    expect_equal(fit$rho, c(1, 2, 0, 0, 1) / 4, tolerance = 1e-12)
    expect_equal(fit$empirical.rho, fit$rho, tolerance = 1e-12)
    expect_equal(fit$accounting$mass, 1, tolerance = 1e-12)
    expect_equal(fit$accounting$min.rho, 0, tolerance = 1e-12)
    expect_equal(fit$subject$n.visits, 4L)
    expect_equal(fit$subject$n.unique.visited, 3L)
    expect_equal(fit$subject$max.multiplicity, 2)
    expect_equal(fit$subject$repeat.fraction, 1 / 3, tolerance = 1e-12)
})

test_that("OD0 empirical density validates inputs and control values", {
    X <- matrix(1:6, ncol = 2L)

    expect_error(
        fit.density(
            X, weights = c(1, -1, 2), method = "empirical"
        ),
        "nonnegative"
    )
    expect_error(
        fit.density(
            X, weights = c(0, 0, 0), method = "empirical"
        ),
        "positive total mass"
    )
    expect_error(
        fit.subject.od(X, subject.index = c(1L, 4L), method = "empirical"),
        "outside 1:nrow"
    )
    expect_error(
        fit.density(
            X,
            weights = c(1, 1, 1),
            method = "empirical",
            density.control = list(mass.tol = -1)
        ),
        "mass.tol"
    )
})

test_that("OD0 correction helper clips and renormalizes raw density fields", {
    X <- matrix(1:4, ncol = 1L)
    out <- geosmooth:::.state.density.finalize(
        method.id = "test",
        X = X,
        fitted.raw = c(0.5, -0.1, 0.6, 0),
        density.control = list(clip.negative = TRUE, renormalize = TRUE)
    )

    expect_identical(out$status, "ok")
    expect_equal(out$rho, c(0.5, 0, 0.6, 0) / 1.1, tolerance = 1e-12)
    expect_equal(out$accounting$raw.mass, 1.0, tolerance = 1e-12)
    expect_equal(out$accounting$neg.mass, 0.1, tolerance = 1e-12)
    expect_equal(out$accounting$clip.mass, 0.1, tolerance = 1e-12)
    expect_equal(out$accounting$normalization.constant, 1.1, tolerance = 1e-12)
})

test_that("OD0 density dispatcher is restricted to density-native methods", {
    X <- matrix(seq(0, 1, length.out = 6), ncol = 1L)
    weights <- c(1, 0, 2, 0, 0, 1)

    fit <- fit.density(X = X, weights = weights, method = "empirical")
    expect_s3_class(fit, "density_fit")
    expect_equal(fit$rho, weights / sum(weights), tolerance = 1e-12)

    expect_error(
        fit.density(X = X, weights = weights, method = "lps_count"),
        "'arg' should be one of"
    )
})

test_that("OD0 density-native methods reject chart-dimension arguments", {
    X <- matrix(seq(0, 1, length.out = 6), ncol = 1L)
    weights <- c(1, 0, 2, 0, 0, 1)
    graph <- list(
        adj.list = list(2L, c(1L, 3L), c(2L, 4L),
                        c(3L, 5L), c(4L, 6L), 5L),
        weight.list = rep(list(1), 6L)
    )

    expect_error(
        fit.density(
            X = X,
            weights = weights,
            method = "empirical",
            chart.dim = "local.auto"
        ),
        "does not use local charts"
    )
    expect_error(
        fit.density(
            X = X,
            weights = weights,
            method = "graph_random_walk",
            graph = graph,
            chart.dim.grid = c(1L, "auto")
        ),
        "does not use local charts"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(1L, 3L, 6L),
            method = "empirical",
            chart.dim = "auto"
        ),
        "does not use local charts"
    )
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(1L, 3L, 6L),
            method = "graph_random_walk",
            graph = graph,
            chart.dim = "local.auto"
        ),
        "does not use local charts"
    )
})

test_that("OD0 normalize.density clips and normalizes numeric and fit objects", {
    X <- matrix(1:4, ncol = 1L)

    numeric.fit <- normalize.density(
        c(0.5, -0.1, 0.6, 0),
        X = X,
        density.control = list(clip.negative = TRUE, renormalize = TRUE)
    )
    expect_s3_class(numeric.fit, "density_fit")
    expect_equal(numeric.fit$rho, c(0.5, 0, 0.6, 0) / 1.1,
                 tolerance = 1e-12)
    expect_equal(numeric.fit$accounting$clip.mass, 0.1, tolerance = 1e-12)
    expect_equal(numeric.fit$empirical.rho, rep(NA_real_, 4L))

    lps.like <- list(
        method.id = "lps",
        X.eval = X,
        fitted.values = c(1, 2, 0, 1)
    )
    class(lps.like) <- c("lps", "list")
    lps.density <- normalize.density(lps.like, keep.source.fit = FALSE)
    expect_s3_class(lps.density, "density_fit")
    expect_identical(lps.density$method.id, "normalized_lps")
    expect_equal(lps.density$rho, c(1, 2, 0, 1) / 4, tolerance = 1e-12)
    expect_equal(lps.density$empirical.rho, rep(NA_real_, 4L))
    expect_identical(lps.density$diagnostics$source.class, "lps")
})

test_that("OD0 normalize.density exposes chart dimensions in a uniform diagnostic", {
    X <- matrix(seq(0, 1, length.out = 4), ncol = 1L)
    make.source <- function(cls, dims) {
        source <- list(
            X.eval = X,
            fitted.values = c(1, 2, 3, 4),
            selected = list(requested.chart.dim = "local.auto"),
            diagnostics = list(
                chart.dim = list(
                    requested = "local.auto",
                    mode = "local.auto",
                    resolved = stats::median(dims),
                    by.anchor = as.integer(dims),
                    summary = list(
                        n.anchor = length(dims),
                        min = min(dims),
                        max = max(dims),
                        median = stats::median(dims),
                        n.unique = length(unique(dims))
                    ),
                    auto = TRUE,
                    local.auto = TRUE,
                    auto.diagnostics = NULL,
                    support.metric = "both",
                    selection.metric = "operator",
                    source.path = "test"
                )
            )
        )
        class(source) <- c(cls, "list")
        source
    }
    sources <- list(
        make.source("chart_kernel", c(1L, 1L, 2L, 2L)),
        make.source("local_likelihood", c(2L, 2L, 1L, 1L)),
        make.source("ps_lps", c(1L, 2L, 1L, 2L)),
        make.source("lps", c(2L, 1L, 2L, 1L))
    )

    for (source in sources) {
        fit <- normalize.density(
            source,
            X = X,
            keep.source.fit = FALSE
        )
        expect_s3_class(fit, "density_fit")
        expect_true(is.integer(fit$diagnostics$chart.dim$by.anchor))
        expect_equal(
            fit$diagnostics$chart.dim$by.anchor,
            source$diagnostics$chart.dim$by.anchor
        )
        expect_equal(length(fit$diagnostics$chart.dim$by.anchor), nrow(X))
        expect_true(is.character(fit$diagnostics$chart.dim$source.path))
    }
})

test_that("OD0 private smoothness helpers have deterministic placeholder behavior", {
    expect_equal(
        geosmooth:::.state.density.local.maxima.count(
            values = c(1, 3, 2),
            adj.list = list(2L, c(1L, 3L), 2L)
        ),
        1L
    )
    expect_equal(
        geosmooth:::.state.density.local.maxima.count(
            values = c(0, 0, 1, 1, 0),
            adj.list = list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
        ),
        0L
    )

    summary <- geosmooth:::.state.density.raw.basin.summary(
        basin.assignment = c("a", "a", "b"),
        rho = c(0.2, 0.3, 0.5)
    )
    expect_equal(summary$raw.basin.size.summary$size, c(2L, 1L))
    expect_equal(summary$raw.basin.mass.summary$mass, c(0.5, 0.5),
                 tolerance = 1e-12)
})

test_that("OD0 finalize wires strict local-maxima smoothness diagnostics", {
    X <- matrix(seq(0, 1, length.out = 5), ncol = 1L)
    fit <- geosmooth:::.state.density.finalize(
        method.id = "test",
        X = X,
        fitted.raw = c(0, 0.2, 0.1, 0.5, 0.1),
        empirical.rho = rep(NA_real_, 5L)
    )
    expect_equal(fit$smoothness$n.local.maxima, 2L)
    expect_identical(fit$smoothness$local.maxima.reason,
                     "computed_from_auto_1d_path")

    plateau <- geosmooth:::.state.density.finalize(
        method.id = "test",
        X = X,
        fitted.raw = c(0, 0, 0, 0, 0),
        empirical.rho = rep(NA_real_, 5L)
    )
    expect_equal(plateau$smoothness$n.local.maxima, 0L)

    multivariate <- geosmooth:::.state.density.finalize(
        method.id = "test",
        X = cbind(seq(0, 1, length.out = 5), 0),
        fitted.raw = c(0, 0.2, 0.1, 0.5, 0.1),
        empirical.rho = rep(NA_real_, 5L)
    )
    expect_identical(multivariate$status, "ok")
    expect_true(is.na(multivariate$smoothness$n.local.maxima))
    expect_identical(
        multivariate$smoothness$local.maxima.reason,
        "not_computed_no_adjacency_for_multivariate_support"
    )
})
