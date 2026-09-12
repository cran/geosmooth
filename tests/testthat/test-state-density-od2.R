make.od2.path.graph <- function(n) {
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

test_that("OD2 LPS count and binary no-repeat workflows agree after normalization", {
    X <- matrix(seq(0, 1, length.out = 24), ncol = 1L)
    subject.index <- c(4L, 8L, 13L, 19L, 23L)
    foldid <- rep(seq_len(4L), length.out = nrow(X))
    common <- list(
        foldid = foldid,
        support.grid = 7L,
        degree.grid = 0L,
        kernel.grid = "gaussian",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )

    count <- do.call(
        fit.subject.od,
        c(list(X = X, subject.index = subject.index, method = "lps_count"),
          common)
    )
    binary <- do.call(
        fit.subject.od,
        c(list(X = X, subject.index = subject.index,
               method = "lps_logistic_binary"),
          common)
    )

    expect_s3_class(count, "density_fit")
    expect_s3_class(binary, "density_fit")
    expect_identical(count$status, "ok")
    expect_identical(binary$status, "ok")
    expect_identical(count$method.id, "lps_count")
    expect_identical(binary$method.id, "lps_logistic_binary")
    expect_equal(sum(count$rho), 1, tolerance = 1e-12)
    expect_equal(sum(binary$rho), 1, tolerance = 1e-12)
    expect_equal(count$rho, binary$rho, tolerance = 1e-10)
    expect_true(is.finite(count$smoothness$n.local.maxima))
    expect_true(is.finite(binary$smoothness$n.local.maxima))

    count.raw <- count$diagnostics$source.fit$fitted.values
    binary.raw <- binary$diagnostics$source.fit$fitted.values
    expect_equal(count.raw * length(subject.index), binary.raw,
                 tolerance = 1e-10)
    expect_identical(count$diagnostics$response.summary$type,
                     "normalized_count_mass")
    expect_identical(binary$diagnostics$response.summary$type,
                     "binary_visit_indicator")
    expect_identical(
        binary$diagnostics$binary.workflow$probability.link,
        "identity_lps_least_squares_clipped"
    )
    expect_match(binary$diagnostics$binary.workflow$note,
                 "not the.*local-logistic IRLS path")
})

test_that("OD2 LPS occupation workflow reports density accounting in dimension > 1", {
    set.seed(401)
    n <- 30L
    X <- cbind(seq(0, 1, length.out = n),
               sin(seq(0, 2 * pi, length.out = n)))
    subject.index <- c(3L, 5L, 8L, 8L, 14L, 22L, 29L)
    graph <- make.od2.path.graph(n)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "lps_count",
        graph = graph,
        support.grid = 9L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = 1e10
    )

    expect_identical(fit$status, "ok")
    expect_identical(fit$method.id, "lps_count")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_true(all(fit$rho >= -1e-12))
    expect_equal(fit$empirical.rho, tabulate(subject.index, nbins = n) /
                     length(subject.index), tolerance = 1e-12)
    expect_true(all(c("raw.mass", "neg.mass", "clip.mass",
                      "normalization.constant") %in% names(fit$accounting)))
    expect_identical(fit$diagnostics$source.method, "fit.lps")
    expect_identical(fit$diagnostics$outcome.family, "gaussian")
    expect_true(is.finite(fit$smoothness$n.local.maxima))
    expect_equal(fit$subject$max.multiplicity, 2)
})

test_that("OD2 PS-LPS occupation workflow reports synchronization telemetry", {
    set.seed(402)
    n <- 26L
    X <- cbind(runif(n), runif(n))
    subject.index <- c(2L, 4L, 7L, 11L, 16L, 20L, 24L)
    foldid <- rep(seq_len(3L), length.out = n)
    graph <- make.od2.path.graph(n)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "ps_lps_count",
        graph = graph,
        foldid = foldid,
        support.grid = 9L,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        chart.dim = 2L,
        lambda.sync.grid = c(0, 0.1),
        lambda.sync.search = "grid",
        local.candidate.search = "full",
        lambda.ridge = 1e-8,
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = 1e10,
        sync.neighbor.size = 3L,
        cv.folds = 3L
    )

    expect_s3_class(fit, "density_fit")
    expect_identical(fit$status, "ok")
    expect_identical(fit$method.id, "ps_lps_count")
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
    expect_identical(fit$diagnostics$source.method, "fit.ps.lps")
    expect_true(is.finite(fit$diagnostics$sync.energy))
    expect_true(is.finite(fit$diagnostics$mean.sync.squared.disagreement))
    expect_true(is.data.frame(fit$diagnostics$lambda.sync.search.telemetry))
    expect_true(nrow(fit$diagnostics$lambda.sync.search.telemetry) > 0L)
    expect_true("lambda.sync" %in% names(fit$diagnostics$selection))
    expect_true(is.finite(fit$smoothness$n.local.maxima))
})

test_that("OD2 LPS defaults deactivate sparse subject charts", {
    X <- matrix(seq(0, 1, length.out = 30), ncol = 1L)
    subject.index <- c(14L, 16L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "lps_count",
        support.grid = 5L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "coordinates",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )

    activation <- fit$diagnostics$source.fit$diagnostics$chart.activation
    rows <- fit$diagnostics$source.fit$chart.activation.diagnostics
    expect_true(isTRUE(activation$enabled))
    expect_equal(activation$n.positive.min, 2L)
    expect_equal(activation$core.weight.rule, "chart_quantile")
    expect_equal(activation$core.weight.quantile, 0.25)
    expect_true(activation$n.inactive > 0L)
    expect_true(any(rows$reason == "no_subject_mass"))
    expect_true(any(rows$reason == "insufficient_positive_support"))
    expect_true(all(fit$diagnostics$source.fit$fitted.values.raw[!rows$active] == 0))
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})

test_that("OD2 LPS visit-CV uses native local-PCA backend under sparse activation", {
    X <- cbind(seq(0, 1, length.out = 24),
               sin(seq(0, 1, length.out = 24)))
    subject.index <- c(10L, 12L, 14L, 16L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "lps_count",
        od.cv = "visit",
        visit.cv.folds = 2L,
        support.grid = 7L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        chart.dim.grid = 1L,
        backend = "cpp.local.pca",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "mean"
    )

    source <- fit$diagnostics$source.fit
    expect_identical(source$backend, "cpp.local.pca")
    expect_identical(source$backend.used, "cpp.local.pca")
    expect_true(isTRUE(source$diagnostics$chart.activation$enabled))
    expect_true(any(!source$chart.activation.diagnostics$active))
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})

test_that("OD2 PS-LPS defaults remove inactive charts from synchronization", {
    X <- matrix(seq(0, 1, length.out = 24), ncol = 1L)
    subject.index <- c(11L, 13L)

    fit <- fit.subject.od(
        X = X,
        subject.index = subject.index,
        method = "ps_lps_count",
        support.size = 5L,
        degree = 1L,
        kernel = "tricube",
        chart.dim = 1L,
        lambda.sync.grid = 0.1,
        lambda.sync.selection = "fixed",
        lambda.sync.search = "grid",
        local.candidate.search = "full",
        lambda.ridge = 1e-8,
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-8),
        ridge.condition.max = 1e10,
        sync.neighbor.size = 2L,
        cv.folds = 3L
    )

    source <- fit$diagnostics$source.fit
    activation <- source$diagnostics$chart.activation
    frame.summary <- source$frame.design.summary
    expect_true(isTRUE(activation$enabled))
    expect_true(activation$n.inactive > 0L)
    expect_true(any(!frame.summary$active))
    expect_true(all(source$fitted.values[!frame.summary$active] == 0))
    if (length(source$sync.rows)) {
        active <- frame.summary$active
        expect_true(all(vapply(source$sync.rows, function(sr) {
            isTRUE(active[[sr$i]]) && isTRUE(active[[sr$j]])
        }, logical(1L))))
    }
    expect_equal(sum(fit$rho), 1, tolerance = 1e-12)
})

test_that("OD2 smoother workflows reject reserved pass-through arguments", {
    X <- matrix(seq(0, 1, length.out = 12), ncol = 1L)
    expect_error(
        fit.subject.od(
            X = X,
            subject.index = c(2L, 5L, 8L),
            method = "lps_count",
            y = rep(0, nrow(X))
        ),
        "reserved argument"
    )
})
