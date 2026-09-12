test_that("PS-LPS normalized product overlap weights have pair mass equal to overlap size", {
    set.seed(11)
    X <- matrix(runif(80), ncol = 2)
    y <- sin(X[, 1])
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 10L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim.by.anchor = rep(2L, nrow(X))
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 3L,
        overlap.weight = "normalized.product"
    )
    expect_gt(length(sync.rows), 0L)
    mass.error <- vapply(sync.rows, function(sr) {
        abs(sum(sr$omega) - length(sr$point))
    }, numeric(1L))
    expect_lt(max(mass.error), 1e-6)
})

test_that("PS-LPS with zero synchronization and zero ridge reproduces ordinary LPS fitted values", {
    set.seed(12)
    X <- matrix(runif(120), ncol = 2)
    y <- sin(4 * X[, 1]) + X[, 2]^2
    foldid <- rep(seq_len(5L), length.out = nrow(X))
    lps <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = 18L,
        degree.grid = 2L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = 2L,
        backend = "R",
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 18L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim = 2L,
        lambda.sync.grid = 0,
        lambda.ridge = 0,
        sync.neighbor.size = 3L,
        design.basis = "monomial",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf
    )
    expect_equal(ps$selected$lambda.ridge[[1L]], 0)
    expect_equal(max(abs(lps$fitted.values - ps$fitted.values)), 0,
                 tolerance = 1e-10)
})

test_that("PS-LPS weighted QR drop and adaptive ridge return telemetry", {
    set.seed(291)
    X <- matrix(runif(36 * 3), 36, 3)
    y <- sin(2 * pi * X[, 1]) + rnorm(nrow(X), sd = 0.03)
    foldid <- rep(1:4, length.out = length(y))

    fit <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 8L,
        degree = 2L,
        kernel = "tricube",
        chart.dim = 3L,
        lambda.sync.grid = c(0, 1e-3),
        lambda.sync.search = "grid",
        lambda.ridge = 0,
        design.basis = "weighted.qr.drop",
        design.drop.tol = 1e-7,
        ridge.multiplier.grid = c(0, 1e-10, 1e-8, 1e-6, 1e-4),
        ridge.condition.max = 1e8
    )

    kept <- fit$frame.design.summary$design.columns.kept
    original <- fit$frame.design.summary$design.columns.original
    expect_true(any(kept < original))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_true("ridge.status" %in% names(fit$cv.table))
    expect_true(all(nzchar(fit$cv.table$ridge.status)))
})

test_that("PS-LPS zero-sync local failures are not weighted-mean rescues", {
    set.seed(292)
    X <- matrix(runif(10 * 3), 10, 3)
    y <- seq_len(nrow(X))
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 5L,
        degree = 2L,
        kernel = "tricube",
        chart.dim.by.anchor = rep(3L, nrow(X)),
        design.basis = "monomial"
    )
    solved <- .ps.lps.solve.independent(
        frames = frames,
        y = y,
        response.weights = rep(1, length(y)),
        lambda.ridge = 0,
        ridge.multiplier.grid = 0,
        ridge.condition.max = 1
    )

    expect_true(any(grepl("^unstable_", solved$ridge.status)))
    expect_true(any(!is.finite(solved$fitted.values)))
    expect_false(any(grepl("fallback_mean", solved$ridge.status)))
})

test_that("PS-LPS orthogonal polynomial frames use transformed chart designs", {
    set.seed(293)
    X <- matrix(runif(40 * 3), 40, 3)
    y <- sin(X[, 1])
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 8L,
        degree = 2L,
        kernel = "tricube",
        chart.dim.by.anchor = rep(3L, nrow(X)),
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-8
    )

    expect_true(all(vapply(frames, `[[`, character(1L),
                           "solver.design.basis") ==
                    "orthogonal.polynomial.transformed"))
    expect_true(any(vapply(frames, `[[`, integer(1L), "q") <
                    vapply(frames, `[[`, integer(1L),
                           "design.columns.original")))
})

test_that("PS-LPS zero synchronization matches LPS with orthogonal basis", {
    set.seed(294)
    X <- matrix(runif(90 * 3), 90, 3)
    y <- cos(3 * X[, 1]) + X[, 2]^2
    foldid <- rep(seq_len(5L), length.out = nrow(X))
    lps <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = 26L,
        degree.grid = 2L,
        kernel.grid = "tricube",
        coordinate.method = "local.pca",
        chart.dim = 3L,
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-8,
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 26L,
        degree = 2L,
        kernel = "tricube",
        chart.dim = 3L,
        lambda.sync.grid = 0,
        lambda.ridge = 0,
        design.basis = "orthogonal.polynomial.drop",
        design.drop.tol = 1e-8,
        sync.neighbor.size = 3L
    )

    expect_equal(ps$fitted.values, lps$fitted.values, tolerance = 1e-8)
    expect_true(is.finite(ps$cv.table$cv.rmse.observed[[1L]]))
})

test_that("PS-LPS zero synchronization and zero ridge reproduces LPS local.auto vector-chart fitted values", {
    set.seed(1)
    X <- matrix(runif(150), ncol = 3)
    X[seq_len(25L), 3L] <- 0.02 * runif(25L)
    y <- sin(4 * X[, 1]) + X[, 2]^2 - 0.5 * X[, 3]
    foldid <- rep(seq_len(5L), length.out = nrow(X))
    lps <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = 18L,
        degree.grid = 2L,
        kernel.grid = "gaussian",
        coordinate.method = "local.pca",
        chart.dim = "local.auto",
        backend = "R"
    )
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = lps$selected$support.size[[1L]],
        degree = lps$selected$degree[[1L]],
        kernel = lps$selected$kernel[[1L]],
        chart.dim = lps$diagnostics$chart.dim$by.anchor,
        lambda.sync.grid = 0,
        lambda.ridge = 0,
        sync.neighbor.size = 3L
    )
    expect_equal(length(lps$diagnostics$chart.dim$by.anchor), nrow(X))
    expect_gt(length(unique(lps$diagnostics$chart.dim$by.anchor)), 1L)
    expect_equal(max(abs(lps$fitted.values - ps$fitted.values)), 0,
                 tolerance = 1e-10)
})

test_that("PS-LPS chart.dim auto matches explicit resolved scalar dimension", {
    set.seed(21)
    X <- matrix(runif(150), ncol = 3)
    y <- sin(3 * X[, 1]) + 0.5 * X[, 2]^2 - X[, 3]
    foldid <- rep(seq_len(5L), length.out = nrow(X))
    fit.auto <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 18L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim = "auto",
        lambda.sync.grid = c(0, 0.1),
        lambda.ridge = 0,
        sync.neighbor.size = 3L,
        local.candidate.search = "full"
    )
    fit.explicit <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 18L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim = fit.auto$chart.dim,
        lambda.sync.grid = c(0, 0.1),
        lambda.ridge = 0,
        sync.neighbor.size = 3L
    )
    expect_identical(fit.auto$requested.chart.dim, "auto")
    expect_identical(fit.auto$chart.dim.mode, "global.auto")
    expect_true(all(fit.auto$diagnostics$chart.dim$by.anchor ==
                    fit.auto$chart.dim))
    expect_equal(fit.auto$fitted.values, fit.explicit$fitted.values,
                 tolerance = 1e-10)
    expect_equal(fit.auto$selected$lambda.sync,
                 fit.explicit$selected$lambda.sync)
})

test_that("PS-LPS chart.dim local.auto matches explicit resolved vector dimensions", {
    set.seed(22)
    X <- matrix(runif(180), ncol = 4)
    X[seq_len(40L), 4L] <- 0.005 * runif(40L)
    y <- cos(3 * X[, 1]) + X[, 2] - 0.2 * X[, 3]^2
    foldid <- rep(seq_len(5L), length.out = nrow(X))
    fit.local <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 20L,
        degree = 2L,
        kernel = "tricube",
        chart.dim = "local.auto",
        lambda.sync.grid = c(0, 0.1),
        lambda.ridge = 0,
        sync.neighbor.size = 3L
    )
    fit.explicit <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 20L,
        degree = 2L,
        kernel = "tricube",
        chart.dim = fit.local$diagnostics$chart.dim$by.anchor,
        lambda.sync.grid = c(0, 0.1),
        lambda.ridge = 0,
        sync.neighbor.size = 3L
    )
    expect_identical(fit.local$requested.chart.dim, "local.auto")
    expect_identical(fit.local$chart.dim.mode, "local.auto")
    expect_equal(length(fit.local$diagnostics$chart.dim$by.anchor), nrow(X))
    expect_gt(length(unique(fit.local$diagnostics$chart.dim$by.anchor)), 1L)
    expect_equal(fit.local$fitted.values, fit.explicit$fitted.values,
                 tolerance = 1e-10)
    expect_equal(fit.local$selected$lambda.sync,
                 fit.explicit$selected$lambda.sync)
})

test_that("PS-LPS local grid selection matches explicit fixed-candidate search", {
    set.seed(23)
    X <- matrix(runif(120), ncol = 3)
    y <- sin(2 * X[, 1]) + X[, 2]^2 - 0.25 * X[, 3]
    foldid <- rep(seq_len(5L), length.out = nrow(X))
    support.grid <- c(14L, 18L)
    kernel.grid <- c("gaussian", "tricube")
    grid.fit <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = support.grid,
        degree.grid = 2L,
        kernel.grid = kernel.grid,
        chart.dim = "auto",
        lambda.sync.grid = c(0, 0.1),
        lambda.ridge = 0,
        sync.neighbor.size = 3L
    )
    explicit <- do.call(rbind, lapply(support.grid, function(k) {
        do.call(rbind, lapply(kernel.grid, function(ker) {
            fit <- fit.ps.lps(
                X = X,
                y = y,
                foldid = foldid,
                support.size = k,
                degree = 2L,
                kernel = ker,
                chart.dim = "auto",
                lambda.sync.grid = c(0, 0.1),
                lambda.ridge = 0,
                sync.neighbor.size = 3L
            )
            data.frame(
                support.size = k,
                kernel = ker,
                cv.rmse.observed = fit$selected$cv.rmse.observed[[1L]],
                lambda.sync = fit$selected$lambda.sync[[1L]],
                stringsAsFactors = FALSE
            )
        }))
    }))
    explicit.best <- explicit[order(explicit$cv.rmse.observed,
                                    explicit$lambda.sync,
                                    explicit$support.size,
                                    explicit$kernel), ][1L, ]
    expect_equal(nrow(grid.fit$local.candidate.table), 4L)
    expect_equal(nrow(grid.fit$lambda.cv.table), 8L)
    expect_equal(grid.fit$support.size, explicit.best$support.size)
    expect_identical(grid.fit$kernel, explicit.best$kernel)
    expect_equal(grid.fit$selected$cv.rmse.observed[[1L]],
                 explicit.best$cv.rmse.observed)
    expect_equal(grid.fit$selected$lambda.sync[[1L]],
                 explicit.best$lambda.sync)
})

test_that("PS-LPS support-grid default uses screened routine search", {
    set.seed(2301)
    X <- matrix(runif(150), ncol = 3)
    y <- cos(3 * X[, 1]) + X[, 2] - 0.25 * X[, 3]^2
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = rep(seq_len(5L), length.out = length(y)),
        support.grid = 10:16,
        degree.grid = 2L,
        kernel.grid = c("gaussian", "tricube"),
        chart.dim = 2L,
        lambda.sync.grid = c(0, 0.1, 1),
        lambda.sync.search = "guarded",
        lambda.sync.search.control = list(max.candidates = 3L,
                                          boundary.expand = FALSE),
        local.candidate.search.control = list(
            top.n = 2L,
            max.candidates = 4L,
            neighbor.radius = 0L,
            guard.support.quantiles = c(0, 1)
        ),
        lambda.ridge = 1e-8,
        sync.neighbor.size = 3L
    )

    evaluated <- ps$local.candidate.table$local.candidate.status ==
        "evaluated"
    expect_equal(ps$local.candidate.search, "screened")
    expect_equal(ps$selection.contract,
                 "screened_lps_cv_then_materialized_fold_cv_with_lambda_sync")
    expect_lt(sum(evaluated), nrow(ps$local.candidate.table))
    expect_lte(sum(evaluated), 4L)
    expect_true(all(ps$lambda.cv.table$local.candidate.id %in%
                    ps$local.candidate.table$local.candidate.id[evaluated]))
})

test_that("PS-LPS with positive ridge and zero synchronization nests ridge-LPS", {
    set.seed(13)
    X <- matrix(runif(90), ncol = 3)
    y <- cos(3 * X[, 1]) + 0.2 * X[, 2]
    lambda.ridge <- 1e-8
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 15L,
        degree = 2L,
        kernel = "tricube",
        chart.dim.by.anchor = rep(2L, nrow(X))
    )
    ridge.lps <- .ps.lps.solve(
        frames = frames,
        y = y,
        response.weights = rep(1, length(y)),
        lambda.sync = 0,
        lambda.ridge = lambda.ridge,
        sync.rows = list()
    )
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = rep(seq_len(5L), length.out = nrow(X)),
        support.size = 15L,
        degree = 2L,
        kernel = "tricube",
        chart.dim = 2L,
        lambda.sync.grid = 0,
        lambda.ridge = lambda.ridge,
        sync.neighbor.size = 3L
    )
    expect_equal(ps$fitted.values, ridge.lps$fitted.values,
                 tolerance = 1e-12)
    expect_gt(ps$ridge.max, 0)
})

test_that("PS-LPS cached solve reproduces direct positive-synchronization solve", {
    set.seed(15)
    X <- matrix(runif(120), ncol = 3)
    y <- sin(3 * X[, 1]) + X[, 2]^2 - 0.25 * X[, 3]
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 16L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim.by.anchor = rep(2L, nrow(X))
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 4L,
        overlap.weight = "normalized.product"
    )
    cache <- .ps.lps.prepare.system.cache(frames, sync.rows)
    response.weights <- rep(1, length(y))
    direct <- .ps.lps.solve(
        frames = frames,
        y = y,
        response.weights = response.weights,
        lambda.sync = 3,
        lambda.ridge = 1e-8,
        sync.rows = sync.rows
    )
    cached <- .ps.lps.solve.cached(
        cache = cache,
        y = y,
        response.weights = response.weights,
        lambda.sync = 3,
        lambda.ridge = 1e-8
    )
    expect_equal(cached$fitted.values, direct$fitted.values,
                 tolerance = 1e-10)
    expect_equal(cached$coefficients, direct$coefficients,
                 tolerance = 1e-10)
    expect_equal(cached$sync.energy, direct$sync.energy,
                 tolerance = 1e-10)
    expect_equal(cached$n.system.rows, direct$n.system.rows)
    expect_equal(cached$n.system.cols, direct$n.system.cols)
})

test_that("PS-LPS cached solve reproduces direct fold-weighted and coefficients-only solves", {
    set.seed(16)
    X <- matrix(runif(160), ncol = 4)
    y <- cos(2 * X[, 1]) + X[, 2] - X[, 3]^2
    chart.dim <- rep(c(2L, 3L), length.out = nrow(X))
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 18L,
        degree = 2L,
        kernel = "tricube",
        chart.dim.by.anchor = chart.dim
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 3L,
        overlap.weight = "normalized.product"
    )
    cache <- .ps.lps.prepare.system.cache(frames, sync.rows)
    foldid <- rep(seq_len(5L), length.out = length(y))
    response.weights <- as.numeric(foldid != 2L)
    direct <- .ps.lps.solve(
        frames = frames,
        y = y,
        response.weights = response.weights,
        lambda.sync = 10,
        lambda.ridge = 1e-8,
        sync.rows = sync.rows
    )
    cached <- .ps.lps.solve.cached(
        cache = cache,
        y = y,
        response.weights = response.weights,
        lambda.sync = 10,
        lambda.ridge = 1e-8
    )
    expect_equal(cached$fitted.values, direct$fitted.values,
                 tolerance = 1e-10)
    direct.diag <- .ps.lps.solve(
        frames = frames,
        y = y,
        response.weights = response.weights,
        lambda.sync = 10,
        lambda.ridge = 1e-8,
        sync.rows = sync.rows,
        coefficients.only = TRUE
    )
    cached.diag <- .ps.lps.solve.cached(
        cache = cache,
        y = y,
        response.weights = response.weights,
        lambda.sync = 10,
        lambda.ridge = 1e-8,
        coefficients.only = TRUE
    )
    expect_equal(cached.diag$total.local.gcv.ps,
                 direct.diag$total.local.gcv.ps,
                 tolerance = 1e-10)
    expect_equal(cached.diag$mean.sync.squared.disagreement,
                 direct.diag$mean.sync.squared.disagreement,
                 tolerance = 1e-10)
})

test_that("PS-LPS normal cache reproduces direct positive-synchronization solve", {
    set.seed(18)
    X <- matrix(runif(140), ncol = 4)
    y <- sin(2 * X[, 1]) - 0.3 * X[, 2] + X[, 3]^2
    chart.dim <- rep(c(2L, 3L), length.out = nrow(X))
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 17L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim.by.anchor = chart.dim
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 4L,
        overlap.weight = "normalized.product"
    )
    cache <- .ps.lps.prepare.system.cache(frames, sync.rows)
    response.weights <- rep(1, length(y))
    normal.cache <- .ps.lps.prepare.normal.cache(
        cache = cache,
        y = y,
        response.weights = response.weights,
        lambda.sync = 4
    )
    direct <- .ps.lps.solve(
        frames = frames,
        y = y,
        response.weights = response.weights,
        lambda.sync = 4,
        lambda.ridge = 1e-7,
        sync.rows = sync.rows
    )
    normal.cached <- .ps.lps.solve.normal.cached(
        normal.cache = normal.cache,
        lambda.ridge = 1e-7
    )
    expect_s3_class(normal.cache, "ps_lps_normal_cache")
    expect_gt(normal.cache$n.system.rows, 0L)
    expect_equal(normal.cached$fitted.values, direct$fitted.values,
                 tolerance = 1e-10)
    expect_equal(normal.cached$coefficients, direct$coefficients,
                 tolerance = 1e-10)
    expect_true(is.finite(
        normal.cached$solve.phase.timings$phase_normal_cache_sec
    ))
})

test_that("PS-LPS normal cache can be reused across ridge values", {
    set.seed(19)
    X <- matrix(runif(180), ncol = 3)
    y <- cos(2 * X[, 1]) + X[, 2] - 0.1 * X[, 3]
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 18L,
        degree = 2L,
        kernel = "tricube",
        chart.dim.by.anchor = rep(2L, nrow(X))
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 3L,
        overlap.weight = "normalized.product"
    )
    cache <- .ps.lps.prepare.system.cache(frames, sync.rows)
    response.weights <- as.numeric(rep(seq_len(5L), length.out = length(y)) != 3L)
    normal.cache <- .ps.lps.prepare.normal.cache(
        cache = cache,
        y = y,
        response.weights = response.weights,
        lambda.sync = 2
    )
    for (ridge in c(0, 1e-10, 1e-6)) {
        direct <- .ps.lps.solve(
            frames = frames,
            y = y,
            response.weights = response.weights,
            lambda.sync = 2,
            lambda.ridge = ridge,
            sync.rows = sync.rows
        )
        cached <- .ps.lps.solve.normal.cached(
            normal.cache = normal.cache,
            lambda.ridge = ridge
        )
        expect_equal(cached$fitted.values, direct$fitted.values,
                     tolerance = 1e-10)
        expect_equal(cached$total.local.gcv.ps, direct$total.local.gcv.ps,
                     tolerance = 1e-10)
    }
})

test_that("PS-LPS component cache reproduces direct solves across lambda.sync values", {
    set.seed(20)
    X <- matrix(runif(210), ncol = 3)
    y <- sin(3 * X[, 1]) + 0.5 * X[, 2]^2 - X[, 3]
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 19L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim.by.anchor = rep(2L, nrow(X))
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 4L,
        overlap.weight = "normalized.product"
    )
    cache <- .ps.lps.prepare.system.cache(frames, sync.rows)
    foldid <- rep(seq_len(5L), length.out = length(y))
    response.weights <- as.numeric(foldid != 4L)
    component.cache <- .ps.lps.prepare.component.cache(
        cache = cache,
        y = y,
        response.weights = response.weights
    )
    expect_s3_class(component.cache, "ps_lps_component_cache")
    expect_gt(component.cache$data.nrow, 0L)
    expect_gt(component.cache$sync.nrow, 0L)
    for (lambda.sync in c(0.25, 1, 8)) {
        direct <- .ps.lps.solve(
            frames = frames,
            y = y,
            response.weights = response.weights,
            lambda.sync = lambda.sync,
            lambda.ridge = 1e-8,
            sync.rows = sync.rows
        )
        cached <- .ps.lps.solve.component.cached(
            component.cache = component.cache,
            lambda.sync = lambda.sync,
            lambda.ridge = 1e-8
        )
        expect_equal(cached$fitted.values, direct$fitted.values,
                     tolerance = 1e-10)
        expect_equal(cached$coefficients, direct$coefficients,
                     tolerance = 1e-10)
        expect_equal(cached$total.local.gcv.ps, direct$total.local.gcv.ps,
                     tolerance = 1e-10)
    }
})

test_that("PS-LPS max-dimension local PCA supports reproduce smaller chart fits", {
    set.seed(221)
    n <- 42L
    t <- seq(-1, 1, length.out = n)
    X <- cbind(t, t^2, sin(2 * t), cos(3 * t), t^3)
    y <- sin(4 * t) + 0.25 * t^2
    support.size <- 13L
    kernel <- "tricube"
    degree <- 2L
    max.dim <- 4L
    max.supports <- rcpp_ps_lps_local_pca_supports(
        X = X,
        support_size = support.size,
        chart_dim_by_anchor = rep(max.dim, n),
        kernel = kernel
    )

    for (chart.dim in c(1L, 2L, 4L)) {
        direct <- .ps.lps.prepare.geometry.cache(
            X = X,
            support.size = support.size,
            degree = degree,
            kernel = kernel,
            chart.dim = chart.dim,
            auto.chart.support.metric = "coordinates",
            auto.chart.selection.metric = "coordinates",
            sync.neighbor.size = 3L,
            overlap.weight = "normalized.product",
            design.basis = "orthogonal.polynomial.drop",
            design.drop.tol = sqrt(.Machine$double.eps)
        )
        reused <- .ps.lps.prepare.geometry.cache(
            X = X,
            support.size = support.size,
            degree = degree,
            kernel = kernel,
            chart.dim = chart.dim,
            auto.chart.support.metric = "coordinates",
            auto.chart.selection.metric = "coordinates",
            sync.neighbor.size = 3L,
            overlap.weight = "normalized.product",
            design.basis = "orthogonal.polynomial.drop",
            design.drop.tol = sqrt(.Machine$double.eps),
            local.pca.supports = max.supports
        )
        expect_equal(
            .ps.lps.frame.design.summary(reused$frames),
            .ps.lps.frame.design.summary(direct$frames)
        )
        for (lambda.sync in c(0, 0.1)) {
            direct.fit <- .ps.lps.solve(
                frames = direct$frames,
                y = y,
                response.weights = rep(1, n),
                lambda.sync = lambda.sync,
                lambda.ridge = 1e-8,
                sync.rows = direct$sync.rows
            )
            reused.fit <- .ps.lps.solve(
                frames = reused$frames,
                y = y,
                response.weights = rep(1, n),
                lambda.sync = lambda.sync,
                lambda.ridge = 1e-8,
                sync.rows = reused$sync.rows
            )
            expect_equal(reused.fit$fitted.values,
                         direct.fit$fitted.values,
                         tolerance = 1e-9)
        }
    }
})

test_that("fit.ps.lps component-cache integration matches direct mixed-grid tuning loop", {
    set.seed(22)
    X <- matrix(runif(150), ncol = 3)
    y <- sin(2 * X[, 1]) + X[, 2]^2 - 0.2 * X[, 3]
    foldid <- rep(seq_len(5L), length.out = length(y))
    lambda.grid <- c(0, 0.2, 1, 5)
    support.size <- 15L
    degree <- 2L
    kernel <- "gaussian"
    chart.dim <- 2L
    lambda.ridge <- 1e-8
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = support.size,
        degree = degree,
        kernel = kernel,
        chart.dim.by.anchor = rep(chart.dim, nrow(X))
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 3L,
        overlap.weight = "normalized.product"
    )
    direct.cv <- data.frame(
        lambda.sync = lambda.grid,
        cv.rmse.observed = NA_real_
    )
    for (ll in seq_along(lambda.grid)) {
        pred <- rep(NA_real_, length(y))
        for (fold in sort(unique(foldid))) {
            fit.fold <- .ps.lps.solve(
                frames = frames,
                y = y,
                response.weights = as.numeric(foldid != fold),
                lambda.sync = lambda.grid[[ll]],
                lambda.ridge = lambda.ridge,
                sync.rows = sync.rows
            )
            pred[foldid == fold] <- fit.fold$fitted.values[foldid == fold]
        }
        direct.cv$cv.rmse.observed[[ll]] <- .klp.rmse(pred, y)
    }
    direct.best <- order(direct.cv$cv.rmse.observed,
                         direct.cv$lambda.sync)[[1L]]
    direct.final <- .ps.lps.solve(
        frames = frames,
        y = y,
        response.weights = rep(1, length(y)),
        lambda.sync = direct.cv$lambda.sync[[direct.best]],
        lambda.ridge = lambda.ridge,
        sync.rows = sync.rows
    )
    fitted.cached <- fit.ps.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.size = support.size,
        degree = degree,
        kernel = kernel,
        chart.dim = chart.dim,
        lambda.sync.grid = lambda.grid,
        lambda.ridge = lambda.ridge,
        sync.neighbor.size = 3L
    )
    expect_equal(fitted.cached$cache.backend, "component")
    expect_equal(fitted.cached$cv.table$cv.rmse.observed,
                 direct.cv$cv.rmse.observed,
                 tolerance = 1e-10)
    expect_equal(fitted.cached$selected$lambda.sync[[1L]],
                 direct.cv$lambda.sync[[direct.best]])
    expect_equal(fitted.cached$fitted.values, direct.final$fitted.values,
                 tolerance = 1e-10)
    expect_true(is.finite(
        fitted.cached$solve.phase.timings$phase_component_cache_sec
    ))
})

test_that("PS-LPS positive-sync caches reject nonpositive lambda.sync", {
    set.seed(21)
    X <- matrix(runif(90), ncol = 3)
    y <- X[, 1] - X[, 2]^2
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = y,
        support.size = 15L,
        degree = 2L,
        kernel = "tricube",
        chart.dim.by.anchor = rep(2L, nrow(X))
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = 3L,
        overlap.weight = "normalized.product"
    )
    cache <- .ps.lps.prepare.system.cache(frames, sync.rows)
    component.cache <- .ps.lps.prepare.component.cache(
        cache = cache,
        y = y,
        response.weights = rep(1, length(y))
    )
    expect_error(
        .ps.lps.prepare.normal.cache(
            cache = cache,
            y = y,
            response.weights = rep(1, length(y)),
            lambda.sync = 0
        ),
        "must be positive"
    )
    expect_error(
        .ps.lps.solve.component.cached(
            component.cache = component.cache,
            lambda.sync = 0,
            lambda.ridge = 1e-8
        ),
        "must be positive"
    )
})

test_that("PS-LPS reports synchronization energy even when lambda.sync is zero", {
    set.seed(14)
    X <- matrix(runif(80), ncol = 2)
    y <- sin(5 * X[, 1]) + rnorm(nrow(X), sd = 0.01)
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = rep(seq_len(4L), length.out = nrow(X)),
        support.size = 12L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim = 2L,
        lambda.sync.grid = 0,
        lambda.ridge = 0,
        sync.neighbor.size = 3L
    )
    expect_gt(ps$sync.energy, 0)
    expect_gt(ps$mean.sync.squared.disagreement, 0)
    expect_equal(ps$mean.sync.disagreement,
                 ps$mean.sync.squared.disagreement)
})

test_that("PS-LPS guarded lambda search recovers interior optimum", {
    evaluator <- function(lambda) {
        data.frame(
            lambda.sync = lambda,
            cv.rmse.observed = (log10(lambda) - log10(27))^2,
            stringsAsFactors = FALSE
        )
    }
    out <- .ps.lps.search.lambda.sync(
        evaluate = evaluator,
        lambda.grid = c(1, 3, 9, 27, 81, 243),
        control = list(rel.tol = 0, boundary.factor = 3)
    )
    expect_equal(out$selected$lambda.sync[[1L]], 27)
    expect_true("refine" %in% out$telemetry$stage)
})

test_that("PS-LPS guarded lambda search expands right boundary", {
    evaluator <- function(lambda) {
        data.frame(
            lambda.sync = lambda,
            cv.rmse.observed = (log10(lambda) - log10(27))^2,
            stringsAsFactors = FALSE
        )
    }
    out <- .ps.lps.search.lambda.sync(
        evaluate = evaluator,
        lambda.grid = c(0, 1, 3, 9),
        control = list(
            rel.tol = 0,
            boundary.factor = 3,
            max.boundary.expansions = 2L
        )
    )
    expect_equal(out$selected$lambda.sync[[1L]], 27)
    expect_true(any(out$telemetry$stage == "coarse" &
                    out$telemetry$lambda.sync == 0))
    expect_true(any(out$telemetry$stage == "boundary_expand_right" &
                    out$telemetry$lambda.sync == 27))
})

test_that("PS-LPS guarded lambda search expands near-best right boundary", {
    evaluator <- function(lambda) {
        cv <- c(
            "0" = 9,
            "1" = 2,
            "3" = 1,
            "9" = 1.005,
            "27" = 0.8,
            "81" = 0.7
        )[[as.character(lambda)]]
        data.frame(
            lambda.sync = lambda,
            cv.rmse.observed = cv,
            stringsAsFactors = FALSE
        )
    }
    out <- .ps.lps.search.lambda.sync(
        evaluate = evaluator,
        lambda.grid = c(0, 1, 3, 9),
        control = list(
            rel.tol = 0,
            boundary.guard.rel.tol = 0.01,
            boundary.factor = 3,
            max.boundary.expansions = 2L
        )
    )
    expect_equal(out$selected$lambda.sync[[1L]], 81)
    expect_equal(
        out$telemetry$lambda.sync[
            out$telemetry$stage == "boundary_expand_right"
        ],
        c(27, 81)
    )
})

test_that("PS-LPS guarded lambda search expands left positive boundary", {
    evaluator <- function(lambda) {
        data.frame(
            lambda.sync = lambda,
            cv.rmse.observed = (log10(lambda) - log10(1))^2,
            stringsAsFactors = FALSE
        )
    }
    out <- .ps.lps.search.lambda.sync(
        evaluate = evaluator,
        lambda.grid = c(0, 9, 27, 81),
        control = list(
            rel.tol = 0,
            boundary.factor = 3,
            max.boundary.expansions = 2L
        )
    )
    expect_equal(out$selected$lambda.sync[[1L]], 1)
    expect_true(any(out$telemetry$stage == "boundary_expand_left" &
                    out$telemetry$lambda.sync == 1))
})

test_that("PS-LPS guarded lambda search keeps zero as diagnostic", {
    evaluator <- function(lambda) {
        data.frame(
            lambda.sync = lambda,
            cv.rmse.observed = if (lambda == 0) 10 else (log10(lambda) + 1)^2,
            stringsAsFactors = FALSE
        )
    }
    out <- .ps.lps.search.lambda.sync(
        evaluate = evaluator,
        lambda.grid = c(0, 0.1, 1, 10),
        control = list(rel.tol = 0, boundary.factor = 10)
    )
    expect_equal(out$selected$lambda.sync[[1L]], 0.1)
    expect_true(0 %in% out$evaluated$lambda.sync)
})

test_that("PS-LPS guarded lambda search enforces max.candidates globally", {
    evaluator <- function(lambda) {
        data.frame(
            lambda.sync = lambda,
            cv.rmse.observed = (log10(lambda) - 1)^2,
            stringsAsFactors = FALSE
        )
    }
    out <- .ps.lps.search.lambda.sync(
        evaluate = evaluator,
        lambda.grid = c(1, 3, 10, 30, 100, 300),
        control = list(max.candidates = 2L, rel.tol = 0)
    )
    expect_lte(nrow(out$evaluated), 2L)
    expect_lte(length(unique(out$telemetry$lambda.sync)), 2L)
})

test_that("fit.ps.lps guarded lambda search returns telemetry", {
    set.seed(23)
    X <- matrix(runif(120), ncol = 3)
    y <- sin(2 * X[, 1]) + X[, 2]^2
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = rep(seq_len(5L), length.out = length(y)),
        support.size = 15L,
        degree = 2L,
        kernel = "gaussian",
        chart.dim = 2L,
        lambda.sync.grid = c(0, 0.1, 1, 10),
        lambda.sync.search = "guarded",
        lambda.sync.search.control = list(rel.tol = 0),
        lambda.ridge = 1e-8,
        sync.neighbor.size = 3L
    )
    expect_equal(ps$lambda.sync.search, "guarded")
    expect_true(nrow(ps$cv.table) <= 6L)
    expect_true(ps$selected$lambda.sync[[1L]] %in% ps$cv.table$lambda.sync)
    expect_true(all(c("stage", "lambda.sync", "boundary", "expansion",
                      "selected.after.stage") %in%
                    names(ps$lambda.sync.search.telemetry)))
})

test_that("fit.ps.lps selected-only lambda diagnostics preserve selected fit", {
    set.seed(231)
    X <- matrix(runif(90), ncol = 3)
    y <- sin(3 * X[, 1]) + 0.25 * X[, 2]
    foldid <- rep(seq_len(3L), length.out = length(y))
    common <- list(
        X = X,
        y = y,
        foldid = foldid,
        support.size = 12L,
        degree = 1L,
        kernel = "tricube",
        chart.dim = 2L,
        lambda.sync.grid = c(0, 0.1),
        lambda.sync.search = "grid",
        lambda.ridge = 1e-8,
        sync.neighbor.size = 3L
    )
    fit.all <- do.call(fit.ps.lps, c(common, list(lambda.diagnostics = "all")))
    fit.selected <- do.call(
        fit.ps.lps,
        c(common, list(lambda.diagnostics = "selected"))
    )

    expect_equal(fit.selected$selected$lambda.sync,
                 fit.all$selected$lambda.sync)
    expect_equal(fit.selected$fitted.values, fit.all$fitted.values,
                 tolerance = 1e-8)
    selected.lambda <- fit.selected$selected$lambda.sync[[1L]]
    selected.row <- fit.selected$cv.table$lambda.sync == selected.lambda
    expect_true(all(is.finite(fit.selected$cv.table$sync.energy[selected.row])))
    if (any(!selected.row)) {
        expect_true(all(is.na(fit.selected$cv.table$sync.energy[!selected.row])))
    }
})

test_that("fit.ps.lps fixed lambda can reuse geometry cache", {
    set.seed(232)
    X <- matrix(runif(90), ncol = 3)
    y <- cos(2 * X[, 1]) - 0.15 * X[, 2]
    common <- list(
        X = X,
        y = y,
        foldid = rep(seq_len(3L), length.out = length(y)),
        support.size = 12L,
        degree = 1L,
        kernel = "tricube",
        chart.dim = 2L,
        lambda.sync.grid = 0.1,
        lambda.sync.selection = "fixed",
        lambda.ridge = 1e-8,
        sync.neighbor.size = 3L
    )
    fit.plain <- do.call(fit.ps.lps, common)
    cache <- .ps.lps.prepare.geometry.cache(
        X = X,
        support.size = 12L,
        degree = 1L,
        kernel = "tricube",
        chart.dim = 2L,
        auto.chart.support.metric = "coordinates",
        auto.chart.selection.metric = "coordinates",
        sync.neighbor.size = 3L,
        overlap.weight = "normalized.product",
        design.basis = "monomial",
        design.drop.tol = sqrt(.Machine$double.eps)
    )
    fit.cached <- do.call(
        fit.ps.lps,
        c(common, list(ps.lps.geometry.cache = cache))
    )

    expect_equal(fit.cached$selected$lambda.sync[[1L]], 0.1)
    expect_true(is.na(fit.cached$selected$cv.rmse.observed[[1L]]))
    expect_equal(fit.cached$lambda.sync.search.telemetry$stage, "fixed")
    expect_true(isTRUE(fit.cached$geometry.cache.used))
    expect_equal(fit.cached$fitted.values, fit.plain$fitted.values,
                 tolerance = 1e-8)
})

test_that("fit.ps.lps screened local search evaluates an auditable subset", {
    set.seed(24)
    X <- matrix(runif(180), ncol = 3)
    y <- sin(2 * X[, 1]) + X[, 2]^2 - 0.5 * X[, 3]
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = rep(seq_len(5L), length.out = length(y)),
        support.grid = 10:16,
        degree.grid = 2L,
        kernel.grid = c("gaussian", "tricube"),
        chart.dim = 2L,
        lambda.sync.grid = c(0, 0.1, 1),
        lambda.sync.search = "guarded",
        lambda.sync.search.control = list(max.candidates = 3L,
                                          boundary.expand = FALSE),
        local.candidate.search = "screened",
        local.candidate.search.control = list(
            top.n = 2L,
            max.candidates = 4L,
            neighbor.radius = 0L,
            guard.support.quantiles = c(0, 1)
        ),
        lambda.ridge = 1e-8,
        sync.neighbor.size = 3L
    )

    expect_equal(ps$local.candidate.search, "screened")
    expect_equal(ps$selection.contract,
                 "screened_lps_cv_then_materialized_fold_cv_with_lambda_sync")
    expect_equal(nrow(ps$local.candidate.table), 14L)
    expect_lte(sum(ps$local.candidate.table$local.candidate.status ==
                       "evaluated"), 4L)
    evaluated <- ps$local.candidate.table[
        ps$local.candidate.table$local.candidate.status == "evaluated",
        ,
        drop = FALSE
    ]
    expect_true(all(c("local.candidate.elapsed.sec",
                      "lambda.search.elapsed.sec",
                      "unique.lambda.count",
                      "system.cache.elapsed.sec") %in%
                    names(ps$local.candidate.table)))
    expect_true(all(is.finite(evaluated$local.candidate.elapsed.sec)))
    expect_true(all(is.finite(evaluated$lambda.search.elapsed.sec)))
    expect_true(all(evaluated$unique.lambda.count > 0L))
    expect_true(is.list(ps$ps.lps.local.grid.timing))
    expect_equal(ps$ps.lps.local.grid.timing$evaluated_local_candidate_count,
                 nrow(evaluated))
    expect_true(any(ps$local.candidate.table$local.candidate.status ==
                        "screened_out"))
    expect_true(all(is.finite(
        ps$local.candidate.table$screening.cv.rmse.observed
    )))
    expect_true(ps$selected$local.candidate.id[[1L]] %in%
                    ps$local.candidate.table$local.candidate.id[
                        ps$local.candidate.table$local.candidate.status ==
                            "evaluated"
                    ])
    expect_true(all(ps$lambda.cv.table$local.candidate.id %in%
                    ps$local.candidate.table$local.candidate.id[
                        ps$local.candidate.table$local.candidate.status ==
                            "evaluated"
                    ]))
})

test_that("fit.ps.lps subgrid local search skips LPS screening pass", {
    set.seed(25)
    X <- matrix(runif(180), ncol = 3)
    y <- cos(2 * X[, 1]) - X[, 2] + 0.25 * X[, 3]^2
    ps <- fit.ps.lps(
        X = X,
        y = y,
        foldid = rep(seq_len(5L), length.out = length(y)),
        support.grid = 10:16,
        degree.grid = 2L,
        kernel.grid = c("gaussian", "tricube"),
        chart.dim = 2L,
        lambda.sync.grid = c(0, 0.1, 1),
        lambda.sync.search = "guarded",
        lambda.sync.search.control = list(max.candidates = 3L,
                                          boundary.expand = FALSE),
        local.candidate.search = "subgrid",
        local.candidate.search.control = list(
            max.candidates = 4L,
            guard.support.quantiles = c(0, 0.5, 1)
        ),
        lambda.ridge = 1e-8,
        sync.neighbor.size = 3L
    )

    expect_equal(ps$local.candidate.search, "subgrid")
    expect_equal(ps$selection.contract,
                 "subgrid_then_materialized_fold_cv_with_lambda_sync")
    expect_null(ps$local.candidate.screen.lps.selected)
    expect_equal(nrow(ps$local.candidate.table), 14L)
    expect_lte(sum(ps$local.candidate.table$local.candidate.status ==
                       "evaluated"), 4L)
    expect_true(all(c("local.candidate.elapsed.sec",
                      "lambda.search.elapsed.sec",
                      "unique.lambda.count") %in%
                    names(ps$local.candidate.table)))
    expect_true(any(ps$local.candidate.table$screening.reason ==
                        "subgrid_guard"))
    expect_true(any(ps$local.candidate.table$local.candidate.status ==
                        "screened_out"))
})
