tier0.poly.column.count <- function(dim, degree) {
    1L + dim + if (degree >= 2L) {
        as.integer(dim * (dim + 1L) / 2L)
    } else {
        0L
    }
}

tier0.polynomial.truth <- function(X, degree) {
    y <- 0.7 + 0.4 * X[, 1L]
    if (ncol(X) >= 2L) {
        y <- y - 0.3 * X[, 2L]
    }
    if (ncol(X) >= 3L) {
        y <- y + 0.2 * X[, 3L]
    }
    if (degree >= 2L) {
        y <- y + 0.5 * X[, 1L]^2
        if (ncol(X) >= 2L) {
            y <- y - 0.25 * X[, 1L] * X[, 2L] + 0.35 * X[, 2L]^2
        }
        if (ncol(X) >= 3L) {
            y <- y - 0.18 * X[, 3L]^2 + 0.12 * X[, 1L] * X[, 3L]
        }
    }
    y
}

tier0.orthonormal.frame <- function(ambient.dim, intrinsic.dim, seed) {
    set.seed(seed)
    qr.Q(qr(matrix(stats::rnorm(ambient.dim * intrinsic.dim),
                  nrow = ambient.dim)))[, seq_len(intrinsic.dim), drop = FALSE]
}

tier0.reproduction.fit <- function(X, y, support.size, degree, kernel,
                                   basis, coordinate.method = "coordinates",
                                   chart.dim = NULL) {
    fit.lps(
        X = X,
        y = y,
        foldid = rep(1:2, length.out = nrow(X)),
        support.grid = support.size,
        degree.grid = degree,
        kernel.grid = kernel,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        local.chart.method = "pca",
        backend = "R",
        design.basis = basis,
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
}

tier0.fixed.lps.fit <- function(X, y, coordinate.method = "coordinates",
                                chart.dim = NULL) {
    fit.lps(
        X = X,
        y = y,
        foldid = rep(1:2, length.out = nrow(X)),
        support.grid = 18L,
        degree.grid = 1L,
        kernel.grid = "tricube",
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        local.chart.method = "pca",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = 0,
        ridge.condition.max = Inf,
        unstable.action = "na"
    )
}

tier0.extract.smoother.matrix <- function(X, coordinate.method = "coordinates",
                                          chart.dim = NULL) {
    n <- nrow(X)
    S <- matrix(NA_real_, nrow = n, ncol = n)
    for (j in seq_len(n)) {
        y <- numeric(n)
        y[[j]] <- 1
        fit <- tier0.fixed.lps.fit(
            X = X,
            y = y,
            coordinate.method = coordinate.method,
            chart.dim = chart.dim
        )
        S[, j] <- fit$fitted.values
    }
    S
}

test_that("E0.1 LPS reproduces ambient polynomials represented by the local design", {
    set.seed(9101)
    kernels <- c("gaussian", "tricube", "epanechnikov", "triangular")
    bases <- c("orthogonal.polynomial.drop", "monomial",
               "weighted.qr", "weighted.qr.drop")
    cases <- expand.grid(
        ambient.dim = c(2L, 3L),
        degree = c(1L, 2L),
        kernel = kernels,
        basis = bases,
        stringsAsFactors = FALSE
    )

    for (ii in seq_len(nrow(cases))) {
        cc <- cases[ii, ]
        n <- 200L
        X <- matrix(
            stats::runif(n * cc$ambient.dim, -1, 1),
            ncol = cc$ambient.dim
        )
        y <- tier0.polynomial.truth(X, cc$degree)
        n.cols <- tier0.poly.column.count(cc$ambient.dim, cc$degree)
        support.size <- min(n - 1L, max(15L, 3L * n.cols))

        fit <- tier0.reproduction.fit(
            X = X,
            y = y,
            support.size = support.size,
            degree = cc$degree,
            kernel = cc$kernel,
            basis = cc$basis
        )
        err <- max(abs(fit$fitted.values - y), na.rm = TRUE)
        tol <- if (identical(cc$basis, "monomial")) 1e-6 else 1e-8

        expect_false(anyNA(fit$fitted.values))
        expect_lt(err, tol)
        expect_equal(nrow(fit$cv.table), 1L)
        expect_equal(
            fit$local.chart.diagnostics.summary$min.design.rank,
            n.cols
        )
    }
})

test_that("E0.1 LPS reproduces intrinsic polynomials on flat embedded subspaces", {
    set.seed(9102)
    kernels <- c("gaussian", "tricube", "epanechnikov", "triangular")
    bases <- c("orthogonal.polynomial.drop", "monomial",
               "weighted.qr", "weighted.qr.drop")
    cases <- expand.grid(
        intrinsic.dim = c(1L, 2L),
        degree = c(1L, 2L),
        kernel = kernels,
        basis = bases,
        stringsAsFactors = FALSE
    )

    for (ii in seq_len(nrow(cases))) {
        cc <- cases[ii, ]
        ambient.dim <- cc$intrinsic.dim + 2L
        n <- 200L
        U <- matrix(
            stats::runif(n * cc$intrinsic.dim, -1, 1),
            ncol = cc$intrinsic.dim
        )
        Q <- tier0.orthonormal.frame(
            ambient.dim = ambient.dim,
            intrinsic.dim = cc$intrinsic.dim,
            seed = 9200L + ii
        )
        X <- U %*% t(Q)
        y <- tier0.polynomial.truth(U, cc$degree)
        n.cols <- tier0.poly.column.count(cc$intrinsic.dim, cc$degree)
        support.size <- min(n - 1L, max(15L, 3L * n.cols))

        fit <- tier0.reproduction.fit(
            X = X,
            y = y,
            support.size = support.size,
            degree = cc$degree,
            kernel = cc$kernel,
            basis = cc$basis,
            coordinate.method = "local.pca",
            chart.dim = cc$intrinsic.dim
        )
        err <- max(abs(fit$fitted.values - y), na.rm = TRUE)
        tol <- if (identical(cc$basis, "monomial")) 1e-6 else 1e-8

        expect_false(anyNA(fit$fitted.values))
        expect_lt(err, tol)
        expect_equal(nrow(fit$cv.table), 1L)
        expect_equal(
            fit$local.chart.diagnostics.summary$min.design.rank,
            n.cols
        )
        expect_identical(
            fit$local.chart.diagnostics.summary$local.chart.method,
            "pca"
        )
        expect_identical(fit$local.chart.diagnostics.summary$fallback.count, 0L)
    }
})

test_that("E0.1 negative control does not reproduce an under-specified polynomial", {
    set.seed(9333)
    n <- 80L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    y <- tier0.polynomial.truth(X, degree = 2L)

    fit <- tier0.reproduction.fit(
        X = X,
        y = y,
        support.size = 28L,
        degree = 1L,
        kernel = "tricube",
        basis = "orthogonal.polynomial.drop"
    )
    err <- max(abs(fit$fitted.values - y), na.rm = TRUE)

    expect_false(anyNA(fit$fitted.values))
    expect_gt(err, 1e-3)
})

test_that("E0.2 fixed-configuration LPS is a linear smoother in ambient coordinates", {
    set.seed(9103)
    n <- 36L
    X <- matrix(stats::runif(2L * n, -1, 1), ncol = 2L)
    S <- tier0.extract.smoother.matrix(X)

    y1 <- stats::rnorm(n)
    y2 <- stats::rnorm(n)
    y3 <- 0.6 * y1 - 1.4 * y2

    f1 <- tier0.fixed.lps.fit(X, y1)$fitted.values
    f2 <- tier0.fixed.lps.fit(X, y2)$fitted.values
    f3 <- tier0.fixed.lps.fit(X, y3)$fitted.values

    expect_equal(drop(S %*% y1), f1, tolerance = 1e-10)
    expect_equal(drop(S %*% y2), f2, tolerance = 1e-10)
    expect_equal(drop(S %*% y3), f3, tolerance = 1e-10)
    expect_equal(f3, 0.6 * f1 - 1.4 * f2, tolerance = 1e-10)

    eps.index <- seq(1L, n, length.out = 6L)
    for (jj in eps.index) {
        y.bump <- y1
        y.bump[[jj]] <- y.bump[[jj]] + 1
        direct.delta <- tier0.fixed.lps.fit(X, y.bump)$fitted.values - f1
        expect_equal(direct.delta, S[, jj], tolerance = 1e-10)
    }

    expect_true(is.finite(sum(diag(S))))
    expect_gt(sum(diag(S)), 0)
    expect_equal(sum(diag(S)), sum(vapply(seq_len(n), function(jj) {
        y.bump <- y1
        y.bump[[jj]] <- y.bump[[jj]] + 1
        tier0.fixed.lps.fit(X, y.bump)$fitted.values[[jj]] - f1[[jj]]
    }, numeric(1L))), tolerance = 1e-10)
})

test_that("E0.2 fixed-configuration LPS is a linear smoother in local PCA charts", {
    set.seed(9104)
    n <- 34L
    U <- matrix(stats::runif(n, -1, 1), ncol = 1L)
    X <- cbind(U[, 1L], 0.5 * U[, 1L])
    S <- tier0.extract.smoother.matrix(
        X,
        coordinate.method = "local.pca",
        chart.dim = 1L
    )

    y1 <- stats::rnorm(n)
    y2 <- stats::rnorm(n)
    y3 <- -0.25 * y1 + 1.8 * y2

    f1 <- tier0.fixed.lps.fit(
        X, y1, coordinate.method = "local.pca", chart.dim = 1L
    )$fitted.values
    f2 <- tier0.fixed.lps.fit(
        X, y2, coordinate.method = "local.pca", chart.dim = 1L
    )$fitted.values
    f3 <- tier0.fixed.lps.fit(
        X, y3, coordinate.method = "local.pca", chart.dim = 1L
    )$fitted.values

    expect_equal(drop(S %*% y1), f1, tolerance = 1e-10)
    expect_equal(drop(S %*% y2), f2, tolerance = 1e-10)
    expect_equal(drop(S %*% y3), f3, tolerance = 1e-10)
    expect_equal(f3, -0.25 * f1 + 1.8 * f2, tolerance = 1e-10)
    expect_true(all(is.finite(diag(S))))
    expect_gt(sum(diag(S)), 0)
    expect_equal(sum(diag(S)), sum(vapply(seq_len(n), function(jj) {
        y.bump <- y1
        y.bump[[jj]] <- y.bump[[jj]] + 1
        tier0.fixed.lps.fit(
            X,
            y.bump,
            coordinate.method = "local.pca",
            chart.dim = 1L
        )$fitted.values[[jj]] - f1[[jj]]
    }, numeric(1L))), tolerance = 1e-10)
})
