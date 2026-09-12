make.plateau.curved.X <- function(n = 30L) {
    t <- seq(-1, 1, length.out = n)
    cbind(t, t^2, sin(2 * t), cos(2 * t))
}

test_that("plateau_kd selects one geometry-only k-d candidate for LPS", {
    X <- make.plateau.curved.X(30L)
    y <- sin(pi * seq(-1, 1, length.out = nrow(X)))

    fit <- fit.lps(
        X = X,
        y = y,
        support.grid = 7:13,
        degree.grid = 1L,
        kernel.grid = "gaussian",
        bandwidth.multiplier.grid = 1,
        coordinate.method = "local.pca",
        chart.dim.grid = 1:4,
        selection.strategy = "plateau_kd",
        backend = "R",
        design.basis = "orthogonal.polynomial.drop",
        ridge.multiplier.grid = c(0, 1e-10),
        ridge.condition.max = Inf,
        cv.folds = 3L,
        cv.seed = 21L
    )

    expect_s3_class(fit, "lps")
    expect_identical(fit$selection.strategy, "plateau_kd")
    expect_equal(nrow(fit$cv.table), 1L)
    expect_true(fit$selected$support.size[[1L]] %in% 7:13)
    expect_true(as.integer(fit$selected$chart.dim[[1L]]) %in% 1:4)
    expect_identical(
        fit$diagnostics$coupled.kd.selection$selection.strategy,
        "plateau_kd"
    )
    expect_true(isTRUE(
        fit$diagnostics$coupled.kd.selection$geometry.only
    ))
    expect_true(is.data.frame(
        fit$diagnostics$coupled.kd.selection$plateau.anchor.diagnostics
    ))
    expect_true(nrow(
        fit$diagnostics$coupled.kd.selection$plateau.anchor.diagnostics
    ) > 0L)
    expect_equal(
        fit$diagnostics$coupled.kd.selection$evaluated.candidates,
        1L
    )
})

test_that("plateau_kd requires an explicit numeric chart-dimension grid", {
    X <- make.plateau.curved.X(24L)
    y <- X[, 1]

    expect_error(
        fit.lps(
            X = X,
            y = y,
            support.grid = 7:9,
            degree.grid = 1L,
            kernel.grid = "gaussian",
            coordinate.method = "local.pca",
            selection.strategy = "plateau_kd",
            backend = "R",
            cv.folds = 3L
        ),
        "requires 'chart.dim.grid'"
    )
})
