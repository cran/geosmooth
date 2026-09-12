local_projector <- function(basis) {
    basis %*% t(basis)
}

projector_error <- function(basis, oracle) {
    sqrt(sum((local_projector(basis) - local_projector(oracle))^2))
}

test_that("Harlim second-order chart agrees with PCA on a flat affine plane", {
    grid <- expand.grid(x = -1:1, y = -1:1)
    X <- cbind(grid$x, grid$y, 0)
    center <- c(0, 0, 0)

    chart <- rcpp_local_second_order_svd_chart(
        X_support = X,
        center = center,
        chart_dim = 2L,
        center_mode = "anchor"
    )
    pca <- rcpp_local_pca_chart(
        X_support = X,
        center = center,
        chart_dim = 2L,
        center_mode = "anchor",
        dim_rule = "fixed"
    )

    expect_false(chart$fallback.used)
    expect_equal(chart$fallback.reason, "none")
    expect_equal(
        local_projector(chart$basis),
        local_projector(pca$basis),
        tolerance = 1e-10
    )
    expect_true(max(abs(chart$curvature.coefficients)) < 1e-10)
    expect_equal(
        chart$curvature.diagnostics$curvature.fitted.frobenius,
        0,
        tolerance = 1e-10
    )
    expect_equal(
        chart$curvature.monomials[, c("a", "b", "multiplier")],
        data.frame(
            a = c(1L, 2L, 1L),
            b = c(1L, 2L, 2L),
            multiplier = c(1, 1, 2)
        )
    )
})

test_that("Harlim second-order chart handles small local coordinate scale", {
    grid <- expand.grid(x = -1:1, y = -1:1)
    X <- 1e-9 * cbind(grid$x, grid$y, 0)

    chart <- rcpp_local_second_order_svd_chart(
        X_support = X,
        center = c(0, 0, 0),
        chart_dim = 2L,
        center_mode = "anchor"
    )

    expect_false(chart$fallback.used)
    expect_equal(chart$curvature.diagnostics$design.rank, 3L)
    expect_equal(chart$curvature.diagnostics$first.rank, 2L)
    expect_equal(chart$curvature.diagnostics$second.rank, 2L)
})

test_that("Harlim second-order chart corrects a noiseless parabola support", {
    t <- seq(-0.4, 0.4, length.out = 9L)
    X <- cbind(t, t^2)
    center <- c(0, 0)
    oracle <- matrix(c(1, 0), ncol = 1)

    chart <- rcpp_local_second_order_svd_chart(
        X_support = X,
        center = center,
        chart_dim = 1L,
        center_mode = "anchor"
    )
    pca <- rcpp_local_pca_chart(
        X_support = X,
        center = center,
        chart_dim = 1L,
        center_mode = "anchor",
        dim_rule = "fixed"
    )

    expect_false(chart$fallback.used)
    expect_lte(projector_error(chart$basis, oracle),
               projector_error(pca$basis, oracle) + 1e-12)
    expect_equal(abs(chart$curvature.coefficients[1, 2]), 2,
                 tolerance = 1e-8)
    expect_true(max(abs(chart$corrected.residual[, 2])) < 1e-10)
})

test_that("Harlim second-order chart returns structured failure for too few rows", {
    X <- matrix(c(0, 0, 0), nrow = 1)

    chart <- rcpp_local_second_order_svd_chart(
        X_support = X,
        center = c(0, 0, 0),
        chart_dim = 2L,
        center_mode = "anchor"
    )

    expect_true(chart$fallback.used)
    expect_equal(chart$fallback.reason, "plain_pca_fallback_not_feasible")
    expect_equal(chart$primary.failure.reason, "too_few_effective_support")
    expect_false(chart$curvature.diagnostics$plain.pca.fallback.feasible)
    expect_true(all(is.na(chart$coordinates)))
    expect_equal(dim(chart$coordinates), c(1L, 2L))
    expect_equal(dim(chart$basis), c(3L, 2L))
})

test_that("Harlim second-order chart returns structured failure for all-zero effective weights", {
    grid <- expand.grid(x = -1:1, y = -1:1)
    X <- cbind(grid$x, grid$y, 0)

    chart <- rcpp_local_second_order_svd_chart(
        X_support = X,
        center = c(0, 0, 0),
        chart_dim = 2L,
        center_mode = "anchor",
        weights = rep(0, nrow(X))
    )

    expect_true(chart$fallback.used)
    expect_equal(chart$fallback.reason, "plain_pca_fallback_not_feasible")
    expect_equal(chart$primary.failure.reason, "too_few_effective_support")
    expect_false(chart$curvature.diagnostics$plain.pca.fallback.feasible)
    expect_true(all(is.na(chart$basis)))
})
