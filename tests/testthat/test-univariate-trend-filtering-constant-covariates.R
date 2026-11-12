test_that("constant covariates are rejected before weighted centering", {
    x <- seq(0, 1, length.out = 6)
    y <- rep(c(0, 1), 3)
    for (fit in list(fit.logistic.trend.filtering, fit.gaussian.trend.filtering)) {
        for (backend in c("dense", "banded")) {
            for (weights in list(rep(1, 6), 1:6, c(1, 3, 7, 11, 17, 23))) {
                for (value in c(0, 1, -3, 1e100)) {
                    z <- cbind(varying = x, constant = rep(value, 6))
                    expect_error(fit(x, y, 1, covariates = z, weights = weights,
                                     backend = backend), "constant columns")
                }
            }
        }
    }
})

test_that("small nonconstant covariates retain their standardized fit", {
    x <- seq(0, 1, length.out = 12)
    y <- rep(c(0, 1), 6)
    z <- matrix(sin(seq_along(x)), ncol = 1)
    for (fit in list(fit.logistic.trend.filtering, fit.gaussian.trend.filtering)) {
        for (backend in c("dense", "banded")) {
            reference <- fit(x, y, 1, covariates = z, nuisance.ridge = 0.1,
                             backend = backend)
            scaled <- fit(x, y, 1, covariates = z * 1e-100, nuisance.ridge = 0.1,
                          backend = backend)
            expect_true(all(reference$diagnostics$converged))
            expect_true(all(scaled$diagnostics$converged))
            expect_equal(fitted(scaled), fitted(reference), tolerance = 1e-8)
        }
    }
})
