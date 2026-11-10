test_that("banded and dense solvers agree for weights, offsets and covariates", {
    x <- seq(-2, 3, length.out = 64)
    z <- cbind(cos(2*x), sin(3*x)); off <- sin(x)/5
    y <- plogis(-.3 + .6*sin(4*x) + .4*z[, 1] + off)
    w <- rep(c(1, 2, 4, 3), 16)
    for (ridge in c(0, .04)) {
        a <- fit.logistic.trend.filtering(x, y, c(Inf, .1, .001), weights = w,
            covariates = z, offset = off, nuisance.ridge = .16,
            curve.ridge = ridge, backend = "dense")
        b <- fit.logistic.trend.filtering(x, y, c(Inf, .1, .001), weights = w,
            covariates = z, offset = off, nuisance.ridge = .16,
            curve.ridge = ridge, backend = "banded")
        expect_true(all(b$diagnostics$converged))
        expect_equal(as.numeric(fitted(a)), as.numeric(fitted(b)), tolerance = 1e-6)
        expect_lt(max(abs(a$diagnostics$objective-b$diagnostics$objective))/sum(w), 1e-8)
        xx <- seq(-3, 4, length.out = 99); zz <- cbind(cos(2*xx), sin(3*xx))
        expect_equal(as.numeric(predict(a, xx, covariates = zz, offset = sin(xx)/5)),
            as.numeric(predict(b, xx, covariates = zz, offset = sin(xx)/5)), tolerance = 2e-6)
        expect_lte(max(b$diagnostics$stationarity), b$control$tolerance)
        expect_lte(max(b$diagnostics$complementarity), b$control$gap.tolerance)
    }
})

test_that("banded fits retain exact ties and preserve arbitrary lambda order", {
    x <- rep(seq(0, 1, length.out = 15), c(rep(4, 14), 5))
    y <- rep(c(0, 1, 0, 1, 1), length.out = length(x))
    a <- fit.logistic.trend.filtering(x, y, c(1, .01, .001), backend = "dense")
    ix <- rev(seq_along(x))
    b <- fit.logistic.trend.filtering(x[ix], y[ix], c(.001, 1, .01), backend = "banded")
    expect_true(all(b$diagnostics$converged))
    expect_equal(as.numeric(fitted(a)), as.numeric(predict(b, x, lambda.index = c(2, 3, 1))), tolerance = 2e-6)
})

test_that("banded geometry and coefficients reconstruct the full objective", {
    x <- sort(c(seq(0, 1, length.out = 28), .413, .413+1e-10))
    y <- .25 + .5*sin(6*x)^2
    a <- fit.logistic.trend.filtering(x, y, .001, backend = "dense")
    b <- fit.logistic.trend.filtering(x, y, .001, backend = "banded")
    expect_true(b$diagnostics$converged)
    expect_equal(as.numeric(fitted(a)), as.numeric(fitted(b)), tolerance = 2e-6)
    expect_equal(b$design$center, a$design$center, tolerance = 1e-9)
    expect_equal(unname(b$design$scales), unname(a$design$scales), tolerance = 1e-9)
    f <- b$fits[[1]]; d <- b$design
    H <- sweep(sweep(geosmooth:::.logistic.tf.basis(x, d$knots), 2, d$center), 2, d$scales, "/")
    A <- cbind(1, H); eta <- as.vector(A %*% f$coefficients)
    expect_equal(eta, as.numeric(predict(b, type = "link")), tolerance = 1e-7)
    grad <- as.vector(crossprod(A, plogis(eta)-y))/length(x)
    penalty <- c(numeric(3), .001/length(x)/d$scales[-c(1, 2)])
    expect_lte(max(abs(grad+penalty*f$dual)), 2e-9)
    objective <- sum(pmax(eta, 0)-y*eta+log1p(exp(-abs(eta)))) +
        .001*sum(abs(f$coefficients[-(1:3)])/d$scales[-c(1, 2)])
    expect_lt(abs(f$objective-objective), 1e-8)
    expect_true(all(abs(f$dual)<=1))
})

test_that("auto chooses structured fitting before allocating a dense basis", {
    x <- seq(0, 1, length.out = 400); y <- .2+.6*sin(3*x)^2
    f <- fit.logistic.trend.filtering(x, y, .1,
        control = list(max.workspace.bytes = 1e6))
    expect_identical(f$backend, "banded")
    expect_true(f$diagnostics$converged)
    expect_lt(f$workspace.estimate.bytes, 1e6)
    expect_error(fit.logistic.trend.filtering(x, y, .1, backend = "dense",
        control = list(max.workspace.bytes = 1e6)), "workspace")
    small <- fit.logistic.trend.filtering(x[1:40], y[1:40], .1,
        control = list(max.workspace.bytes = 1e5))
    expect_identical(small$backend, "banded")
    expect_true(small$diagnostics$converged)
})

test_that("banded CV isolates training and failures remain explicit", {
    x <- seq(0, 1, length.out = 36); y <- .3+.4*sin(4*x)^2
    folds <- rep(1:3, 12); lambda <- c(.1, .001)
    f <- fit.logistic.trend.filtering(x, y, lambda, backend = "banded",
        lambda.selection = "cv", foldid = folds)
    expect_true(all(f$cv.table$eligible))
    for (k in 1:3) {
        tr <- folds != k
        b <- fit.logistic.trend.filtering(x[tr], y[tr], lambda * mean(tr), backend = "banded")
        eta <- predict(b, x[!tr], type = "link")
        loss <- colSums(pmax(eta, 0)-y[!tr]*eta+log1p(exp(-abs(eta))))
        expect_equal(as.numeric(loss), as.numeric(f$cv.fold.loss.sums[k, ]), tolerance = 1e-8)
    }
    fail <- fit.logistic.trend.filtering(x, y, .001, backend = "banded",
        control = list(max.iterations = 1), warn = FALSE)
    expect_false(fail$diagnostics$converged)
    expect_true(all(is.na(predict(fail))))
    separated <- fit.logistic.trend.filtering(x, as.numeric(x>.5), .1,
        backend = "banded", warn = FALSE)
    expect_false(separated$diagnostics$converged)
    expect_identical(separated$diagnostics$status, "separation_detected")
    expect_true(any(separated$fits[[1]]$separation.direction != 0))
    saturated <- fit.logistic.trend.filtering(x, rep(c(0, 1), 18), 0,
        backend = "banded", warn = FALSE)
    expect_false(saturated$diagnostics$converged)
    expect_identical(saturated$diagnostics$status, "separation_detected")
    expect_true(any(saturated$fits[[1]]$separation.curve.direction != 0))
})

test_that("tied binary data do not stall when likelihood changes reach roundoff", {
    set.seed(35009)
    n <- 4000L; x <- round(sort(runif(n))*99)/99
    z <- cbind(rnorm(n), rbinom(n, 1, .4)); offset <- -.15*cos(3*pi*x)
    truth <- plogis(-1+.8*sin(6*pi*x)+.25*z[, 1]-.4*z[, 2]+offset)
    y <- rbinom(n, 1, truth)
    m <- fit.logistic.trend.filtering(x, y, c(1, .1, .01, .001),
        covariates = z, offset = offset, nuisance.ridge = .16,
        curve.ridge = .04, backend = "banded")
    expect_true(all(m$diagnostics$converged))
    expect_lte(max(m$diagnostics$stationarity), m$control$tolerance)
    expect_true(all(m$diagnostics$iterations < m$control$max.iterations))
    expect_true(all(is.finite(predict(m))))
    # Independently reconstruct the global, not just working-set, gradient.
    f <- m$fits[[4]]; d <- m$design
    Z <- sweep(sweep(z, 2, d$zcenter), 2, d$zscale, "/")
    H <- sweep(sweep(geosmooth:::.logistic.tf.basis(x, d$knots), 2, d$center), 2, d$scales, "/")
    A <- cbind(1, Z, H); eta <- as.vector(A %*% f$coefficients + offset)
    curve <- as.vector(H %*% tail(f$coefficients, ncol(H)))
    grad <- as.vector(crossprod(A, (plogis(eta)-y)/n))
    grad[2:3] <- grad[2:3] + .16/n*f$coefficients[2:3]
    grad[-(1:3)] <- grad[-(1:3)] + .04/n^2*as.vector(crossprod(H, curve))
    penalty <- c(numeric(5), .001/n/d$scales[-c(1, 2)])
    expect_lt(max(abs(grad+penalty*f$dual)), 1e-8)
})


test_that("clustered terminal knots use a stable compact projection", {
    # The old normal-equation projection threw before returning this path.
    # These are synthetic fractional outcomes; no private AGP data are used.
    x <- c(seq(0, .8, length.out = 24), .9, .90001, .90002, 1)
    y <- .2 + .5*sin(4*x)^2
    lambda <- c(Inf, .001, 0)
    dense <- fit.logistic.trend.filtering(x, y, lambda, backend = "dense")
    banded <- fit.logistic.trend.filtering(x, y, lambda, backend = "banded")
    expect_true(all(banded$diagnostics$converged))
    expect_lte(max(banded$diagnostics$stationarity), banded$control$tolerance)
    expect_lte(max(banded$diagnostics$complementarity), banded$control$gap.tolerance)
    expect_equal(as.numeric(fitted(banded)), as.numeric(fitted(dense)), tolerance = 2e-6)
    expect_lt(max(abs(banded$diagnostics$objective-dense$diagnostics$objective)), 1e-7)
    # Exercise nonuniform weights, row order and exact predictor ties as well.
    ix <- rev(rep(seq_along(x), 2))
    weighted <- fit.logistic.trend.filtering(x[ix], y[ix], lambda,
        weights = rep(c(1, 3), each = length(x)), backend = "banded")
    reference <- fit.logistic.trend.filtering(x[ix], y[ix], lambda,
        weights = rep(c(1, 3), each = length(x)), backend = "dense")
    expect_true(all(weighted$diagnostics$converged))
    expect_equal(as.numeric(fitted(weighted)), as.numeric(fitted(reference)), tolerance = 2e-6)
})
