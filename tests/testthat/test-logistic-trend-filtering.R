test_that("quadratic limit matches an independent weighted logistic GLM", {
    x <- seq(-2, 3, length.out = 35)
    y <- plogis(-.4 + .3*x - .2*x^2 + .15*sin(3*x))
    w <- 1 + seq_along(x) %% 4
    z <- matrix(cos(2*x), ncol = 1L)
    off <- sin(x)/10
    f <- fit.logistic.trend.filtering(x, y, Inf, weights = w,
        covariates = z, offset = off)
    g <- suppressWarnings(glm(y ~ x + I(x^2) + z, weights = w,
        offset = off, family = binomial()))
    expect_true(f$diagnostics$converged)
    expect_equal(as.numeric(fitted(f)), as.numeric(fitted(g)), tolerance = 1e-7)
    expect_equal(as.numeric(predict(f, type = "link")), as.numeric(predict(g)), tolerance = 1e-6)
})

test_that("trial weights, repeated binary records, ordering and exact ties agree", {
    x <- seq(0, 1, length.out = 14)
    trials <- rep(c(5, 7), 7); successes <- rep(c(2, 4, 3, 1, 5, 2, 3), 2)
    xx <- rep(x, trials)
    yy <- unlist(Map(function(s, n) c(rep(1, s), rep(0, n-s)), successes, trials))
    a <- fit.logistic.trend.filtering(x, successes/trials, c(1, .01), weights = trials)
    b <- fit.logistic.trend.filtering(xx, yy, c(1, .01))
    expect_true(all(a$diagnostics$converged)); expect_true(all(b$diagnostics$converged))
    expect_equal(as.numeric(predict(a, xx)), as.numeric(fitted(b)), tolerance = 2e-6)
    expect_equal(a$diagnostics$objective, b$diagnostics$objective, tolerance = 1e-7)
    ix <- rev(seq_along(xx))
    c <- fit.logistic.trend.filtering(xx[ix], yy[ix], c(.01, 1))
    expect_equal(as.numeric(fitted(b)), as.numeric(predict(c, xx, lambda.index = c(2, 1))), tolerance = 2e-6)
})

test_that("reported KKT checks and objective can be reconstructed independently", {
    x <- seq(0, 1, length.out = 22)
    y <- .2 + .6*(sin(5*x)^2)
    z <- matrix(cos(9*x), ncol = 1)
    w <- rep(c(1, 3), 11)
    m <- fit.logistic.trend.filtering(x, y, .003, weights = w, covariates = z,
        nuisance.ridge = .16, curve.ridge = .04)
    f <- m$fits[[1]]; d <- m$design
    H <- sweep(sweep(geosmooth:::.logistic.tf.basis(x, d$knots), 2, d$center), 2, d$scales, "/")
    Z <- sweep(sweep(z, 2, d$zcenter), 2, d$zscale, "/")
    A <- cbind(1, Z, H); b <- f$coefficients
    eta <- as.vector(A %*% b); curve <- as.vector(H %*% b[-c(1, 2)])
    penalty <- c(numeric(4), 1/d$scales[-c(1, 2)])
    objective <- sum(w*(pmax(eta, 0) - y*eta + log1p(exp(-abs(eta))))) +
        .003*sum(penalty*abs(b)) + .16/2*b[2]^2 + .04/(2*sum(w))*sum(w*curve^2)
    gradient <- as.vector(crossprod(A, w*(plogis(eta)-y)))
    gradient[2] <- gradient[2] + .16*b[2]
    gradient[-c(1, 2)] <- gradient[-c(1, 2)] + .04/sum(w)*as.vector(crossprod(H, w*curve))
    gradient <- (gradient + .003*penalty*f$dual)/sum(w)
    complementarity <- .003*sum(penalty*(abs(b)-f$dual*b))/sum(w)
    expect_true(f$converged)
    expect_equal(f$objective, objective, tolerance = 1e-9)
    expect_lt(abs(f$stationarity - max(abs(gradient))), 1e-12)
    expect_lt(abs(f$complementarity - complementarity), 1e-12)
    expect_lte(f$stationarity, m$control$tolerance)
    expect_lte(f$complementarity, m$control$gap.tolerance)
    expect_true(all(abs(f$dual) <= 1))
})

test_that("failed fits stay explicit and cannot silently supply predictions", {
    x <- seq(0, 1, length.out = 20); y <- .2 + .6*sin(5*x)^2
    m <- fit.logistic.trend.filtering(x, y, .01,
        control = list(max.iterations = 1L), warn = FALSE)
    expect_false(m$diagnostics$converged)
    expect_identical(m$diagnostics$status, "iteration_limit")
    expect_true(all(is.na(fitted(m))))
    expect_true(all(is.finite(predict(m, allow.unconverged = TRUE))))
    expect_error(fit.logistic.trend.filtering(x, rep(0, 20), 1), "finite unpenalized intercept")
    separated <- fit.logistic.trend.filtering(x, as.numeric(x > .5), 1, warn = FALSE)
    expect_false(separated$diagnostics$converged)
    expect_identical(separated$diagnostics$status, "separation_detected")
    # An explicit curve ridge removes that unpenalized recession direction.
    stable <- fit.logistic.trend.filtering(x, as.numeric(x > .5), 1, curve.ridge = .1)
    expect_true(stable$diagnostics$converged)
})

test_that("prediction uses saved transformations and flags boundary values", {
    x <- seq(-2, 2, length.out = 20); y <- .2 + .6*sin(4*x)^2
    z <- matrix(cos(x), ncol = 1)
    m <- fit.logistic.trend.filtering(x, y, c(.1, .01), covariates = z,
        offset = x/10, nuisance.ridge = .1)
    expect_error(predict(m, c(-3, 3), covariates = matrix(0, 2, 1)), "Supply offset")
    outside <- predict(m, c(-3, 3), covariates = matrix(0, 2, 1), offset = c(0, 0))
    edge <- predict(m, range(x), covariates = matrix(0, 2, 1), offset = c(0, 0))
    expect_equal(as.numeric(outside), as.numeric(edge), tolerance = 1e-12)
    expect_identical(attr(outside, "outside.range"), c(TRUE, TRUE))
    expect_equal(as.numeric(fitted(m)), as.numeric(predict(m, x, covariates = z, offset = x/10)))
})

test_that("CV uses independent training transformations and rejects failed folds", {
    x <- seq(-1, 1, length.out = 24); y <- .3 + .4*cos(3*x)^2
    folds <- rep(1:3, 8); lambda <- c(1, .02)
    m <- fit.logistic.trend.filtering(x, y, lambda,
        lambda.selection = "cv", foldid = folds)
    expect_true(all(m$cv.table$eligible))
    expect_identical(m$selected.index, which.min(m$cv.table$log.loss))
    for (k in 1:3) {
        take <- folds != k
        f <- fit.logistic.trend.filtering(x[take], y[take], lambda * mean(take))
        eta <- predict(f, x[!take], type = "link")
        loss <- colSums(pmax(eta, 0) - y[!take]*eta + log1p(exp(-abs(eta))))
        expect_equal(unname(m$cv.fold.loss.sums[k, ]), unname(loss), tolerance = 1e-8)
    }
    f <- fit.logistic.trend.filtering(x, y, lambda, lambda.selection = "cv",
        foldid = folds, control = list(max.iterations = 1), warn = FALSE)
    expect_true(all(!f$cv.table$eligible)); expect_null(f$selected.index)
})

test_that("invalid inputs and excessive dense workspace are rejected", {
    x <- seq(0, 1, length.out = 12); y <- rep(c(0, 1), 6)
    expect_error(fit.logistic.trend.filtering(x, y, -1), "lambda")
    expect_error(fit.logistic.trend.filtering(x, y, 1, weights = rep(0, 12)), "positive")
    expect_error(fit.logistic.trend.filtering(x, y, 1, control = list(typo = 3)), "control")
    expect_error(fit.logistic.trend.filtering(x, y, 1,
        control = list(max.workspace.bytes = 1)), "workspace")
    expect_error(fit.logistic.trend.filtering(x, y, 1,
        covariates = matrix(1, 12, 1)), "constant")
    expect_error(fit.logistic.trend.filtering(x, y, 1,
        covariates = matrix(rep(c(-1e200, 1e200), 6), 12, 1)), "Nonfinite covariate")
    expect_error(fit.logistic.trend.filtering(x, y, 1,
        lambda.selection = "cv"), "foldid")
})

test_that("near-coincident locations preserve the operator and weak-penalty fit", {
    x <- sort(c(seq(0, 1, length.out = 16), .413, .413 + 1e-10))
    y <- .25 + .5*sin(6*x)^2
    f <- fit.logistic.trend.filtering(x, y, .001)
    expect_true(f$diagnostics$converged)
    # Reconstruct third divided differences without using the solver basis.
    D <- diff(diag(length(x)))
    for (k in 1:2) D <- diff(diag(k/(x[(k+1):length(x)]-x[1:(length(x)-k)])) %*% D)
    B <- cbind(1, geosmooth:::.logistic.tf.basis(x, x))
    expect_lt(max(abs(D %*% B[, 1:3])), 1e-3)
    expect_equal(unname(D %*% B[, -(1:3)]), diag(length(x)-3), tolerance = 1e-3)
    expect_true(all(is.finite(fitted(f))))
})

test_that("tied predictors retain distinct covariates and offsets", {
    x <- rep(seq(0, 1, length.out = 10), each = 3)
    z <- matrix(rep(c(-1, 0, 1), 10), ncol = 1)
    off <- sin(seq_along(x))/8
    y <- plogis(-.3 + x - .4*x*x + .6*z[, 1] + off)
    f <- fit.logistic.trend.filtering(x, y, Inf, covariates = z, offset = off)
    g <- suppressWarnings(glm(y ~ x + I(x^2) + z, offset = off, family = binomial()))
    expect_true(f$diagnostics$converged)
    expect_equal(as.numeric(fitted(f)), as.numeric(fitted(g)), tolerance = 1e-7)
    expect_gt(diff(range(fitted(f)[1:3])), .1)
})


test_that("large Bernoulli sums do not stall the quadratic initializer", {
    # Three predictor groups have known finite fitted probabilities. Ordinary
    # objective summation previously exhausted 1,000 Newton iterations here.
    x <- rep(c(-1, 0, 1), each = 2000)
    y <- c(rep(0, 1885), rep(1, 115), rep(0, 1700), rep(1, 300),
           rep(0, 1400), rep(1, 600))
    A <- cbind(1, scale(cbind(x, x*x)))
    f <- geosmooth:::logistic_tf_native(A, matrix(0, 3, 3), y,
        rep(1, length(y)), numeric(length(y)), numeric(3), Inf,
        1e-9, 1e-10, 1000L)[[1]]
    expect_true(f$converged)
    expect_lt(f$iterations, 20L)
    expect_lte(f$stationarity, 1e-9)
    expect_equal(as.numeric(tapply(f$fitted.values, x, mean)),
        c(115, 300, 600)/2000, tolerance = 1e-9)
})
