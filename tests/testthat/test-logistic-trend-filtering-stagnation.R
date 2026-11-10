# Frozen generators from the practical-validation examples (60 observations,
# seeds 7102 and 7104). No external research files are needed for these tests.
# Saved failures: 2 legacy banded fold candidates, 11 training-mass dense fold
# candidates, and 4 dense full-data candidates. Replay complete warm-start paths.
test_that("binary paths pass after centered-dual and stable-objective repairs", {
    lambda <- c(Inf, 10^seq(-.5, -4, length.out = 12))
    for (shape in c("wave", "hinge")) {
        set.seed(if (shape == "wave") 7102 else 7104)
        x <- seq(0, 1, length.out = 60)
        signal <- if (shape == "wave") sin(2*pi*x) else 4*pmax(x-.45, 0)^2 - 1.2*x
        y <- rbinom(60, 1, plogis(2*signal-.4))
        folds <- sample(rep(1:3, length.out = 60))
        for (fold in 0:3) for (convention in c("training_mass", "none")) {
            if (fold == 0 && convention == "none") next
            train <- if (fold == 0) rep(TRUE, 60) else folds != fold
            factor <- if (convention == "training_mass") mean(train) else 1
            fits <- lapply(c("dense", "banded"), function(backend)
                fit.logistic.trend.filtering(x[train], y[train], lambda*factor,
                    curve.ridge = .02*factor, backend = backend, warn = FALSE))
            label <- paste(shape, fold, convention)
            for (f in fits) {
                expect_true(all(f$diagnostics$converged), info = label)
                expect_lt(max(f$diagnostics$iterations), 200)
                expect_lte(max(f$diagnostics$stationarity), 1e-9)
                expect_lte(max(f$diagnostics$complementarity), 1e-10)
                # Independently reconstruct the original objective and full
                # standardized KKT equations, including all inactive knots.
                d <- f$design; u <- (x[train]-d$lo)/d$span; n <- sum(train)
                H <- sweep(sweep(geosmooth:::.logistic.tf.basis(u, d$knots),
                    2, d$center), 2, d$scales, "/")
                A <- cbind(1, H)
                for (j in seq_along(lambda)) {
                    b <- f$fits[[j]]$coefficients; dual <- f$fits[[j]]$dual
                    eta <- as.vector(A %*% b); curve <- as.vector(H %*% b[-1])
                    pen <- if (is.finite(lambda[j]))
                        c(numeric(3), f$lambda[j]/n/d$scales[-(1:2)]) else numeric(ncol(A))
                    gradient <- as.vector(crossprod(A, (plogis(eta)-y[train])/n))
                    gradient[-1] <- gradient[-1] + f$curve.ridge/n * as.vector(crossprod(H, curve/n))
                    residual <- gradient + pen*dual
                    if (is.infinite(lambda[j])) residual <- residual[1:3]
                    expect_lt(max(abs(residual)), 2e-9)
                    objective <- sum(pmax(eta, 0)-y[train]*eta+log1p(exp(-abs(eta)))) +
                        f$curve.ridge*mean(curve^2)/2 + n*sum(pen*abs(b))
                    expect_equal(f$fits[[j]]$objective, objective, tolerance = 1e-9)
                    expect_lte(max(abs(dual)), 1)
                    expect_lte(sum(pen*(abs(b)-dual*b)), 1.01e-10)
                }
            }
            a <- fits[[1]]; b <- fits[[2]]
            expect_lt(max(abs(a$diagnostics$objective-b$diagnostics$objective)), 1e-7)
            expect_lt(max(abs(predict(a, lambda.index=seq_along(lambda))-
                predict(b, lambda.index=seq_along(lambda)))), 1e-5)
            if (convention == "none" && fold == 2) {
                # Previously a feasible inactive knot was added and removed
                # about 480 times; a larger iteration budget did not help.
                j <- if (shape == "wave") 10 else 13
                expect_lt(b$fits[[j]]$working.set.additions, 30)
                expect_lt(b$fits[[j]]$working.set.removals, 30)
            }
        }
    }
})

test_that("the shared dense solver accepts a rectangular ridge residual factor", {
    x <- seq(-1, 1, length.out = 21)
    A <- cbind(1, x, x^2)
    y <- sin(x) + .3*cos(4*x); w <- seq(.5, 2, length.out = length(x))
    off <- .1*x
    ridge.factor <- rbind(c(0, 2, -1), c(0, 0, .7))
    f <- geosmooth:::gaussian_tf_native(A, ridge.factor, y, w, off, numeric(3),
        0, 1e-9, 1e-10, 1000L)[[1]]
    expected <- solve(crossprod(A, A*w)+crossprod(ridge.factor), crossprod(A, w*(y-off)))
    expect_true(f$converged)
    expect_equal(f$coefficients, as.vector(expected), tolerance = 1e-9)
    objective <- sum(w*(A %*% expected+off-y)^2)/2 + sum((ridge.factor %*% expected)^2)/2
    expect_equal(f$objective, as.numeric(objective), tolerance = 1e-10)
})
