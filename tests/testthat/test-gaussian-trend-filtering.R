gaussian_tf_difference <- function(x) {
    n <- length(x)
    D <- diff(diag(n))
    for (k in 1:2) D <- diff(diag(k/(x[(k+1):n]-x[1:(n-k)])) %*% D)
    D
}

test_that("order-2 limits match independent weighted least squares and tied means", {
    x <- c(.03, .1, .17, .25, .39, .51, .58, .72, .81, .93, 1)
    y <- sin(2*pi*x) + seq_along(x)/30
    w <- seq(0.5, 2, length.out = length(x))
    expected <- stats::lm.wfit(cbind(1,x,x^2), y, w)$fitted.values
    for (backend in c("dense", "banded")) {
        fit <- fit.gaussian.trend.filtering(x, y, c(Inf, 0), weights = w, backend = backend)
        expect_true(all(fit$diagnostics$converged))
        expect_identical(fit$order, 2L)
        expect_equal(as.numeric(fitted(fit, lambda.index = 1)), expected, tolerance = 1e-8)
        expect_equal(as.numeric(fitted(fit, lambda.index = 2)), y, tolerance = 1e-7)
        expect_equal(as.numeric(predict(fit, c(-1,2), lambda.index=1)),
                     expected[c(1,length(x))], tolerance=1e-8)
        expect_identical(attr(predict(fit, c(-1,2)), "outside.range"), c(TRUE,TRUE))
        expect_output(print(fit), "Order-2 Gaussian")
        expect_invisible(print(fit))

        xx <- rep(x, each=2); yy <- rep(y, each=2)+rep(c(-.2,.4),length(x))
        ww <- rep(c(1,3),length(x)); ord <- rev(seq_along(xx))
        tied <- fit.gaussian.trend.filtering(xx[ord],yy[ord],0,weights=ww[ord],backend=backend)
        expect_true(tied$diagnostics$converged)
        expect_equal(as.numeric(fitted(tied)), rep(y+.25,each=2)[ord], tolerance=1e-7)
        expect_equal(as.numeric(residuals(tied)), yy[ord]-rep(y+.25,each=2)[ord], tolerance=1e-7)
    }
})

test_that("penalized curves satisfy an independently constructed divided-difference KKT system", {
    x <- c(0,.05,.13,.21,.38,.44,.57,.66,.8,.91,1)
    y <- sin(2*pi*x)+.2*cos(seq_along(x))
    w <- seq(.7,1.3,length.out=length(x)); lambda <- c(.01,.001)
    D <- gaussian_tf_difference(x)
    fits <- lapply(c("dense","banded"), function(backend)
        fit.gaussian.trend.filtering(x,y,lambda,weights=w,backend=backend))
    for (fit in fits) for (j in seq_along(lambda)) {
        expect_true(fit$diagnostics$converged[j])
        beta <- as.numeric(fitted(fit,lambda.index=j)); differences <- as.vector(D %*% beta)
        dual <- as.vector(solve(tcrossprod(D), D %*% (w*(y-beta)))/lambda[j])
        expect_lte(max(abs(w*(beta-y)+lambda[j]*crossprod(D,dual))), 1e-7)
        expect_lte(max(abs(dual)), 1+1e-5)
        expect_lte(abs(lambda[j]*sum(abs(differences)-dual*differences)), 1e-7)
        objective <- sum(w*(y-beta)^2)/2+lambda[j]*sum(abs(differences))
        expect_equal(fit$diagnostics$objective[j],objective,tolerance=1e-7)
    }
    expect_equal(as.numeric(fitted(fits[[1]])),as.numeric(fitted(fits[[2]])),tolerance=2e-6)
})

test_that("Gaussian native fits agree with an independent generalized-lasso solver", {
    skip_if_not_installed("genlasso")
    x <- c(0,.05,.13,.21,.38,.44,.57,.66,.8,.91,1)
    y <- sin(2*pi*x)+.2*cos(seq_along(x)); lambda <- c(.01,.001)
    reference <- genlasso::trendfilter(y, pos=x, ord=2, minlam=0)
    expected <- stats::coef(reference, lambda=lambda)$beta
    for (backend in c("dense","banded")) {
        fit <- fit.gaussian.trend.filtering(x,y,lambda,backend=backend)
        expect_true(all(fit$diagnostics$converged))
        expect_equal(unname(as.matrix(fitted(fit))), unname(expected),
                     tolerance=2e-6, ignore_attr=TRUE)
    }
})

test_that("adjusted Gaussian quadratic fits reproduce the penalized normal equations", {
    x <- seq(0,1,length.out=30); z <- sin(seq_along(x)); off <- cos(seq_along(x))/5
    y <- sin(2*pi*x)+z+off; w <- seq(.5,2,length.out=length(x)); wn <- w/sum(w)
    Z <- (z-sum(wn*z))/sqrt(sum(wn*(z-sum(wn*z))^2))
    H <- cbind(x,x^2); centered <- sweep(H,2,colSums(H*wn))
    A <- cbind(1,Z,H); Q <- matrix(0,4,4)
    Q[2,2] <- .2; Q[3:4,3:4] <- .3*crossprod(centered,centered*wn)
    coefficients <- solve(crossprod(A,A*w)+Q,crossprod(A,w*(y-off)))
    expected <- as.vector(A %*% coefficients+off)
    for (backend in c("dense","banded")) {
        fit <- fit.gaussian.trend.filtering(x,y,Inf,weights=w,covariates=matrix(z),
            offset=off,nuisance.ridge=.2,curve.ridge=.3,backend=backend)
        expect_true(fit$diagnostics$converged)
        expect_equal(as.numeric(fitted(fit)),expected,tolerance=1e-8)
        expect_error(predict(fit,x,covariates=matrix(z)),"Supply offset")
    }
})

test_that("CV scores use every held-out observation and training-only geometry", {
    x <- seq(-2,3,length.out=32); y <- sin(x)+.1*cos(seq_along(x))
    w <- seq(.5,2,length.out=length(x)); folds <- rep(1:4,each=8)
    fit <- fit.gaussian.trend.filtering(x,y,c(Inf,.01),weights=w,
        lambda.selection="cv",foldid=folds,backend="banded")
    scores <- matrix(NA_real_,4,2)
    for (k in 1:4) {
        train <- folds!=k
        ff <- fit.gaussian.trend.filtering(x[train],y[train],
            fit$lambda * sum(w[train])/sum(w),
            weights=w[train],backend="banded")
        pred <- predict(ff,x[!train])
        scores[k,] <- colSums((y[!train]-pred)^2*w[!train])
    }
    expect_equal(fit$cv.fold.loss.sums,scores)
    expect_equal(fit$cv.table$mse,colSums(scores)/sum(w))
    expect_identical(fit$selected.index,which.min(fit$cv.table$mse))
    expect_length(fitted(fit),length(x))
})

test_that("Gaussian controls preserve failures and reject malformed input", {
    x <- seq(0,1,length.out=20); y <- sin(2*pi*x)
    for (backend in c("dense","banded")) {
        expect_warning(fit <- fit.gaussian.trend.filtering(x,y,.01,backend=backend,
            control=list(max.iterations=1L)),"failed convergence")
        expect_false(fit$diagnostics$converged)
        expect_true(all(is.na(fitted(fit))))
        expect_true(all(is.finite(fitted(fit,allow.unconverged=TRUE))))
    }
    for (value in list(NA_real_,c(1,Inf),"y",matrix(y)))
        expect_error(fit.gaussian.trend.filtering(x,value,.1),"y must")
    expect_error(fit.gaussian.trend.filtering(x,y,-1),"lambda")
    expect_error(fit.gaussian.trend.filtering(rep(1:3,length.out=20),y,.1),"four distinct")
    expect_error(fit.gaussian.trend.filtering(x,y,.1,control=list(typo=1)),"Unknown")
    expect_error(fit.gaussian.trend.filtering(x,y,.1,control=list(max.workspace.bytes=1)),"workspace")
    expect_error(fit.gaussian.trend.filtering(x,y,.1,lambda.selection="cv"),"foldid")
    fit <- fit.gaussian.trend.filtering(x,y,.1)
    expect_error(predict(fit,lambda.index=1.5),"lambda.index")
})

test_that("automatic Gaussian dispatch honors workspace and predictor scaling", {
    x <- seq(0,1,length.out=300); y <- 2-x+x^2+.01*cos(seq_along(x))
    fit <- fit.gaussian.trend.filtering(x,y,Inf,control=list(max.workspace.bytes=1e6))
    expect_identical(fit$backend,"banded")
    expect_true(fit$diagnostics$converged)
    expect_lte(fit$workspace.estimate.bytes,1e6)
    expect_error(fit.gaussian.trend.filtering(x,y,Inf,backend="dense",
        control=list(max.workspace.bytes=1e6)),"workspace")
    shifted <- fit.gaussian.trend.filtering(1e8+3*x,y,Inf,backend="banded")
    expect_equal(as.numeric(fitted(fit)),as.numeric(fitted(shifted)),tolerance=1e-8)
})
