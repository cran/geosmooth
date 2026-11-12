test_that("a fixed projection reproduces both native families, offsets and ties", {
    x <- rep(seq(-2, 3, length.out = 12), each = 2)
    X <- cbind(x = x, constant = 1); w <- rep(c(.7, 1.3), 12)
    off <- .1*cos(seq_along(x))
    for (family in c("gaussian", "binomial")) for (backend in c("dense", "banded")) {
        y <- if (family == "gaussian") sin(x)+off else .3+.35*sin(x)^2
        for (rho in c(0, .03)) {
            f <- fit.projection.trend.filtering(X, y, .0002, family = family,
                directions = matrix(c(1, 0), 1), direction.search = "fixed",
                weights = w, offset = off, curve.ridge = rho, backend = backend)
            native <- if (family == "gaussian") fit.gaussian.trend.filtering else
                fit.logistic.trend.filtering
            reference <- native(x, y, .0002*sum(w), weights = w, offset = off,
                                curve.ridge = rho*sum(w), backend = backend)
            expect_true(f$converged)
            expect_equal(as.numeric(fitted(f)), as.numeric(fitted(reference)), tolerance = 1e-6)
            expect_equal(f$checks$objective, reference$diagnostics$objective/sum(w), tolerance = 1e-8)
            expect_equal(colSums(f$terms*w)/sum(w), 0, tolerance = 1e-12)
            q <- c(-4, .2, 5); Q <- cbind(x = q, constant = 1)
            expect_equal(as.numeric(predict(f, Q, offset = rep(.1, 3))),
                as.numeric(predict(reference, q, offset = rep(.1, 3))), tolerance = 1e-6)
            expect_identical(as.vector(attr(predict(f, Q, offset = rep(.1, 3)), "outside.range")),
                             c(TRUE, FALSE, TRUE))
            expect_error(predict(f, Q), "Supply offset")
        }
    }
})

test_that("joint quadratic fits agree with independent penalized regression", {
    X <- as.matrix(expand.grid(a = seq(-1, 1, length.out = 5), b = seq(-.8, 1.2, length.out = 5)))
    w <- seq(.8, 1.2, length.out = nrow(X)); wn <- w/sum(w); rho <- .15
    off <- .05*sin(seq_len(nrow(X)))
    blocks <- lapply(1:2, function(j) {
        B <- cbind(X[, j], X[, j]^2)
        sweep(B, 2, colSums(B*wn))
    })
    A <- cbind(1, blocks[[1]], blocks[[2]])
    Q <- matrix(0, 5, 5)
    Q[2:3, 2:3] <- rho*crossprod(blocks[[1]], blocks[[1]]*wn)
    Q[4:5, 4:5] <- rho*crossprod(blocks[[2]], blocks[[2]]*wn)
    for (family in c("gaussian", "binomial")) {
        y <- .4 + .5*X[, 1]^2 - .6*X[, 2] + .1*sin(seq_len(nrow(X)))
        if (family == "binomial") y <- plogis(y)
        loss <- function(beta) {
            eta <- as.vector(A %*% beta)+off
            value <- if (family == "gaussian") (y-eta)^2/2 else
                pmax(eta, 0)-y*eta+log1p(exp(-abs(eta)))
            sum(wn*value) + sum(beta*(Q %*% beta))/2
        }
        grad <- function(beta) {
            eta <- as.vector(A %*% beta)+off
            r <- if (family == "gaussian") eta-y else plogis(eta)-y
            as.vector(crossprod(A, wn*r)+Q %*% beta)
        }
        ref <- optim(numeric(5), loss, grad, method = "BFGS",
                     control = list(reltol = 1e-14, maxit = 2000))
        f <- fit.projection.trend.filtering(X, y, Inf, family = family,
            direction.search = "fixed", weights = w, offset = off, curve.ridge = rho,
            control = list(tolerance = 1e-9))
        expect_true(f$converged)
        expect_equal(f$checks$objective, ref$value, tolerance = 1e-8)
        expect_equal(as.numeric(predict(f, type = "link")),
                     as.vector(A %*% ref$par)+off, tolerance = 1e-6)
        expect_lte(max(diff(f$objective.history)), 1e-10)
    }
})

ptf_test_difference <- function(x) {
    n <- length(x); D <- diff(diag(n))
    for (k in 1:2) D <- diff(diag(k/(x[(k+1):n]-x[1:(n-k)])) %*% D)
    D
}

test_that("finite-penalty joint fits satisfy independently constructed KKT equations", {
    X <- cbind(c(0,.12,.27,.43,.57,.72,.87,1), c(.43,1,.12,.72,0,.87,.27,.57))
    w <- seq(.8, 1.2, length.out = 8); wn <- w/sum(w)
    lambda <- .003; rho <- .2
    for (family in c("gaussian", "binomial")) {
        y <- sin(4*X[, 1])+.5*cos(5*X[, 2])
        if (family == "binomial") y <- plogis(y)
        f <- fit.projection.trend.filtering(X, y, lambda, family = family,
            direction.search = "fixed", weights = w, curve.ridge = rho,
            control = list(max.sweeps = 600, tolerance = 1e-9))
        expect_true(f$converged)
        eta <- as.numeric(predict(f, type = "link"))
        residual <- wn*(if (family == "gaussian") eta-y else plogis(eta)-y)
        roughness <- 0
        for (j in 1:2) {
            ix <- order(X[, j]); D <- ptf_test_difference(X[ix, j])
            h <- f$terms[ix, j]; g <- (residual+rho*wn*f$terms[, j])[ix]
            dual <- as.vector(solve(tcrossprod(D), -D %*% g)/lambda)
            differences <- as.vector(D %*% h)
            expect_lte(max(abs(g+lambda*crossprod(D, dual))), 1e-7)
            expect_lte(max(abs(dual)), 1+1e-5)
            expect_lte(abs(lambda*sum(abs(differences)-dual*differences)), 1e-7)
            roughness <- roughness + sum(abs(differences))
        }
        loss <- if (family == "gaussian") (y-eta)^2/2 else pmax(eta,0)-y*eta+log1p(exp(-abs(eta)))
        expected <- sum(wn*loss)+lambda*roughness+rho/2*sum(f$terms^2*wn)
        expect_equal(f$checks$objective, expected, tolerance = 1e-8)
        p <- geosmooth:::.ptf.problem(X,y,w,NULL,family,"scale",lambda,rho,
                                    f$control,"dense")
        bad <- f; bad$terms[, 1] <- bad$terms[, 1]+seq(-.1,.1,length.out=8)
        expect_false(geosmooth:::.ptf.diagnose(p,bad)$converged)
        bad <- f; bad$components[[1]]$fit$fits[[1]]$dual[4] <- 2
        expect_gt(geosmooth:::.ptf.diagnose(p,bad)$dual.feasibility, 0)
    }
})

test_that("weight normalization, sign normalization and serialization are stable", {
    X <- cbind(a=seq(-1,1,length.out=20), b=sin(seq_len(20)))
    y <- cos(2*X[,1]); w <- seq(.5,2,length.out=20)
    args <- list(X=X,y=y,lambda=.002,directions=matrix(c(1,.2),1),
                 direction.search="fixed",curve.ridge=.02)
    f <- do.call(fit.projection.trend.filtering,c(args,list(weights=w)))
    g <- do.call(fit.projection.trend.filtering,c(args,list(weights=7*w)))
    expect_equal(as.numeric(fitted(f)),as.numeric(fitted(g)),tolerance=1e-7)
    expect_equal(f$checks$objective,g$checks$objective,tolerance=1e-9)
    args$directions <- -3*args$directions
    h <- do.call(fit.projection.trend.filtering,c(args,list(weights=w)))
    expect_equal(f$directions,h$directions)
    expect_equal(as.numeric(fitted(f)),as.numeric(fitted(h)),tolerance=1e-8)
    restored <- unserialize(serialize(f,NULL))
    expect_equal(predict(f,X[,2:1]),predict(restored,X))
    terms <- predict(f,type="terms")
    expect_equal(as.numeric(predict(f,type="link")),rowSums(terms)+attr(terms,"intercept"))
    expect_output(print(f),"Order-2 projection")
    expect_equal(summary(f)$nterms,1)
    expect_equal(as.numeric(residuals(f)),y-as.numeric(fitted(f)))
})

test_that("constant models and numerical failures remain explicit", {
    X <- matrix(seq(-1,1,length.out=20),ncol=1)
    y <- sin(4*X[,1])
    null <- fit.projection.trend.filtering(X,y,Inf,direction.search="fixed",nterms=0)
    expect_true(null$converged)
    expect_equal(as.numeric(fitted(null)),rep(mean(y),20))
    quadratic <- fit.projection.trend.filtering(X,y,Inf,direction.search="fixed",curve.ridge=0)
    expect_gt(sd(as.numeric(fitted(quadratic))),.01)
    f <- fit.projection.trend.filtering(X,y,.001,direction.search="fixed",
        control=list(inner=list(max.iterations=1)),warn=FALSE)
    expect_false(f$converged)
    expect_match(f$status,"component_solver")
    expect_true(all(is.na(fitted(f))))
    expect_error(fit.projection.trend.filtering(X,rep(0,20),.1,family="binomial"),"finite unpenalized")
    separated <- fit.projection.trend.filtering(X,as.numeric(X[,1]>0),.01,
        family="binomial",direction.search="fixed",curve.ridge=0,warn=FALSE)
    expect_false(separated$converged)
    stable <- fit.projection.trend.filtering(X,as.numeric(X[,1]>0),.01,
        family="binomial",direction.search="fixed",curve.ridge=.05)
    expect_true(stable$converged)
})

test_that("invalid geometry and controls are rejected before fitting", {
    X <- cbind(a=seq_len(12),b=rep(1:3,4),constant=1); y <- sin(seq_len(12))
    fit <- function(...) fit.projection.trend.filtering(X,y,.01,direction.search="fixed",...)
    expect_error(fit(directions=matrix(c(0,1,0),1)),"four distinct")
    expect_error(fit(directions=matrix(c(0,0,1),1)),"constant predictor")
    expect_error(fit(directions=rbind(c(1,0,0),c(-2,0,0))),"Duplicate")
    expect_error(fit(directions=matrix(c(1,0,0),1),nterms=2),"nterms")
    expect_error(fit(weights=rep(0,12)),"positive")
    expect_error(fit(control=list(typo=1)),"Unknown")
    expect_error(fit(control=list(inner=list(typo=1)),nterms=0),"inner control")
})

test_that("empty column names from cbind are aligned by position safely", {
    x <- seq(-1,1,length.out=20); X <- cbind(x,cos(3*x))
    expect_identical(colnames(X),c("x",""))
    f <- fit.projection.trend.filtering(X,sin(x),.001,
        directions=matrix(c(1,0),1),direction.search="fixed")
    expect_equal(predict(f,X),predict(f))
    expect_equal(predict(f,X[,2:1,drop=FALSE]),predict(f))
    expect_equal(as.numeric(predict(f,unname(X))),as.numeric(predict(f)))
    expect_length(predict(f,X[1:3,,drop=FALSE]),3)
})

test_that("unrepresentable normalized spacings remain explicit failures", {
    X <- matrix(c(-1,0,1e-20,.3,1),ncol=1)
    f <- fit.projection.trend.filtering(X,c(0,1,0,1,0),.01,
        direction.search="fixed",preprocessing="none",warn=FALSE)
    expect_false(f$converged)
    expect_match(f$status,"collapsed during normalization")
    expect_true(all(is.na(predict(f))))
})
