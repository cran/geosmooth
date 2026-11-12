ptf_fd_hessian <- function(f, p, h = 1e-4) {
    H <- matrix(0, p, p); zero <- numeric(p)
    for (j in seq_len(p)) for (k in seq_len(p)) {
        ej <- ek <- zero; ej[j] <- h; ek[k] <- h
        H[j, k] <- (f(ej+ek)-f(ej-ek)-f(-ej+ek)+f(-ej-ek))/(4*h*h)
    }
    H
}

test_that("stored TF2 slopes match independent one-sided and interior differences", {
    x <- rep(c(-1.7,-1.2,-.6,-.25,.12,.4,.83,1.7,2.4), each = 2)
    for (family in c("gaussian", "binomial")) for (backend in c("dense", "banded")) {
        y <- sin(2.3*x)+.15*cos(seq_along(x))
        if (family == "binomial") y <- plogis(y)
        fit <- geosmooth:::.fit.univariate.trend.filtering(x, y, .0001*length(x),
            curve.ridge = .01*length(x), family = family, backend = backend)
        expect_true(fit$diagnostics$converged)
        value <- function(z) as.numeric(predict(fit, z, type = "link"))
        knots <- sort(unique(x)); h <- min(diff(knots))*1e-5
        interior <- (head(knots, -1)+tail(knots, -1))/2
        expect_equal(geosmooth:::.ptf.curve.slopes(fit, interior),
            (value(interior+h)-value(interior-h))/(2*h), tolerance = 2e-5)
        t <- knots[-c(1, length(knots))]
        left <- (value(t)-value(t-h))/h
        right <- (value(t+h)-value(t))/h
        expect_equal(geosmooth:::.ptf.curve.slopes(fit, t), left, tolerance = 2e-4)
        expect_equal(geosmooth:::.ptf.curve.slopes(fit, knots[c(1, length(knots))]),
            c((value(knots[1]+h)-value(knots[1]))/h,
              (value(tail(knots, 1))-value(tail(knots, 1)-h))/h), tolerance = 2e-4)
        d <- fit$design; b <- fit$fits[[1]]$coefficients[-1]/d$scales
        jumps <- numeric(length(t))
        for (j in 3:length(b)) jumps[j-1] <- b[j]*(d$knots[j]-d$knots[j-1])/(2*d$span)
        expect_equal(right-left, jumps, tolerance = 4e-4)
        outside <- c(min(x)-2, max(x)+2)
        expect_equal(geosmooth:::.ptf.curve.slopes(fit, outside), c(0, 0))
        expect_equal(value(outside), value(range(x)), tolerance = 1e-10)
        diagnostic <- geosmooth:::.ptf.curve.diagnostics(list(fit = fit))
        expect_gt(diagnostic$active.knots, 0)
        expect_gte(diagnostic$largest.jump, max(abs(right-left))-4e-4)
    }
})

test_that("quadratic-limit slopes are informative and constant slopes are exactly zero", {
    x <- c(-2,-1.3,-.9,-.1,.2,.8,1.1,2)
    for (family in c("gaussian", "binomial")) {
        eta <- .2+.3*x+.15*x^2
        y <- if (family == "gaussian") eta else plogis(eta)
        fit <- geosmooth:::.fit.univariate.trend.filtering(x, y, Inf, family = family,
                                                         backend = "banded")
        expect_equal(geosmooth:::.ptf.curve.slopes(fit), .3+.3*x, tolerance = 1e-6)
        expect_equal(as.numeric(predict(fit, c(-3, 3), type = "link")),
                     eta[c(1, length(eta))], tolerance = 1e-6)
        fit$fits[[1]]$coefficients[-1] <- 0
        expect_identical(geosmooth:::.ptf.curve.slopes(fit), rep(0, length(x)))
    }
})

test_that("family GN curvature matches the linearized loss, with floor checked separately", {
    Z <- cbind(seq(-1.2, 1.5, length.out = 9), cos(seq_len(9)))
    g <- seq(-.7, 1.1, length.out = 9); eta <- sin(seq_len(9))
    w <- seq(1, 3, length.out = 9); w <- w/sum(w)
    for (family in c("gaussian", "binomial")) {
        y <- if (family == "gaussian") cos(seq_len(9)) else (seq_len(9) %% 4)/3
        eq <- geosmooth:::.ptf.gn.equations(Z, y, eta, g, w, family)
        loss <- function(delta) {
            v <- eta+as.vector((Z*g) %*% delta)
            sum(w*if (family == "gaussian") (v-y)^2/2 else
                pmax(v, 0)-y*v+log1p(exp(-abs(v))))
        }
        expect_equal(eq$unfloored, ptf_fd_hessian(loss, 2), tolerance = 1e-7)
        h <- 1e-5
        gradient <- vapply(1:2, function(j) {
            d <- numeric(2); d[j] <- h
            (loss(d)-loss(-d))/(2*h)
        }, numeric(1))
        expect_equal(eq$b, -gradient, tolerance = 1e-8)
        expect_equal(eq$A, eq$unfloored, tolerance = 1e-14)
    }
    eta <- rep(c(-25,25,0), 3); y <- rep(c(1,0,.5), 3)
    eq <- geosmooth:::.ptf.gn.equations(Z, y, eta, g, w, "binomial")
    correction <- crossprod(eq$J, eq$J*(w*(eq$stabilized.variance-eq$variance)))
    expect_lt(max(abs(eq$A-eq$unfloored-correction)), 1e-14)
    expect_equal(eq$b, as.vector(crossprod(eq$J, w*(y-plogis(eta)))), tolerance = 1e-14)
    expect_equal(eq$b, as.vector(crossprod(eq$J, eq$working.weights*eq$working.residual)))
    expect_equal(eq$floored, 6)
    expect_gte(min(eigen(correction, symmetric = TRUE)$values), -1e-14)
    nonlinear <- function(d) ((1+d)^2-2)^2/2
    expect_equal(as.numeric(ptf_fd_hessian(nonlinear, 1)), 2, tolerance = 1e-6)
    expect_equal(geosmooth:::.ptf.gn.equations(matrix(1), 2, 1, 2, 1, "gaussian")$A[1], 4)
})

test_that("Jacobian rows use the actual isotropically preprocessed scores", {
    X <- cbind(seq(-1,1,length.out=12), cos(seq_len(12)))*70+30
    p <- geosmooth:::.ptf.problem(X, sin(X[,1]), NULL, NULL, "gaussian", "isotropic",
        .01, .01, geosmooth:::.ptf.control(list()), "banded")
    u <- geosmooth:::.ptf.canonical(c(1,.3))
    x <- geosmooth:::.ptf.projection(p,u)
    f <- fit.gaussian.trend.filtering(x, sin(x), .001, backend="banded")
    q <- (head(sort(unique(x)),-1)+tail(sort(unique(x)),-1))/2
    Z <- p$Z[seq_along(q),,drop=FALSE]
    g <- geosmooth:::.ptf.curve.slopes(f,q)
    J <- Z*g; h <- 1e-5
    for (j in 1:2) {
        numerical <- (as.numeric(predict(f,q+h*Z[,j]))-as.numeric(predict(f,q-h*Z[,j])))/(2*h)
        expect_equal(J[,j], numerical, tolerance=1e-5)
    }
    expect_gt(p$preprocessing$scale[1], 10)
})

test_that("tangent damping avoids the normalize-after-solving counterexample", {
    u <- c(1,0); mu <- .1; A <- matrix(c(2,1,1,1),2)-diag(mu,2); b <- c(3,1)
    unconstrained <- solve(A+diag(mu,2), b)
    projected <- as.vector((diag(2)-tcrossprod(u)) %*% unconstrained)
    expect_equal(unconstrained,c(2,-1))
    expect_gt(-sum(b*projected),0)
    step <- geosmooth:::.ptf.gn.increment(A,b,u,mu)
    expect_equal(step$delta,c(0,1),tolerance=1e-12)
    expect_lt(-sum(b*step$delta),0)
    expect_equal(sum(u*step$delta),0,tolerance=1e-14)
    expect_equal(sum(geosmooth:::.ptf.retract(u,step$delta)^2),1,tolerance=1e-14)
    zero <- geosmooth:::.ptf.gn.increment(matrix(0,2,2),c(0,0),u)
    expect_equal(zero$delta,c(0,0)); expect_equal(zero$mu,1e-8)
    expect_identical(zero$status,"zero_tangent_score")
    expect_identical(geosmooth:::.ptf.gn.increment(matrix(1),1,1)$status,"rank_one")
})

test_that("PCA tied blocks have deterministic bases and explicit realized ranks", {
    U <- qr.Q(qr(matrix(c(1,2,3,2,-1,1),3,2)))
    R <- matrix(c(cos(.8),sin(.8),-sin(.8),cos(.8)),2)
    basis <- geosmooth:::.ptf.projector.basis
    expect_equal(basis(tcrossprod(U),2),basis(tcrossprod(U%*%R),2),tolerance=1e-12)
    expect_equal(crossprod(basis(tcrossprod(U),2)),diag(2),tolerance=1e-12)
    X <- rbind(diag(3),-diag(3))
    pca <- geosmooth:::.ptf.pca
    f <- pca(X,theta=.5)
    expect_equal(f$rank,3)
    expect_equal(f$basis,diag(3),tolerance=1e-12)
    g <- pca(X,rank=1)
    expect_true(g$rank.expanded); expect_equal(g$rank,3)
    expect_false(pca(X,rank=4)$available)
    expect_equal(pca(matrix(1,6,3))$rank,0)
    expect_equal(pca(X,theta=1)$rank,3)
    expect_error(basis(matrix(0,2,2),1),"recover")
    X <- X %*% diag(sqrt(c(1,1-5e-9,.2)))
    expect_equal(pca(X,theta=.6)$rank,2)
    expect_equal(pca(X,weights=rep(7,6))$basis,pca(X)$basis,tolerance=1e-12)
})

test_that("seed intercepts use binomial curvature weights and SIR has its own availability", {
    X <- cbind(seq(-2,2,length.out=60),sin(seq_len(60)))
    y <- (seq_len(60)%%5)/4; w <- seq(1,4,length.out=60)
    p <- geosmooth:::.ptf.problem(X,y,w,X[,1]/2,"binomial","isotropic",.001,.01,
        geosmooth:::.ptf.control(list()),"banded")
    state <- geosmooth:::.ptf.backfit(p,geosmooth:::.ptf.initial(p,matrix(0,0,2)))
    s <- geosmooth:::.ptf.seed.data(p,state)
    eta <- p$offset+state$intercept; v <- pmax(plogis(eta)*plogis(-eta),1e-6)
    expect_equal(s$w,p$w*v)
    expect_equal(sum(s$w*s$r),0,tolerance=1e-14)
    expect_equal(colSums(s$Z*s$w),c(0,0),tolerance=1e-14)
    Z <- matrix(rep(c(-1,2,0,-2,1),each=15),ncol=1)
    response <- rep(-2:2,each=15); weights <- rep(1/75,75)
    expect_equal(sum(Z[,1]*response*weights),0,tolerance=1e-14)
    sir <- geosmooth:::.ptf.sir(Z,response,weights)
    expect_identical(sir$status,"available"); expect_equal(sir$eigenvalue,1,tolerance=1e-12)
    problem <- geosmooth:::.ptf.problem(Z,response,weights,NULL,"gaussian","isotropic",.01,.01,
        geosmooth:::.ptf.control(list()),"banded")
    null <- geosmooth:::.ptf.backfit(problem,geosmooth:::.ptf.initial(problem,matrix(0,0,1)))
    seeds <- geosmooth:::.ptf.supervised.seeds(problem,null)
    expect_true(seeds$zero.seed)
    expect_null(seeds$directions$pls)
    expect_identical(seeds$sir$status,"available")
    scaled <- geosmooth:::.ptf.sir(Z,response,7*weights)
    expect_equal(scaled$direction,sir$direction); expect_equal(scaled$eigenvalue,sir$eigenvalue)
    expect_identical(geosmooth:::.ptf.sir(Z,rep(0,75),weights)$status,"sir_unavailable")
    r <- rep(0:2,each=20); ww <- c(50,rep(1,59))
    zz <- matrix(seq_len(60),ncol=1); zz <- zz-sum(zz[,1]*ww)/sum(ww)
    tied <- geosmooth:::.ptf.sir(zz,r,ww)
    expect_equal(tied$slices,rep(1:3,each=20))
    small <- geosmooth:::.ptf.sir(matrix(c(-1,1),ncol=1),c(0,1),c(1,1),classes=c(0,1))
    expect_identical(small$status,"available")
    p <- geosmooth:::.ptf.problem(X[1:24,],rep(c(0,1),12),NULL,NULL,"binomial",
        "isotropic",.01,.01,geosmooth:::.ptf.control(list()),"banded")
    state <- geosmooth:::.ptf.backfit(p,geosmooth:::.ptf.initial(p,matrix(0,0,2)))
    first <- geosmooth:::.ptf.supervised.seeds(p,state)
    expect_true(first$first.bernoulli); expect_identical(first$sir$status,"available")
    p$weights[] <- 16
    proportions <- geosmooth:::.ptf.supervised.seeds(p,state)
    expect_false(proportions$first.bernoulli)
    expect_identical(proportions$sir$status,"sir_unavailable")
})

test_that("keyed start banks preserve RNG, ordering, duplicate origins and context", {
    X <- cbind(seq(-1,1,length.out=60),sin(seq_len(60)))
    p <- geosmooth:::.ptf.problem(X,X[,1],NULL,NULL,"gaussian","isotropic",.01,.01,
        geosmooth:::.ptf.control(list()),"banded")
    state <- geosmooth:::.ptf.backfit(p,geosmooth:::.ptf.initial(p,matrix(0,0,2)))
    set.seed(731); previous <- .Random.seed
    a <- geosmooth:::.ptf.gn.starts(p,state,1,83,key=c("dataset","A","fold","2"))
    b <- geosmooth:::.ptf.gn.starts(p,state,1,83,key=c("dataset","A","fold","2"))
    c <- geosmooth:::.ptf.gn.starts(p,state,1,83,key=c("dataset","A","fold","3"))
    expect_identical(.Random.seed,previous); expect_identical(a,b)
    expect_false(identical(a$seed,c$seed))
    expect_lte(nrow(a$directions),9)
    expect_equal(a$directions[1:2,],diag(2),ignore_attr=TRUE)
    expect_identical(a$origins[[1]][1],"axis_1")
    expect_equal(nrow(geosmooth:::.ptf.unique.directions(a$directions)),nrow(a$directions))
    pls <- geosmooth:::.ptf.gn.starts(p,state,1,83,method="pls")
    expect_equal(pls$random.draws,0)
    expect_identical(pls$supervised$sir$status,"not_requested")
    p$y[] <- 1; state$intercept <- 1
    fallback <- geosmooth:::.ptf.gn.starts(p,state,1,83,method="pls")
    expect_equal(nrow(fallback$directions),1)
    expect_identical(fallback$origins[[1]],"random_fallback")
    expect_equal(fallback$random.draws,1)
})
