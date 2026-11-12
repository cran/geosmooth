test_that("direction refinement improves a rotated ridge for both losses", {
    X <- as.matrix(expand.grid(a=seq(-1,1,length.out=7),b=seq(-1,1,length.out=7)))
    truth <- 2*sin(2*(.6*X[,1]+.8*X[,2]))
    for (family in c("gaussian","binomial")) {
        y <- if (family=="gaussian") truth else plogis(truth)
        control <- list(starts=1,sparse.directions=0,dense.directions=0,
                        refinement.passes=3,refinement.proposals=24)
        base.control <- control; base.control$refinement.passes <- 0
        base <- fit.projection.trend.filtering(X,y,.0001,family=family,nterms=1,
            curve.ridge=.01,preprocessing="isotropic",control=base.control,seed=42)
        fit <- fit.projection.trend.filtering(X,y,.0001,family=family,nterms=1,
            curve.ridge=.01,preprocessing="isotropic",control=control,seed=42)
        expect_true(base$converged); expect_true(fit$converged)
        expect_lt(fit$checks$objective,base$checks$objective-.01)
        expect_lt(mean((as.numeric(fitted(fit))-y)^2),mean((as.numeric(fitted(base))-y)^2)/2)
        expect_lte(max(diff(fit$objective.history)),0)
        expect_gt(sum(fit$search.history$accepted & fit$search.history$stage=="refinement"),0)
        expect_equal(sqrt(rowSums(fit$directions^2)),1,tolerance=1e-12)
        expect_identical(fit$search.status,"component_limit")
    }
})

test_that("learned fits and search traces reproduce without changing the caller RNG", {
    X <- cbind(seq(-1,1,length.out=18),sin(seq_len(18)),cos(seq_len(18)))
    y <- .5+sin(X[,1]+X[,2])
    control <- list(starts=2,sparse.directions=2,dense.directions=2,
                    refinement.passes=1,refinement.proposals=4)
    set.seed(345); rng <- .Random.seed
    f <- fit.projection.trend.filtering(X,y,.002,nterms=1,seed=7,control=control)
    expect_identical(.Random.seed,rng)
    g <- fit.projection.trend.filtering(X,y,.002,nterms=1,seed=7,control=control)
    expect_equal(f$directions,g$directions)
    expect_equal(f$search.history,g$search.history)
    expect_equal(fitted(f),fitted(g))
    expect_equal(nrow(f$search.starts),2)
    expect_gte(f$solver.calls,nrow(f$search.history))
})

test_that("multiple learned terms are jointly certified and null signals remain constant", {
    X <- as.matrix(expand.grid(a=seq(-1,1,length.out=5),b=seq(-1,1,length.out=5)))
    control <- list(starts=1,sparse.directions=0,dense.directions=0,refinement.passes=0)
    y <- X[,1]*X[,2]
    f <- fit.projection.trend.filtering(X,y,Inf,nterms=2,curve.ridge=.05,
        directions=rbind(c(1,1),c(1,-1)),control=control,preprocessing="isotropic")
    expect_true(f$converged)
    expect_equal(f$nterms,2)
    expect_equal(length(unique(f$components[[1]]$fit$x)),9)
    expect_lt(mean((as.numeric(fitted(f))-y)^2),mean(y^2)/10)
    null <- fit.projection.trend.filtering(X,rep(.3,nrow(X)),.01,
        family="binomial",nterms=2,control=control)
    expect_true(null$converged)
    expect_equal(null$nterms,0)
    expect_identical(null$search.status,"no_improving_addition")
    expect_equal(as.numeric(predict(null)),rep(.3,nrow(X)),tolerance=1e-10)
})

test_that("learned binomial fits handle observed zeros and ones", {
    X <- cbind(a=seq(-1,1,length.out=40),b=sin(seq_len(40)))
    y <- as.numeric(sin(5*X[,1]+X[,2])>0)
    f <- fit.projection.trend.filtering(X,y,.001,family="binomial",nterms=1,
        curve.ridge=.05,seed=72,
        control=list(starts=1,sparse.directions=2,dense.directions=2,
                     refinement.passes=1,refinement.proposals=4))
    expect_true(f$converged)
    expect_true(all(is.finite(predict(f,type="link"))))
    expect_true(all(predict(f)>0 & predict(f)<1))
})

test_that("search failures cannot masquerade as a successful intercept-only model", {
    X <- cbind(seq(-1,1,length.out=20),sin(seq_len(20)))
    f <- fit.projection.trend.filtering(X,sin(X[,1]),.001,nterms=1,warn=FALSE,
        control=list(starts=1,sparse.directions=0,dense.directions=0,
                     inner=list(max.iterations=1)))
    expect_false(f$converged)
    expect_identical(f$search.status,"all_candidates_failed")
    expect_true(all(is.na(predict(f))))
    expect_true(any(grepl("component_solver",f$search.history$status)))
})
