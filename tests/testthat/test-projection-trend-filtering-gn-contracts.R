test_that("native budgets are enforced inside backfitting and before failed native calls", {
    X <- cbind(seq(-1,1,length.out=12),cos(seq_len(12)))
    for (family in c("gaussian","binomial")) {
        y <- sin(2*X[,1])+X[,2]/3
        if (family=="binomial") y <- plogis(y)
        p <- geosmooth:::.ptf.problem(X,y,NULL,NULL,family,"isotropic",.001,.1,
            geosmooth:::.ptf.control(list()),"banded")
        budget <- geosmooth:::.ptf.native.budget(2,2,"gn",per.start=1)
        p$before.fit <- function() geosmooth:::.ptf.reserve.native(budget)
        p$counter <- budget
        result <- geosmooth:::.ptf.backfit(p,geosmooth:::.ptf.initial(p,diag(2)))
        expect_identical(result$status,"budget_exhausted")
        expect_equal(budget$native,1)
        expect_equal(budget$sweeps,1)
        expect_false(result$converged)
        expect_identical(budget$blocked,"start")
    }
    budget <- geosmooth:::.ptf.native.budget(1,1,"gn",per.start=1)
    called <- 0L
    testthat::local_mocked_bindings(gaussian_tf_banded_native=function(...) {
        called <<- called+1L
        stop("injected native failure")
    },.package="geosmooth")
    expect_error(geosmooth:::.fit.univariate.trend.filtering(seq_len(8),sin(seq_len(8)),.01,
        family="gaussian",backend="banded",before.fit=function()
            geosmooth:::.ptf.reserve.native(budget)),"injected native failure")
    expect_equal(called,1); expect_equal(budget$native,1)
    expect_error(geosmooth:::.fit.univariate.trend.filtering(seq_len(8),sin(seq_len(8)),.01,
        family="gaussian",backend="banded",before.fit=function()
            geosmooth:::.ptf.reserve.native(budget)),"budget_exhausted")
    expect_equal(called,1); expect_equal(budget$native,1)
})

test_that("eligibility policy distinguishes certificates, caps, geometry and numerical failure", {
    decide <- geosmooth:::.ptf.addition.outcome
    expect_true(decide(completed=1,admissible=1)$eligible)
    expect_identical(decide(completed=1,admissible=1)$status,"no_improving_addition")
    expect_identical(decide()$status,"no_admissible_direction")
    expect_false(decide(admissible=1,failures=1)$eligible)
    expect_identical(decide(admissible=1,failures=1)$status,"all_candidates_failed")
    expect_identical(decide(admissible=2,failures=1,interrupted=1)$status,"budget_exhausted")
    expect_identical(decide(completed=1,admissible=2,interrupted=1,improving=TRUE)$status,"accepted")
    expect_true(decide(completed=1,admissible=2,failures=1)$eligible)
    expect_false(decide(completed=1,improving=TRUE,global.exhausted=TRUE)$eligible)
    expect_identical(decide(completed=1,selected.failed=TRUE)$status,"addition_refit_failed")
    b <- geosmooth:::.ptf.native.budget(3,5,"gn")
    expect_equal(b$total.limit,40*3*18+400)
    expect_equal(geosmooth:::.ptf.native.budget(3,5,"pls")$total.limit,40*3+400)
    b$shared.limit <- 1
    geosmooth:::.ptf.budget.enter(b,"shared")
    geosmooth:::.ptf.reserve.native(b)
    expect_error(geosmooth:::.ptf.reserve.native(b),"budget_exhausted")
    expect_true(b$global.exhausted); expect_identical(b$blocked,"shared")
    expect_false(geosmooth:::.ptf.eligible(list(converged=TRUE,cv.eligible=FALSE)))
})

test_that("Gaussian and binomial GN fits use original responses and accept certified improvements", {
    X <- cbind(seq(-1,1,length.out=24),cos(seq_len(24)))
    for (family in c("gaussian","binomial")) {
        eta <- .2+1.1*(.6*X[,1]+.8*X[,2])^2
        y <- if (family=="gaussian") eta else plogis(eta)
        w <- rep(c(4,8,16),8); off <- .1*sin(seq_len(24))
        f <- fit.projection.trend.filtering(X,y,Inf,family=family,nterms=1,
            direction.search="gn",weights=w,offset=off,curve.ridge=.03,
            preprocessing="isotropic",backend="banded",seed=71)
        expect_true(f$converged); expect_true(f$cv.eligible); expect_true(f$numerical.certified)
        expect_equal(f$counters$completed.candidate.fits,1)
        expect_lte(max(diff(f$objective.history)),0)
        expect_lte(max(f$search.starts$native.fits),40)
        expect_lte(f$counters$native.fits,f$counters$candidate.cap)
        expect_equal(f$counters$native.fits,sum(f$search.starts$native.fits))
        expect_equal(sqrt(rowSums(f$directions^2)),rep(1,f$nterms),tolerance=1e-12)
        expect_true(any(f$search.history$status=="accepted" & f$search.history$stage=="refinement"))
        for (component in f$components) {
            expect_equal(component$fit$y,y); expect_equal(component$fit$weights,w)
        }
        expect_true(all(is.finite(predict(f,type="link"))))
        restored <- unserialize(serialize(f,NULL))
        expect_equal(predict(restored),predict(f))
    }
})

test_that("PLS multi-component backfits retain offsets and original binomial proportions", {
    X <- cbind(seq(-1,1,length.out=24),sin(seq_len(24)))
    original <- geosmooth:::.fit.univariate.trend.filtering
    seen <- list()
    testthat::local_mocked_bindings(.fit.univariate.trend.filtering=function(...) {
        args <- list(...); seen[[length(seen)+1L]] <<- args
        original(...)
    },.package="geosmooth")
    for (family in c("gaussian","binomial")) {
        eta <- X[,1]+.7*X[,2]^2; y <- if (family=="gaussian") eta else plogis(eta)
        weights <- rep(c(3,8,16),8); offset <- .1*X[,2]
        seen <- list()
        f <- fit.projection.trend.filtering(X,y,Inf,family=family,weights=weights,
            offset=offset,curve.ridge=.2,nterms=2,direction.search="pls",backend="banded",
            preprocessing="isotropic",control=list(gn.max.steps=0))
        expect_true(f$converged)
        expect_gte(length(seen),1)
        for (args in seen) {
            expect_equal(args[[2]],y)
            expect_equal(args$weights,weights)
            expect_length(args$offset,nrow(X))
        }
        expect_equal(colSums(f$terms*weights)/sum(weights),rep(0,f$nterms),tolerance=1e-12)
        expect_lte(f$counters$native.fits,f$counters$candidate.cap)
    }
})

test_that("cap-limited starts may nominate certificates but global caps cannot", {
    X <- cbind(seq(-1,1,length.out=18),cos(seq_len(18)))
    for (family in c("gaussian","binomial")) {
        y <- (.6*X[,1]+.8*X[,2])^2
        if (family=="binomial") y <- plogis(y)
        f <- fit.projection.trend.filtering(X,y,Inf,family=family,nterms=1,direction.search="gn",
            curve.ridge=.05,preprocessing="isotropic",backend="banded",
            control=list(gn.max.fits=1))
        expect_true(f$cv.eligible)
        expect_true(any(f$search.starts$status=="budget_exhausted"))
        expect_lte(max(f$search.starts$native.fits),1)
        p <- geosmooth:::.ptf.problem(X,y,NULL,NULL,family,"isotropic",Inf,.05,
            geosmooth:::.ptf.control(list()),"banded")
        b <- geosmooth:::.ptf.native.budget(1,2,"gn",total=1)
        stopped <- geosmooth:::.ptf.gn.search(p,1,12,budget=b)
        expect_identical(stopped$search$status,"budget_exhausted")
        expect_false(stopped$state$cv.eligible)
        expect_true(stopped$state$checks$converged)
        expect_equal(stopped$state$counters$native.fits,1)
    }
    failed <- fit.projection.trend.filtering(X,sin(X[,1]),.001,nterms=1,
        direction.search="gn",backend="banded",warn=FALSE,control=list(gn.max.fits=0))
    expect_false(failed$cv.eligible); expect_true(failed$checks$converged)
    expect_equal(failed$counters$candidate.fits,1)
    expect_equal(failed$counters$completed.candidate.fits,0)
    expect_true(all(is.na(predict(failed))))
})

test_that("geometric exclusions, flat curves and rank-one search remain distinct", {
    X <- matrix(rep(1:3,4),ncol=1)
    f <- fit.projection.trend.filtering(X,sin(X[,1]),.001,nterms=1,direction.search="gn")
    expect_identical(f$search.status,"no_admissible_direction")
    expect_true(f$converged); expect_equal(f$counters$native.fits,0)
    X <- matrix(c(-1,0,1e-20,.3,1),ncol=1)
    f <- fit.projection.trend.filtering(X,c(0,1,0,1,0),.01,nterms=1,
        direction.search="gn",preprocessing="none")
    expect_identical(f$search.status,"no_admissible_direction")
    expect_equal(f$counters$native.fits,0)
    X <- matrix(seq(-1,1,length.out=12),ncol=1)
    f <- fit.projection.trend.filtering(X,X[,1]^2,Inf,nterms=1,direction.search="gn")
    expect_identical(f$preprocessing$method,"isotropic")
    expect_true(f$converged); expect_true("rank_one" %in% f$search.history$status)
    expect_equal(f$counters$native.fits,1)
    for (family in c("gaussian","binomial")) {
        f <- fit.projection.trend.filtering(X,rep(.3,12),Inf,family=family,nterms=1,
            direction.search="gn",backend="banded")
        expect_true(f$converged); expect_identical(f$search.status,"no_improving_addition")
        expect_equal(f$nterms,0)
    }
    expect_error(fit.projection.trend.filtering(X,rep(0,12),Inf,
        family="binomial",direction.search="gn"),"finite unpenalized")
})

test_that("weighted anchors preserve mean-loss scaling for Gaussian and binomial fits", {
    X <- cbind(seq(-1,1,length.out=16),cos(seq_len(16)))
    for (family in c("gaussian","binomial")) {
        y <- sin(2*X[,1])+.1*X[,2]
        if (family=="binomial") y <- plogis(y)
        make <- function(w) geosmooth:::.ptf.problem(X,y,w,NULL,family,"isotropic",.001,.03,
            geosmooth:::.ptf.control(list()),"banded")
        w <- rep(c(3,8,16,5),4)
        a <- geosmooth:::.ptf.weighted.anchor(make(w),diag(2))
        b <- geosmooth:::.ptf.weighted.anchor(make(7*w),diag(2))
        c <- geosmooth:::.ptf.weighted.anchor(make(rep(1,16)),diag(2))
        d <- geosmooth:::.ptf.weighted.anchor(make(rep(16,16)),diag(2))
        expect_equal(a$value,b$value,tolerance=1e-7)
        expect_equal(c$value,d$value,tolerance=1e-7)
        expect_true(all(a$records$status=="certified"))
        expect_gte(a$native.fits,4)
    }
})

test_that("finite-penalty GN acceptance agrees with an independent divided-difference objective", {
    X <- cbind(seq(-1,1,length.out=20),cos(seq_len(20)))
    for (family in c("gaussian","binomial")) {
        y <- sin(3*(.6*X[,1]+.8*X[,2]))
        if (family=="binomial") y <- plogis(y)
        args <- list(X=X,y=y,lambda=.0003,family=family,nterms=1,
            direction.search="gn",curve.ridge=.03,preprocessing="isotropic",backend="banded")
        base <- do.call(fit.projection.trend.filtering,c(args,list(control=list(gn.max.steps=0))))
        fit <- do.call(fit.projection.trend.filtering,c(args,list(control=list(gn.max.steps=2))))
        expect_true(base$converged); expect_true(fit$converged)
        expect_lte(fit$checks$objective,base$checks$objective+1e-9)
        roughness <- 0
        for (j in seq_len(fit$nterms)) {
            x <- fit$components[[j]]$fit$x
            ix <- order(x); t <- (x[ix]-min(x))/diff(range(x)); n <- length(t)
            D <- diff(diag(n))
            for (k in 1:2) D <- diff(diag(k/(t[(k+1):n]-t[1:(n-k)])) %*% D)
            roughness <- roughness+sum(abs(D %*% fit$terms[ix,j]))
        }
        eta <- as.numeric(predict(fit,type="link"))
        loss <- if (family=="gaussian") (y-eta)^2/2 else
            pmax(eta,0)-y*eta+log1p(exp(-abs(eta)))
        objective <- mean(loss)+.0003*roughness+.03/2*mean(rowSums(fit$terms^2))
        expect_equal(fit$checks$objective,objective,tolerance=1e-8)
    }
})

test_that("failed refinement trials retain preceding certified states", {
    X <- cbind(seq(-1,1,length.out=18),sin(seq_len(18)))
    y <- sin(3*(.6*X[,1]+.8*X[,2]))
    args <- list(X=X,y=y,lambda=.001,nterms=1,direction.search="gn",
        preprocessing="isotropic",backend="banded",curve.ridge=.03)
    base <- do.call(fit.projection.trend.filtering,c(args,list(control=list(gn.max.steps=0))))
    testthat::local_mocked_bindings(.ptf.gn.trial=function(...)
        list(state=NULL,status="injected_trial_failure"),.package="geosmooth")
    fit <- do.call(fit.projection.trend.filtering,args)
    expect_true(fit$converged)
    expect_equal(fit$checks$objective,base$checks$objective,tolerance=1e-10)
    expect_true("injected_trial_failure" %in% fit$search.history$status)
})

test_that("CV requires search eligibility and new modes use deterministic fold streams", {
    X <- cbind(seq(-1,1,length.out=18),cos(seq_len(18)))
    for (family in c("gaussian","binomial")) {
        y <- X[,1]^2+X[,2]/5
        if (family=="binomial") y <- plogis(y)
        args <- list(X=X,y=y,foldid=rep(1:3,6),lambda.grid=Inf,nterms.grid=c(0,1),
            direction.search="pls",family=family,backend="banded",keep.fits=TRUE,
            control=list(gn.max.steps=0),preprocessing="isotropic")
        a <- do.call(cv.projection.trend.filtering,args)
        b <- do.call(cv.projection.trend.filtering,args)
        expect_identical(a$status,"complete")
        expect_equal(a$predictions,b$predictions)
        expect_false(identical(a$fold.fits[[2]][[1]]$search.banks[[1]]$seed,
                               a$fold.fits[[2]][[2]]$search.banks[[1]]$seed))
        args$control$gn.max.fits <- 0
        failed <- do.call(cv.projection.trend.filtering,args)
        expect_false(failed$cv.table$eligible[2])
        expect_true(failed$cv.table$eligible[1])
        expect_equal(failed$selected.index,1)
        stopped <- failed$fold.results[failed$fold.results$candidate==2,]
        expect_true(all(stopped$numerical.certified))
        expect_false(any(stopped$search.eligible | stopped$cv.eligible))
        expect_true(all(stopped$effective.nterms==0 & stopped$requested.nterms==1))
        expect_true(all(stopped$native.fits==0))
    }
    fake <- fit.projection.trend.filtering(X,X[,1],Inf,nterms=0,direction.search="fixed")
    fake$cv.eligible <- FALSE
    testthat::local_mocked_bindings(fit.projection.trend.filtering=function(...) fake,
                                   .package="geosmooth")
    f <- cv.projection.trend.filtering(X,X[,1],rep(1:3,6),Inf,nterms.grid=0,warn=FALSE)
    expect_identical(f$status,"no_eligible_candidate")
})
