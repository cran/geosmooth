test_that("practical patience counts accepted gains without redefining acceptance", {
    init <- geosmooth:::.ptf.practical
    step <- geosmooth:::.ptf.practical.update
    r <- init(1, .01, 2)
    r <- step(r, TRUE, 1, .995)
    expect_equal(r$streak, 1); expect_false(r$stopped)
    r <- step(r, FALSE, .995, NA_real_)
    expect_equal(r$streak, 1)
    r <- step(r, TRUE, .995, .98)
    expect_equal(r$streak, 0)
    r <- step(r, TRUE, .98, .975)
    r <- step(r, TRUE, .975, .97)
    expect_true(r$stopped)
    expect_false(step(init(0, 1), TRUE, 2, 1)$stopped)
    expect_equal(step(init(1, 0), TRUE, 2, 1)$streak, 0)
    expect_true(step(init(1, .5), TRUE, 1, .5)$small)
    expect_equal(init(1, .01)$streak, 0)
    p <- list(control = list(improvement.tolerance = 1e-9))
    b <- 1; boundary <- b-p$control$improvement.tolerance*max(1,abs(b))
    for (a in c(boundary-.Machine$double.eps, boundary, boundary+.Machine$double.eps))
        expect_identical(geosmooth:::.ptf.improves(p,a,b),
            is.finite(a) && a < b-p$control$improvement.tolerance*max(1,abs(b)))
    expect_false(geosmooth:::.ptf.improves(p, boundary, b))
})

test_that("optional observation preserves Gaussian and binomial joint fitting", {
    X <- cbind(seq(-1,1,length.out=18), cos(seq_len(18)))
    weights <- rep(c(2,3,5),6); offset <- .1*sin(seq_len(18))
    for (family in c("gaussian", "binomial")) {
        eta <- .3+X[,1]+.4*X[,2]^2
        y <- if (family == "gaussian") eta else plogis(eta)
        p <- geosmooth:::.ptf.problem(X,y,weights,offset,family,"isotropic",Inf,.1,
            geosmooth:::.ptf.control(list()),"banded")
        initial <- geosmooth:::.ptf.initial(p,diag(2))
        plain <- geosmooth:::.ptf.backfit(p, initial)
        events <- list()
        p$observer <- geosmooth:::.ptf.observer(function(event) events[[length(events)+1L]] <<- event)
        observed <- geosmooth:::.ptf.backfit(p, initial)
        expect_identical(observed, plain); expect_true(observed$converged)
        expect_length(p$observer$stack,0)
        expect_gte(p$observer$serialization.seconds, 0)
        expect_lte(p$observer$serialization.seconds, p$observer$diagnostic.seconds+1e-12)
        expect_true(all(vapply(Filter(function(e) e$event == "end", events),
            function(e) is.finite(e$timing$serialization.elapsed) &&
                e$timing$serialization.elapsed >= 0, logical(1))))
        native <- Filter(function(e) e$kind == "native" && e$event == "end",events)
        expect_equal(length(native),p$observer$native)
        expect_true(all(vapply(native,function(e)e$native.entries==1,logical(1))))
        sweeps <- Filter(function(e)e$kind == "sweep" && e$event == "end",events)
        expect_equal(length(sweeps),observed$sweeps)
        expect_true(all(vapply(sweeps,function(e)e$data$completed.updates==2,logical(1))))
        geosmooth:::.ptf.observe.snapshot(p,observed,"fixture")
        snap <- tail(events,1)[[1]]$data
        terms <- vapply(snap$state$components,function(c)
            as.numeric(predict(c$fit,c$fit$x,offset=numeric(nrow(X)),type="link"))-c$center,
            numeric(nrow(X)))
        expect_equal(terms,observed$terms,tolerance=1e-12)
        expect_equal(snap$link,offset+observed$intercept+rowSums(terms),tolerance=1e-12)
        null <- geosmooth:::.ptf.backfit(p,geosmooth:::.ptf.initial(p,matrix(0,0,2)))
        expected <- sum(p$w*geosmooth:::.ptf.loss(y,offset+null$intercept,family))
        expect_equal(null$checks$loss,expected,tolerance=1e-14)
        p2 <- geosmooth:::.ptf.problem(X,y,weights*7,offset,family,"isotropic",Inf,.1,
            geosmooth:::.ptf.control(list()),"banded")
        null2 <- geosmooth:::.ptf.backfit(p2,geosmooth:::.ptf.initial(p2,matrix(0,0,2)))
        expect_equal(null2$checks$loss,expected,tolerance=1e-12)
    }
})

test_that("observer failures escape joint rollback and numerical error handlers", {
    X <- cbind(seq(-1,1,length.out=12),sin(seq_len(12)))
    for (family in c("gaussian","binomial")) {
        y <- if(family=="gaussian") X[,1]^2+X[,2] else plogis(X[,1]^2+X[,2])
        for (where in c("native_begin","native_end","sweep_end","snapshot_end")) {
            p <- geosmooth:::.ptf.problem(X,y,NULL,NULL,family,"isotropic",Inf,.1,
                geosmooth:::.ptf.control(list()),"banded")
            p$observer <- geosmooth:::.ptf.observer(function(e) {
                if (paste(e$kind,e$event,sep="_")==where) stop("injected observer failure")
            })
            initial <- geosmooth:::.ptf.initial(p,diag(2))
            if(where=="snapshot_end") {
                fit <- geosmooth:::.ptf.backfit(p,initial)
                expect_error(geosmooth:::.ptf.observe.snapshot(p,fit,"fixture"),
                    class="projection_observer_error")
            } else expect_error(geosmooth:::.ptf.backfit(p,initial),class="projection_observer_error")
            if(where=="native_begin") expect_equal(p$observer$native,0)
            if(where=="native_end") expect_equal(p$observer$native,1)
        }
    }
    e <- geosmooth:::.ptf.infrastructure.error("quota",resource=TRUE)
    expect_error(geosmooth:::.ptf.numerical.error(e),class="projection_resource_error")
    ordinary <- simpleError("ordinary failure")
    expect_identical(geosmooth:::.ptf.numerical.error(ordinary),ordinary)
})

test_that("evaluated-trial visits retain the second accepted state and eighth-slot flags", {
    counter <- function() { e<-new.env(); e$visits<-e$events<-list(); e }
    g <- rep(1,6)
    testthat::local_mocked_bindings(
        .ptf.curve.slopes=function(...)g,
        .ptf.gn.equations=function(...)list(A=diag(2),b=c(0,1),floored=0),
        .ptf.gn.increment=function(...)list(status="proposal",mu=.01,delta=c(0,.1),tangent.score=1),
        .ptf.backfit=function(problem,state) { state$solver.calls<-0L; state },
        .package="geosmooth")
    base <- list(directions=matrix(c(1,0),1),terms=matrix(1,6,1),
        components=list(list(fit=list(design=list(span=1)))),
        checks=list(objective=1),converged=TRUE,status="converged")
    for (reject in c(0L,6L)) {
        p <- list(observer=NULL,control=list(improvement.tolerance=1e-9),
            practical=geosmooth:::.ptf.practical(1,.0001),counter=counter(),
            offset=rep(0,6),Z=matrix(0,6,2),y=rep(0,6),w=rep(1/6,6),family="gaussian")
        i <- 0L
        evaluate <- function(state,u,...) {
            i <<- i+1L
            if(i<=reject)return(NULL)
            state$directions[1,] <- u
            state$checks$objective <- state$checks$objective-1e-5
            list(state=state,row=i)
        }
        got <- geosmooth:::.ptf.gn.visit(p,base,1,8,1,1,1,evaluate,function(...)NULL,
            function(a,b)is.finite(a)&&a<b-1e-9*max(1,abs(b)))
        visit <- p$counter$visits[[1]]
        expect_identical(visit$stopping,"practical_gain_stop")
        expect_equal(visit$evaluated,reject+2)
        expect_equal(got$state$checks$objective,1-2e-5)
        expect_identical(visit$proposal.limit.reached,reject==6L)
        expect_equal(visit$unused.slots,6L-reject)
        expect_length(got$accepted.rows,2)
    }
})
