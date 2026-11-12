test_that("CV pools complete held-out losses and preserves mean-loss penalties", {
    X <- matrix(seq(-1,1,length.out=24),ncol=1); w <- seq(.5,2,length.out=24)
    folds <- cbind(rep(1:3,8),rep(3:1,each=8))
    for (family in c("gaussian","binomial")) {
        y <- sin(2*X[,1]); if (family=="binomial") y <- plogis(y)
        cv <- cv.projection.trend.filtering(X,y,folds,c(Inf,.002),nterms.grid=c(0,1),
            family=family,weights=w,direction.search="fixed",keep.fits=TRUE)
        expect_identical(cv$status,"complete")
        expect_true(all(cv$cv.table$eligible))
        for (j in seq_len(nrow(cv$cv.table))) {
            losses <- vapply(cv$predictions,function(P) {
                eta <- P[,j]
                loss <- if (family=="gaussian") (y-eta)^2 else pmax(eta,0)-y*eta+log1p(exp(-abs(eta)))
                sum(w*loss)/sum(w)
            },numeric(1))
            expect_equal(cv$cv.table$cv.loss[j],mean(losses))
            for (f in cv$fold.fits[[j]]) {
                expect_equal(f$lambda,cv$cv.table$lambda[j])
                if (f$nterms) expect_equal(f$components[[1]]$fit$curve.ridge,
                                          sum(f$weights)*cv$cv.table$curve.ridge[j])
            }
        }
        expect_identical(cv$selected.index,which.min(cv$cv.table$cv.loss))
        expect_equal(predict(cv),predict(cv$fit))
    }
})

test_that("held-out responses cannot change training preprocessing or learned directions", {
    X <- cbind(a=seq(-1,1,length.out=24),b=sin(seq_len(24)))
    y <- cos(2*X[,1]+X[,2]); folds <- rep(1:3,8)
    control <- list(starts=1,sparse.directions=1,dense.directions=1,
                    refinement.passes=1,refinement.proposals=4)
    args <- list(X=X,foldid=folds,lambda.grid=.002,nterms.grid=1,
                 keep.fits=TRUE,control=control,seed=17)
    a <- do.call(cv.projection.trend.filtering,c(args,list(y=y)))
    altered <- y; altered[folds==1] <- 3+altered[folds==1]
    b <- do.call(cv.projection.trend.filtering,c(args,list(y=altered)))
    f <- a$fold.fits[[1]][["1:1"]]; g <- b$fold.fits[[1]][["1:1"]]
    expect_identical(f$preprocessing,g$preprocessing)
    expect_identical(f$directions,g$directions)
    expect_equal(f$terms,g$terms)
    expect_equal(f$X,X[folds!=1,,drop=FALSE])
    expect_false(identical(a$predictions[[1]],b$predictions[[1]]))
})

test_that("group leakage, malformed folds and one-class training failures are explicit", {
    X <- matrix(seq(-1,1,length.out=12),ncol=1); y <- rep(c(0,1),each=6)
    groups <- rep(seq_len(6),each=2)
    expect_error(cv.projection.trend.filtering(X,y,rep(1:3,4),.01,groups=groups),"group")
    expect_error(cv.projection.trend.filtering(X,y,rep(1,12),.01),"two folds")
    expect_error(cv.projection.trend.filtering(X,y,1:4,.01),"foldid")
    cv <- cv.projection.trend.filtering(X,y,rep(1:2,each=6),.01,nterms.grid=1,
        family="binomial",direction.search="fixed",warn=FALSE)
    expect_identical(cv$status,"no_eligible_candidate")
    expect_null(cv$selected.index)
    expect_true(all(!cv$fold.results$converged))
    expect_true(all(grepl("finite unpenalized intercept",cv$fold.results$status)))
    expect_error(predict(cv),"No validated")
})
