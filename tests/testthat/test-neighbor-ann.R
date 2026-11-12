# Large enough to exercise ANN retrieval rather than the exhaustive small-set path.
neighbor.compare.exact <- function(a, b) {
    fields <- c("probability", "neighbors", "cases", "radius", "boundary.ties",
                "split.tie", "available.people", "available.cases", "status")
    expect_equal(a[fields], b[fields], tolerance = 0)
}

test_that("exact ANN matches brute force on normalized and repeated-participant profiles", {
    set.seed(716)
    X <- matrix(rexp(480 * 12), 480)
    X[100:110, ] <- X[rep(1, 11), ]
    y <- sample(c(0,1,NA), 480, TRUE, prob = c(.6,.3,.1))
    id <- sprintf("record-%04d", seq_len(480))
    group <- rep(seq_len(240), 2)
    query <- rbind(X[1:12,], matrix(rexp(12 * 12), 12))
    query.group <- c(group[1:12], rep(0,12))
    for (norm in c("none", "sum", "max")) for (m in c("knn","case_raw","case_corrected")) {
        args <- list(X=X,y=y,id=id,group=group,method=m,k=2,normalization=norm)
        a <- do.call(fit.neighbor.probability,c(args,list(search="ann")))
        b <- do.call(fit.neighbor.probability,c(args,list(search="brute")))
        sizes <- c(2:8,20,60,1000)
        # Separate a local path from unsupported sizes so ANN must use its tree.
        for (ks in list(2:8, c(2,20,60), sizes)) {
            pa <- predict(a,query,query.group=query.group,k=ks,type="diagnostics")
            pb <- predict(b,query,query.group=query.group,k=ks,type="diagnostics")
            neighbor.compare.exact(pa,pb)
            expect_identical(predict(a,query,query.group=query.group,k=ks),pa$probability)
        }
        expect_identical(a$search,"ann")
        expect_identical(a$eps,0)
    }
})

test_that("exact ANN completes huge tied shells and diagnostic near-ties", {
    # More records at a tied boundary than the initial search retrieves.
    X <- matrix(c(rep(0,120),rep(1,120),1+5e-13,rep(2,60)),ncol=1)
    y <- rep(c(1,0,0),length.out=nrow(X))
    id <- sprintf("tie-%04d",seq_len(nrow(X)))
    group <- rep(seq_len(151),length.out=nrow(X))
    query <- matrix(c(0,.5,1,3),ncol=1)
    for (m in c("knn","case_raw","case_corrected")) {
        a <- fit.neighbor.probability(X,y,method=m,k=3,id=id,group=group,search="ann")
        b <- fit.neighbor.probability(X,y,method=m,k=3,id=id,group=group,search="brute")
        neighbor.compare.exact(predict(a,query,k=2:20,type="diagnostics"),
                               predict(b,query,k=2:20,type="diagnostics"))
        set.seed(95); p <- sample(seq_len(nrow(X)))
        perm <- fit.neighbor.probability(X[p,,drop=FALSE],y[p],method=m,k=3,
                                         id=id[p],group=group[p],search="ann")
        expect_identical(predict(a,query,k=2:20),predict(perm,query,k=2:20))
    }
})

test_that("ANN handles rare, absent and excluded cases", {
    X <- matrix(seq_len(300),ncol=1)
    for (y in list(rep(0,300),rep(1,300),as.numeric(seq_len(300)>298))) {
        for (m in c("knn","case_raw","case_corrected")) {
            a <- fit.neighbor.probability(X,y,method=m,k=2,search="ann")
            b <- fit.neighbor.probability(X,y,method=m,k=2,search="brute")
            neighbor.compare.exact(predict(a,X[c(1,299,300),,drop=FALSE],
                                           query.group=c("1","299","300"),k=c(2,20),type="diagnostics"),
                                   predict(b,X[c(1,299,300),,drop=FALSE],
                                           query.group=c("1","299","300"),k=c(2,20),type="diagnostics"))
        }
    }
})

test_that("grouped repeated CV and serialization preserve exact ANN results", {
    set.seed(294)
    X <- matrix(runif(360*4),360)
    y <- sample(c(0,1,NA),360,TRUE,prob=c(.5,.45,.05))
    group <- rep(1:180,2)
    fold <- cbind(rep(rep(1:3,60),2),rep(sample(rep(1:3,60)),2))
    for (m in c("knn","case_raw","case_corrected")) {
        a <- fit.neighbor.probability(X,y,group=group,foldid=fold,k.grid=c(2,4,8),
                                      method=m,search="ann")
        b <- fit.neighbor.probability(X,y,group=group,foldid=fold,k.grid=c(2,4,8),
                                      method=m,search="brute")
        expect_identical(a$cv.table,b$cv.table)
        expect_identical(a$k,b$k)
        expect_identical(predict(unserialize(serialize(a,NULL)),X[1:3,]),predict(a,X[1:3,]))
    }
})

test_that("positive eps is explicit, reproducible and preserves closest answering visits", {
    set.seed(310)
    X <- matrix(rnorm(800*8),800)
    y <- rep(c(0,1,NA,1),200)
    id <- sprintf("approx-%04d",1:800)
    group <- rep(1:400,2)
    q <- X[1:10,]
    for (m in c("knn","case_raw","case_corrected")) {
        a <- fit.neighbor.probability(X,y,k=3,method=m,id=id,group=group,eps=2)
        d <- predict(a,q,query.group=group[1:10],k=2:10,type="diagnostics")
        b <- fit.neighbor.probability(X,y,k=3,method=m,id=id,group=group,search="brute")
        exact <- predict(b,q,query.group=group[1:10],k=2:10,type="diagnostics")
        expect_true(all(is.finite(d$probability) & d$probability>=0 & d$probability<=1))
        expect_identical(d$available.people,exact$available.people)
        expect_identical(d$available.cases,exact$available.cases)
        expect_identical(d$boundary.complete,FALSE)
        expect_identical(d$eps,2)
        expect_identical(predict(a,q,k=2:10),predict(a,q,k=2:10))
        p <- sample(1:800)
        perm <- fit.neighbor.probability(X[p,],y[p],k=3,method=m,id=id[p],group=group[p],eps=2)
        expect_identical(predict(a,q,k=2:10),predict(perm,q,k=2:10))
        # When every participant is requested, approximation reduces to exhaustive ranks.
        expect_identical(predict(a,q,k=400),predict(b,q,k=400))
    }
    # Closest answering visit must win even if a different visit enters candidates.
    X <- matrix(c(seq(100,499),seq(0,.399,length.out=400)),ncol=1)
    y <- c(rep(1,400),rep(0,400))
    a <- fit.neighbor.probability(X,y,k=10,group=rep(1:400,2),eps=2)
    expect_equal(as.numeric(predict(a,matrix(0,ncol=1))),0)
})

test_that("search routing and invalid epsilon values are checked", {
    X <- matrix(1:100,ncol=1); y <- rep(0:1,50)
    expect_identical(fit.neighbor.probability(X,y,k=2)$search,"ann")
    expect_identical(fit.neighbor.probability(X,y,k=2,metric="chebyshev")$search,"brute")
    expect_identical(fit.neighbor.probability(y=y,distances=as.matrix(dist(X)),k=2)$search,"brute")
    expect_error(fit.neighbor.probability(X,y,k=2,search="ann",metric="chebyshev"),"Euclidean")
    expect_error(fit.neighbor.probability(y=y,distances=as.matrix(dist(X)),k=2,search="ann"),"Euclidean")
    expect_error(fit.neighbor.probability(X,y,k=2,search="brute",eps=.1),"requires ANN")
    expect_error(fit.neighbor.probability(X,y,k=2,metric="chebyshev",eps=.1),"requires ANN")
    for (eps in list(-1,Inf,NaN,NA_real_,numeric(),c(0,1),"0",1e200))
        expect_error(fit.neighbor.probability(X,y,k=2,eps=eps),"eps")
})

test_that("exact ANN recovers high-dimensional sparse ranks and large excluded groups", {
    set.seed(819)
    X <- matrix(rexp(360*340),360); X[X < 1.5] <- 0
    y <- sample(c(0,1,NA),360,TRUE,prob=c(.6,.35,.05))
    group <- c(rep("repeat",120),paste0("p",1:240))
    X[1:120,] <- X[rep(121,120),]
    q <- X[c(1,121,220,360),]
    for (m in c("knn","case_raw","case_corrected")) for (scale in c(1e-100,1,1e100)) {
        a <- fit.neighbor.probability(X*scale,y,k=3,group=group,method=m,search="ann")
        b <- fit.neighbor.probability(X*scale,y,k=3,group=group,method=m,search="brute")
        neighbor.compare.exact(predict(a,q*scale,query.group=rep("repeat",4),k=c(3,10),type="diagnostics"),
                               predict(b,q*scale,query.group=rep("repeat",4),k=c(3,10),type="diagnostics"))
    }
})
