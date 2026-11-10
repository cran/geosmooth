intrinsic.compare <- function(a, b) {
    for (field in c("probability", "neighbors", "cases", "boundary.ties",
                    "split.tie", "available.people", "available.cases", "status"))
        expect_identical(a[[field]], b[[field]], info=field)
    expect_equal(a$radius, b$radius, tolerance=1e-14)
}

# Explicit exhaustive query-to-training matrices, using the public dgraphs API.
intrinsic.exhaustive <- function(X, Q) {
    Z <- rbind(Q,X)
    dgraphs::linf.simplex.distances(Z,sources=seq_len(nrow(Q)),
                                   targets=nrow(Q)+seq_len(nrow(X)))
}

test_that("streamed intrinsic paths match exhaustive surface distances for all methods", {
    set.seed(11903)
    for (d in c(1L,3L,12L,340L)) {
        X <- matrix(runif(72*d),72,d); X[X<.5] <- 0; X[,1] <- pmax(X[,1],.01)
        X[2:4,] <- X[rep(1,3),]
        y <- sample(c(0,1,NA),72,TRUE,prob=c(.5,.4,.1))
        id <- sprintf("intrinsic-%04d",1:72); group <- rep(1:36,2)
        Q <- rbind(X[1:5,,drop=FALSE],matrix(runif(5*d),5,d))
        qgroup <- c(group[1:5],rep(0,5))
        D <- dgraphs::linf.simplex.distances(X)
        Dq <- intrinsic.exhaustive(X,Q)
        for (method in c("knn","case_raw","case_corrected")) {
            a <- fit.neighbor.probability(X,y,method=method,k=3,id=id,group=group,
                                          metric="linf.simplex",normalization="max")
            b <- fit.neighbor.probability(y=y,distances=D,method=method,k=3,id=id,group=group)
            for (ks in list(2:8,c(2,20,100))) {
                pa <- predict(a,Q,query.group=qgroup,k=ks,type="diagnostics")
                pb <- predict(b,distances=Dq,query.group=qgroup,k=ks,type="diagnostics")
                intrinsic.compare(pa,pb)
                expect_identical(predict(a,Q,query.group=qgroup,k=ks),pa$probability)
            }
            expect_identical(a$search,"brute")
            expect_identical(a$distance.backend,"dgraphs")
            expect_null(a$distances)
            # Serialization retains coordinates, not a native pointer/cache.
            expect_identical(predict(unserialize(serialize(a,NULL)),Q),predict(a,Q))
        }
    }
})

test_that("intrinsic metric changes neighbors and admits intermediate faces", {
    # Same query and endpoints as the independently checked surface example.
    X <- rbind(c(0,1,0,0),c(1,1,1,1)); q <- rbind(c(1,0,0,0))
    a <- fit.neighbor.probability(X,c(1,0),k=1,normalization="max",metric="linf.simplex")
    b <- fit.neighbor.probability(X,c(1,0),k=1,normalization="max",metric="euclidean")
    expect_equal(unname(predict(a,q)),matrix(0,1,1))
    expect_equal(unname(predict(b,q)),matrix(1,1,1))
    expect_equal(predict(a,q,type="diagnostics")$radius[1,1],sqrt(3))
    # The third coordinate supplies a shorter path than endpoint-face unfolding.
    X <- rbind(c(0,1,.99),c(1,1,0)); q <- rbind(c(1,0,.99))
    for (method in c("knn","case_raw","case_corrected")) {
        a <- fit.neighbor.probability(X,c(1,1),k=2,method=method,
                                      normalization="max",metric="linf.simplex")
        diag <- predict(a,q,k=if(method=="case_corrected")2 else 1:2,type="diagnostics")
        expect_equal(min(diag$radius),if(method=="case_corrected")1.01*sqrt(2) else sqrt(1+.99^2))
        # Full path ordered radii: same-face endpoint is closer than the cross-face point.
        expect_equal(max(diag$radius),1.01*sqrt(2),tolerance=1e-14)
    }
})

test_that("same-face intrinsic paths equal Euclidean chords", {
    set.seed(431)
    X <- cbind(1,matrix(runif(160*6),160)); y <- rep(c(0,1,1,NA),40)
    Q <- cbind(1,matrix(runif(8*6),8))
    for (method in c("knn","case_raw","case_corrected")) {
        a <- fit.neighbor.probability(X,y,k=3,method=method,metric="linf.simplex")
        b <- fit.neighbor.probability(X,y,k=3,method=method,search="brute")
        intrinsic.compare(predict(a,Q,k=3:20,type="diagnostics"),
                           predict(b,Q,k=3:20,type="diagnostics"))
    }
})

test_that("intrinsic CV, participant exclusions and tied ordering match exhaustive references", {
    set.seed(610)
    X <- matrix(sample(0:3,96*6,TRUE),96); X[,1] <- pmax(1,X[,1])
    X[1:8,] <- X[rep(9,8),]
    y <- sample(c(0,1,NA),96,TRUE,prob=c(.45,.5,.05))
    group <- rep(1:48,2); id <- paste0("s",1:96)
    folds <- cbind(rep(rep(1:3,16),2),rep(sample(rep(1:3,16)),2))
    D <- dgraphs::linf.simplex.distances(X)
    for (method in c("knn","case_raw","case_corrected")) {
        a <- fit.neighbor.probability(X,y,method=method,id=id,group=group,
            metric="linf.simplex",normalization="max",foldid=folds,k.grid=c(2,3,6))
        b <- fit.neighbor.probability(y=y,distances=D,method=method,id=id,group=group,
                                     foldid=folds,k.grid=c(2,3,6))
        expect_identical(a$cv.table,b$cv.table)
        expect_identical(a$k,b$k)
        intrinsic.compare(predict(a,type="diagnostics"),predict(b,type="diagnostics"))
        p <- sample(1:96)
        c <- fit.neighbor.probability(X[p,],y[p],method=method,k=a$k,id=id[p],group=group[p],
                                      metric="linf.simplex",normalization="max")
        expect_identical(predict(a,X[1:8,],query.group=group[1:8]),
                         predict(c,X[1:8,],query.group=group[1:8]))
    }
})

test_that("intrinsic missing, discordant and unsupported responses remain explicit", {
    X <- rbind(c(1,0),c(1,.1),c(0,1),c(1,.5),c(1,1))
    group <- c("a","a","b","c","d"); y <- c(0,1,1,NA,0)
    Q <- rbind(c(1,0),c(1,.1),c(0,1))
    D <- dgraphs::linf.simplex.distances(X); Dq <- intrinsic.exhaustive(X,Q)
    for (method in c("knn","case_raw","case_corrected")) {
        for (labels in list(y,rep(0,5),rep(1,5),c(NA,NA,1,NA,NA))) {
            a <- fit.neighbor.probability(X,labels,method=method,k=2,group=group,metric="linf.simplex")
            b <- fit.neighbor.probability(y=labels,distances=D,method=method,k=2,group=group)
            intrinsic.compare(predict(a,Q,query.group=c("a","a","b"),k=2:5,type="diagnostics"),
                               predict(b,distances=Dq,query.group=c("a","a","b"),k=2:5,type="diagnostics"))
        }
    }
    a <- fit.neighbor.probability(X,y,k=1,group=group,metric="linf.simplex")
    expect_equal(as.numeric(predict(a,Q[2,,drop=FALSE])),1)
})

test_that("intrinsic coordinates and search options are validated", {
    X <- rbind(c(2,1),c(1,2),c(2,2)); y <- c(1,0,1)
    a <- fit.neighbor.probability(X,y,k=1,metric="linf.simplex",normalization="max")
    expect_error(fit.neighbor.probability(X,y,k=1,metric="linf.simplex"),"max-normalized")
    expect_error(fit.neighbor.probability(X,y,k=1,metric="linf.simplex",normalization="sum"),"Intrinsic")
    expect_error(fit.neighbor.probability(X,y,k=1,metric="linf.simplex",search="ann"),"Euclidean")
    expect_error(fit.neighbor.probability(X,y,k=1,metric="linf.simplex",eps=.1),"requires ANN")
    expect_error(fit.neighbor.probability(-X,y,k=1,metric="linf.simplex",normalization="max"),"nonnegative")
    expect_error(fit.neighbor.probability(X*0,y,k=1,metric="linf.simplex",normalization="max"),"positive")
    expect_error(predict(a,-X),"nonnegative")
    expect_error(predict(a,X*0),"positive")
    b <- fit.neighbor.probability(X/2*(1+1e-10),y,k=1,metric="linf.simplex")
    expect_equal(predict(a,X),predict(b,X/2),tolerance=1e-14)
    # Scaling each raw profile cannot change an intrinsic estimate.
    c <- fit.neighbor.probability(X*c(2,3,4),y,k=1,metric="linf.simplex",normalization="max")
    expect_identical(predict(a,X),predict(c,X*7))
})

test_that("prediction and CV release native caches on success and error", {
    set.seed(272)
    X <- matrix(runif(90*4),90); X <- X/apply(X,1,max); y <- rep(0:1,45)
    original <- geosmooth:::.neighbor.index
    retained <- list()
    testthat::local_mocked_bindings(.neighbor.index=function(object,train) {
        ptr <- original(object,train)
        if(!is.null(ptr)) retained[[length(retained)+1L]] <<- list(ptr=ptr,train=train)
        ptr
    }, .package="geosmooth")
    for(metric in c("euclidean","linf.simplex")) {
        retained <- list()
        f <- fit.neighbor.probability(X,y,k=3,metric=metric)
        predict(f,X[1:3,])
        expect_length(retained,1L)
        expect_error(geosmooth:::.neighbor.at.rows(f,1L,retained[[1]]$train,3L,FALSE,
                     retained[[1]]$ptr),"Invalid (ANN|intrinsic) index")
        retained <- list()
        fit.neighbor.probability(X,y,metric=metric,foldid=rep(1:3,30),k.grid=2:4)
        expect_length(retained,3L)
        for(cache in retained)
            expect_error(geosmooth:::.neighbor.at.rows(f,1L,cache$train,3L,FALSE,cache$ptr),
                          "Invalid (ANN|intrinsic) index")
        retained <- list()
        expect_error(geosmooth:::.neighbor.call(f,matrix(1,1,3),NULL,0L,1:90,3L,FALSE),
                      "coordinate dimensions")
        expect_error(geosmooth:::.neighbor.at.rows(f,1L,retained[[1]]$train,3L,FALSE,
                     retained[[1]]$ptr),"Invalid (ANN|intrinsic) index")
    }
})
