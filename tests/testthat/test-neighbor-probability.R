test_that("all three C++ paths match the original sPTB Python predictions and tuning", {
    fixture <- jsonlite::fromJSON(testthat::test_path("fixtures", "neighbor-python.json"),
                                  simplifyVector = FALSE)
    mat <- function(x) do.call(rbind, lapply(x, unlist))
    for (a in fixture$fixtures) {
        y <- unlist(a$y); id <- unlist(a$ids)
        D <- mat(a$D); X <- mat(a$X)
        train <- unlist(a$train); query <- unlist(a$query)
        norm <- if (a$metric == "l1_euclidean") "sum" else "max"
        metric <- if (a$metric == "max_linf") "chebyshev" else "euclidean"
        for (m in c("knn", "case_raw", "case_corrected")) {
            fixed <- fit.neighbor.probability(y = y[train], distances = D[train,train],
                                              method = m, k = 2, id = id[train])
            expected <- mat(a$predictions[[m]])
            sizes <- if (m == "case_corrected") seq_len(ncol(expected))+1L else
                seq_len(ncol(expected))
            actual <- predict(fixed, distances = D[query,train], k = sizes)
            expect_equal(unname(actual), expected, tolerance = 1e-15)
            coords <- fit.neighbor.probability(X[train,], y[train], method = m, k = 2,
                                              metric = metric, normalization = norm,
                                              id = id[train])
            expect_equal(unname(predict(coords, X[query,], k = sizes)), expected,
                         tolerance = 1e-15)
            got <- predict(fixed, distances = D[query,train], type = "diagnostics")
            actual.diag <- lapply(seq_along(query), function(j) list(
                neighbors=got$neighbors[j,1], cases=got$cases[j,1],
                radius=got$radius[j,1], boundary.ties=got$boundary.ties[j,1],
                split.tie=got$split.tie[j,1]))
            expected.diag <- lapply(a$diagnostics[[m]], function(d) list(
                neighbors=d$neighbors, cases=d$cases, radius=d$radius,
                boundary.ties=d$boundary_ties, split.tie=d$split_tie))
            expect_equal(actual.diag, expected.diag)
            tuned <- fit.neighbor.probability(y = y[train], distances = D[train,train],
                        method = m, id = id[train], foldid = mat(a$folds))
            expect_equal(tuned$k, a$selected[[m]])
            rows <- Filter(function(s) s$method == m, a$scores)
            expect_equal(tuned$cv.table$brier, vapply(rows, function(s) s$brier, 0),
                         tolerance = 1e-14)
            expect_equal(tuned$cv.table$k, unlist(a$grids[[m]]))
            expect_identical(tuned$backend, "C++")
        }
    }
})

test_that("the kth case is included and corrected estimates obey the exact formula", {
    X <- matrix(1:10, ncol=1); y <- as.numeric(1:10 %in% c(2,5,9))
    query <- matrix(0, ncol=1)
    raw <- fit.neighbor.probability(X,y,method="case_raw",k=3)
    corr <- fit.neighbor.probability(X,y,method="case_corrected",k=3)
    nn <- fit.neighbor.probability(X,y,k=5)
    expect_equal(as.numeric(predict(raw,query)), 3/9)
    expect_equal(as.numeric(predict(corr,query)), 2/8)
    expect_equal(as.numeric(predict(nn,query)), 2/5)
    expect_lte(as.numeric(predict(corr,query)), as.numeric(predict(raw,query)))
    expect_true(is.na(predict(raw,query,k=4)[1]))
    expect_true(is.na(predict(nn,query,k=11)[1]))
    expect_error(fit.neighbor.probability(X,y,method="case_corrected",k=1), ">= 2")
    expect_error(fit.neighbor.probability(X,y,k=1.5), "integers")
})

test_that("participants are excluded and closest answering visits are label-blind", {
    X <- matrix(c(1,2,3,4,5),ncol=1)
    y <- c(0,1,1,NA,0); group <- c("a","a","b","c","d")
    f <- fit.neighbor.probability(X,y,k=2,group=group)
    expect_equal(as.numeric(predict(f,matrix(0,ncol=1))), 1/2) # a first visit, b
    expect_equal(as.numeric(predict(f,matrix(0,ncol=1),query.group="a")), 1/2) # b,d
    expect_equal(as.numeric(predict(f,matrix(2.1,ncol=1),k=1)), 1) # a second visit
    expect_equal(as.numeric(predict(f)[1:2,1]), c(.5,.5)) # exclude both a visits
    y2 <- y; y2[1:2] <- 1-y2[1:2]
    f2 <- fit.neighbor.probability(X,y2,k=2,group=group)
    expect_equal(predict(f)[1:2,,drop=FALSE],predict(f2)[1:2,,drop=FALSE])
    expect_error(fit.neighbor.probability(X,y,group=group,foldid=1:5), "cross folds")
    missing <- fit.neighbor.probability(X,c(NA,NA,1,NA,0),k=2)
    expect_equal(as.numeric(predict(missing,matrix(0,ncol=1))),.5)
})

test_that("IDs determine exact ties reproducibly under record permutation", {
    X <- matrix(rep(1,6),ncol=1); y <- c(1,0,1,0,1,0)
    id <- c("a","b","é","d","e","f")
    p <- c(6,4,2,5,3,1)
    a <- fit.neighbor.probability(X,y,k=2,id=id)
    b <- fit.neighbor.probability(X[p,,drop=FALSE],y[p],k=2,id=id[p])
    expect_equal(predict(a,matrix(1,ncol=1),k=1:6), predict(b,matrix(1,ncol=1),k=1:6))
    d <- predict(a,matrix(1,ncol=1),type="diagnostics")
    expect_equal(d$boundary.ties[1,1],6)
    expect_true(d$split.tie[1,1])
})

test_that("missing cases, unreachable vertices and one-class data are explicit", {
    X <- matrix(1:4,ncol=1)
    no <- fit.neighbor.probability(X,rep(0,4),method="case_raw",k=1)
    expect_true(all(is.na(predict(no))))
    yes <- fit.neighbor.probability(X,rep(1,4),method="case_corrected",k=2)
    expect_equal(as.numeric(predict(yes)),rep(1,4))
    D <- as.matrix(dist(1:4))
    f <- fit.neighbor.probability(y=c(1,1,0,0),distances=D,k=2)
    d <- predict(f,distances=matrix(c(1,Inf,Inf,Inf),nrow=1),type="diagnostics")
    expect_equal(d$available.people,1L)
    expect_identical(unname(d$status[1,1]),"insufficient_people")
    expect_error(fit.neighbor.probability(X,rep(NA_real_,4),k=1), "observed")
    expect_error(fit.neighbor.probability(X,c(1,0,2,0),k=1), "binary")
    expect_error(fit.neighbor.probability(X,c(1,0,Inf,0),k=1), "binary")
    expect_error(fit.neighbor.probability(X,c(1,0,NaN,0),k=1), "binary")
    expect_error(fit.neighbor.probability(X,c(1,0,1,0),k=1,id=rep("a",4)), "unique")
    expect_error(fit.neighbor.probability(X,c(1,0,1,0),k=1,group=c("a",NA,"b","c")), "nonmissing")
    expect_error(predict(yes,matrix(1,ncol=1),query.group=NA), "nonmissing")
    expect_error(predict(yes,matrix(c(1,2),nrow=1)), "feature")
    expect_error(predict(yes,distances=matrix(1,1,4)), "newdata")
    expect_error(predict(f,matrix(1,ncol=1)), "distances")
    expect_error(predict(f,distances=matrix(NA_real_,1,4)), "nonnegative")
    expect_error(fit.neighbor.probability(matrix(0,4,1),c(1,0,1,0),k=1,
                                          normalization="sum"), "positive mass")
})

test_that("tuning is participant-balanced with repeated records and respects selection rules", {
    X <- matrix(c(1,1.1,2,3,4,5,6,7,8),ncol=1)
    y <- c(0,1,1,0,1,0,1,0,1)
    group <- c("a","a",letters[2:8]); fold <- c(1,1,2,3,1,2,3,1,2)
    f <- fit.neighbor.probability(X,y,group=group,foldid=fold,k.grid=1:3)
    loss <- matrix(NA_real_,length(y),3)
    for (v in unique(fold)) {
        tr <- which(fold!=v); te <- which(fold==v)
        g <- fit.neighbor.probability(X[tr,,drop=FALSE],y[tr],id=as.character(tr),
                                      group=group[tr],k=1)
        loss[te,] <- (predict(g,X[te,,drop=FALSE],query.group=group[te],k=1:3)-y[te])^2
    }
    person <- do.call(rbind,lapply(unique(group),function(h)colMeans(loss[group==h,,drop=FALSE])))
    expect_equal(f$cv.table$brier,unname(colMeans(person)))
    expect_equal(f$cv.table$se,unname(apply(person,2,sd)/sqrt(nrow(person))))
    minimum <- max(which(f$cv.table$brier<=min(f$cv.table$brier)+1e-14))
    good <- which(f$cv.table$brier<=f$cv.table$brier[minimum]+f$cv.table$se[minimum]+1e-14)
    large <- fit.neighbor.probability(X,y,group=group,foldid=fold,k.grid=1:3,selection="one.se")
    small <- fit.neighbor.probability(X,y,group=group,foldid=fold,k.grid=1:3,selection="locality.one.se")
    expect_equal(large$k,max(good))
    expect_equal(small$k,min(good))
    expect_error(fit.neighbor.probability(X,y,foldid=fold,k.grid=100), "No candidate")
    expect_error(fit.neighbor.probability(X,y,k=2,foldid=fold), "not both")
    expect_error(fit.neighbor.probability(X,y), "Supply fixed")
})

test_that("tiny label configurations agree with independent rank arithmetic", {
    X <- matrix(c(0,1,2,4,8,16),ncol=1)
    q <- matrix(0.1,ncol=1)
    predictions <- references <- list()
    for (bits in 0:63) {
        y <- as.numeric(intToBits(bits))[1:6]
        ord <- order(abs(X[,1]-.1))
        stop <- which(y[ord]==1)
        for (m in c("knn","case_raw","case_corrected")) {
            sizes <- if (m=="case_corrected") 2:6 else 1:6
            fit <- fit.neighbor.probability(X,y,method=m,k=min(sizes))
            actual <- as.numeric(predict(fit,q,k=sizes))
            expected <- vapply(sizes,function(k) {
                if (m=="knn") return(mean(y[ord[1:k]]))
                if (k>length(stop)) return(NA_real_)
                if (m=="case_raw") k/stop[k] else (k-1)/(stop[k]-1)
            },0)
            key <- paste(bits,m,sep="/")
            predictions[[key]] <- actual
            references[[key]] <- expected
        }
    }
    expect_equal(predictions,references)
})
