# Independent dense heat-penalty reference and public policy checks.
test_that('independent dense heat-penalty equations reproduce ordinary accepted paths', {
    set.seed(9172026)
    for (case in seq_len(30L)) {
        n <- sample(4:10,1)
        W <- matrix(0,n,n)
        for(i in seq_len(n-1L)) W[i,i+1L] <- runif(1,.4,2)
        for(i in seq_len(n-2L)) for(j in seq.int(i+2L,n)) {
            if(runif(1)<.3) W[i,j] <- runif(1,.4,2)
        }
        W <- W+t(W)
        a <- lapply(seq_len(n),function(i) which(W[i,]>0))
        w <- lapply(seq_len(n),function(i) W[i,a[[i]]])
        b <- metric.graph.lowpass.basis(a,w,weight.type='conductance',spectrum='full',
                                       numerical.policy='verified.fallback')
        y <- rnorm(n)
        O <- sort(sample.int(n,sample(3:n,1)))
        y[-O] <- NA_real_
        L <- diag(rowSums(W))-W
        ev <- eigen(L,symmetric=TRUE)
        times <- c(.01,.1,.5)/max(ev$values)
        path <- apply.metric.graph.lowpass.path(b,y,times)
        R <- diag(n)[O,,drop=FALSE]
        for(j in seq_along(times)) {
            penalty <- ev$vectors %*% diag(expm1(times[j]*pmax(0,ev$values))) %*% t(ev$vectors)
            T <- solve(crossprod(R)+penalty,t(R))
            pred <- as.vector(T %*% y[O])
            df <- sum(diag(R %*% T))
            gcv <- length(O)*sum((y[O]-pred[O])^2)/(length(O)-df)^2
            expect_equal(unname(path$fitted.values[,j]),pred,tolerance=1e-9)
            expect_equal(path$effective.df[j],df,tolerance=1e-9)
            expect_equal(path$gcv.scores[j],gcv,tolerance=1e-8)
        }
    }
})

test_that('legacy defaults and public policy overrides do not silently enable fallback', {
    a <- list(2L,c(1L,3L),c(2L,4L),4L)
    # Correct symmetric four-node chain.
    a[[4L]] <- 3L
    w <- list(1,c(1,2),c(2,1),1)
    y <- c(0,NA,2,3)
    b <- metric.graph.lowpass.basis(a,w,weight.type='conductance',spectrum='full')
    v <- apply.metric.graph.lowpass.path(b,y,c(.01,.1,1),numerical.policy='verified.fallback')
    expect_identical(b$parameters$numerical.policy,'conservative')
    expect_identical(v$numerical$policy,'verified.fallback')
    bv <- b; bv$parameters$numerical.policy <- 'verified.fallback'
    c <- apply.metric.graph.lowpass.path(bv,y,c(.01,.1,1),numerical.policy='conservative')
    expect_null(c$numerical)
    expect_equal(c$fitted.values,v$fitted.values,tolerance=1e-10)
    expect_error(apply.metric.graph.lowpass.path(b,c(0,NaN,2,3),.1,
                 numerical.policy='verified.fallback'),'NaN')
    expect_error(apply.metric.graph.lowpass.path(b,c(0,Inf,2,3),.1,
                 numerical.policy='verified.fallback'),'finite|Inf')
    path <- apply.metric.graph.lowpass.path(b,y,c(1,.01,.1,.01),numerical.policy='verified.fallback')
    expect_equal(path$eta.grid,c(1,.01,.1,.01))
    one <- apply.metric.graph.lowpass.path(b,c(5,NA,NA,NA),c(.1,1),numerical.policy='verified.fallback')
    expect_equal(unname(one$fitted.values),matrix(5,4,2),tolerance=1e-10)
    expect_identical(one$tuning[[1]]$status,'no_valid_candidate')
    expect_true(all(is.na(one$gcv.scores)))
    z <- apply.metric.graph.lowpass.path(b,y,0,numerical.policy='verified.fallback')
    expect_equal(as.vector(z$fitted.values),y)
    expect_equal(nrow(z$numerical$diagnostics),0L)
})
