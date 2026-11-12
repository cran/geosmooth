# Response-mask, permutation and extreme-scale invariance checks.
test_that('component centering survives node permutations and mixed response masks', {
    set.seed(6172026)
    W <- matrix(0,9,9)
    W[1,2] <- 1; W[2,3] <- 2; W[3,4] <- .5
    W[5,6] <- 3; W[6,7] <- 1; W[7,8] <- 2
    W <- W+t(W)
    Y <- cbind(c(0,NA,1,2,NA,3,-1,2,7),c(NA,1,NA,3,2,NA,0,1,4),
               c(5,NA,5,5,NA,-2,-2,-2,1))
    shift <- c(rep(2^30,4),rep(-2^28,4),2^35)
    makebasis <- function(M) {
        a <- lapply(seq_len(nrow(M)),function(i) which(M[i,]>0))
        w <- lapply(seq_len(nrow(M)),function(i) M[i,a[[i]]])
        metric.graph.lowpass.basis(a,w,weight.type='conductance',spectrum='full',
                                  numerical.policy='verified.fallback')
    }
    times <- c(.001,.1,1,100)
    b <- makebasis(W)
    p <- apply.metric.graph.lowpass.path(b,Y,times)
    for (iteration in 1:6) {
        perm <- sample.int(9)
        r <- apply.metric.graph.lowpass.path(makebasis(W[perm,perm]),(Y+shift)[perm,],times)
        expect_equal(r$gcv.scores,p$gcv.scores,tolerance=1e-10)
        expect_equal(r$effective.df,p$effective.df,tolerance=1e-10)
        # Absolute fit error after removing the baseline allows rounding at 2^35.
        returned <- sweep(r$fitted.values,1,shift[perm],'-')
        expect_lt(max(abs(returned-p$fitted.values[perm,,,drop=FALSE])),1e-5)
        expect_identical(r$tuning[[3]]$status,'unresolved_flat_gcv')
        expect_true(all(r$gcv.scores[,3]==0))
        expect_true(all(r$numerical$diagnostics$gcv.final.scaled.error<=1))
        expect_identical(vapply(r$tuning,`[[`,'','status'),vapply(p$tuning,`[[`,'','status'))
    }
})

test_that('unit restoration preserves checked scores including very small and large finite scales', {
    a <- list(2L,c(1L,3L),c(2L,4L),3L)
    w <- list(1,c(1,2),c(2,3),3)
    b <- metric.graph.lowpass.basis(a,w,weight.type='conductance',spectrum='full',
                                  numerical.policy='verified.fallback')
    y <- c(-1,NA,2,4)
    p <- apply.metric.graph.lowpass.path(b,y,c(.01,.1,1))
    for (scale in c(1e-120,1e-20,1e20,1e120)) {
        r <- apply.metric.graph.lowpass.path(b,y*scale,p$eta.grid)
        expect_equal(r$gcv.scores/scale/scale,p$gcv.scores,tolerance=1e-10)
        expect_equal(r$fitted.values/scale,p$fitted.values,tolerance=1e-10)
        expect_true(all(r$numerical$diagnostics$gcv.final.scaled.error<=1))
    }
})
