verified_graph <- function(w = 1:4) {
    n <- length(w) + 1L
    a <- lapply(seq_len(n), function(i) as.integer(c(i-1L,i+1L)[c(i>1L,i<n)]))
    list(a=a,w=lapply(seq_len(n),function(i) w[pmin(i,a[[i]])]))
}
verified_basis <- function(w = 1:4) {
    g <- verified_graph(w)
    metric.graph.lowpass.basis(g$a,g$w,weight.type="conductance",spectrum="full",
                              numerical.policy="verified.fallback")
}
verified_fit <- function(y, grid=c(.01,.1,1), w=1:4, ...) {
    g <- verified_graph(w)
    fit.metric.graph.lowpass(g$a,g$w,y,weight.type="conductance",spectrum="full",
                            numerical.policy="verified.fallback",eta.grid=grid,...)
}
verified_penalty_reference <- function(L,y,time) {
    ev <- eigen(L,symmetric=TRUE)
    penalty <- tcrossprod(sweep(ev$vectors,2,expm1(time*pmax(ev$values,0)),`*`),ev$vectors)
    O <- which(!is.na(y)); R <- diag(length(y))[O,,drop=FALSE]
    T <- solve(crossprod(R)+penalty,t(R)); fitted <- as.vector(T%*%y[O])
    df <- sum(diag(T[O,,drop=FALSE]))
    list(fitted=fitted,df=df,gcv=length(O)*sum((y[O]-fitted[O])^2)/(length(O)-df)^2)
}

test_that("every declared agreement criterion is enforced including solve residuals", {
    limits <- geosmooth:::.metric.graph.verified.limits()
    zero <- as.list(limits*0)
    expect_silent(geosmooth:::.metric.graph.verify.agreement(zero))
    for (name in names(limits)) {
        bad <- zero; bad[[name]] <- 2*limits[[name]]
        expect_error(geosmooth:::.metric.graph.verify.agreement(bad),name,fixed=TRUE)
        bad[[name]] <- NaN
        expect_error(geosmooth:::.metric.graph.verify.agreement(bad),"failed numerical agreement")
    }
})

test_that("verified fits agree with direct heat penalty, with mask grouping", {
    b <- verified_basis()
    Y <- cbind(c(2,NA,4,NA,-1), c(1,NA,3,NA,8), 1:5, c(NA,3,NA,4,1))
    times <- c(.01,.1,.3)
    path <- apply.metric.graph.lowpass.path(b,Y,times)
    expect_equal(path$numerical$response.mask,c(1L,1L,2L,3L))
    expect_equal(nrow(path$numerical$diagnostics),9L)
    for (j in seq_along(times)) for (k in seq_len(ncol(Y))) {
        r <- verified_penalty_reference(as.matrix(b$operator$laplacian$matrix),Y[,k],times[j])
        expect_equal(path$fitted.values[,j,k],r$fitted,tolerance=1e-9)
        expect_equal(path$effective.df[j,k],r$df,tolerance=1e-9)
        expect_equal(path$gcv.scores[j,k],r$gcv,tolerance=1e-8)
    }
    # Perturb the independent calculation. A successful Cholesky solve alone
    # must not make the candidate acceptable.
    old <- geosmooth:::.metric.graph.verified.reference
    testthat::local_mocked_bindings(.metric.graph.verified.reference=function(basis) {
        ref <- old(basis); ref$eigenvalues <- ref$eigenvalues*1.1; ref
    }, .package="geosmooth")
    expect_error(apply.metric.graph.lowpass.path(b,Y,times),"failed numerical agreement")
})

test_that("identity, tiny times, equilibrium and global scale remain well defined", {
    b <- verified_basis()
    Y <- cbind(c(2,NA,4,NA,-1),rep(7,5))
    grid <- c(0,1e-200,1e-20,.1,1e6)
    path <- apply.metric.graph.lowpass.path(b,Y,grid)
    expect_equal(path$fitted.values[,1,],Y)
    expect_true(all(is.na(path$gcv.scores[1,])))
    expect_true(all(path$status[1,]=="identity_undefined_gcv"))
    expect_equal(path$gcv.scores[2,],path$gcv.scores[3,],tolerance=1e-12)
    expect_equal(path$fitted.values[,5,1],rep(5/3,5),tolerance=1e-10)
    expect_equal(path$fitted.values[,2:5,2],matrix(7,5,4),ignore_attr=TRUE,tolerance=1e-10)
    for (factor in c(1e-100,1e100)) {
        scaled <- apply.metric.graph.lowpass.path(verified_basis(factor*(1:4)),Y,grid/factor)
        expect_equal(scaled$fitted.values,path$fitted.values,ignore_attr=TRUE,tolerance=1e-9)
        expect_equal(scaled$gcv.scores,path$gcv.scores,tolerance=1e-9)
    }
})

test_that("flat tuning is unresolved and endpoints survive fit and refit", {
    b <- verified_basis()
    y <- c(2,NA,NA,NA,5)
    path <- apply.metric.graph.lowpass.path(b,y,c(.01,.1,1))
    expect_equal(path$tuning[[1]]$status,"unresolved_flat_gcv")
    expect_true(is.na(path$tuning[[1]]$best.idx))
    expect_error(verified_fit(y),"GCV is flat")
    fit <- verified_fit(0:4,c(0,.0001,.1,1))
    expect_true(fit$gcv$tuning$lower.endpoint)
    expect_false(fit$gcv$tuning$upper.endpoint)
    expect_false(fit$gcv$tuning$bracketed)
    high <- verified_fit(c(1,-1,1,-1,1),c(.0001,.001,.01))
    expect_true(high$gcv$tuning$upper.endpoint)
    rerun <- refit(fit,0:4,per.column.gcv=TRUE,eta.grid=c(0,.0001,.1,1))
    expect_equal(rerun$tuning[[1]],fit$gcv$tuning)
    expect_error(refit(fit,y,per.column.gcv=TRUE,eta.grid=c(.01,.1,1)),"GCV is flat")
    fixed <- refit(fit,y)
    expect_true(all(is.finite(fixed$fitted.values)))
    expect_equal(fixed$tuning[[1]]$status,"single_candidate")
    expect_true(fixed$tuning[[1]]$lower.endpoint && fixed$tuning[[1]]$upper.endpoint)
    expect_false(fixed$tuning[[1]]$bracketed)
    expect_error(verified_fit(0:4,eta.search="guarded.gcv"),"requires eta.search")
})

test_that("weak bridges retain spectral gate and must match analytic answers or refuse", {
    for (e in c(1e-14,1e-16)) {
        expect_error(verified_basis(c(1,e,1)),"Unresolved|unresolved")
    }
    b <- verified_basis(c(1,1e-10,1))
    expect_error(apply.metric.graph.lowpass.path(b,c(0,0,1,1),1e10,
                                               numerical.policy="conservative"),"unresolved")
    result <- try(apply.metric.graph.lowpass.path(b,c(0,0,1,1),1e10),silent=TRUE)
    if (!inherits(result,"try-error")) {
        e <- 1e-10; slow <- e-e^2/(sqrt(1+e^2)+1); fast <- 2+e+e^2/(sqrt(1+e^2)+1)
        B <- matrix(c(1,-1,-1,1+2*e),2)
        H <- exp(-1e10*slow)*(fast*diag(2)-B)/(fast-slow)
        u <- as.vector(H%*%c(-.5,-.5))
        expect_equal(result$fitted.values[,1],c(.5+u,.5-rev(u)),tolerance=1e-7)
    } else expect_match(as.character(result),"failed numerical agreement|unresolved|Unresolved")
    b <- verified_basis(c(1,1e-12,1))
    result <- try(apply.metric.graph.lowpass.path(b,c(7,7,NA,NA),.1),silent=TRUE)
    if (!inherits(result,"try-error")) {
        expect_equal(result$fitted.values[,1],rep(7,4),tolerance=1e-7)
    } else expect_match(as.character(result),"failed numerical agreement|unresolved|Unresolved")
})

test_that("disconnected components keep their observed-data identification rule", {
    adj <- list(2L,1L,4L,3L,integer())
    w <- list(1,1,2,2,numeric())
    b <- metric.graph.lowpass.basis(adj,w,weight.type="conductance",spectrum="full",
                                   numerical.policy="verified.fallback")
    path <- apply.metric.graph.lowpass.path(b,c(1,NA,3,NA,4),c(.1,1))
    expect_equal(path$fitted.values,matrix(c(1,1,3,3,4),5,2),ignore_attr=TRUE,tolerance=1e-10)
    expect_true(all(is.na(path$gcv.scores)))
    expect_error(apply.metric.graph.lowpass.path(b,c(1,NA,NA,NA,4),.1),"component needs")
    # Partial or normalized bases cannot silently enter the bounded fallback.
    g <- verified_graph()
    expect_error(metric.graph.lowpass.basis(g$a,g$w,numerical.policy="verified.fallback"),"requires spectrum")
    expect_error(metric.graph.lowpass.basis(g$a,g$w,spectrum="full",laplacian.type="symmetric.normalized",
                                           numerical.policy="verified.fallback"),"unnormalized")
    expect_error(apply.metric.graph.lowpass.path(verified_basis(),1:5,.1,filter.type="tikhonov"),"only heat_kernel")
})
