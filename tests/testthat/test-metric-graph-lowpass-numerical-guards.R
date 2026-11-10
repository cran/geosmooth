weighted_path <- function(w) {
    n <- length(w) + 1L
    a <- lapply(seq_len(n), function(i) as.integer(c(i-1L,i+1L)[c(i>1L,i<n)]))
    list(a = a, w = lapply(seq_len(n), function(i) w[pmin(i,a[[i]])]))
}
weighted_basis <- function(w, ...) {
    g <- weighted_path(w)
    metric.graph.lowpass.basis(g$a,g$w,weight.type="conductance",spectrum="full",...)
}
weighted_fit <- function(w, y, ...) {
    g <- weighted_path(w)
    fit.metric.graph.lowpass(g$a,g$w,y,weight.type="conductance",spectrum="full",...)
}

test_that("weak bridges cannot return falsely verified heat fits", {
    for (e in c(1e-14,1e-16)) {
        expect_error(weighted_fit(c(1,e,1), rep(7,4),eta.grid=1/e),
                     "Unresolved|unresolved")
        expect_error(weighted_fit(c(1,e,1), c(0,0,1,1),eta.grid=1/e),
                     "Unresolved|unresolved")
    }
    # A resolved eigenvalue can still be inadequate at a long requested time.
    b <- weighted_basis(c(1,1e-10,1))
    expect_error(apply.metric.graph.lowpass.path(b,c(0,0,1,1),1e10),"unresolved")
    expect_error(weighted_fit(c(1,1e-12,1),c(7,7,NA,NA),eta.grid=.1),
                 "ill-conditioned|unresolved")
})

test_that("known component nullspace is preserved and analytic heat is accurate", {
    e <- 1e-3
    b <- weighted_basis(c(1,e,1))
    expect_equal(b$spectral$eigenvectors[,1], rep(.5,4), tolerance=1e-15)
    path <- apply.metric.graph.lowpass.path(b,rep(7,4),c(.1,100,1e5))
    expect_equal(path$fitted.values,matrix(7,4,3),ignore_attr=TRUE,tolerance=1e-12)
    # Symmetry reduction on antisymmetric pairs, evaluated with a stable root.
    slow <- e - e^2/(sqrt(1+e^2)+1)
    fast <- 2 + e + e^2/(sqrt(1+e^2)+1)
    B <- matrix(c(1,-1,-1,1+2*e),2)
    y <- c(0,0,1,1); z <- c(-.5,-.5); t <- 1/e
    h <- exp(-t*slow)*(fast*diag(2)-B)/(fast-slow) +
         exp(-t*fast)*(B-slow*diag(2))/(fast-slow)
    u <- as.vector(h%*%z)
    reference <- c(.5+u,.5-rev(u))
    actual <- apply.metric.graph.lowpass.path(b,y,t)$fitted.values[,1]
    expect_equal(actual,reference,tolerance=1e-10)
    expect_equal(weighted_fit(c(1,.01,1),c(7,7,NA,NA),eta.grid=.1)$fitted.values,
                 rep(7,4),tolerance=1e-10)
    # No response-dependent constant special case: all fits use the same operator.
    expect_equal(b$spectral$nullity,1L)
    expect_true(all(b$spectral$eigenvalues[-1]>0))
})

test_that("automatic grids and metadata respect conductance scaling", {
    y <- c(2,NA,4,NA,-1)
    base <- weighted_fit(1:4,y)
    b0 <- weighted_basis(1:4)
    g0 <- metric.graph.heat.eta.grid(b0)
    for (factor in c(1e-100,1e-20,1e20,1e100)) {
        fit <- weighted_fit(factor*(1:4),y)
        b <- weighted_basis(factor*(1:4))
        expect_equal(fit$fitted.values,base$fitted.values,tolerance=1e-10)
        expect_equal(fit$gcv$gcv.scores,base$gcv$gcv.scores,tolerance=1e-10)
        expect_equal(fit$gcv$eta.grid*factor,base$gcv$eta.grid,tolerance=1e-10)
        expect_equal(as.numeric(metric.graph.heat.eta.grid(b))*factor,as.numeric(g0),tolerance=1e-10)
        expect_equal(b$spectral$smallest.positive.eigenvalue/factor,
                     b0$spectral$smallest.positive.eigenvalue,tolerance=1e-12)
    }
})

test_that("boundary flags use valid competitive positive times on both ends", {
    low <- weighted_fit(1:4,0:4,eta.grid=c(0,.0001,.1,1))
    expect_true(low$gcv$search$rounds$lower.endpoint.selected)
    expect_false(low$gcv$search$rounds$upper.endpoint.selected)
    high <- weighted_fit(1:4,c(1,-1,1,-1,1),eta.grid=c(.0001,.001,.01),eta.search="guarded.gcv")
    expect_equal(high$gcv$search$status,"upper_endpoint_selected")
    expect_true(tail(high$gcv$search$rounds$upper.endpoint.selected,1))
    expect_false(tail(high$gcv$search$rounds$lower.endpoint.selected,1))
})
