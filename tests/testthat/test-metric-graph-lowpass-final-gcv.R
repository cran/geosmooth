# Analytic invariance tests: no independent eigensolver is needed to know that
# componentwise constants lie in the unnormalized Laplacian's nullspace.
final_gcv_graph <- function(bridge = 1e-8) {
    list(adj=list(2L,c(1L,3L),c(2L,4L),3L),
         w=list(1,c(1,bridge),c(bridge,1),1))
}
final_gcv_basis <- function(g) {
    metric.graph.lowpass.basis(g$adj,g$w,weight.type="conductance",spectrum="full",
                              numerical.policy="verified.fallback")
}
final_gcv_fit <- function(g,y,grid) {
    fit.metric.graph.lowpass(g$adj,g$w,y,weight.type="conductance",spectrum="full",
                            numerical.policy="verified.fallback",eta.grid=grid)
}

test_that("final scores and tuning are invariant to additive shifts on weak bridges", {
    grid <- c(.01,.1,1,1000)
    for (bridge in c(1e-3,1e-6,1e-8)) {
        g <- final_gcv_graph(bridge); b <- final_gcv_basis(g)
        baseline <- apply.metric.graph.lowpass.path(b,c(0,NA,NA,1),grid)
        # A fixed-time seed allows testing automatically tuned refit independently
        # of whether the original response can identify a smoothing time.
        seed <- final_gcv_fit(g,c(0,NA,NA,1),grid[1])
        for (shift in c(-1e6,0,10,1000,1e6)) {
            y <- c(shift,NA,NA,shift+1)
            path <- apply.metric.graph.lowpass.path(b,y,grid)
            # With two observations on a connected graph, GCV=(y1-y2)^2.
            expect_equal(path$gcv.scores,rep(1,4),tolerance=1.001e-7)
            expect_equal(path$gcv.scores,baseline$gcv.scores,tolerance=1e-12)
            expect_identical(path$tuning[[1]]$status,"unresolved_flat_gcv")
            expect_false(path$tuning[[1]]$bracketed)
            expect_true(is.na(path$tuning[[1]]$best.idx))
            expect_error(final_gcv_fit(g,y,grid),"GCV is flat")
            expect_error(refit(seed,y,per.column.gcv=TRUE,eta.grid=grid),"GCV is flat")
            fixed <- refit(seed,y)
            expect_equal(fixed$gcv.scores,1,tolerance=1.001e-7)
            expect_equal(fixed$fitted.values-shift,baseline$fitted.values[,1],tolerance=2e-10)
            expect_true(all(path$numerical$diagnostics$gcv.final.scaled.error<=1))
        }
    }
})

test_that("componentwise centering handles disconnected and constant responses", {
    g <- list(adj=list(2L,c(1L,3L),2L,5L,c(4L,6L),5L,integer()),
              w=list(1,c(1,2),2,3,c(3,4),4,numeric()))
    b <- final_gcv_basis(g); grid <- c(.01,.1,1)
    y <- c(1,NA,3,4,-2,NA,7)
    shifts <- c(rep(1e6,3),rep(-1000,3),10)
    base <- apply.metric.graph.lowpass.path(b,y,grid)
    actual <- apply.metric.graph.lowpass.path(b,y+shifts,grid)
    expect_equal(actual$gcv.scores,base$gcv.scores,tolerance=1e-12)
    expect_equal(actual$fitted.values-shifts,base$fitted.values,tolerance=2e-10)
    expect_equal(actual$tuning,base$tuning)
    seed <- final_gcv_fit(g,y,grid)
    tuned <- refit(seed,y+shifts,per.column.gcv=TRUE,eta.grid=grid)
    expect_equal(tuned$eta.optimal,seed$gcv$eta.optimal)
    expect_equal(tuned$gcv.scores,seed$gcv$gcv.optimal,tolerance=1e-12)
    for (constant in list(rep(1000,7),c(1e6,NA,1e6,-1000,-1000,NA,10))) {
        p <- apply.metric.graph.lowpass.path(b,constant,grid)
        expect_identical(p$gcv.scores,rep(0,3))
        expect_identical(p$tuning[[1]]$status,"unresolved_flat_gcv")
        expect_error(final_gcv_fit(g,constant,grid),"GCV is flat")
        expect_error(refit(seed,constant,per.column.gcv=TRUE,eta.grid=grid),"GCV is flat")
        expected <- if (anyNA(constant)) c(rep(1e6,3),rep(-1000,3),10) else constant
        expect_equal(refit(seed,constant)$fitted.values,expected,tolerance=0)
    }
})

test_that("final score checks cover restored units for nonflat multivariate responses", {
    g <- final_gcv_graph(1);b <- final_gcv_basis(g);grid <- c(.02,.2,2)
    Y <- cbind(c(2,NA,-1,4),c(1,-3,5,7))
    base <- apply.metric.graph.lowpass.path(b,Y,grid)
    shifted <- apply.metric.graph.lowpass.path(b,sweep(Y,2,c(1000,-1e6),`+`),grid)
    expect_equal(shifted$gcv.scores,base$gcv.scores,tolerance=1e-12)
    expect_equal(shifted$tuning,base$tuning)
    for (scale in c(.01,1e4)) {
        p <- apply.metric.graph.lowpass.path(b,Y*scale,grid)
        expect_equal(p$gcv.scores/scale^2,base$gcv.scores,tolerance=1e-10)
        expect_true(all(p$numerical$diagnostics$gcv.final.scaled.error<=1))
    }
    # Fault-inject only the final-unit check: a passing standardized check is
    # insufficient, even if every other diagnostic is zero.
    limits <- geosmooth:::.metric.graph.verified.limits()
    checks <- as.list(limits*0); checks$gcv.final.scaled.error <- 1.01
    expect_error(geosmooth:::.metric.graph.verify.agreement(checks),"gcv.final.scaled.error",fixed=TRUE)
})
