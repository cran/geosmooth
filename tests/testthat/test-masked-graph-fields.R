field_operator <- function() {
    metric.graph.lowpass.operator(list(2L, c(1L,3L), 2L, integer()),
        list(1, c(1,1), 1, numeric()))
}
test_that("masked solve agrees with analytic solution and leaves unsupported components NA", {
    op <- field_operator()
    fit <- fit.masked.graph.field(op, c(1,3), c(2,8), lambda = 0.5)
    L <- as.matrix(op$laplacian$matrix)[1:3,1:3]
    reference <- solve(diag(c(.5,0,.5)) + .5*L, c(1,0,4))
    expect_equal(fit$fitted[1:3], as.numeric(reference), tolerance=1e-10)
    expect_equal(fit$fitted[2], 5, tolerance=1e-10)
    expect_true(is.na(fit$fitted[4]))
    expect_lt(fit$relative.residual, 1e-10)
    constant <- fit.masked.graph.field(op, 2, 7, 1)
    expect_equal(constant$fitted[1:3], rep(7,3), tolerance=1e-10)
})
test_that("repeated vertex observations retain weighted evidence", {
    fit <- fit.masked.graph.field(field_operator(), c(1,1), c(2,8), 1, c(1,3))
    expect_equal(fit$fitted[1:3], rep(6.5,3), tolerance=1e-10)
})
test_that("CV excludes full participants and uses training-only constant predictions", {
    group <- rep(letters[1:12], each=2)
    vertex <- rep(c(1,3),12); y <- seq_along(vertex)
    set.seed(91); before <- .Random.seed
    cv <- cv.masked.graph.field(field_operator(), vertex, y, group,
        lambdas=c(.01,1), outer.folds=3, inner.folds=2, seed=12)
    expect_identical(.Random.seed, before)
    expect_true(all(vapply(split(cv$fold.id, group), function(x) length(unique(x))==1, logical(1))))
    for (fold in unique(cv$fold.id)) {
        train <- cv$fold.id != fold
        expect_equal(unique(cv$baseline[!train]), mean(y[train]))
    }
    expect_equal(cv$coverage, 1)
    expect_true(all(is.na(cv$final$fitted[4])))
    changed <- y; held <- cv$fold.id == 1; changed[held] <- changed[held] + 10000
    # Fix seed as well as folds to ensure identical inner splits.
    cv3 <- cv.masked.graph.field(field_operator(), vertex, changed, group,
        lambdas=c(.01,1), inner.folds=2, seed=12, fold.id=cv$fold.id)
    expect_equal(cv$predicted[held], cv3$predicted[held])
})
test_that("unsupported held-out records remain visible in coverage", {
    group <- letters[1:12]; vertex <- c(rep(1,11),4)
    cv <- cv.masked.graph.field(field_operator(), vertex, seq_len(12), group,
        lambdas=c(.01,1), outer.folds=3, inner.folds=2)
    expect_true(is.na(cv$predicted[12]))
    expect_equal(cv$coverage, 11/12)
})
test_that("graph parallelism reproduces serial predictions", {
    skip_on_os("windows")
    skip_if(identical(Sys.getenv("GEOSMOOTH_TEST_CORES"), "1"),
        "The parallelism comparison requires two worker processes")
    studies <- list(one=list(vertex=rep(c(1,3),6), y=1:12, group=letters[1:12]))
    ops <- list(a=field_operator(), b=field_operator())
    a <- compare.masked.graph.fields(ops, studies, workers=1,
        lambdas=c(.01,1), outer.folds=3, inner.folds=2)
    b <- compare.masked.graph.fields(ops, studies, workers=2,
        lambdas=c(.01,1), outer.folds=3, inner.folds=2)
    expect_equal(a,b)
})
test_that("invalid observations and leakage folds fail explicitly", {
    expect_error(fit.masked.graph.field(field_operator(),1,NA,1),"finite")
    expect_error(fit.masked.graph.field(field_operator(),1,1,0),"lambda")
    expect_error(cv.masked.graph.field(field_operator(),c(1,1),c(1,2),c("a","a"),
        fold.id=c(1,2)),"fold.id")
})

test_that("nested graph selection uses training-fold evidence", {
    studies <- list(HMP=list(vertex=rep(c(1,3),6), y=1:12, group=letters[1:12]))
    fits <- compare.masked.graph.fields(list(a=field_operator(),b=field_operator()),
        studies, lambdas=c(.01,1), outer.folds=3, inner.folds=2)
    selected <- select.masked.graph.field(fits, "HMP")
    expect_equal(selected$predicted, fits$a$HMP$predicted)
    expect_identical(selected$final.graph, "a")
    broken <- fits; broken$b$HMP$fold.id <- rev(broken$b$HMP$fold.id)
    expect_error(select.masked.graph.field(broken,"HMP"),"identical")
})

test_that("complete observations recover the full-spectrum Tikhonov solution", {
    op <- metric.graph.lowpass.operator(list(2L,c(1L,3L),2L),list(1,c(1,1),1))
    L <- as.matrix(op$laplacian$matrix); eig <- eigen(L,symmetric=TRUE)
    y <- c(1,9,3); lambda <- .2
    reference <- eig$vectors %*% ((crossprod(eig$vectors,y))/(1+3*lambda*eig$values))
    fit <- fit.masked.graph.field(op,1:3,y,lambda)
    expect_equal(fit$fitted,as.numeric(reference),tolerance=1e-10)
    normalized <- metric.graph.lowpass.operator(list(2L,c(1L,3L),2L),
        list(1,c(1,1),1),laplacian.type="symmetric.normalized")
    expect_error(fit.masked.graph.field(normalized,1,2,1),"unnormalized")
})

test_that("support summaries count people rather than repeated specimens", {
    support <- support.masked.graph.field(field_operator(),c(1,1,3),c('a','a','b'))
    expect_equal(support$nearest.observation.hops,c(0,1,0,Inf))
    expect_equal(support$local.participants,c(1L,2L,1L,0L))
    expect_equal(support$component.supported,c(TRUE,TRUE,TRUE,FALSE))
})
