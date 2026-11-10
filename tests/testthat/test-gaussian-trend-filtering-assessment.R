gaussian_assessment_fixture <- function(backend = "banded") {
    x <- seq(0,1,length.out=36)
    fit.gaussian.trend.filtering(x, sin(2*pi*x)+.15*cos(seq_along(x)),
        c(Inf,.01),backend=backend)
}

test_that("Gaussian permutation statistics, optional tests and Holm adjustment reconstruct", {
    fit <- gaussian_assessment_fixture()
    a <- assess.gaussian.trend.filtering(fit,permutations=9,inner.folds=3,
        outer.folds=3,seed=42,additional.tests=c("total.variation","absolute.deviation"))
    expect_identical(a$status,"complete")
    expect_equal(a$null.mean,mean(fit$y))
    expect_identical(a$settings$order,2L)
    meanfit <- as.numeric(fitted(a$fit))
    expected <- max(0,sum((fit$y-mean(fit$y))^2)-sum((fit$y-meanfit)^2))
    expect_equal(a$statistic,expected,tolerance=1e-10)
    expect_equal(unname(a$curve.summaries["absolute.deviation"]),mean(abs(meanfit-mean(fit$y))))
    grid <- seq(0,1,length.out=10001)
    expect_equal(unname(a$curve.summaries["total.variation"]),
                 sum(abs(diff(as.numeric(predict(a$fit,grid))))),tolerance=1e-6)
    for (b in seq_len(9)) {
        set.seed(a$seeds$permutations[b]); yy <- fit$y[sample.int(length(fit$y))]
        refit <- fit.gaussian.trend.filtering(fit$x,yy,fit$lambda,
            backend=fit$backend,lambda.selection="cv",foldid=a$folds$tuning)
        pred <- as.numeric(fitted(refit))
        expect_equal(a$permutations$statistic[b],
            max(0,sum((yy-mean(yy))^2)-sum((yy-pred)^2)),tolerance=1e-10)
        expect_equal(a$permutation.statistics$absolute.deviation[b],mean(abs(pred-mean(yy))))
        expect_equal(a$permutation.statistics$total.variation[b],
                     sum(abs(diff(as.numeric(predict(refit,grid))))),tolerance=2e-6)
        expect_equal(a$permutation.diagnostics[[b]]$cv.table,refit$cv.table)
    }
    for (j in seq_len(nrow(a$tests))) {
        tt <- a$tests[j,]; values <- a$permutation.statistics[[tt$test]]
        p <- (1+sum(values >= tt$statistic-tt$comparison.tolerance))/10
        expect_equal(tt$p.value,p)
        expect_equal(c(tt$p.value.lower,tt$p.value.upper),rep(p,2))
    }
    expect_equal(a$tests$p.value.holm,p.adjust(a$tests$p.value,"holm"))
    expect_equal(a$p.value,a$tests$p.value[1])
    expect_output(print(a),"Residual sum-of-squares reduction")
    expect_invisible(print(a))
})

test_that("Gaussian nested CV preserves outer outcome isolation and scores all rows", {
    fit <- gaussian_assessment_fixture()
    a <- assess.gaussian.trend.filtering(fit,permutations=1,inner.folds=3,outer.folds=3,seed=42)
    p <- a$prediction$predictions
    expect_identical(a$prediction$status,"complete")
    expect_equal(a$prediction$scores$mse[2],mean((fit$y-p$mean)^2))
    expect_equal(a$prediction$scores$mae[2],mean(abs(fit$y-p$mean)))
    expect_equal(unname(a$prediction$improvement),
        c(mean(p$constant.squared.error-p$squared.error),
          mean(p$constant.absolute.error-p$absolute.error)))
    for (k in 1:3) {
        train <- a$folds$outer!=k
        expect_equal(unique(p$constant.mean[!train]),mean(fit$y[train]))
        expect_identical(as.integer(names(a$folds$inner[[k]])),which(train))
        ff <- fit.gaussian.trend.filtering(fit$x[train],fit$y[train],fit$lambda * mean(train),
            lambda.selection="cv",foldid=unname(a$folds$inner[[k]]),backend=fit$backend)
        pred <- predict(ff,fit$x[!train])
        expect_equal(p$mean[!train],as.numeric(pred))
        expect_identical(p$outside.range[!train],attr(pred,"outside.range"))
        expect_equal(a$prediction$diagnostics[[k]]$cv.table,ff$cv.table)
    }
    held <- a$folds$outer==1
    changed <- fit; changed$y[held] <- changed$y[held]+5
    changed$foldid <- rank(changed$y); changed$selected.index <- 1L
    b <- assess.gaussian.trend.filtering(changed,permutations=1,inner.folds=3,outer.folds=3,seed=42)
    expect_identical(a$folds,b$folds)
    expect_identical(p$mean[held],b$prediction$predictions$mean[held])
    expect_identical(a$prediction$diagnostics[[1]],b$prediction$diagnostics[[1]])
})

test_that("Gaussian summaries include interior turns and empirical tie frequencies", {
    x <- c(0,.1,.15,.75,1)
    fit <- fit.gaussian.trend.filtering(x,1-20*(x-.4)^2,Inf,backend="banded")
    expect_equal(.gaussian.tf.assessment.variation(fit),2-(-2.2)-(-6.2),tolerance=1e-10)
    expect_gt(.gaussian.tf.assessment.variation(fit),sum(abs(diff(as.numeric(fitted(fit))))))
    for (y in list(rep(3,5),2+4*x)) {
        ff <- fit.gaussian.trend.filtering(x,y,Inf,backend="banded")
        expect_equal(.gaussian.tf.assessment.variation(ff),abs(tail(y,1)-head(y,1)),tolerance=1e-10)
    }
    x <- rep(seq(0,1,length.out=12),c(8,rep(2,11)))
    fit <- fit.gaussian.trend.filtering(x,x^2+.1*cos(seq_along(x)),c(Inf,.01),backend="banded")
    a <- assess.gaussian.trend.filtering(fit,permutations=1,inner.folds=3,seed=11,predictive=FALSE)
    expected <- mean(abs(as.numeric(fitted(a$fit))-mean(fit$y)))
    uniform <- mean(abs(as.numeric(predict(a$fit,unique(x)))-mean(fit$y)))
    expect_equal(unname(a$curve.summaries["absolute.deviation"]),expected)
    expect_gt(abs(expected-uniform),1e-4)
})

test_that("RNG replay, optional tests and prediction preserve Gaussian randomizations", {
    fit <- gaussian_assessment_fixture()
    set.seed(729); previous <- .Random.seed
    a <- assess.gaussian.trend.filtering(fit,permutations=3,inner.folds=3,outer.folds=3,seed=42)
    expect_identical(.Random.seed,previous)
    b <- assess.gaussian.trend.filtering(fit,permutations=3,inner.folds=3,seed=42,
        predictive=FALSE,additional.tests=c("absolute.deviation","total.variation"))
    expect_identical(.Random.seed,previous)
    expect_identical(a$permutations,b$permutations)
    expect_identical(a$fit,b$fit)
    expect_identical(a$p.value,b$p.value)
    expect_identical(a$curve.summaries,b$curve.summaries)
    expect_null(b$prediction)
    expect_identical(a,assess.gaussian.trend.filtering(fit,permutations=3,
        inner.folds=3,outer.folds=3,seed=42))
    local({
        on.exit(assign(".Random.seed",previous,envir=.GlobalEnv))
        rm(".Random.seed",envir=.GlobalEnv)
        invisible(assess.gaussian.trend.filtering(fit,permutations=1,inner.folds=3,predictive=FALSE))
        expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
    })
})

test_that("constant outcomes give exact ties and complete p-values of one", {
    for (backend in c("dense","banded")) {
        fit <- fit.gaussian.trend.filtering(seq(0,1,length.out=24),rep(3,24),c(Inf,.01),backend=backend)
        a <- assess.gaussian.trend.filtering(fit,permutations=3,inner.folds=3,outer.folds=3,
            additional.tests=c("absolute.deviation","total.variation"))
        expect_identical(a$status,"complete")
        expect_equal(a$statistic,0)
        expect_equal(a$tests$p.value,rep(1,3))
        expect_equal(a$tests$p.value.holm,rep(1,3))
        expect_equal(unname(a$curve.summaries),c(0,0),tolerance=1e-10)
        expect_equal(unname(a$prediction$improvement),c(0,0),tolerance=1e-10)
    }
})

test_that("Gaussian failures retain missing permutations and suppress incomplete prediction scores", {
    fit <- gaussian_assessment_fixture()
    original <- .gaussian.tf.assessment.refit; calls <- 0L
    local_mocked_bindings(.gaussian.tf.assessment.refit=function(object,y,rows,folds,summarize=FALSE) {
        calls <<- calls+1L
        if (calls %in% c(2L,6L)) object$control$max.iterations <- 1L
        original(object,y,rows,folds,summarize)
    },.package="geosmooth")
    expect_warning(a <- assess.gaussian.trend.filtering(fit,permutations=4,
        inner.folds=3,outer.folds=3,seed=42,
        additional.tests=c("absolute.deviation","total.variation")),"failed fits")
    expect_identical(a$status,"incomplete_permutations")
    expect_equal(a$permutation.summary$failed,1)
    expect_equal(nrow(a$permutations),4)
    expect_true(is.na(a$p.value))
    r <- sum(a$permutations$statistic >= a$statistic-a$permutation.summary$comparison.tolerance,na.rm=TRUE)
    expect_equal(unname(a$p.value.bounds),c(1+r,2+r)/5)
    expect_true(all(is.na(a$tests$p.value.holm)))
    expect_identical(a$prediction$status,"incomplete")
    expect_true(all(is.na(a$prediction$improvement)))
    expect_true(all(is.na(a$prediction$scores[2,c("mse","mae")])))
    expect_true(any(is.finite(a$prediction$predictions$mean)))
    fit$control$max.iterations <- 1L
    b <- assess.gaussian.trend.filtering(fit,permutations=1,inner.folds=3,predictive=FALSE,warn=FALSE)
    expect_identical(b$status,"observed_fit_failed")
    expect_true(all(is.na(b$p.value.bounds)))
})

test_that("summary failures do not erase the primary test or shrink the Holm family", {
    local_mocked_bindings(.gaussian.tf.assessment.variation=function(fit) stop("Deliberate variation failure"),
                         .package="geosmooth")
    a <- assess.gaussian.trend.filtering(gaussian_assessment_fixture(),permutations=2,
        inner.folds=3,predictive=FALSE,warn=FALSE,additional.tests="total.variation")
    expect_identical(a$status,"complete")
    expect_true(is.finite(a$p.value))
    expect_true(is.na(a$tests$p.value[2]))
    expect_true(all(is.na(a$tests$p.value.holm)))
    expect_match(a$curve.summary.details$message,"Deliberate variation failure")
})

test_that("Gaussian assessment rejects unsupported nulls and invalid arguments", {
    fit <- gaussian_assessment_fixture()
    for (field in c("weights","covariates","offset","offset.used")) {
        changed <- fit
        changed[[field]] <- switch(field,weights=rep(2,36),covariates=matrix(fit$x),
                                  offset=rep(.1,36),offset.used=TRUE)
        expect_error(assess.gaussian.trend.filtering(changed,permutations=1),"unit weights, no covariates and no offset")
    }
    changed <- fit; changed$order <- 1L
    expect_error(assess.gaussian.trend.filtering(changed),"order=2L")
    changed <- fit; changed$y[1] <- NA
    expect_error(assess.gaussian.trend.filtering(changed),"finite aligned")
    expect_error(assess.gaussian.trend.filtering(NULL),"object must")
    for (value in list(0,-1,NA,Inf,1.5,"2",c(1,2)))
        expect_error(assess.gaussian.trend.filtering(fit,permutations=value),"permutations")
    expect_error(assess.gaussian.trend.filtering(fit,inner.folds=1),"inner.folds")
    expect_error(assess.gaussian.trend.filtering(fit,outer.folds=100),"Fold counts")
    expect_error(assess.gaussian.trend.filtering(fit,predictive=NA),"predictive")
    expect_error(assess.gaussian.trend.filtering(fit,seed=-1),"seed")
    expect_error(assess.gaussian.trend.filtering(fit,additional.tests="likelihood"),"additional.tests")
})

test_that("dense and banded Gaussian assessments agree", {
    assess <- function(backend) assess.gaussian.trend.filtering(gaussian_assessment_fixture(backend),
        permutations=3,inner.folds=3,outer.folds=3,seed=42,
        additional.tests=c("absolute.deviation","total.variation"))
    a <- assess("banded"); b <- assess("dense")
    expect_identical(a$status,"complete"); expect_identical(b$status,"complete")
    expect_identical(a$folds,b$folds)
    expect_equal(a$permutation.statistics,b$permutation.statistics,tolerance=2e-6)
    expect_equal(a$curve.summaries,b$curve.summaries,tolerance=2e-6)
    expect_equal(a$tests$p.value,b$tests$p.value)
    expect_equal(a$prediction$predictions$mean,b$prediction$predictions$mean,tolerance=2e-6)
})

test_that("a finite exchangeable null has conservative exact rejection rates", {
    # All 90 distinct assignments of {-1,-1,0,0,1,1} have equal probability.
    # This checks a small conditional null, not broad calibration or power.
    x <- seq(0,1,length.out=6); arrangements <- list()
    for (minus in combn(6,2,simplify=FALSE)) for (zero in combn(setdiff(1:6,minus),2,simplify=FALSE)) {
        y <- rep(1,6); y[minus] <- -1; y[zero] <- 0
        arrangements[[length(arrangements)+1L]] <- y
    }
    statistics <- vapply(arrangements,function(y) {
        fit <- fit.gaussian.trend.filtering(x,y,c(Inf,.01),backend="banded")
        a <- assess.gaussian.trend.filtering(fit,permutations=1,inner.folds=6,
            predictive=FALSE,seed=17,warn=FALSE)
        expect_identical(a$status,"complete")
        c(squared.error=a$statistic,a$curve.summaries)
    },numeric(3))
    exact <- t(apply(statistics,1,function(values) vapply(values,function(v)
        mean(values>=v-1e-10*max(1,abs(v))),numeric(1))))
    holm <- apply(exact,2,p.adjust,method="holm")
    for (alpha in c(.01,.05,.1,.2)) {
        expect_true(all(rowMeans(exact<=alpha)<=alpha))
        expect_lte(mean(colSums(holm<=alpha)>0),alpha)
    }
})
