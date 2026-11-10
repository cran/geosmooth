assessment_fixture <- function() {
    fit.logistic.trend.filtering(seq(0, 1, length.out = 60),
        rep(c(0, 1, 0, 1, 1), 12), c(Inf, .01),
        backend = "banded", curve.ridge = .1)
}

test_that("the statistic and inclusive randomization p-value reproduce independently", {
    fit <- assessment_fixture()
    a <- assess.logistic.trend.filtering(fit, permutations = 9,
        inner.folds = 3, outer.folds = 3, seed = 42)
    expect_identical(a$status, "complete")
    loglik <- sum(dbinom(fit$y, 1, fitted(a$fit), log = TRUE))
    null <- sum(dbinom(fit$y, 1, mean(fit$y), log = TRUE))
    expect_equal(a$statistic, max(0, 2*(loglik-null)), tolerance = 1e-10)
    statistics <- vapply(a$seeds$permutations, function(seed) {
        set.seed(seed)
        y <- fit$y[sample.int(length(fit$y))]
        refit <- fit.logistic.trend.filtering(fit$x, y, fit$lambda,
            backend = fit$backend, curve.ridge = fit$curve.ridge,
            lambda.selection = "cv", foldid = a$folds$tuning, control = fit$control)
        max(0, 2*sum(dbinom(y, 1, fitted(refit), log = TRUE) -
                     dbinom(y, 1, mean(y), log = TRUE)))
    }, numeric(1))
    expect_equal(a$permutations$statistic, statistics, tolerance = 1e-10)
    tolerance <- a$permutation.summary$comparison.tolerance
    expected <- (1 + sum(statistics >= a$statistic-tolerance))/10
    expect_equal(a$p.value, expected)
    expect_equal(unname(a$p.value.bounds), rep(expected, 2))
    expect_gte(a$p.value, .1)
    expect_true(all(vapply(a$permutation.diagnostics, function(x)
        nrow(x$cv.table) == length(fit$lambda) && length(x$cv.fold.status) == 3,
        logical(1))))
    expect_output(print(a), "Permutation p-value")
    expect_invisible(print(a))
})

test_that("nested CV uses training-only baselines, tuning and transformations", {
    fit <- assessment_fixture()
    a <- assess.logistic.trend.filtering(fit, permutations = 1,
        inner.folds = 3, outer.folds = 3, seed = 42)
    p <- a$prediction$predictions
    expect_identical(a$prediction$status, "complete")
    expect_equal(a$prediction$scores$log.loss[2],
                 -mean(dbinom(fit$y, 1, p$probability, log = TRUE)), tolerance = 1e-10)
    expect_equal(a$prediction$scores$brier[2], mean((fit$y-p$probability)^2))
    expect_equal(unname(a$prediction$improvement),
        c(mean(p$constant.log.loss-p$log.loss), mean(p$constant.brier-p$brier)))
    for (k in 1:3) {
        train <- a$folds$outer != k
        expect_equal(unique(p$constant.probability[!train]), mean(fit$y[train]))
        inner <- a$folds$inner[[k]]
        expect_identical(as.integer(names(inner)), which(train))
        ff <- fit.logistic.trend.filtering(fit$x[train], fit$y[train], fit$lambda * mean(train),
            lambda.selection = "cv", foldid = unname(inner), backend = fit$backend,
            curve.ridge = fit$curve.ridge * mean(train), control = fit$control)
        expect_equal(as.numeric(predict(ff, fit$x[!train])), p$probability[!train], tolerance = 1e-10)
        expect_equal(ff$cv.table, a$prediction$diagnostics[[k]]$cv.table)
        expect_equal(range(fit$x[train]), a$prediction$diagnostics[[k]]$training.range)
    }
    held <- a$folds$outer == 1
    changed <- fit
    changed$y[held] <- 1-changed$y[held]
    # Even a stored, outcome-dependent fit fold assignment is ignored.
    changed$foldid <- changed$y + 1L
    b <- assess.logistic.trend.filtering(changed, permutations = 1,
        inner.folds = 3, outer.folds = 3, seed = 42)
    expect_identical(a$folds, b$folds)
    expect_equal(p$eta[held], b$prediction$predictions$eta[held], tolerance = 1e-12)
    expect_equal(a$prediction$diagnostics[[1]], b$prediction$diagnostics[[1]])
})

test_that("random state is preserved and predictive assessment does not alter the test", {
    fit <- assessment_fixture()
    set.seed(912)
    previous <- .Random.seed
    a <- assess.logistic.trend.filtering(fit, permutations = 3,
        inner.folds = 3, outer.folds = 3, seed = 42)
    expect_identical(.Random.seed, previous)
    b <- assess.logistic.trend.filtering(fit, permutations = 3,
        inner.folds = 3, seed = 42, predictive = FALSE)
    expect_identical(.Random.seed, previous)
    expect_identical(a$permutations, b$permutations)
    expect_identical(a$p.value, b$p.value)
    expect_null(b$prediction)
    expect_null(b$folds$outer)
    expect_identical(a, assess.logistic.trend.filtering(fit, permutations = 3,
        inner.folds = 3, outer.folds = 3, seed = 42))
    local({
        on.exit(assign(".Random.seed", previous, envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv)
        invisible(assess.logistic.trend.filtering(fit, permutations = 1,
            inner.folds = 3, seed = 42, predictive = FALSE))
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    })
})

test_that("failed randomizations and outer fits are retained without selective scoring", {
    fit <- fit.logistic.trend.filtering(seq(0, 1, length.out = 12),
        rep(c(0, 0, 0, 1), 3), c(Inf, .01), backend = "banded", curve.ridge = .1)
    expect_warning(a <- assess.logistic.trend.filtering(fit, permutations = 39,
        inner.folds = 3, outer.folds = 3, seed = 42), "failed fits")
    expect_identical(a$status, "incomplete_permutations")
    expect_true(is.na(a$p.value))
    expect_gt(a$permutation.summary$failed, 0)
    expect_equal(nrow(a$permutations), 39)
    r <- sum(a$permutations$statistic >= a$statistic -
             a$permutation.summary$comparison.tolerance, na.rm = TRUE)
    missing <- sum(!is.finite(a$permutations$statistic))
    expect_equal(unname(a$p.value.bounds), c(1+r, 1+r+missing)/40)
    expect_identical(a$prediction$status, "incomplete")
    expect_true(all(is.na(a$prediction$improvement)))
    expect_true(all(is.na(a$prediction$scores[2, c("log.loss", "brier")])))
    expect_equal(nrow(a$prediction$predictions), 12)
    expect_true(any(vapply(a$permutation.diagnostics, function(x)
        any(x$cv.fold.status != "complete"), logical(1))))
    broken <- assessment_fixture()
    broken$control$max.iterations <- 1L
    a <- assess.logistic.trend.filtering(broken, permutations = 2,
        inner.folds = 3, seed = 42, predictive = FALSE, warn = FALSE)
    expect_identical(a$status, "observed_fit_failed")
    expect_true(is.na(a$p.value))
    expect_true(all(is.na(a$p.value.bounds)))
})

test_that("unsupported nulls and malformed controls are rejected", {
    fit <- assessment_fixture()
    for (field in c("weights", "y", "covariates", "offset")) {
        changed <- fit
        changed[[field]] <- switch(field, weights = rep(2, 60), y = rep(.5, 60),
            covariates = matrix(fit$x), offset = rep(.1, 60))
        expect_error(assess.logistic.trend.filtering(changed, permutations = 1), "independent binary")
    }
    changed <- fit; changed$offset.used <- TRUE
    expect_error(assess.logistic.trend.filtering(changed, permutations = 1), "no offset")
    expect_error(assess.logistic.trend.filtering(NULL), "object")
    for (value in list(0, -1, NA_real_, Inf, 1.5, c(1,2), "2"))
        expect_error(assess.logistic.trend.filtering(fit, permutations = value), "permutations")
    expect_error(assess.logistic.trend.filtering(fit, inner.folds = 1), "inner.folds")
    expect_error(assess.logistic.trend.filtering(fit, outer.folds = 100), "Fold counts")
    expect_error(assess.logistic.trend.filtering(fit, predictive = NA), "predictive")
    expect_error(assess.logistic.trend.filtering(fit, seed = -1), "seed")
    changed <- fit; changed$x <- c(0, 1, 2, rep(3, 57))
    expect_error(assess.logistic.trend.filtering(changed, permutations = 1,
        inner.folds = 2, outer.folds = 60), "four distinct")
})

test_that("enumeration of a small conditional null gives conservative exact rejection rates", {
    # Under the constant Bernoulli null, all 70 labelings with four events among
    # eight independent observations are equally probable. This enumerates the
    # complete conditional null, rather than relying on a noisy simulation.
    x <- seq(0, 1, length.out = 8)
    labels <- combn(8, 4)
    statistics <- apply(labels, 2, function(events) {
        y <- numeric(8); y[events] <- 1
        fit <- fit.logistic.trend.filtering(x, y, c(Inf, .01),
            backend = "banded", curve.ridge = .5)
        a <- assess.logistic.trend.filtering(fit, permutations = 1,
            inner.folds = 8, seed = 17, predictive = FALSE, warn = FALSE)
        expect_identical(a$status, "complete")
        c(likelihood = a$statistic, a$curve.summaries)
    })
    expect_true(all(is.finite(statistics)))
    exact.p <- t(apply(statistics, 1, function(values) vapply(values, function(observed)
        mean(values >= observed-1e-10*max(1, abs(observed))), numeric(1))))
    holm <- apply(exact.p, 2, stats::p.adjust, method = "holm")
    for (alpha in c(.01, .05, .1, .2)) {
        expect_true(all(rowMeans(exact.p <= alpha) <= alpha))
        expect_lte(mean(colSums(holm <= alpha) > 0), alpha)
    }
})

test_that("assessment agrees between dense and banded objectives", {
    banded <- assessment_fixture()
    dense <- fit.logistic.trend.filtering(banded$x, banded$y, banded$lambda,
        backend = "dense", curve.ridge = banded$curve.ridge)
    a <- assess.logistic.trend.filtering(banded, permutations = 9,
        inner.folds = 3, outer.folds = 3, seed = 42,
        additional.tests = c("absolute.deviation", "total.variation"))
    b <- assess.logistic.trend.filtering(dense, permutations = 9,
        inner.folds = 3, outer.folds = 3, seed = 42,
        additional.tests = c("absolute.deviation", "total.variation"))
    expect_identical(a$status, "complete")
    expect_identical(b$status, "complete")
    expect_identical(a$folds, b$folds)
    expect_equal(a$statistic, b$statistic, tolerance = 2e-6)
    expect_equal(a$permutations$statistic, b$permutations$statistic, tolerance = 2e-6)
    expect_equal(a$p.value, b$p.value)
    expect_equal(a$prediction$predictions$probability,
                 b$prediction$predictions$probability, tolerance = 2e-6)
    expect_equal(a$prediction$improvement, b$prediction$improvement, tolerance = 2e-6)
    expect_equal(a$curve.summaries, b$curve.summaries, tolerance = 2e-6)
    expect_equal(a$permutation.statistics, b$permutation.statistics, tolerance = 2e-6)
    expect_equal(a$tests$p.value, b$tests$p.value)
})

test_that("curve summaries recover analytic probability variation including interior extrema", {
    # Construct known functions in the documented unscaled falling-factorial
    # basis; no fitting or test-statistic helper is used to obtain expectations.
    model <- function(coef) {
        x <- c(0, .1, .15, .75, 1)
        structure(list(x = x, backend = "dense", selected.index = 1L, lambda = 0,
            offset.used = FALSE,
            fits = list(list(converged = TRUE, coefficients = coef)),
            design = list(lo = 0, span = 1, knots = x, center = numeric(4),
                scales = rep(1,4), zcenter = numeric(), zscale = numeric())),
            class = "logistic.trend.filtering.fit")
    }
    tv <- .logistic.tf.assessment.variation
    expect_equal(tv(model(c(qlogis(.2),0,0,0,0))), 0)
    expect_equal(tv(model(c(-2,4,0,0,0))), plogis(2)-plogis(-2), tolerance=1e-12)
    # eta = 1 - 20*(x-.4)^2, maximum at .4 inside the long [.15,.75] piece.
    quadratic <- model(c(-2.2,16,-20,0,0))
    expected <- 2*plogis(1)-plogis(-2.2)-plogis(-6.2)
    expect_equal(tv(quadratic), expected, tolerance=1e-12)
    expect_gt(expected, sum(abs(diff(as.numeric(predict(quadratic, quadratic$x))))))
    # A knot changes curvature and derivative; use an independent dense grid
    # to check a piecewise curve with several turns, not just a quadratic.
    piecewise <- model(c(-1,8,-30,150,-300))
    grid <- sort(unique(c(seq(0,1,length.out=20001),piecewise$x)))
    expect_equal(tv(piecewise), sum(abs(diff(as.numeric(predict(piecewise,grid))))),
                 tolerance=1e-6)
    # Translation/scaling of x preserves variation and handles large origins.
    shifted <- quadratic; shifted$x <- 1e8 + 3*quadratic$x
    shifted$design$lo <- 1e8; shifted$design$span <- 3
    expect_equal(tv(shifted), expected, tolerance=1e-7)
})

test_that("optional tests share existing fits and independently reconstructed statistics", {
    fit <- assessment_fixture()
    old <- assess.logistic.trend.filtering(fit, permutations=9, inner.folds=3,
        outer.folds=3, seed=42)
    a <- assess.logistic.trend.filtering(fit, permutations=9, inner.folds=3,
        outer.folds=3, seed=42,
        additional.tests=c("absolute.deviation","total.variation"))
    expect_identical(a$folds, old$folds)
    expect_identical(a$seeds, old$seeds)
    expect_identical(a$permutations, old$permutations)
    expect_identical(a$prediction, old$prediction)
    expect_identical(a$fit, old$fit)
    expect_identical(a$p.value, old$p.value)
    expect_identical(a$statistic, old$statistic)
    expect_identical(a$p.value.bounds, old$p.value.bounds)
    expect_equal(old$tests$test, "likelihood")
    expect_equal(old$tests$p.value.holm, old$p.value)
    expect_equal(a$curve.summaries, old$curve.summaries)
    expect_equal(unname(a$curve.summaries["absolute.deviation"]),
                 mean(abs(as.numeric(fitted(a$fit))-mean(fit$y))))
    grid <- seq(min(fit$x),max(fit$x),length.out=10001)
    expect_equal(unname(a$curve.summaries["total.variation"]),
                 sum(abs(diff(as.numeric(predict(a$fit,grid))))), tolerance=1e-6)
    for (b in seq_len(9)) {
        set.seed(a$seeds$permutations[b])
        yy <- fit$y[sample.int(length(fit$y))]
        f <- fit.logistic.trend.filtering(fit$x,yy,fit$lambda,
            lambda.selection="cv", foldid=a$folds$tuning, backend="banded",
            curve.ridge=fit$curve.ridge)
        expect_equal(a$permutation.statistics$absolute.deviation[b],
                     mean(abs(as.numeric(fitted(f))-mean(yy))))
        expect_equal(a$permutation.statistics$total.variation[b],
                     sum(abs(diff(as.numeric(predict(f,grid))))), tolerance=1e-6)
    }
    for (j in seq_len(nrow(a$tests))) {
        tt <- a$tests[j,]
        values <- a$permutation.statistics[[tt$test]]
        p <- (1+sum(values >= tt$statistic-tt$comparison.tolerance))/10
        expect_equal(tt$p.value, p)
        expect_equal(c(tt$p.value.lower,tt$p.value.upper),rep(p,2))
        expect_equal(tt$completed,9)
    }
    expect_equal(a$tests$p.value.holm, p.adjust(a$tests$p.value,"holm"))
    expect_output(print(a),"Holm adjustment")
    b <- assess.logistic.trend.filtering(fit, permutations=1, inner.folds=3,
        seed=42, predictive=FALSE, additional.tests="total.variation")
    expect_equal(b$tests$test,c("likelihood","total.variation"))
    expect_equal(names(b$permutation.statistics),c("replicate","likelihood","total.variation"))
})

test_that("absolute deviation counts participants at tied predictors", {
    x <- rep(seq(0,1,length.out=12), c(10,rep(2,11)))
    y <- c(rep(0,8),rep(1,2),rep(c(0,1),11))
    fit <- fit.logistic.trend.filtering(x,y,c(Inf,.01),curve.ridge=.5,backend="banded")
    a <- assess.logistic.trend.filtering(fit,permutations=1,inner.folds=3,
        seed=11,predictive=FALSE,additional.tests="absolute.deviation")
    p <- as.numeric(predict(a$fit))
    expect_equal(unname(a$curve.summaries["absolute.deviation"]),mean(abs(p-mean(y))))
    expect_gt(abs(mean(abs(p-mean(y))) -
        mean(abs(as.numeric(predict(a$fit,unique(x)))-mean(y)))),1e-6)
})

test_that("optional tests retain failures, bounds and the complete adjustment family", {
    fit <- fit.logistic.trend.filtering(seq(0,1,length.out=12),
        rep(c(0,0,0,1),3),c(Inf,.01),backend="banded",curve.ridge=.1)
    a <- assess.logistic.trend.filtering(fit,permutations=39,inner.folds=3,
        seed=42,predictive=FALSE,warn=FALSE,
        additional.tests=c("absolute.deviation","total.variation"))
    expect_true(all(a$tests$status == "incomplete_permutations"))
    expect_true(all(is.na(a$tests$p.value)))
    expect_true(all(is.na(a$tests$p.value.holm)))
    for (j in seq_len(nrow(a$tests))) {
        tt <- a$tests[j,]; v <- a$permutation.statistics[[tt$test]]
        r <- sum(v >= tt$statistic-tt$comparison.tolerance,na.rm=TRUE)
        m <- sum(!is.finite(v))
        expect_equal(c(tt$p.value.lower,tt$p.value.upper),c(1+r,1+r+m)/40)
    }
    fit$control$max.iterations <- 1L
    b <- assess.logistic.trend.filtering(fit,permutations=1,inner.folds=3,
        seed=42,predictive=FALSE,warn=FALSE,
        additional.tests=c("absolute.deviation","total.variation"))
    expect_true(all(is.na(b$curve.summaries)))
    expect_true(all(is.na(b$tests$p.value.lower)))
    expect_true(all(is.na(b$tests$p.value.upper)))
    # Exact ties (including a flat curve) must all count as exceedances.
    tt <- .logistic.tf.assessment.test("total.variation",0,rep(0,9))
    expect_equal(tt$p.value,1)
})

test_that("additional tests reject malformed specifications before fitting", {
    fit <- assessment_fixture()
    for (value in list(NULL,NA_character_,1,"all","likelihood",
                       c("total.variation","total.variation")))
        expect_error(assess.logistic.trend.filtering(fit,permutations=1,
            additional.tests=value),"additional.tests")
})

test_that("summary arithmetic failures do not erase other tests or reduce the family", {
    fit <- assessment_fixture()
    local_mocked_bindings(.logistic.tf.assessment.variation = function(fit)
        stop("Deliberate summary failure"), .package="geosmooth")
    a <- assess.logistic.trend.filtering(fit,permutations=2,inner.folds=3,
        seed=42,predictive=FALSE,warn=FALSE,
        additional.tests=c("absolute.deviation","total.variation"))
    expect_identical(a$status,"complete")
    expect_true(is.finite(a$p.value))
    expect_true(is.finite(a$tests$p.value[2]))
    expect_true(is.na(a$tests$p.value[3]))
    expect_true(all(is.na(a$tests$p.value.holm)))
    expect_true(is.na(a$tests$p.value.lower[3]))
    expect_match(a$curve.summary.details$message,"Deliberate summary failure")
    expect_true(all(is.na(a$permutation.statistics$total.variation)))
    expect_true(all(is.finite(a$permutation.statistics$likelihood)))
})
