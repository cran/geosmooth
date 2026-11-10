## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", fig.width=7.4,
                    fig.height=4.8, dpi=120)
options(knitr.kable.NA="Unavailable")
library(geosmooth)

## ----folds-and-scores---------------------------------------------------------
make_folds <- function(y, k=5L, seed=1L) {
    set.seed(seed)
    fold <- integer(length(y))
    for (value in c(0, 1)) {
        at <- which(y == value)
        fold[at] <- sample(rep(seq_len(k), length.out=length(at)))
    }
    fold
}
log_loss <- function(y, eta) {
    mean((1-y)*pmax(eta, 0) + y*pmax(-eta, 0) +
         log1p(exp(-abs(eta))))
}
probability_scores <- function(y, eta, truth=NULL) {
    p <- plogis(eta)
    c(`Log loss`=log_loss(y, eta), `Brier score`=mean((p-y)^2),
      `Probability RMSE`=if (is.null(truth)) NA_real_ else sqrt(mean((p-truth)^2)))
}
lambda_path <- c(Inf, 10^seq(0, -4, by=-0.5))
spline_df <- c(3L, 5L, 7L)

## ----analysis-helpers, include=FALSE------------------------------------------
method_names <- c("Constant probability", "Linear logit", "Quadratic logit",
                  "Natural spline", "Trend filtering")
curve_colors <- c("Linear logit"="#0072B2", "Quadratic logit"="#009E73",
                  "Natural spline"="#8B5AA5", "Trend filtering"="#D55E00")

fit_reference <- function(x, y, kind, df=NULL, covariates=NULL) {
    lo <- min(x); span <- diff(range(x))
    if (!is.finite(span) || span <= 0) stop("Training predictor needs a positive range.")
    d <- data.frame(y=y, u=(x-lo)/span)
    zcenter <- zscale <- numeric()
    if (!is.null(covariates)) {
        zcenter <- colMeans(covariates)
        zz <- sweep(covariates, 2, zcenter)
        zscale <- sqrt(colMeans(zz^2))
        if (any(zscale <= 0)) stop("Constant training covariate.")
        zz <- sweep(zz, 2, zscale, "/")
        for (j in seq_len(ncol(zz))) d[[paste0("z", j)]] <- zz[, j]
    }
    rhs <- switch(kind, constant="1", linear="u", quadratic="u + I(u^2)",
                  spline=sprintf("splines::ns(u, df=%d)", df))
    if (length(zcenter)) rhs <- paste(rhs, paste0("z", seq_along(zcenter), collapse=" + "), sep=" + ")
    messages <- character()
    fit <- withCallingHandlers(
        stats::glm(stats::as.formula(paste("y ~", rhs)), data=d,
                   family=stats::binomial(), control=stats::glm.control(maxit=100)),
        warning=function(w) {
            messages <<- c(messages, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    if (!fit$converged || any(!is.finite(stats::coef(fit))) || length(messages))
        stop(paste(c("Reference fit did not pass convergence checks.", messages), collapse=" "))
    list(model=fit, lo=lo, span=span, zcenter=zcenter, zscale=zscale)
}
predict_reference <- function(object, x, covariates=NULL) {
    d <- data.frame(u=(x-object$lo)/object$span)
    if (length(object$zcenter)) {
        zz <- sweep(sweep(covariates, 2, object$zcenter), 2, object$zscale, "/")
        for (j in seq_len(ncol(zz))) d[[paste0("z", j)]] <- zz[, j]
    }
    eta <- as.numeric(stats::predict(object$model, newdata=d, type="link"))
    if (any(!is.finite(eta))) stop("Nonfinite reference prediction.")
    eta
}
fit_comparison <- function(x, y, folds, covariates=NULL) {
    basic <- lapply(c("constant", "linear", "quadratic"), function(kind)
        fit_reference(x, y, kind, covariates=covariates))
    names(basic) <- method_names[1:3]
    losses <- matrix(NA_real_, max(folds), length(spline_df))
    spline_status <- matrix("complete", max(folds), length(spline_df))
    for (k in sort(unique(folds))) for (j in seq_along(spline_df)) {
        train <- folds != k; test <- !train
        ztrain <- if (is.null(covariates)) NULL else covariates[train, , drop=FALSE]
        ztest <- if (is.null(covariates)) NULL else covariates[test, , drop=FALSE]
        result <- tryCatch({
            fit <- fit_reference(x[train], y[train], "spline", spline_df[j], ztrain)
            log_loss(y[test], predict_reference(fit, x[test], ztest))*sum(test)
        }, error=function(e) e)
        if (inherits(result, "error")) spline_status[k, j] <- conditionMessage(result)
        else losses[k, j] <- result
    }
    cv_loss <- colSums(losses)/length(y)
    if (!any(is.finite(cv_loss))) stop("No spline candidate passed all training folds.")
    selected_df <- spline_df[which.min(replace(cv_loss, !is.finite(cv_loss), Inf))]
    spline <- fit_reference(x, y, "spline", selected_df, covariates)
    # Candidate failures are displayed in the diagnostic tables below.
    tf <- fit.logistic.trend.filtering(x, y, lambda=lambda_path,
        covariates=covariates, lambda.selection="cv", foldid=folds,
        backend="banded", warn=FALSE)
    if (is.null(tf$selected.index)) stop("No trend-filtering candidate passed every fold and the full fit.")
    list(basic=basic, spline=spline, tf=tf,
         spline_cv=data.frame(df=spline_df, log.loss=cv_loss,
                              eligible=is.finite(cv_loss)),
         spline_fold_status=spline_status)
}
predict_comparison <- function(object, x, covariates=NULL) {
    ans <- vapply(object$basic, function(f) predict_reference(f, x, covariates), numeric(length(x)))
    ans <- cbind(ans, predict_reference(object$spline, x, covariates),
                 as.numeric(stats::predict(object$tf, x, covariates=covariates, type="link")))
    colnames(ans) <- method_names
    stopifnot(all(is.finite(ans)))
    ans
}
score_comparison <- function(object, x, y, truth=NULL, covariates=NULL) {
    eta <- predict_comparison(object, x, covariates)
    out <- t(vapply(seq_len(ncol(eta)), function(j)
        probability_scores(y, eta[, j], truth), numeric(3)))
    data.frame(Model=colnames(eta), out, row.names=NULL, check.names=FALSE)
}
penalty_label <- function(value) {
    if (is.infinite(value)) "Quadratic limit" else format(value, scientific=TRUE, digits=3)
}
tuning_summary <- function(object) {
    tf <- object$tf; j <- tf$selected.index
    data.frame(`Selected lambda`=penalty_label(tf$lambda[j]),
        `Adaptive knots`=tf$diagnostics$working.knots[j],
        `Eligible lambda candidates`=sum(tf$cv.table$eligible),
        `Lambda candidates`=length(tf$lambda),
        `Selected spline df`=object$spline_cv$df[which.min(replace(
            object$spline_cv$log.loss, !object$spline_cv$eligible, Inf))], check.names=FALSE)
}
cv_plot <- function(object, title) {
    tab <- object$tf$cv.table
    finite <- is.finite(tab$lambda) & is.finite(tab$log.loss)
    ylim <- range(tab$log.loss[is.finite(tab$log.loss)])
    plot(log10(tab$lambda[finite]), tab$log.loss[finite], type="b", pch=16,
         col="#0072B2", xlab=expression(log[10](lambda)), ylab="Training CV log loss",
         main=title, ylim=ylim)
    quadratic <- tab$log.loss[is.infinite(tab$lambda)]
    if (is.finite(quadratic)) abline(h=quadratic, lty=2, col="#777777")
    selected <- object$tf$lambda[object$tf$selected.index]
    if (is.finite(selected)) abline(v=log10(selected), col="#D55E00", lty=3, lwd=2)
}
curve_plot <- function(x, y, grid, eta, title, truth=NULL, xlab="Predictor x") {
    plot(grid, rep(NA_real_, length(grid)), type="n", ylim=c(0,1),
         xlab=xlab, ylab="", yaxt="n", main=title)
    axis(2, at=seq(0,1,by=.2), las=1)
    mtext("Outcome probability", side=2, line=2.9)
    if (!is.null(truth)) lines(grid, truth, col="#222222", lty=2, lwd=2)
    for (name in names(curve_colors)) lines(grid, plogis(eta[, name]),
        col=curve_colors[name], lwd=if (name=="Trend filtering") 2.3 else 1.2)
    rug(x[y==0], side=1, col=grDevices::adjustcolor("#555555", .25), ticksize=.025)
    rug(x[y==1], side=3, col=grDevices::adjustcolor("#555555", .25), ticksize=.025)
    labs <- names(curve_colors); cols <- unname(curve_colors); types <- rep(1, length(labs))
    if (!is.null(truth)) {labs <- c("Known truth", labs); cols <- c("#222222", cols); types <- c(2, types)}
    legend("topleft", labs, col=cols, lty=types, lwd=1.7, bty="n", cex=.63, ncol=2)
}
bin_summary <- function(x, y, bins=8L) {
    breaks <- unique(as.numeric(stats::quantile(x, seq(0,1,length.out=bins+1))))
    group <- if (length(breaks)<2) rep(1L, length(x)) else
        cut(x, breaks=breaks, include.lowest=TRUE, labels=FALSE)
    out <- do.call(rbind, lapply(split(seq_along(x), group), function(i) {
        n <- length(i); p <- mean(y[i]); z <- stats::qnorm(.975)
        center <- (p+z^2/(2*n))/(1+z^2/n)
        half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n)
        data.frame(x=mean(x[i]), observed=p, n=n, lower=center-half, upper=center+half)
    }))
    rownames(out) <- NULL; out
}

## ----synthetic-data-----------------------------------------------------------
truth_eta <- list(
    `Linear logit`=function(x) -2+4*x,
    `Local changes`=function(x) -1.4+2.8*exp(-((x-.35)/.08)^2)-
        1.5*exp(-((x-.75)/.15)^2))
synthetic <- vector("list", length(truth_eta)); names(synthetic) <- names(truth_eta)
for (j in seq_along(truth_eta)) {
    set.seed(9100L+j)
    train_x <- runif(600); train_p <- plogis(truth_eta[[j]](train_x))
    train_y <- rbinom(600, 1, train_p)
    test_x <- runif(2000); test_p <- plogis(truth_eta[[j]](test_x))
    test_y <- rbinom(2000, 1, test_p)
    folds <- make_folds(train_y, seed=110L+j)
    fit <- fit_comparison(train_x, train_y, folds)
    synthetic[[j]] <- list(x=train_x, y=train_y, test_x=test_x,
        test_y=test_y, test_p=test_p, folds=folds, fit=fit)
}
knitr::kable(do.call(rbind, lapply(names(synthetic), function(name) {
    d <- synthetic[[name]]
    data.frame(Example=name, `Training observations`=length(d$y),
        `Training events`=sum(d$y), `Test observations`=length(d$test_y),
        `Test events`=sum(d$test_y), check.names=FALSE)
})), row.names=FALSE)

## ----direct-fit, eval=FALSE---------------------------------------------------
# d <- synthetic[["Local changes"]]
# fit <- fit.logistic.trend.filtering(d$x, d$y, lambda=lambda_path,
#     lambda.selection="cv", foldid=d$folds, backend="banded")
# fit$cv.table
# fit$diagnostics
# predict(fit, seq(0, 1, length.out=300))

## ----synthetic-curves, fig.height=4.5, fig.cap="Known probabilities and training-data fits. Bottom rugs show controls and top rugs show events. Both panels use the same probability scale; separate test observations assess prediction below.", fig.alt="Two panels compare known probability curves with linear, quadratic, spline and trend-filtering fits for a linear-logit signal and a localized peak/trough signal."----
grid <- seq(0,1,length.out=300)
old_par <- par(mfrow=c(1,2), mar=c(4,4,2.5,.5))
for (name in names(synthetic)) {
    d <- synthetic[[name]]
    curve_plot(d$x, d$y, grid, predict_comparison(d$fit, grid),
               name, truth=plogis(truth_eta[[name]](grid)))
}
par(old_par)

## ----synthetic-cv, fig.height=4, fig.cap="Training-only penalty selection on the same five folds used by the spline comparator. Blue points show finite penalties, the gray dashed line the quadratic limit, and an orange dotted line the selected finite penalty. The quadratic limit cannot be placed on a log10(lambda) axis.", fig.alt="Two cross-validation plots show log loss against finite log10 penalty values, with quadratic-limit reference lines and selected-penalty markers."----
old_par <- par(mfrow=c(1,2), mar=c(4,4,2.5,.5))
for (name in names(synthetic)) cv_plot(synthetic[[name]]$fit, name)
par(old_par)
knitr::kable(do.call(rbind, lapply(names(synthetic), function(name)
    cbind(Example=name, tuning_summary(synthetic[[name]]$fit)))), row.names=FALSE)

## ----smoothing-path, fig.cap="Three members of the same local-changes path: the quadratic limit, the CV-selected candidate, and the weakest supplied penalty. Coincident curves are retained when the selected candidate is an endpoint.", fig.alt="The known probability curve is overlaid with quadratic-limit, selected and weak-penalty trend-filtering curves."----
d <- synthetic[["Local changes"]]; tf <- d$fit$tf
indices <- c(which(is.infinite(tf$lambda)), tf$selected.index, which.min(tf$lambda))
path_p <- plogis(predict(tf, grid, lambda.index=indices, type="link"))
old_par <- par(mar=c(4,4,1,.5))
plot(grid, plogis(truth_eta[["Local changes"]](grid)), type="l", ylim=c(0,1),
     xlab="Predictor x", ylab="Outcome probability", lwd=2, lty=2, col="#222222")
matlines(grid, path_p, col=c("#009E73", "#D55E00", "#0072B2"), lty=c(1,1,3), lwd=2)
legend("topright", c("Known truth", "Quadratic limit", "Selected penalty", "Weakest penalty"),
       col=c("#222222", "#009E73", "#D55E00", "#0072B2"), lty=c(2,1,1,3), bty="n", cex=.8)
par(old_par)

## ----synthetic-test-scores----------------------------------------------------
synthetic_scores <- do.call(rbind, lapply(names(synthetic), function(name) {
    d <- synthetic[[name]]
    cbind(Example=name, score_comparison(d$fit, d$test_x, d$test_y, d$test_p))
}))
knitr::kable(synthetic_scores, digits=4, row.names=FALSE)

## ----binomial-ties------------------------------------------------------------
set.seed(9140)
x_group <- seq(0,1,length.out=12)
trials <- rep(c(8L, 12L, 16L), 4)
successes <- rbinom(length(x_group), trials, plogis(-1+2*x_group))
x_individual <- rep(x_group, trials)
y_individual <- unlist(Map(function(s, n) c(rep(1, s), rep(0, n-s)), successes, trials))
individual <- fit.logistic.trend.filtering(x_individual, y_individual,
    lambda=.01, backend="banded")
grouped <- fit.logistic.trend.filtering(x_group, successes/trials,
    weights=trials, lambda=.01, backend="banded")
maximum_difference <- max(abs(as.numeric(predict(individual, grid))-
                              as.numeric(predict(grouped, grid))))
stopifnot(maximum_difference < 1e-7,
          abs(individual$fits[[1]]$objective-grouped$fits[[1]]$objective) < 1e-6)
knitr::kable(data.frame(`Binary records`=length(y_individual),
    `Distinct locations`=length(x_group), `Maximum probability difference`=maximum_difference,
    check.names=FALSE), digits=9)

## ----pima-training------------------------------------------------------------
pima_train <- MASS::Pima.tr
pima_y <- as.numeric(pima_train$type == "Yes")
pima_folds <- make_folds(pima_y, seed=9201)
pima <- fit_comparison(pima_train$glu, pima_y, pima_folds)
# Freeze the predeclared optional adjusted comparison before opening the test set.
z_train <- as.matrix(pima_train[, c("age", "bmi")])
pima_adjusted <- fit_comparison(pima_train$glu, pima_y, pima_folds,
                               covariates=z_train)
knitr::kable(tuning_summary(pima), row.names=FALSE)

## ----pima-tuning-diagnostics--------------------------------------------------
pima_cv <- pima$tf$cv.table
pima_cv$lambda <- vapply(pima_cv$lambda, penalty_label, character(1))
pima_cv$selected <- seq_len(nrow(pima_cv)) == pima$tf$selected.index
pima_cv$full.fit.status <- pima$tf$diagnostics$status
knitr::kable(pima_cv, digits=4, row.names=FALSE)
knitr::kable(pima$spline_cv, digits=4, row.names=FALSE)

## ----pima-test----------------------------------------------------------------
pima_test <- MASS::Pima.te
pima_test_y <- as.numeric(pima_test$type == "Yes")
pima_eta <- predict_comparison(pima, pima_test$glu)
pima_scores <- score_comparison(pima, pima_test$glu, pima_test_y)
knitr::kable(data.frame(Split=c("Training", "Test"),
    Women=c(nrow(pima_train), nrow(pima_test)),
    Events=c(sum(pima_y), sum(pima_test_y)),
    `Glucose minimum`=c(min(pima_train$glu), min(pima_test$glu)),
    `Glucose maximum`=c(max(pima_train$glu), max(pima_test$glu)), check.names=FALSE))
knitr::kable(pima_scores[, c("Model", "Log loss", "Brier score")], digits=4, row.names=FALSE)
pima_outside <- sum(pima_test$glu < min(pima_train$glu) |
                    pima_test$glu > max(pima_train$glu))

## ----pima-curves, fig.cap="Training-data fits and descriptive training proportions. Points and 95% Wilson intervals summarize glucose groups; they are not confidence bands for a fitted curve. Bottom/top rugs mark training controls/events. Models used every original training record, without binning.", fig.alt="Glucose probability curves for four models, with binned training proportions, descriptive binomial intervals and case/control rugs."----
pima_grid <- seq(min(pima_train$glu), max(pima_train$glu), length.out=250)
old_par <- par(mar=c(4,4,2.5,.5))
curve_plot(pima_train$glu, pima_y, pima_grid,
           predict_comparison(pima, pima_grid), "Diabetes status and glucose",
           xlab="Plasma glucose (recorded dataset scale)")
bins <- bin_summary(pima_train$glu, pima_y)
segments(bins$x, bins$lower, bins$x, bins$upper, col="#555555")
points(bins$x, bins$observed, pch=21, bg="white", cex=.9)
par(old_par)

## ----pima-calibration, fig.width=5.4, fig.height=5, fig.cap="Calibration on the untouched 332-woman test set. Each model forms its own groups from predicted probabilities. The dashed diagonal is perfect agreement; intervals are descriptive 95% Wilson intervals.", fig.alt="Held-out observed event proportions versus mean predicted probability for linear logistic regression and trend filtering, with binomial intervals and a diagonal agreement line."----
old_par <- par(pty="s")
plot(c(0,1), c(0,1), type="n", xlab="Mean predicted probability",
     ylab="Observed test event proportion", xlim=c(0,1), ylim=c(0,1),
     xaxs="i", yaxs="i")
abline(0,1,lty=2,col="#777777")
for (name in c("Linear logit", "Trend filtering")) {
    b <- bin_summary(plogis(pima_eta[, name]), pima_test_y, bins=5)
    color <- curve_colors[name]
    segments(b$x, b$lower, b$x, b$upper, col=color)
    points(b$x, b$observed, col=color, pch=if (name=="Linear logit") 16 else 17)
}
legend("topleft", c("Linear logit", "Trend filtering"),
       col=curve_colors[c("Linear logit", "Trend filtering")], pch=c(16,17), bty="n")
par(old_par)
calibration <- do.call(rbind, lapply(c("Linear logit", "Trend filtering"), function(name) {
    eta <- pima_eta[, name]
    fit <- stats::glm(pima_test_y ~ eta, family=stats::binomial())
    stopifnot(fit$converged, all(is.finite(stats::coef(fit))))
    data.frame(Model=name, Intercept=unname(stats::coef(fit)[1]),
               Slope=unname(stats::coef(fit)[2]))
}))
knitr::kable(calibration, digits=3, row.names=FALSE)

## ----pima-adjustment----------------------------------------------------------
z_test <- as.matrix(pima_test[, c("age", "bmi")])
adjusted_scores <- score_comparison(pima_adjusted, pima_test$glu, pima_test_y,
                                   covariates=z_test)
show_adjusted <- adjusted_scores$Model %in% c("Linear logit", "Quadratic logit", "Trend filtering")
knitr::kable(adjusted_scores[show_adjusted, c("Model", "Log loss", "Brier score")],
             digits=4, row.names=FALSE)
knitr::kable(tuning_summary(pima_adjusted), row.names=FALSE)
profile <- apply(z_train, 2, stats::median)

## ----pima-adjusted-curves, fig.cap="Conditional probability curves with age and BMI fixed at their training medians. The adjustment changes the question; these curves are neither population averages nor causal effects. Sparse support for particular glucose/age/BMI combinations remains possible.", fig.alt="Two panels show unadjusted and age/BMI-adjusted glucose curves, with adjustment evaluated at the training-median covariate profile."----
z_grid <- matrix(rep(profile, each=length(pima_grid)), ncol=2,
                 dimnames=list(NULL, colnames(z_train)))
old_par <- par(mfrow=c(1,2), mar=c(4,4,2.5,.5))
curve_plot(pima_train$glu, pima_y, pima_grid, predict_comparison(pima, pima_grid),
           "Glucose only", xlab="Glucose (recorded scale)")
curve_plot(pima_train$glu, pima_y, pima_grid,
           predict_comparison(pima_adjusted, pima_grid, z_grid),
           "Adjusted, median profile", xlab="Glucose (recorded scale)")
par(old_par)

## ----selected-diagnostics-----------------------------------------------------
selected_diagnostics <- do.call(rbind, lapply(
    list(`Synthetic linear`=synthetic[[1]]$fit,
         `Synthetic local changes`=synthetic[[2]]$fit,
         `Pima glucose`=pima, `Pima adjusted`=pima_adjusted), function(object) {
        object$tf$diagnostics[object$tf$selected.index,
            c("lambda", "converged", "status", "stationarity", "complementarity", "backend", "working.knots")]
    }))
diagnostic_table <- data.frame(
    Example=rownames(selected_diagnostics),
    Lambda=vapply(selected_diagnostics$lambda, penalty_label, character(1)),
    Converged=selected_diagnostics$converged,
    Stationarity=formatC(selected_diagnostics$stationarity, format="e", digits=2),
    Complementarity=formatC(selected_diagnostics$complementarity, format="e", digits=2),
    Knots=selected_diagnostics$working.knots)
knitr::kable(diagnostic_table, row.names=FALSE)
stopifnot(all(selected_diagnostics$converged))

## ----explicit-failure---------------------------------------------------------
failed <- fit.logistic.trend.filtering(synthetic[[2]]$x, synthetic[[2]]$y,
    lambda=.001, backend="banded", control=list(max.iterations=1L), warn=FALSE)
failed$diagnostics[, c("converged", "status")]
stopifnot(!failed$diagnostics$converged, all(is.na(predict(failed))))

## ----optional-scaling, eval=FALSE---------------------------------------------
# set.seed(9301)
# n <- 24000L
# x <- runif(n)
# y <- rbinom(n, 1, plogis(-1+sin(6*pi*x)))
# system.time(large_fit <- fit.logistic.trend.filtering(x, y,
#     lambda=c(1, .1, .01, .001), backend="banded"))
# large_fit$diagnostics

## ----global-assessment--------------------------------------------------------
assessment_examples <- list()
for (j in 1:2) {
    set.seed(9400L+j)
    x <- seq(0, 1, length.out=160)
    probability <- if (j == 1) rep(.35, length(x)) else
        plogis(-2+4.8*exp(-((x-.4)/.15)^2))
    y <- rbinom(length(x), 1, probability)
    fit <- fit.logistic.trend.filtering(x, y, lambda=c(Inf, .1, .01, .001),
                                       backend="banded")
    assessment_examples[[j]] <- assess.logistic.trend.filtering(fit,
        permutations=99, inner.folds=3, outer.folds=5, seed=9500L+j,
        additional.tests=c("absolute.deviation", "total.variation"))
}
names(assessment_examples) <- c("Constant truth", "Localized peak")
assessment_scores <- do.call(rbind, lapply(names(assessment_examples), function(name) {
    a <- assessment_examples[[name]]
    stopifnot(a$status == "complete", a$prediction$status == "complete")
    data.frame(Example=name, `Permutation p-value`=a$p.value,
        `Log-loss improvement`=unname(a$prediction$improvement["log.loss"]),
        `Brier improvement`=unname(a$prediction$improvement["brier"]), check.names=FALSE)
}))
knitr::kable(assessment_scores, row.names=FALSE, digits=4)

## ----probability-shape-summaries----------------------------------------------
shape_scores <- do.call(rbind, lapply(names(assessment_examples), function(name) {
    a <- assessment_examples[[name]]
    data.frame(Example=name,
        `Mean absolute deviation`=unname(a$curve.summaries["absolute.deviation"]),
        `Total variation`=unname(a$curve.summaries["total.variation"]),
        check.names=FALSE)
}))
knitr::kable(shape_scores, row.names=FALSE, digits=4)

## ----probability-shape-tests--------------------------------------------------
shape_tests <- do.call(rbind, lapply(names(assessment_examples), function(name) {
    tt <- assessment_examples[[name]]$tests
    stopifnot(all(tt$status == "complete"))
    data.frame(Example=name, Statistic=tt$test,
        `Individual p-value`=tt$p.value, `Holm p-value`=tt$p.value.holm,
        check.names=FALSE)
}))
knitr::kable(shape_tests, row.names=FALSE, digits=3)

## ----permutation-comparison, fig.height=3.8, fig.cap="Each histogram shows likelihood improvements after randomly shuffling outcomes while retaining predictor locations and the event count. The orange line marks the observed statistic. The horizontal scales include the observed statistic, so an unusually large improvement can lie far beyond the shuffled results.", fig.alt="Two histograms compare shuffled-outcome likelihood improvements with the observed improvement for constant truth and a localized peak."----
old_par <- par(mfrow=c(1, 2), mar=c(4, 4, 2, .5))
for (name in names(assessment_examples)) {
    a <- assessment_examples[[name]]
    hist(a$permutations$statistic, breaks=10, col="grey85", border="white",
         xlim=range(c(0, a$statistic, a$permutations$statistic)), main=name,
         xlab="Likelihood improvement statistic", ylab="Permutation count")
    abline(v=a$statistic, col="#D55E00", lwd=2)
}
par(old_par)

## ----helper-source, ref.label='analysis-helpers', eval=FALSE------------------
# method_names <- c("Constant probability", "Linear logit", "Quadratic logit",
#                   "Natural spline", "Trend filtering")
# curve_colors <- c("Linear logit"="#0072B2", "Quadratic logit"="#009E73",
#                   "Natural spline"="#8B5AA5", "Trend filtering"="#D55E00")
# 
# fit_reference <- function(x, y, kind, df=NULL, covariates=NULL) {
#     lo <- min(x); span <- diff(range(x))
#     if (!is.finite(span) || span <= 0) stop("Training predictor needs a positive range.")
#     d <- data.frame(y=y, u=(x-lo)/span)
#     zcenter <- zscale <- numeric()
#     if (!is.null(covariates)) {
#         zcenter <- colMeans(covariates)
#         zz <- sweep(covariates, 2, zcenter)
#         zscale <- sqrt(colMeans(zz^2))
#         if (any(zscale <= 0)) stop("Constant training covariate.")
#         zz <- sweep(zz, 2, zscale, "/")
#         for (j in seq_len(ncol(zz))) d[[paste0("z", j)]] <- zz[, j]
#     }
#     rhs <- switch(kind, constant="1", linear="u", quadratic="u + I(u^2)",
#                   spline=sprintf("splines::ns(u, df=%d)", df))
#     if (length(zcenter)) rhs <- paste(rhs, paste0("z", seq_along(zcenter), collapse=" + "), sep=" + ")
#     messages <- character()
#     fit <- withCallingHandlers(
#         stats::glm(stats::as.formula(paste("y ~", rhs)), data=d,
#                    family=stats::binomial(), control=stats::glm.control(maxit=100)),
#         warning=function(w) {
#             messages <<- c(messages, conditionMessage(w))
#             invokeRestart("muffleWarning")
#         })
#     if (!fit$converged || any(!is.finite(stats::coef(fit))) || length(messages))
#         stop(paste(c("Reference fit did not pass convergence checks.", messages), collapse=" "))
#     list(model=fit, lo=lo, span=span, zcenter=zcenter, zscale=zscale)
# }
# predict_reference <- function(object, x, covariates=NULL) {
#     d <- data.frame(u=(x-object$lo)/object$span)
#     if (length(object$zcenter)) {
#         zz <- sweep(sweep(covariates, 2, object$zcenter), 2, object$zscale, "/")
#         for (j in seq_len(ncol(zz))) d[[paste0("z", j)]] <- zz[, j]
#     }
#     eta <- as.numeric(stats::predict(object$model, newdata=d, type="link"))
#     if (any(!is.finite(eta))) stop("Nonfinite reference prediction.")
#     eta
# }
# fit_comparison <- function(x, y, folds, covariates=NULL) {
#     basic <- lapply(c("constant", "linear", "quadratic"), function(kind)
#         fit_reference(x, y, kind, covariates=covariates))
#     names(basic) <- method_names[1:3]
#     losses <- matrix(NA_real_, max(folds), length(spline_df))
#     spline_status <- matrix("complete", max(folds), length(spline_df))
#     for (k in sort(unique(folds))) for (j in seq_along(spline_df)) {
#         train <- folds != k; test <- !train
#         ztrain <- if (is.null(covariates)) NULL else covariates[train, , drop=FALSE]
#         ztest <- if (is.null(covariates)) NULL else covariates[test, , drop=FALSE]
#         result <- tryCatch({
#             fit <- fit_reference(x[train], y[train], "spline", spline_df[j], ztrain)
#             log_loss(y[test], predict_reference(fit, x[test], ztest))*sum(test)
#         }, error=function(e) e)
#         if (inherits(result, "error")) spline_status[k, j] <- conditionMessage(result)
#         else losses[k, j] <- result
#     }
#     cv_loss <- colSums(losses)/length(y)
#     if (!any(is.finite(cv_loss))) stop("No spline candidate passed all training folds.")
#     selected_df <- spline_df[which.min(replace(cv_loss, !is.finite(cv_loss), Inf))]
#     spline <- fit_reference(x, y, "spline", selected_df, covariates)
#     # Candidate failures are displayed in the diagnostic tables below.
#     tf <- fit.logistic.trend.filtering(x, y, lambda=lambda_path,
#         covariates=covariates, lambda.selection="cv", foldid=folds,
#         backend="banded", warn=FALSE)
#     if (is.null(tf$selected.index)) stop("No trend-filtering candidate passed every fold and the full fit.")
#     list(basic=basic, spline=spline, tf=tf,
#          spline_cv=data.frame(df=spline_df, log.loss=cv_loss,
#                               eligible=is.finite(cv_loss)),
#          spline_fold_status=spline_status)
# }
# predict_comparison <- function(object, x, covariates=NULL) {
#     ans <- vapply(object$basic, function(f) predict_reference(f, x, covariates), numeric(length(x)))
#     ans <- cbind(ans, predict_reference(object$spline, x, covariates),
#                  as.numeric(stats::predict(object$tf, x, covariates=covariates, type="link")))
#     colnames(ans) <- method_names
#     stopifnot(all(is.finite(ans)))
#     ans
# }
# score_comparison <- function(object, x, y, truth=NULL, covariates=NULL) {
#     eta <- predict_comparison(object, x, covariates)
#     out <- t(vapply(seq_len(ncol(eta)), function(j)
#         probability_scores(y, eta[, j], truth), numeric(3)))
#     data.frame(Model=colnames(eta), out, row.names=NULL, check.names=FALSE)
# }
# penalty_label <- function(value) {
#     if (is.infinite(value)) "Quadratic limit" else format(value, scientific=TRUE, digits=3)
# }
# tuning_summary <- function(object) {
#     tf <- object$tf; j <- tf$selected.index
#     data.frame(`Selected lambda`=penalty_label(tf$lambda[j]),
#         `Adaptive knots`=tf$diagnostics$working.knots[j],
#         `Eligible lambda candidates`=sum(tf$cv.table$eligible),
#         `Lambda candidates`=length(tf$lambda),
#         `Selected spline df`=object$spline_cv$df[which.min(replace(
#             object$spline_cv$log.loss, !object$spline_cv$eligible, Inf))], check.names=FALSE)
# }
# cv_plot <- function(object, title) {
#     tab <- object$tf$cv.table
#     finite <- is.finite(tab$lambda) & is.finite(tab$log.loss)
#     ylim <- range(tab$log.loss[is.finite(tab$log.loss)])
#     plot(log10(tab$lambda[finite]), tab$log.loss[finite], type="b", pch=16,
#          col="#0072B2", xlab=expression(log[10](lambda)), ylab="Training CV log loss",
#          main=title, ylim=ylim)
#     quadratic <- tab$log.loss[is.infinite(tab$lambda)]
#     if (is.finite(quadratic)) abline(h=quadratic, lty=2, col="#777777")
#     selected <- object$tf$lambda[object$tf$selected.index]
#     if (is.finite(selected)) abline(v=log10(selected), col="#D55E00", lty=3, lwd=2)
# }
# curve_plot <- function(x, y, grid, eta, title, truth=NULL, xlab="Predictor x") {
#     plot(grid, rep(NA_real_, length(grid)), type="n", ylim=c(0,1),
#          xlab=xlab, ylab="", yaxt="n", main=title)
#     axis(2, at=seq(0,1,by=.2), las=1)
#     mtext("Outcome probability", side=2, line=2.9)
#     if (!is.null(truth)) lines(grid, truth, col="#222222", lty=2, lwd=2)
#     for (name in names(curve_colors)) lines(grid, plogis(eta[, name]),
#         col=curve_colors[name], lwd=if (name=="Trend filtering") 2.3 else 1.2)
#     rug(x[y==0], side=1, col=grDevices::adjustcolor("#555555", .25), ticksize=.025)
#     rug(x[y==1], side=3, col=grDevices::adjustcolor("#555555", .25), ticksize=.025)
#     labs <- names(curve_colors); cols <- unname(curve_colors); types <- rep(1, length(labs))
#     if (!is.null(truth)) {labs <- c("Known truth", labs); cols <- c("#222222", cols); types <- c(2, types)}
#     legend("topleft", labs, col=cols, lty=types, lwd=1.7, bty="n", cex=.63, ncol=2)
# }
# bin_summary <- function(x, y, bins=8L) {
#     breaks <- unique(as.numeric(stats::quantile(x, seq(0,1,length.out=bins+1))))
#     group <- if (length(breaks)<2) rep(1L, length(x)) else
#         cut(x, breaks=breaks, include.lowest=TRUE, labels=FALSE)
#     out <- do.call(rbind, lapply(split(seq_along(x), group), function(i) {
#         n <- length(i); p <- mean(y[i]); z <- stats::qnorm(.975)
#         center <- (p+z^2/(2*n))/(1+z^2/n)
#         half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n)
#         data.frame(x=mean(x[i]), observed=p, n=n, lower=center-half, upper=center+half)
#     }))
#     rownames(out) <- NULL; out
# }

