#' Test Constancy and Assess Gaussian Trend-Filtering Prediction
#'
#' Assess a continuous outcome using order-2 trend filtering, retuned outcome
#' permutations and a separate nested-CV comparison with a constant mean.
#'
#' @inheritParams assess.logistic.trend.filtering
#' @param object A `gaussian.trend.filtering.fit` with independent, unit-weight
#'   observations, no covariates and no offset. Predictor choice, penalty path,
#'   ridge penalties and numerical settings must be specified independently of
#'   outcomes. The assessment uses order 2 and reuses these explicit settings.
#' @param additional.tests Unique optional names `"absolute.deviation"` and/or
#'   `"total.variation"`. Default `character()` tests squared-error improvement
#'   only. Both observed curve summaries are always returned.
#'
#' @details The randomization null is an outcome distribution that is unchanged
#'   across the fixed predictor locations. Independent outcomes must be
#'   exchangeable under the null; a constant mean plus identically distributed
#'   independent errors is a sufficient model. Normality is unnecessary.
#'   This is not a test of constant mean under arbitrary heteroskedasticity.
#'   Repeated participants, serial dependence, covariate adjustment and
#'   unequal weights require other null procedures and are rejected where
#'   detectable. Independence and exchangeability cannot be checked from a fit.
#'
#'   Fresh balanced folds are generated independently of outcomes; stored fit
#'   selections and fold labels are ignored. The explicit lambda path is tuned
#'   by CV on the full data and again on every shuffled outcome vector, always
#'   using the same tuning folds. The primary statistic is
#'   `max(0, sum((y-mean(y))^2) - sum((y-f.hat(x))^2))`, the nonnegative
#'   reduction in residual sum of squares, in squared response units. It is
#'   not a penalized-objective difference or a chi-squared statistic. Constant
#'   outcomes have zero statistics and permutation p-values of one.
#'
#'   With B requested permutations the inclusive plus-one p-value is
#'   `(1 + sum(permuted.statistic >= observed.statistic))/(B + 1)`, allowing
#'   a conservative numerical tie tolerance of
#'   `1e-10 * max(1, abs(observed.statistic))`. Retuning is part of each
#'   randomization. Failed candidates are ineligible, and a replicate fails
#'   if none is eligible. Failed replicates are retained without retries or
#'   denominator changes. Any failure suppresses the point p-value; lower
#'   and upper bounds count missing statistics as nonexceedances or exceedances.
#'   Observed-fit failure leaves both bounds unavailable. These are failure
#'   bounds, not confidence intervals.
#'
#'   Absolute deviation is `mean(abs(f.hat(x)-mean(y)))`, averaging over
#'   observations with repeated predictors counted at their observed frequency.
#'   Total variation sums absolute changes along the continuous fitted curve
#'   over `[min(x),max(x)]`, including every knot and interior quadratic turning
#'   point. Both are in response units. Total variation is not the third-
#'   difference penalty. These describe the newly tuned observed fit; they are
#'   not bias-corrected population effects or uncertainty intervals. Optional
#'   tests reuse the same fitted permutations, with per-statistic failure
#'   bounds and Holm adjustment across squared error and the requested tests.
#'   An incomplete family has no Holm p-values. This adjustment does not cover
#'   additional predictors or outcomes. Choose tests and penalties in advance.
#'
#'   Nested CV inherits the input fit's `cv.penalty.scale`. Under
#'   `"training_mass"`, outer training fits scale the reference lambda grid
#'   and ridge penalties by their sample-size fraction; inner folds apply
#'   their own training fraction. Full-data observed and permutation refits
#'   retain the reference grid. Older saved objects without this field use
#'   `"none"` to preserve their earlier convention.
#'   Nested CV tunes using only each outer training set. Its constant baseline
#'   is that set's sample mean. Pooled mean squared error (`mse`) and mean
#'   absolute error (`mae`) use the same held-out rows; improvements are
#'   constant minus smoother loss. Positive values favor trend filtering.
#'   The baseline remains the training mean for both scores. Boundary
#'   predictions use the fitter's constant-boundary policy and are flagged.
#'   Incomplete predictions leave overall smoother scores unavailable.
#'   No standard errors are computed from overlapping CV folds. Prediction
#'   scores assess training at the outer-training sample size, and need not
#'   improve when the global association test rejects.
#'
#'   Runs are sequential, with the same fitting cost and failure policy as
#'   \code{\link{assess.logistic.trend.filtering}}. The test does not establish
#'   nonlinearity, causality or a useful predictive effect.
#'
#' @return A `gaussian.trend.filtering.assessment` with the same evidence
#'   structure as \code{\link{assess.logistic.trend.filtering}}, replacing
#'   `null.probability` with `null.mean` and the primary `likelihood` test with
#'   `squared.error`. Scalar `statistic`, `p.value`, `p.value.bounds`, `status`
#'   and the `permutations` table refer to squared-error improvement.
#'   `fit` is the newly tuned observed fit. `curve.summaries`, `tests`, and
#'   `permutation.statistics` report curve summaries and requested tests.
#'   Fold assignments, seeds, candidate diagnostics and failure messages are
#'   retained. Optional `prediction` includes row-level `mean` and
#'   `constant.mean` predictions, losses, fold diagnostics, model scores and
#'   baseline-minus-smoother improvements. No confidence intervals are supplied.
#' @references
#' Phipson and Smyth (2010). Permutation p-values should never be zero.
#' \doi{10.2202/1544-6115.1585}
#' @seealso \code{\link{fit.gaussian.trend.filtering}}, \code{\link{assess.logistic.trend.filtering}}
#' @examples
#' x <- seq(0, 1, length.out = 40)
#' y <- sin(2*pi*x) + 0.2*cos(seq_along(x))
#' fit <- fit.gaussian.trend.filtering(x, y, c(Inf, 0.01), backend = "banded")
#' # Nine permutations illustrate the API; use more for inference.
#' assessment <- assess.gaussian.trend.filtering(fit, permutations = 9,
#'     inner.folds = 3, outer.folds = 3, seed = 42,
#'     additional.tests = c("absolute.deviation", "total.variation"))
#' assessment$tests
#' assessment$prediction$improvement
#' @export
assess.gaussian.trend.filtering <- function(object, permutations = 999L,
    inner.folds = 5L, outer.folds = 5L, seed = 1L, predictive = TRUE, warn = TRUE,
    additional.tests = character()) {
    .assess.univariate.trend.filtering(object, permutations, inner.folds, outer.folds,
        seed, predictive, warn, additional.tests, family = "gaussian", call = match.call())
}

.gaussian.tf.assessment.variation <- function(fit) {
    .univariate.tf.assessment.variation(fit, identity)
}

.gaussian.tf.assessment.refit <- function(object, y, rows, folds, summarize = FALSE) {
    .univariate.tf.assessment.refit(object, y, rows, folds, summarize, gaussian = TRUE)
}

.gaussian.tf.assessment.prediction <- function(object, outer, inner) {
    n <- length(object$y); k <- length(inner)
    predictions <- data.frame(observation = seq_len(n), fold = outer, y = object$y,
        constant.mean = NA_real_, mean = NA_real_, outside.range = NA)
    records <- vector("list", k)
    for (j in seq_len(k)) {
        train <- which(outer != j); test <- which(outer == j)
        predictions$constant.mean[test] <- mean(object$y[train])
        result <- .gaussian.tf.assessment.refit(object, object$y, train, unname(inner[[j]]))
        record <- result$record
        record$training.range <- range(object$x[train])
        if (record$status == "complete") {
            estimate <- stats::predict(result$fit, object$x[test], lambda.index = record$selected.index)
            if (all(is.finite(estimate)) && all(is.finite((object$y[test]-estimate)^2))) {
                predictions$mean[test] <- as.numeric(estimate)
                predictions$outside.range[test] <- attr(estimate, "outside.range")
            } else record$status <- "nonfinite_prediction_or_loss"
        }
        records[[j]] <- record
    }
    predictions$constant.squared.error <- (object$y-predictions$constant.mean)^2
    predictions$squared.error <- (object$y-predictions$mean)^2
    predictions$constant.absolute.error <- abs(object$y-predictions$constant.mean)
    predictions$absolute.error <- abs(object$y-predictions$mean)
    complete <- all(is.finite(predictions$squared.error)) &&
        all(is.finite(predictions$constant.squared.error))
    scores <- data.frame(model = c("Constant mean", "Trend filtering"),
        mse = c(mean(predictions$constant.squared.error),
                if (complete) mean(predictions$squared.error) else NA_real_),
        mae = c(mean(predictions$constant.absolute.error),
                if (complete) mean(predictions$absolute.error) else NA_real_))
    list(status = if (complete) "complete" else "incomplete", scores = scores,
        improvement = c(mse = scores$mse[1L]-scores$mse[2L], mae = scores$mae[1L]-scores$mae[2L]),
        folds = data.frame(fold = seq_len(k), training = n-tabulate(outer, k),
            test = tabulate(outer, k),
            selected.lambda = vapply(records, `[[`, numeric(1), "selected.lambda"),
            status = vapply(records, `[[`, character(1), "status")),
        predictions = predictions, diagnostics = records)
}

#' @rdname assess.gaussian.trend.filtering
#' @param x A `gaussian.trend.filtering.assessment` to print.
#' @param ... Reserved for generic compatibility.
#' @method print gaussian.trend.filtering.assessment
#' @export
print.gaussian.trend.filtering.assessment <- function(x, ...) {
    cat("Order-2 Gaussian trend-filtering assessment\n")
    cat("Null: a common outcome distribution; independent, exchangeable observations\n")
    cat("Global test:", x$status, "\n")
    cat("Residual sum-of-squares reduction:", format(x$statistic), "\n")
    cat("Permutation p-value:", format(x$p.value), "\n")
    cat("Permutations:", x$permutation.summary$completed, "of",
        x$permutation.summary$requested, "completed\n")
    if (x$status == "incomplete_permutations")
        cat("Possible p-value range:", paste(format(x$p.value.bounds), collapse = " to "), "\n")
    cat("Observed curve summaries (response units):\n")
    print(x$curve.summaries)
    if (nrow(x$tests) > 1L) {
        cat("Requested tests (Holm adjustment across these statistics only):\n")
        print(x$tests[, c("test", "statistic", "status", "p.value", "p.value.holm")], row.names = FALSE)
    }
    if (!is.null(x$prediction)) {
        cat("Nested-CV prediction:", x$prediction$status, "\n")
        print(x$prediction$scores, row.names = FALSE)
        cat("Improvement (constant minus smoother):\n")
        print(x$prediction$improvement)
    }
    invisible(x)
}
