#' Test Constancy and Assess Logistic Trend-Filtering Prediction
#'
#' Assess whether an order-2 logistic trend filter detects a nonconstant binary
#' probability, and separately estimate predictive improvement over a constant.
#' This is an additional assessment, not an option that changes the fitted
#' objective or supplies confidence bands.
#'
#' @param object A `logistic.trend.filtering.fit` with independent binary
#'   observations, unit weights, no covariates and no offset. Predictor choice,
#'   penalty path, ridge penalties and solver settings must have been specified
#'   independently of the outcomes. Backend and explicit penalties are reused.
#' @param permutations Number of random outcome permutations, at least one.
#'   The default 999 gives a minimum possible p-value of 0.001. Small values in
#'   examples illustrate usage only.
#' @param inner.folds Number of balanced, outcome-independent tuning folds.
#'   Used for the global test and within each outer training set.
#' @param outer.folds Number of outer folds for predictive assessment.
#' @param seed Nonnegative integer random seed. The caller's RNG state is
#'   restored. Realized folds, permutation seeds and RNG kind are returned.
#' @param predictive Whether to run nested CV in addition to the global test.
#' @param warn Warn if the observed fit, any permutation, or an outer fit fails.
#' @param additional.tests Character vector selecting optional permutation tests:
#'   `"absolute.deviation"`, `"total.variation"`, or both. Default `character()`
#'   tests likelihood improvement only. Both curve summaries are always returned
#'   for the observed fit. Specify the tests before inspecting outcomes.
#'
#' @details The null is a common unknown Bernoulli probability at all observed
#'   predictor locations. The constant estimate is the sample event rate, not
#'   a uniformly weighted integral over the continuous predictor domain.
#'   Independent observations must be exchangeable under this null. Repeated
#'   participants, serial dependence, confounding adjustment and grouped
#'   binomial counts require other null procedures and are not supported here.
#'   Independence cannot be established from a fitted object.
#'
#'   The input object's selected penalty and fold labels are not reused: those
#'   labels may depend on outcomes. A new outcome-independent fold assignment
#'   selects a penalty from the supplied path, including a full-data refit.
#'   Every permutation repeats this complete selection procedure on the same
#'   folds with the same solver settings. The statistic is the nonnegative
#'   part of twice the ordinary binomial log-likelihood improvement over the
#'   constant. It is not a penalized-objective difference and is not compared
#'   with a chi-squared distribution. Lambda order determines selection ties.
#'   Changing the path after inspecting the assessment invalidates the stated
#'   test unless that additional selection is included in the randomization.
#'
#'   Absolute deviation is `mean(abs(p.hat(x) - mean(y)))`, evaluated at all
#'   observed predictor locations, including repeats. It averages over sampled
#'   observations, not uniformly over the predictor domain. A value of 0.04
#'   means an average absolute departure of four percentage points. Total
#'   variation is the sum of absolute probability changes along the fitted
#'   curve over `[min(x), max(x)]`. It includes every knot and interior turning
#'   point of the piecewise quadratic logit, up to floating-point precision;
#'   it is not a coarse-grid approximation or the trend-filtering penalty.
#'   For a monotone curve it equals the endpoint probability difference in
#'   absolute value; an oscillating curve can have total variation above one.
#'   These are descriptive summaries of the newly tuned fit, not uncertainty
#'   intervals or bias-corrected population effect estimates.
#'
#'   Optional tests use the same selected fits, folds and permutations as the
#'   likelihood test, with large values indicating departure from constancy.
#'   No additional model fits or random draws are needed. Each statistic has
#'   its own inclusive comparison tolerance and failure bounds. `tests`
#'   reports individual p-values and Holm-adjusted p-values across likelihood
#'   and the requested additional tests. Holm adjustment remains valid under
#'   dependence; it does not adjust across other predictors or outcomes.
#'   If any requested test is incomplete, all Holm values are unavailable.
#'   Choosing the smallest unadjusted p-value after inspecting several tests
#'   does not preserve the stated test level. The legacy scalar fields always
#'   describe the unadjusted likelihood test, regardless of `additional.tests`.
#'
#'   Outcomes are shuffled uniformly with the predictor fixed, preserving the
#'   total event count. With B successful permutations the Monte Carlo p-value
#'   is `(1 + sum(permuted.statistic >= observed.statistic))/(B + 1)`.
#'   Under exchangeability this is a valid conditional randomization test,
#'   subject to Monte Carlo discreteness. For numerical ties, comparisons allow
#'   a conservative tolerance of `1e-10 * max(1, abs(observed.statistic))`.
#'   A failed candidate is ineligible under the fitter's usual CV rule. If no
#'   candidate is eligible, that replicate fails. Failed replicates are never
#'   dropped or retried: the point p-value is NA and `p.value.bounds` gives the
#'   smallest/largest possible Monte Carlo p-value counting each missing
#'   statistic as a nonexceedance/exceedance. These bounds are not a confidence
#'   interval. An observed-fit failure leaves both bounds unavailable.
#'
#'   Predictive assessment uses separate outer folds. Each outer training set
#'   inherits the input fit's `cv.penalty.scale`. Under `"training_mass"`,
#'   its lambda grid and ridge penalties are multiplied by the outer/full
#'   sample-size ratio, and its inner folds apply their own training ratio.
#'   Observed and permuted full-data refits use the same unscaled reference
#'   grid. Older saved fits without this field retain `"none"` for reproduction.
#'   Each outer training set
#'   tunes the path using only its inner folds and predicts its outer test set.
#'   The constant predictor is that outer training set's event rate. Both use
#'   identical test observations. Pooled mean natural-log loss and Brier score
#'   are reported; improvements are constant loss minus trend-filtering loss,
#'   so positive values favor the smoother. These assess training on the outer
#'   training-set size, not the full sample. No standard error is computed from
#'   overlapping CV folds. Incomplete smoother predictions leave its overall
#'   scores and improvements unavailable rather than scoring a selected subset.
#'
#'   Computation is sequential. A complete run uses approximately
#'   `(permutations + 1) * (inner.folds + 1)` native path fits for the test,
#'   plus `outer.folds * (inner.folds + 1)` for prediction. Per-replicate
#'   diagnostics are retained, but full permuted fitted objects are not.
#'   Rejecting constancy is evidence of association, not necessarily
#'   nonlinearity, predictive usefulness or a causal effect.
#'
#' @return A `logistic.trend.filtering.assessment` list containing `statistic`,
#'   `p.value`, `p.value.bounds`, `status`, `null.probability`, the newly tuned
#'   observed `fit`, permutation statistics and diagnostic records, optional
#'   nested-CV `prediction`, realized `folds`, seeds and solver settings.
#'   `curve.summaries` is a named vector with `absolute.deviation` and
#'   `total.variation`; `curve.summary.details` records their conventions and
#'   any calculation failure. `tests` has one row per requested test, including
#'   likelihood, with statistics, status, completed/failed permutation counts,
#'   exceedances, comparison tolerances, raw p-values, lower/upper failure bounds
#'   and `p.value.holm`. `permutation.statistics` contains a `replicate` column
#'   and one column per tested statistic. The original `permutations` table
#'   retains its likelihood-only schema.
#'   `prediction$scores` contains model losses; `prediction$improvement`
#'   contains baseline-minus-smoother gains; `prediction$predictions` and
#'   `prediction$diagnostics` retain row-level and fold-level evidence.
#' @references
#' Ojala and Garriga (2010). Permutation Tests for Studying Classifier Performance.
#' \url{https://www.jmlr.org/papers/v11/ojala10a.html}
#'
#' Phipson and Smyth (2010). Permutation p-values should never be zero.
#' \doi{10.2202/1544-6115.1585}
#' @seealso \code{\link{fit.logistic.trend.filtering}}
#' @examples
#' x <- seq(0, 1, length.out = 60)
#' y <- rep(c(0, 1, 0, 1, 1), 12)
#' fit <- fit.logistic.trend.filtering(x, y, lambda = c(Inf, 0.01),
#'     backend = "banded", curve.ridge = 0.1)
#' # Only nine permutations for a quick illustration; use more for inference.
#' assessment <- assess.logistic.trend.filtering(fit, permutations = 9,
#'     inner.folds = 3, outer.folds = 3, seed = 42)
#' print(assessment)
#' assessment$prediction$improvement
#' assessment$curve.summaries
#' # Optional tests reuse the fits within an assessment run:
#' shape <- assess.logistic.trend.filtering(fit, permutations = 9,
#'     inner.folds = 3, seed = 42, predictive = FALSE,
#'     additional.tests = c("absolute.deviation", "total.variation"))
#' shape$tests
#' @export
assess.logistic.trend.filtering <- function(object, permutations = 999L,
    inner.folds = 5L, outer.folds = 5L, seed = 1L, predictive = TRUE, warn = TRUE,
    additional.tests = character()) {
    .assess.univariate.trend.filtering(object, permutations, inner.folds, outer.folds,
        seed, predictive, warn, additional.tests, family = "binomial", call = match.call())
}

# Keep the same failure policy and conservative tie rule for every statistic.
.logistic.tf.assessment.test <- function(name, reference, statistics) {
    valid <- is.finite(statistics)
    tolerance <- if (is.finite(reference)) 1e-10 * max(1, abs(reference)) else NA_real_
    r <- if (is.finite(reference)) sum(statistics[valid] >= reference-tolerance) else NA_integer_
    missing <- sum(!valid); B <- length(statistics)
    status <- if (!is.finite(reference)) "observed_statistic_failed" else
        if (missing) "incomplete_permutations" else "complete"
    data.frame(test = name, statistic = reference, status = status,
        requested = B, completed = sum(valid), failed = missing, exceedances = r,
        comparison.tolerance = tolerance, resolution = 1/(B+1),
        p.value = if (status == "complete") (1+r)/(B+1) else NA_real_,
        p.value.lower = (1+r)/(B+1), p.value.upper = (1+r+missing)/(B+1))
}

# A continuous quadratic logit is monotone on either side of its vertex.
# Interpolate in interval-relative coordinates to avoid unscaling potentially
# ill-conditioned global falling-factorial coefficients. Summing probability
# changes at endpoints and vertices gives total variation without a grid.
.logistic.tf.assessment.variation <- function(fit) {
    .univariate.tf.assessment.variation(fit, stats::plogis)
}

.univariate.tf.assessment.variation <- function(fit, inverse.link) {
    knots <- sort(unique(fit$x))
    left <- head(knots, -1L); right <- utils::tail(knots, -1L)
    midpoint <- left + (right-left)/2
    eta <- as.numeric(stats::predict(fit, c(knots, midpoint), type = "link"))
    if (any(!is.finite(eta))) stop("Nonfinite curve in total-variation calculation.")
    m <- length(knots)
    y0 <- eta[seq_len(m-1L)]; y1 <- eta[seq_len(m-1L)+1L]
    ym <- eta[m+seq_len(m-1L)]
    fraction <- (midpoint-left)/(right-left)
    # When no representable midpoint exists, the endpoints are the only
    # distinct evaluable locations. No division by a zero interval fraction.
    usable <- fraction > 0 & fraction < 1
    a <- b <- numeric(m-1L)
    a[usable] <- ((ym[usable]-y0[usable])/fraction[usable] -
                     (y1[usable]-y0[usable]))/(fraction[usable]-1)
    b[usable] <- y1[usable]-y0[usable]-a[usable]
    vertex <- rep(NA_real_, m-1L)
    curved <- usable & a != 0
    vertex[curved] <- -b[curved]/(2*a[curved])
    inside <- which(is.finite(vertex) & vertex > 0 & vertex < 1)
    p0 <- inverse.link(y0); p1 <- inverse.link(y1)
    change <- abs(p1-p0)
    if (length(inside)) {
        v <- vertex[inside]
        pv <- inverse.link(y0[inside] + v*(b[inside]+v*a[inside]))
        change[inside] <- abs(pv-p0[inside]) + abs(p1[inside]-pv)
    }
    if (any(!is.finite(a)) || any(!is.finite(b)) || any(!is.finite(change)))
        stop("Nonfinite arithmetic in total-variation calculation.")
    sum(change)
}

.logistic.tf.assessment.folds <- function(x, k) {
    if (k > length(x)) stop("Too many inner folds for an outer training set.", call. = FALSE)
    labels <- sample(rep(seq_len(k), length.out = length(x)))
    if (any(vapply(seq_len(k), function(j) length(unique(x[labels != j])) < 4L, logical(1))))
        stop("Each training split needs at least four distinct predictor locations; change fold counts or use more data.",
             call. = FALSE)
    labels
}

.logistic.tf.binary.loss <- function(y, eta) {
    ifelse(y == 1, pmax(-eta, 0), pmax(eta, 0)) + log1p(exp(-abs(eta)))
}

.logistic.tf.assessment.refit <- function(object, y, rows, folds, summarize = FALSE) {
    .univariate.tf.assessment.refit(object, y, rows, folds, summarize, gaussian = FALSE)
}

.univariate.tf.assessment.refit <- function(object, y, rows, folds, summarize, gaussian) {
    fitter <- if (gaussian) fit.gaussian.trend.filtering else fit.logistic.trend.filtering
    convention <- object$cv.penalty.scale %||% "none"
    factor <- .tf.cv.penalty.factor(length(rows), length(object$y), convention)
    fit <- tryCatch(fitter(object$x[rows], y[rows], object$lambda * factor,
        lambda.selection = "cv", foldid = folds, backend = object$backend,
        nuisance.ridge = object$nuisance.ridge * factor,
        curve.ridge = object$curve.ridge * factor, cv.penalty.scale = convention,
        control = object$control, warn = FALSE), error = function(e) e)
    record <- list(status = "fit_error", statistic = NA_real_, selected.lambda = NA_real_,
        selected.index = NA_integer_, eligible.candidates = 0L, message = NULL,
        diagnostics = NULL, cv.table = NULL, cv.fold.status = NULL,
        cv.fold.loss.sums = NULL, cv.penalty.scale = convention,
        reference.penalty.factor = factor, cv.fold.lambda = NULL,
        curve.summaries = c(absolute.deviation = NA_real_, total.variation = NA_real_),
        curve.summary.message = NULL)
    if (inherits(fit, "error")) {
        record$message <- conditionMessage(fit)
        return(list(fit = NULL, record = record))
    }
    fit$call[[1L]] <- as.name(if (gaussian) "fit.gaussian.trend.filtering" else "fit.logistic.trend.filtering")
    record$diagnostics <- fit$diagnostics
    record$cv.table <- fit$cv.table
    record$cv.fold.status <- fit$cv.fold.status
    record$cv.fold.loss.sums <- fit$cv.fold.loss.sums
    record$cv.fold.lambda <- fit$cv.fold.lambda
    record$eligible.candidates <- as.integer(sum(fit$cv.table$eligible))
    record$status <- "no_eligible_candidate"
    if (!is.null(fit$selected.index)) {
        eta <- as.numeric(stats::predict(fit, type = "link", lambda.index = fit$selected.index))
        if (all(is.finite(eta))) {
            yy <- y[rows]
            null.loss <- if (gaussian) sum((yy-mean(yy))^2)/2 else
                sum(.logistic.tf.binary.loss(yy, stats::qlogis(mean(yy))))
            fitted.loss <- if (gaussian) sum((yy-eta)^2)/2 else
                sum(.logistic.tf.binary.loss(yy, eta))
            if (!is.finite(null.loss) || !is.finite(fitted.loss)) {
                record$status <- "nonfinite_loss"
                return(list(fit = fit, record = record))
            }
            record$statistic <- max(0, 2 * (null.loss - fitted.loss))
            record$selected.lambda <- fit$lambda[fit$selected.index]
            record$selected.index <- fit$selected.index
            record$status <- "complete"
            if (summarize) {
                record$curve.summaries["absolute.deviation"] <-
                    mean(abs((if (gaussian) eta else stats::plogis(eta))-mean(yy)))
                variation.function <- if (gaussian) .gaussian.tf.assessment.variation else
                    .logistic.tf.assessment.variation
                variation <- tryCatch(variation.function(fit),
                                      error = function(e) e)
                if (inherits(variation, "error"))
                    record$curve.summary.message <- conditionMessage(variation)
                else record$curve.summaries["total.variation"] <- variation
            }
        } else record$status <- "nonfinite_prediction"
    }
    list(fit = fit, record = record)
}

.logistic.tf.assessment.prediction <- function(object, outer, inner) {
    n <- length(object$y); k <- length(inner)
    predictions <- data.frame(observation = seq_len(n), fold = outer, y = object$y,
        constant.probability = NA_real_, probability = NA_real_, eta = NA_real_,
        outside.range = NA)
    records <- vector("list", k)
    for (j in seq_len(k)) {
        train <- which(outer != j); test <- which(outer == j)
        predictions$constant.probability[test] <- mean(object$y[train])
        result <- .logistic.tf.assessment.refit(object, object$y, train, unname(inner[[j]]))
        record <- result$record
        record$training.range <- range(object$x[train])
        if (record$status == "complete") {
            eta <- stats::predict(result$fit, object$x[test], type = "link",
                                  lambda.index = record$selected.index)
            if (all(is.finite(eta))) {
                predictions$eta[test] <- as.numeric(eta)
                predictions$probability[test] <- stats::plogis(eta)
                predictions$outside.range[test] <- attr(eta, "outside.range")
            } else record$status <- "nonfinite_prediction"
        }
        records[[j]] <- record
    }
    predictions$constant.log.loss <- .logistic.tf.binary.loss(object$y,
        stats::qlogis(predictions$constant.probability))
    predictions$log.loss <- .logistic.tf.binary.loss(object$y, predictions$eta)
    predictions$constant.brier <- (object$y - predictions$constant.probability)^2
    predictions$brier <- (object$y - predictions$probability)^2
    complete <- all(is.finite(predictions$eta))
    scores <- data.frame(model = c("Constant probability", "Trend filtering"),
        log.loss = c(mean(predictions$constant.log.loss),
                     if (complete) mean(predictions$log.loss) else NA_real_),
        brier = c(mean(predictions$constant.brier),
                  if (complete) mean(predictions$brier) else NA_real_))
    list(status = if (complete) "complete" else "incomplete", scores = scores,
        improvement = c(log.loss = scores$log.loss[1L] - scores$log.loss[2L],
                        brier = scores$brier[1L] - scores$brier[2L]),
        folds = data.frame(fold = seq_len(k), training = tabulate(outer, k) * -1L + n,
            test = tabulate(outer, k),
            selected.lambda = vapply(records, `[[`, numeric(1), "selected.lambda"),
            status = vapply(records, `[[`, character(1), "status")),
        predictions = predictions, diagnostics = records)
}

#' @rdname assess.logistic.trend.filtering
#' @param x A `logistic.trend.filtering.assessment` to print.
#' @param ... Reserved for generic compatibility.
#' @method print logistic.trend.filtering.assessment
#' @export
print.logistic.trend.filtering.assessment <- function(x, ...) {
    cat("Order-2 logistic trend-filtering assessment\n")
    cat("Null: constant binary probability; independent, exchangeable observations\n")
    cat("Global test:", x$status, "\n")
    cat("Likelihood improvement statistic:", format(x$statistic), "\n")
    cat("Permutation p-value:", format(x$p.value), "\n")
    cat("Permutations:", x$permutation.summary$completed, "of",
        x$permutation.summary$requested, "completed\n")
    if (x$status == "incomplete_permutations")
        cat("Possible p-value range:", paste(format(x$p.value.bounds), collapse = " to "), "\n")
    if (!is.null(x$curve.summaries)) {
        cat("Observed probability-curve summaries:\n")
        print(x$curve.summaries)
    }
    if (!is.null(x$tests) && nrow(x$tests) > 1L) {
        cat("Requested tests (Holm adjustment across these statistics only):\n")
        print(x$tests[, c("test", "statistic", "status", "p.value", "p.value.holm")],
              row.names = FALSE)
    }
    if (!is.null(x$prediction)) {
        cat("Nested-CV prediction:", x$prediction$status, "\n")
        print(x$prediction$scores, row.names = FALSE)
        cat("Improvement (constant minus smoother):\n")
        print(x$prediction$improvement)
    }
    invisible(x)
}
