# Common randomization, failure accounting and nested-CV orchestration.

.assess.univariate.trend.filtering <- function(object, permutations = 999L,
    inner.folds = 5L, outer.folds = 5L, seed = 1L, predictive = TRUE, warn = TRUE,
    additional.tests = character(), family = "binomial", call = NULL) {
    gaussian <- identical(family, "gaussian")
    expected.class <- if (gaussian) "gaussian.trend.filtering.fit" else "logistic.trend.filtering.fit"
    refit <- if (gaussian) .gaussian.tf.assessment.refit else .logistic.tf.assessment.refit
    primary <- if (gaussian) "squared.error" else "likelihood"
    if (!inherits(object, expected.class))
        stop("object must be a ", expected.class, ".", call. = FALSE)
    if (gaussian && (!identical(object$order, 2L) || !is.numeric(object$x) ||
        !is.numeric(object$y) || !is.null(dim(object$x)) || !is.null(dim(object$y)) ||
        length(object$x) != length(object$y) || length(object$y) < 4L ||
        any(!is.finite(object$x)) || any(!is.finite(object$y))))
        stop("object must contain finite aligned numeric x and y and order=2L.", call. = FALSE)
    if (gaussian && (!is.numeric(object$weights) ||
        length(object$weights) != length(object$y) || any(!is.finite(object$weights)) ||
        !is.matrix(object$covariates) || nrow(object$covariates) != length(object$y) ||
        !is.numeric(object$offset) || length(object$offset) != length(object$y) ||
        any(!is.finite(object$offset))))
        stop("object must contain aligned finite weights, covariates and offset.", call. = FALSE)
    if ((!gaussian && any(!object$y %in% c(0, 1))) || any(object$weights != 1) ||
        ncol(object$covariates) != 0L || isTRUE(object$offset.used) || any(object$offset != 0))
        stop(if (gaussian) "Assessment requires independent exchangeable continuous observations with unit weights, no covariates and no offset." else
             "Assessment requires independent binary observations with unit weights, no covariates and no offset.",
             call. = FALSE)
    integer.scalar <- function(value, name, minimum) {
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value != floor(value) || value < minimum || value >= .Machine$integer.max)
            stop(name, " must be an integer >= ", minimum, " and < .Machine$integer.max.", call. = FALSE)
        as.integer(value)
    }
    permutations <- integer.scalar(permutations, "permutations", 1L)
    inner.folds <- integer.scalar(inner.folds, "inner.folds", 2L)
    outer.folds <- integer.scalar(outer.folds, "outer.folds", 2L)
    seed <- integer.scalar(seed, "seed", 0L)
    choices <- c("absolute.deviation", "total.variation")
    if (!is.character(additional.tests) || anyNA(additional.tests) ||
        any(!additional.tests %in% choices) || anyDuplicated(additional.tests))
        stop("additional.tests must contain unique names from 'absolute.deviation' and 'total.variation'.",
             call. = FALSE)
    # Canonical ordering makes family membership independent of argument order.
    additional.tests <- choices[choices %in% additional.tests]
    for (name in c("predictive", "warn")) {
        value <- get(name)
        if (!is.logical(value) || length(value) != 1L || is.na(value))
            stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    n <- length(object$y)
    if (inner.folds > n || (predictive && outer.folds > n))
        stop("Fold counts must not exceed the number of observations.", call. = FALSE)
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (existed) previous <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (existed) assign(".Random.seed", previous, envir = .GlobalEnv) else
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv))
    set.seed(seed)
    streams <- sample.int(.Machine$integer.max, 3L)
    set.seed(streams[1L])
    tuning <- .logistic.tf.assessment.folds(object$x, inner.folds)
    outer <- inner <- NULL
    if (predictive) {
        set.seed(streams[2L])
        outer <- .logistic.tf.assessment.folds(object$x, outer.folds)
        inner <- lapply(seq_len(outer.folds), function(k) {
            rows <- which(outer != k)
            labels <- .logistic.tf.assessment.folds(object$x[rows], inner.folds)
            stats::setNames(labels, rows)
        })
    }
    set.seed(streams[3L])
    permutation.seeds <- sample.int(.Machine$integer.max, permutations, replace = TRUE)
    observed <- refit(object, object$y, seq_len(n), tuning,
                                            summarize = TRUE)
    reference <- observed$record$statistic
    records <- vector("list", permutations)
    for (b in seq_len(permutations)) {
        set.seed(permutation.seeds[b])
        shuffled <- object$y[sample.int(n)]
        result <- refit(object, shuffled, seq_len(n), tuning,
                                               summarize = length(additional.tests) > 0L)
        records[[b]] <- result$record
    }
    statistics <- vapply(records, `[[`, numeric(1), "statistic")
    valid <- is.finite(statistics)
    tolerance <- if (is.finite(reference)) 1e-10 * max(1, abs(reference)) else NA_real_
    exceedances <- if (is.finite(reference)) sum(statistics[valid] >= reference - tolerance) else NA_integer_
    missing <- sum(!valid)
    bounds <- if (is.finite(reference))
        c(lower = (1 + exceedances)/(permutations + 1),
          upper = (1 + exceedances + missing)/(permutations + 1)) else c(lower = NA_real_, upper = NA_real_)
    status <- if (!is.finite(reference)) "observed_fit_failed" else
        if (missing) "incomplete_permutations" else "complete"
    permutation.statistics <- data.frame(replicate = seq_len(permutations),
                                         likelihood = statistics)
    names(permutation.statistics)[2L] <- primary
    references <- c(stats::setNames(reference, primary), observed$record$curve.summaries[additional.tests])
    for (name in additional.tests)
        permutation.statistics[[name]] <- vapply(records, function(r)
            unname(r$curve.summaries[name]), numeric(1))
    tests <- do.call(rbind, lapply(names(references), function(name)
        .logistic.tf.assessment.test(name, unname(references[name]),
                                    permutation.statistics[[name]])))
    tests$p.value.holm <- if (all(tests$status == "complete"))
        stats::p.adjust(tests$p.value, method = "holm") else rep(NA_real_, nrow(tests))
    prediction <- if (!predictive) NULL else if (gaussian)
        .gaussian.tf.assessment.prediction(object, outer, inner) else
        .logistic.tf.assessment.prediction(object, outer, inner)
    if (warn && (any(tests$status != "complete") ||
        any(!is.finite(observed$record$curve.summaries)) ||
        (predictive && prediction$status != "complete")))
        warning("Assessment has failed fits or curve summaries; inspect tests, curve.summary.details, permutation diagnostics and prediction diagnostics. No failed replicate is discarded.",
                call. = FALSE)
    out <- list(call = call, statistic = reference,
        p.value = if (status == "complete") unname(bounds["upper"]) else NA_real_,
        p.value.bounds = bounds, status = status, null.probability = mean(object$y),
        fit = observed$fit, observed = observed$record,
        curve.summaries = observed$record$curve.summaries,
        curve.summary.details = list(domain = range(object$x),
            absolute.deviation.measure = "empirical observations, including repeated predictors",
            total.variation.method = "piecewise quadratic turning points on the probability scale",
            message = observed$record$curve.summary.message),
        tests = tests, permutation.statistics = permutation.statistics,
        permutations = data.frame(replicate = seq_len(permutations), statistic = statistics,
            selected.lambda = vapply(records, `[[`, numeric(1), "selected.lambda"),
            eligible.candidates = vapply(records, `[[`, integer(1), "eligible.candidates"),
            status = vapply(records, `[[`, character(1), "status")),
        permutation.diagnostics = records,
        permutation.summary = list(requested = permutations, completed = sum(valid),
            failed = missing, exceedances = exceedances, resolution = 1/(permutations + 1),
            comparison.tolerance = tolerance),
        prediction = prediction, folds = list(tuning = tuning, outer = outer, inner = inner),
        seeds = list(seed = seed, streams = streams, permutations = permutation.seeds,
                     rng.kind = RNGkind()),
        settings = list(lambda = object$lambda, backend = object$backend,
            cv.penalty.scale = object$cv.penalty.scale %||% "none",
            control = object$control, nuisance.ridge = object$nuisance.ridge,
            curve.ridge = object$curve.ridge, inner.folds = inner.folds,
            outer.folds = if (predictive) outer.folds else NULL,
            additional.tests = additional.tests, adjustment = "holm",
            null = "Independent Bernoulli observations with a common probability"))
    if (gaussian) {
        names(out)[names(out) == "null.probability"] <- "null.mean"
        out$curve.summary.details$total.variation.method <-
            "piecewise quadratic turning points on the response scale"
        out$settings$null <- "Independent exchangeable outcomes with a common distribution"
        out$settings$order <- 2L
    }
    class(out) <- if (gaussian) "gaussian.trend.filtering.assessment" else "logistic.trend.filtering.assessment"
    out
}
