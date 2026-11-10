#' Summarize a Geometric Smoother or Probability-Mass Fit
#'
#' Summaries expose stored settings and diagnostics without recomputing geometry,
#' fitting, or selection. Console output is bounded; the summary's structured
#' fields retain selected controls, criterion tables, local diagnostic counts,
#' solver information, and density accounting for closer inspection.
#'
#' A finite result is not a convergence certificate. Where the algorithm records
#' no stopping certificate, convergence is explicitly unavailable. Local fallback
#' counts and failed evaluations are reported separately. Selection criteria are
#' conditional tuning scores; they are not an independent assessment of prediction
#' error. GCV, CV MSE, RMSE, and likelihood losses have different interpretations.
#'
#' @param object,x A fitted object, or a summary for its print method.
#' @param ... Reserved; unused arguments are rejected.
#' @return `summary()` returns a `summary.geosmooth_fit` list with `model`,
#'   `n.training`, `n.evaluation`, `n.responses`, `n.observed`, `nonfinite`,
#'   `status`, `convergence`, `prediction.domain`, `selected`, `criteria`,
#'   `diagnostics`, `solver`, `mass`, `accounting`, and `warnings` fields.
#'   Print methods return their input invisibly.
#' @examples
#' X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
#' fit <- fit.chart.kernel(X, X[, 1]^2, support.size = 8L)
#' fit
#' summary(fit)
#' mass <- normalize.density(c(1, 2, 3))
#' mass
#' summary(mass)$accounting
#' @name model-summary
NULL

.geosmooth.fit.summary <- function(object) {
    density <- inherits(object, "density_fit")
    values <- if (density) object$rho else object$fitted.values
    labels <- c(lps = "Local polynomial smoothing", malps = "Model-averaged local polynomial smoothing",
        ps_lps = "Prediction-synchronized local polynomial smoothing",
        chart_kernel = "Chart-kernel smoothing", local_likelihood = "Local likelihood",
        lpl_tf = "Local polynomial lifting trend filtering", slpl_tf = "Synchronized lifting trend filtering",
        density_fit = "Probability mass", metric.graph.lowpass.fit = "Metric graph low-pass",
        metric.graph.lowpass.refit = "Metric graph low-pass refit",
        ssrhe.hessian.fit = "Quadratic Hessian regression", ssrhe.hessian.refit = "Quadratic Hessian refit",
        ssrhe.hessian.l1.fit = "Hessian L1 regression", graph.trend.filtering.fit = "Graph trend filtering",
        pttf.trend.filtering.fit = "Parallel-transport trend filtering")
    cls <- intersect(class(object), names(labels))[1L]
    solver <- object$solver
    convergence <- if (is.list(solver) && !is.null(solver$converged) &&
                       length(solver$converged) == 1L && !is.na(solver$converged)) {
        if (isTRUE(solver$converged)) "converged at reported tolerances" else "NOT CONVERGED"
    } else "unavailable (no stored stopping certificate)"
    selected <- object$selected %||% list()
    for (key in c("lambda", "lambda1", "lambda2", "degree", "support.size", "kernel")) {
        if (is.null(selected[[key]]) && !is.null(object[[key]])) selected[[key]] <- object[[key]]
    }
    if (!is.null(object$gcv$eta.optimal)) selected$eta <- object$gcv$eta.optimal
    criteria <- list(cv = object$cv, cv.table = object$cv.table, gcv = object$gcv)
    domain <- if (density) "mass on the stored support" else if (inherits(object, "lps")) {
        "stored evaluation rows; predict() also accepts new coordinates"
    } else if (inherits(object, "malps") && !identical(object$support.metric, "graph.geodesic")) {
        "training rows; predict() also accepts new coordinates"
    } else if (inherits(object, c("chart_kernel", "local_likelihood"))) {
        "stored evaluation rows; supply X.eval when fitting"
    } else "stored training vertices only"
    structure(list(
        model = unname(labels[cls]), method = object$method.id %||% cls,
        n.training = if (!is.null(object$X)) nrow(object$X) else if (!is.null(object$y)) NROW(object$y) else NA_integer_,
        n.evaluation = NROW(values), n.responses = NCOL(values),
        n.observed = if (!is.null(object$observed)) sum(object$observed) else if (!is.null(object$y)) sum(is.finite(object$y)) else NA_integer_,
        nonfinite = sum(!is.finite(values)),
        status = object$status %||% if (any(!is.finite(values))) "nonfinite fitted values" else "finite output",
        convergence = convergence, prediction.domain = domain, selected = selected,
        criteria = criteria, diagnostics = object$diagnostics, solver = solver,
        mass = if (density) list(total = sum(values), nonzero = sum(values > 0), minimum = min(values)) else NULL,
        accounting = object$accounting, warnings = object$warnings %||% character()
    ), class = "summary.geosmooth_fit")
}

#' @rdname model-summary
#' @method print summary.geosmooth_fit
#' @export
print.summary.geosmooth_fit <- function(x, ...) {
    .geosmooth.check.dots(...)
    cat(x$model, "\n")
    cat("  Training rows:", x$n.training, "; evaluation rows:", x$n.evaluation,
        "; responses:", x$n.responses, "\n")
    cat("  Observed responses:", x$n.observed, "; nonfinite fitted entries:", x$nonfinite, "\n")
    cat("  Output:", x$status, "\n  Convergence:", x$convergence, "\n")
    keys <- intersect(c("support.size", "degree", "kernel", "lambda", "lambda1", "lambda2",
                        "lambda.sync", "lambda.ridge", "eta", "likelihood.family",
                        "cv.rmse.observed", "cv.brier.observed", "cv.logloss.observed"), names(x$selected))
    parts <- vapply(keys, function(k) {
        value <- x$selected[[k]]
        if (length(value) == 1L && is.atomic(value)) paste0(k, "=", format(value, digits = 5)) else ""
    }, character(1))
    if (any(nzchar(parts))) cat(paste(strwrap(paste("Selected:",
        paste(parts[nzchar(parts)], collapse = ", ")), width = 76, indent = 2, exdent = 4), collapse = "\n"), "\n")
    if (any(!vapply(x$criteria, is.null, logical(1)))) {
        cat("  Selection scores: stored in summary(fit)$criteria; conditional tuning scores.\n")
    }
    if (!is.null(x$diagnostics$fallback.count)) cat("  Local fallbacks:", x$diagnostics$fallback.count, "\n")
    counts <- x$diagnostics$status.counts
    if (!is.null(counts)) {
        counts <- unlist(counts)
        cat("  Local status:", paste(names(counts)[counts > 0], counts[counts > 0], sep = "=", collapse = ", "), "\n")
    }
    if (!is.null(x$mass)) {
        cat("  Total mass:", format(x$mass$total, digits = 6), "; positive support:", x$mass$nonzero, "\n")
        accounting <- unlist(x$accounting)
        if (length(accounting)) cat(paste(strwrap(paste("Accounting:",
            paste(names(accounting), format(accounting, digits=4), sep="=", collapse=", ")),
            width=76, indent=2, exdent=4), collapse="\n"), "\n")
    }
    cat("  Prediction scope:", x$prediction.domain, "\n")
    if (length(x$warnings)) cat("  Recorded warnings:", length(x$warnings), "(see summary(fit)$warnings)\n")
    invisible(x)
}

#' @rdname model-summary
#' @method summary lps
#' @export
summary.lps <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary malps
#' @export
summary.malps <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary ps_lps
#' @export
summary.ps_lps <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary chart_kernel
#' @export
summary.chart_kernel <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary local_likelihood
#' @export
summary.local_likelihood <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary lpl_tf
#' @export
summary.lpl_tf <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary slpl_tf
#' @export
summary.slpl_tf <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary metric.graph.lowpass.fit
#' @export
summary.metric.graph.lowpass.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary metric.graph.lowpass.refit
#' @export
summary.metric.graph.lowpass.refit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary ssrhe.hessian.fit
#' @export
summary.ssrhe.hessian.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary ssrhe.hessian.refit
#' @export
summary.ssrhe.hessian.refit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary ssrhe.hessian.l1.fit
#' @export
summary.ssrhe.hessian.l1.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary graph.trend.filtering.fit
#' @export
summary.graph.trend.filtering.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary pttf.trend.filtering.fit
#' @export
summary.pttf.trend.filtering.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method summary density_fit
#' @export
summary.density_fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .geosmooth.fit.summary(object)
}

#' @rdname model-summary
#' @method print lpl_tf
#' @export
print.lpl_tf <- function(x, ...) {
    .geosmooth.check.dots(...)
    print(.geosmooth.fit.summary(x))
    invisible(x)
}

#' @rdname model-summary
#' @method print slpl_tf
#' @export
print.slpl_tf <- function(x, ...) {
    .geosmooth.check.dots(...)
    print(.geosmooth.fit.summary(x))
    invisible(x)
}

#' @rdname model-summary
#' @method print density_fit
#' @export
print.density_fit <- function(x, ...) {
    .geosmooth.check.dots(...)
    print(.geosmooth.fit.summary(x))
    invisible(x)
}

#' @rdname model-summary
#' @method print ps_lps
#' @export
print.ps_lps <- function(x, ...) {
    .geosmooth.check.dots(...)
    print(.geosmooth.fit.summary(x))
    invisible(x)
}

#' @rdname model-summary
#' @method print chart_kernel
#' @export
print.chart_kernel <- function(x, ...) {
    .geosmooth.check.dots(...)
    print(.geosmooth.fit.summary(x))
    invisible(x)
}

#' @rdname model-summary
#' @method print local_likelihood
#' @export
print.local_likelihood <- function(x, ...) {
    .geosmooth.check.dots(...)
    print(.geosmooth.fit.summary(x))
    invisible(x)
}
