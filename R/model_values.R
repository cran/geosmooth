#' Extract Stored Fitted Values and Response Residuals
#'
#' Use `fitted(object)` for the stored response-scale fitted values, preserving
#' vector or matrix shape, row order, names, and missing values. It does not
#' recompute a fit or imply prediction at new points. For `harmonic.smoother()`
#' results it returns `harmonic_predictions`. The basic
#' `perform.harmonic.smoothing()` result remains an unclassed list, accessed
#' through `$harmonic_predictions`; density probability masses remain `$rho`.
#'
#' `residuals(object)` returns stored response residuals where available,
#' including missing-label masks. Otherwise it computes `y - fitted.values`
#' only when the evaluation coordinates equal the training coordinates.
#' Residuals at different evaluation points are undefined and cause an error,
#' even when training and evaluation have the same number of rows. Binomial
#' residuals are on the response (probability) scale, not Pearson or deviance
#' residuals. These accessors do not change which refit results are reusable.
#'
#' @param object A supported fitted model; see Methods below.
#' @param ... Reserved; unused arguments are rejected.
#' @return The stored vector or matrix, or response residuals of the same shape.
#' @seealso \code{\link{refit}}, \code{\link{smoother.matrix}}
#' @examples
#' X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
#' model <- fit.malps(X, X[, 1]^2, degree = 1L,
#'                    support.type = "knn", support.size = 8L)
#' head(fitted(model))
#' head(residuals(model))
#' @name model-values
NULL

.model.response.residuals <- function(object) {
    if (!is.null(object$residuals)) return(object$residuals)
    if (!is.null(object$X.eval) && !identical(object$X.eval, object$X)) {
        stop("Residuals require evaluation coordinates identical to training coordinates.", call. = FALSE)
    }
    if (is.null(object$y) || !identical(dim(object$y), dim(object$fitted.values)) ||
        length(object$y) != length(object$fitted.values)) {
        stop("No aligned training response is stored for this fit.", call. = FALSE)
    }
    object$y - object$fitted.values
}

#' @rdname model-values
#' @method fitted lps
#' @export
fitted.lps <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals lps
#' @export
residuals.lps <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted malps
#' @export
fitted.malps <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals malps
#' @export
residuals.malps <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted ps_lps
#' @export
fitted.ps_lps <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals ps_lps
#' @export
residuals.ps_lps <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted chart_kernel
#' @export
fitted.chart_kernel <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals chart_kernel
#' @export
residuals.chart_kernel <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted local_likelihood
#' @export
fitted.local_likelihood <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals local_likelihood
#' @export
residuals.local_likelihood <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted lpl_tf
#' @export
fitted.lpl_tf <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals lpl_tf
#' @export
residuals.lpl_tf <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted slpl_tf
#' @export
fitted.slpl_tf <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals slpl_tf
#' @export
residuals.slpl_tf <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted metric.graph.lowpass.fit
#' @export
fitted.metric.graph.lowpass.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals metric.graph.lowpass.fit
#' @export
residuals.metric.graph.lowpass.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted metric.graph.lowpass.refit
#' @export
fitted.metric.graph.lowpass.refit <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals metric.graph.lowpass.refit
#' @export
residuals.metric.graph.lowpass.refit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted ssrhe.hessian.fit
#' @export
fitted.ssrhe.hessian.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals ssrhe.hessian.fit
#' @export
residuals.ssrhe.hessian.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted ssrhe.hessian.refit
#' @export
fitted.ssrhe.hessian.refit <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals ssrhe.hessian.refit
#' @export
residuals.ssrhe.hessian.refit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted ssrhe.hessian.l1.fit
#' @export
fitted.ssrhe.hessian.l1.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals ssrhe.hessian.l1.fit
#' @export
residuals.ssrhe.hessian.l1.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted graph.trend.filtering.fit
#' @export
fitted.graph.trend.filtering.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals graph.trend.filtering.fit
#' @export
residuals.graph.trend.filtering.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted pttf.trend.filtering.fit
#' @export
fitted.pttf.trend.filtering.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$fitted.values
}

#' @rdname model-values
#' @method residuals pttf.trend.filtering.fit
#' @export
residuals.pttf.trend.filtering.fit <- function(object, ...) {
    .geosmooth.check.dots(...)
    .model.response.residuals(object)
}

#' @rdname model-values
#' @method fitted harmonic_smoother
#' @export
fitted.harmonic_smoother <- function(object, ...) {
    .geosmooth.check.dots(...)
    object$harmonic_predictions
}
