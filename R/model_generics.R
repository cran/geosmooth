#' Refit A Model Using Its Stored Geometry
#'
#' Fit a new response while reusing the geometry or operators in a supported
#' fitted object. R selects the method from the class of \code{object}.
#'
#' @param object A fitted MALPS, local polynomial lifting trend filter,
#'   synchronized lifting trend filter, metric graph low-pass, quadratic
#'   Hessian-energy, or Hessian L1 model.
#' @param y New response. MALPS and the two lifting trend filters accept a
#'   vector; Hessian L1 accepts one response; metric graph low-pass and
#'   quadratic Hessian methods also accept a matrix with one response per
#'   column. MALPS and both Hessian methods reuse the stored response when
#'   \code{y} is omitted or \code{NULL}. The other methods require \code{y}.
#' @param ... Method-specific controls, described in the linked method help.
#'   Unrecognized arguments are errors.
#'
#' @details
#' Refitting does not automatically repeat the original model selection:
#' \describe{
#'   \item{\code{malps}}{\code{\link{refit.malps}} reuses supports and
#'     averaging weights; accepts observation \code{weights}. Robust local
#'     residual weights, if used, are recomputed for the new response.}
#'   \item{\code{lpl_tf}}{\code{\link{refit.lpl_tf}} reuses the residual
#'     operator and penalty, or accepts a new fixed \code{lambda}.}
#'   \item{\code{slpl_tf}}{\code{\link{refit.slpl_tf}} reuses both operators
#'     and penalties, or accepts fixed \code{lambda1} and \code{lambda2}.}
#'   \item{\code{metric.graph.lowpass.fit}}{
#'     \code{\link{refit.metric.graph.lowpass.fit}} reuses the spectral basis
#'     and smoothing parameter. Set \code{per.column.gcv = TRUE} to select a
#'     parameter for each response column.}
#'   \item{\code{ssrhe.hessian.fit}}{
#'     \code{\link{refit.ssrhe.hessian.fit}} reuses the Hessian operator and
#'     fixed penalties, with optional \code{lambda1}, \code{lambda2},
#'     \code{weights}, and \code{ridge}. This also handles CV- and GCV-selected
#'     fits, retaining their selected penalties by default.}
#'   \item{\code{ssrhe.hessian.l1.fit}}{
#'     \code{\link{refit.ssrhe.hessian.l1.fit}} reuses the derivative operator
#'     with fixed or CV-selected penalties; accepts \code{lambda.grid},
#'     \code{lambda.selection}, and solver controls.}
#' }
#' Use the original fitted object for repeated metric graph low-pass or
#' quadratic Hessian refits: their refit results are summaries, not reusable
#' input fits. Other supported families return reusable fitted objects.
#'
#' @section Migration:
#' The six former exported \code{refit.*()} functions are replaced by
#' \code{refit(object, y, ...)}. Replace named \code{fitted.model} and
#' \code{y.new} arguments with \code{object} and \code{y}, respectively.
#' Method-specific controls and return classes are preserved. Call the
#' generic rather than an individual S3 method; the former help names remain
#' available for navigation.
#'
#' @return A refitted object or refit summary, as documented by the method.
#' @seealso \code{\link{smoother.matrix}}
#' @examples
#' X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
#' fit <- fit.malps(X, X[, 1]^2, degree = 1L,
#'                  support.type = "knn", support.size = 8L)
#' updated <- refit(fit, y = X[, 1]^3)
#' head(predict(updated))
#' @export
refit <- function(object, y, ...) {
    UseMethod("refit")
}

#' @rdname refit
#' @export
refit.default <- function(object, y, ...) {
    stop("No refit() method for class: ", paste(class(object), collapse = ", "),
         ". See help('refit', package = 'geosmooth') for supported fits.",
         call. = FALSE)
}

#' Extract A Conditional Smoother Matrix
#'
#' Extract the linear map from responses to predictions for a supported
#' fitted model. R selects the method from the class of \code{object}.
#'
#' @param object An \code{"lps"} or \code{"malps"} fitted object.
#' @param ... Method-specific controls: \code{check.tol} for LPS;
#'   \code{max.n} and \code{allow.robust} for MALPS. Unrecognized arguments
#'   are errors.
#'
#' @details
#' \code{\link{smoother.matrix.lps}} requires a fixed-configuration Gaussian
#' LPS fit using the R backend, the orthogonal polynomial drop basis, singleton
#' tuning grids, and fixed chart dimension. It returns one row per evaluation
#' point and checks that the matrix reproduces the stored raw fitted values.
#'
#' \code{\link{smoother.matrix.malps}} returns a training-point matrix
#' conditional on stored supports, charts, and weights. Its default size
#' limit is \code{max.n = 1000L}. Robust fits are rejected unless
#' \code{allow.robust = TRUE}; that option freezes the final robust weights
#' and does not represent how those weights change with the response.
#'
#' Neither method represents uncertainty from repeating model selection.
#'
#' @section Migration:
#' Replace \code{lps.smoother.matrix(fit, ...)} and
#' \code{malps.smoother.matrix(fit, ...)} with
#' \code{smoother.matrix(fit, ...)}. Existing method-specific arguments and
#' restrictions are preserved. The former functions are no longer exported;
#' their help aliases lead to the corresponding methods.
#'
#' @return A numeric matrix mapping training responses to predictions.
#'   See the method help for dimensions and attributes.
#' @seealso \code{\link{refit}}, \code{\link{lps.pointwise.band}},
#'   \code{\link{malps.gcv}}
#' @examples
#' X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
#' fit <- fit.malps(X, X[, 1]^2, degree = 1L,
#'                  support.type = "knn", support.size = 8L)
#' S <- smoother.matrix(fit)
#' max(abs(as.vector(S %*% fit$y) - predict(fit)))
#' @export
smoother.matrix <- function(object, ...) {
    UseMethod("smoother.matrix")
}

#' @rdname smoother.matrix
#' @export
smoother.matrix.default <- function(object, ...) {
    stop("No smoother.matrix() method for class: ",
         paste(class(object), collapse = ", "),
         ". See help('smoother.matrix', package = 'geosmooth') for supported fits.",
         call. = FALSE)
}

.geosmooth.check.dots <- function(...) {
    dots <- list(...)
    if (length(dots)) {
        labels <- names(dots)
        if (is.null(labels)) labels <- rep("", length(dots))
        labels[!nzchar(labels)] <- "<unnamed>"
        stop("Unused arguments: ", paste(labels, collapse = ", "), call. = FALSE)
    }
    invisible(NULL)
}
