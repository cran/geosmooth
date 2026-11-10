#' Fit Order-2 Gaussian Trend Filtering
#'
#' Fit a continuous response with a piecewise-quadratic curve and an L1
#' penalty on third divided differences. Dense and banded native C++ solvers
#' share the order-2 geometry of \code{\link{fit.logistic.trend.filtering}}.
#'
#' @inheritParams fit.logistic.trend.filtering
#' @param y Finite numeric continuous response, with one value per observation.
#'   Constant responses are supported.
#' @param lambda.selection `"fixed"` returns the supplied penalty path;
#'   `"cv"` selects minimum held-out weighted mean squared error (MSE).
#'
#' @details With training predictor locations normalized to `[0,1]`, minimize
#'   \deqn{\frac{1}{2}\sum_i w_i(y_i-\eta_i)^2
#'   +\lambda\|D^{(3)}f\|_1
#'   +\frac{r_z}{2}\|\gamma\|^2
#'   +\frac{r_f}{2\sum_i w_i}\sum_i w_i(f_i-\bar f_w)^2.}
#'   Here eta is the curve plus standardized covariates and offset. The
#'   intercept is unpenalized. Lambda is on the summed-loss scale; the default
#'   ridge penalties are zero. `lambda=Inf` restricts the curve to a global
#'   quadratic, not a constant. At `lambda=0` with no covariates, offsets or
#'   ridge, predictions at distinct training locations are their weighted
#'   response means. Exact ties retain their individual loss contributions.
#'   Order is fixed at 2; no order-selection step is performed.
#'
#'   Gaussian refers to squared-error fitting. Fitting and CV do not estimate
#'   a noise variance, supply confidence intervals, or assume normality for
#'   the permutation assessment. Additional covariates and offsets are
#'   supported for fitting; \code{\link{assess.gaussian.trend.filtering}} supports only
#'   independent, unit-weight observations without adjustment.
#'
#'   Every CV split recomputes its training-only predictor range, knots and
#'   transformations. The default `cv.penalty.scale="training_mass"` scales
#'   lambda and ridge penalties by the training/full weight-sum ratio; `"none"`
#'   retains the earlier raw-penalty convention. This holds the mean-loss
#'   coefficients constant, not a common raw-coordinate operator across splits.
#'   Selection requires convergence on the full data and all
#'   folds. Equal CV scores select the first supplied lambda. Predictions
#'   outside the training range hold the curve at the nearest boundary;
#'   covariates and offsets continue to contribute. The banded backend scans
#'   all knots without thinning locations. Both solvers require finite
#'   arithmetic, stationarity, dual feasibility and L1 complementarity checks;
#'   failed fits remain explicit. Tolerances apply to mean weighted loss.
#'
#' @return A `gaussian.trend.filtering.fit` containing `order=2L`, the original
#'   data, penalty path, `fits`, convergence `diagnostics`, training `design`,
#'   and, when requested, `cv.table` with MSE and eligibility, fold-level loss
#'   sums and status, and `selected.index`. A fixed path has no selected index.
#'   `cv.penalty.scale`, `cv.fold.lambda`, `cv.fold.penalty.factors`, and
#'   `cv.fold.ridge` record the actual convention and training penalties.
#'   Standard prediction returns NA for unconverged candidates.
#' @references
#' Wang, Smola and Tibshirani (2014). The Falling Factorial Basis and Its
#' Statistical Applications. \url{https://arxiv.org/abs/1405.0558}
#' @seealso \code{\link{assess.gaussian.trend.filtering}}, \code{\link{fit.logistic.trend.filtering}}
#' @examples
#' x <- seq(0, 1, length.out = 40)
#' y <- sin(2*pi*x) + 0.1*cos(seq_along(x))
#' fit <- fit.gaussian.trend.filtering(x, y, lambda = c(Inf, 0.01),
#'     backend = "banded")
#' predict(fit, c(0.25, 0.75), lambda.index = 2)
#' @export
fit.gaussian.trend.filtering <- function(
    x, y, lambda, weights = NULL, covariates = NULL, offset = NULL,
    nuisance.ridge = 0, curve.ridge = 0,
    lambda.selection = c("fixed", "cv"), foldid = NULL,
    control = list(), warn = TRUE, backend = c("auto", "dense", "banded"),
    cv.penalty.scale = c("training_mass", "none")) {
    .fit.univariate.trend.filtering(x, y, lambda, weights, covariates, offset,
        nuisance.ridge, curve.ridge, lambda.selection, foldid, control, warn, backend,
        family = "gaussian", call = match.call(), cv.penalty.scale = cv.penalty.scale)
}

#' Predict from Gaussian Trend Filtering
#' @inheritParams predict.logistic.trend.filtering.fit
#' @param object A `gaussian.trend.filtering.fit` object.
#' @param type `"response"` or `"link"`; both return fitted means because the
#'   Gaussian model uses the identity link.
#' @return A numeric vector for one path index, otherwise a matrix. Attributes
#'   `outside.range` and `lambda` record boundary predictions and penalties.
#' @examples
#' x <- seq(0, 1, length.out = 20)
#' fit <- fit.gaussian.trend.filtering(x, 1+x-x^2, lambda = Inf)
#' predict(fit, c(0.25, 0.75))
#' head(fitted(fit))
#' head(residuals(fit))
#' print(fit)
#' @export
predict.gaussian.trend.filtering.fit <- function(
    object, newdata = NULL, covariates = NULL, offset = NULL,
    type = c("response", "link"), lambda.index = NULL,
    allow.unconverged = FALSE, ...) {
    .predict.univariate.trend.filtering(object, newdata, covariates, offset,
        type, lambda.index, allow.unconverged, family = "gaussian")
}

#' @rdname predict.gaussian.trend.filtering.fit
#' @export
fitted.gaussian.trend.filtering.fit <- function(object, ...) stats::predict(object, ...)

#' @rdname predict.gaussian.trend.filtering.fit
#' @export
residuals.gaussian.trend.filtering.fit <- function(object, ...) object$y - stats::fitted(object, ...)

#' @param x A fitted Gaussian trend-filtering object.
#' @rdname predict.gaussian.trend.filtering.fit
#' @method print gaussian.trend.filtering.fit
#' @export
print.gaussian.trend.filtering.fit <- function(x, ...) {
    cat("Order-2 Gaussian trend filtering (native C++; ", x$backend, ")\n", sep = "")
    cat(length(x$x), "observations;", length(x$design$knots), "distinct locations\n")
    print(x$diagnostics[, c("lambda", "converged", "status", "stationarity", "complementarity")],
          row.names = FALSE)
    if (!is.null(x$selected.index)) cat("Selected lambda:", x$lambda[x$selected.index], "\n")
    invisible(x)
}
