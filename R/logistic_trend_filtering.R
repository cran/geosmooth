#' Fit Order-2 Logistic Trend Filtering
#'
#' Fits a piecewise-quadratic logit curve with an L1 penalty on its third
#' divided differences. An R interface calls native C++ dense or banded
#' solvers. No Python or external optimizer is required.
#'
#' @param x Finite numeric predictor, with at least four distinct values.
#' @param y Numeric response in `[0, 1]`. Binary observations use unit weights;
#'   binomial proportions use trial counts as weights.
#' @param lambda Nonnegative penalties, or `Inf` for the quadratic limit.
#'   Required explicitly: no data-dependent default penalty is assumed.
#' @param weights Positive observation weights, default one.
#' @param covariates Optional finite numeric matrix of additional predictors,
#'   without an intercept. Columns are centered and scaled on training data.
#' @param offset Optional finite vector added to the linear predictor.
#' @param nuisance.ridge,curve.ridge Nonnegative explicit ridge penalties,
#'   default zero. See Details for their precise scaling.
#' @param lambda.selection `"fixed"` returns the supplied path without selecting
#'   a penalty; `"cv"` selects minimum held-out weighted binomial log loss.
#' @param foldid Integer fold labels, required for `lambda.selection="cv"`.
#'   Supply the same label for observations that must be held out together.
#' @param cv.penalty.scale `"training_mass"` (default) multiplies each fold's
#'   lambda and ridge penalties by its training/full weight-sum ratio. This
#'   preserves their coefficients relative to mean loss. `"none"` retains
#'   the earlier same-numerical-penalty convention. Fixed-path fits are unchanged.
#' @param control Named list overriding `tolerance=1e-9`,
#'   `gap.tolerance=1e-10`, `max.iterations=1000L`, or
#'   `max.workspace.bytes=1024^3`. Tolerances apply to mean weighted loss.
#' @param warn Logical; warn when a path fit fails convergence checks.
#' @param backend Native solver: `"auto"` uses `"banded"` above 256 distinct
#'   locations or when the dense workspace estimate exceeds its limit, and
#'   `"dense"` otherwise. The banded backend uses compact working
#'   sets and checks every knot; it does not thin the predictor locations.
#'
#' @details Predictor locations are normalized to `[0,1]` using the training
#'   range. Lambda refers to third divided differences on that normalized
#'   scale, not on the original units of x. The minimized objective is
#'   \deqn{\sum_i w_i[\log(1+\exp(\eta_i))-y_i\eta_i]
#'   +\lambda\|D^{(3)}f\|_1
#'   +\frac{r_z}{2}\|\gamma\|^2
#'   +\frac{r_f}{2\sum_i w_i}\sum_i w_i(f_i-\bar f_w)^2.}
#'   Here eta is the curve plus standardized covariates and offset, gamma
#'   contains standardized covariate coefficients, and the curve includes an
#'   unpenalized intercept. Third differences annihilate quadratic polynomials.
#'   Both backends preserve the falling-factorial penalty exactly. The dense
#'   solver uses a scaled global basis and a logarithmic barrier. The banded
#'   solver uses locally supported linear combinations of the same functions,
#'   a fully corrective working set, banded curve Hessians, and a small Schur
#'   complement for covariates/centering. Every stopping check scans all knots;
#'   omitted working-set coefficients remain exactly zero and violating knots
#'   are admitted until the global KKT checks pass. It does not approximate the
#'   penalty or restrict the permitted knot locations. Exact ties share a curve value but retain their own likelihood
#'   contributions, covariates, and offsets. No jittering or approximate thinning
#'   is performed. Extra covariates can be confounded with the curve; without
#'   identifying constraints, their individual coefficients are not effects
#'   identifiable separately from the curve.
#'
#'   Convergence requires finite arithmetic, small stationarity residual,
#'   dual feasibility and small L1 complementarity. The reported dual vector
#'   lies in `[-1,1]`; complementarity checks whether it is also an approximate
#'   subgradient at the returned coefficients. This is an approximate KKT check,
#'   not a separately computed primal-dual objective gap. Iteration limits and
#'   numerical failures remain explicit. Separation can prevent a finite
#'   minimizer, including in the unpenalized quadratic component. The classical
#'   Gaussian trend-filtering minimax result is not a guarantee for this model.
#'
#'   The dense backend is intended for moderate knot counts. The banded backend
#'   avoids observation-by-full-knot and full-knot-square arrays during fitting
#'   and prediction. Its storage is linear in observations and knots for a fixed
#'   number of covariates and penalties; computation also depends on working-set
#'   changes and Newton convergence. The conservative workspace estimate is
#'   checked before allocation. It is a solver workspace guard, not a process-
#'   wide memory limit. The quadratic initialization uses the dense solver only
#'   on the three polynomial columns and optional covariates.
#'
#'   CV recomputes training normalization, centering, and scaling in each fold.
#'   With the default `cv.penalty.scale="training_mass"`, the penalty
#'   coefficients relative to weighted mean loss are held constant. The
#'   predictor range and standardized covariates still belong to each training
#'   split; scaling does not hold a common raw-coordinate operator fixed.
#'   A candidate with any failed training fit is ineligible for selection.
#'   All fits are refitted on the complete training data. Prediction outside
#'   its range uses the nearest boundary curve value; covariates and offsets
#'   continue to contribute normally. This boundary policy is not extrapolation
#'   of a quadratic polynomial.
#'
#' @return An object of class `logistic.trend.filtering.fit`, with `fits`,
#'   `lambda`, `diagnostics`, training transformations, optional `cv.table` and
#'   `cv.fold.lambda`, `cv.fold.penalty.factors`, `cv.fold.ridge`, the selected
#'   `cv.penalty.scale`, and
#'   `selected.index`. Failed fits are preserved for diagnosis but standard
#'   fitted/prediction methods return NA for them unless explicitly permitted.
#' @seealso \code{\link{assess.logistic.trend.filtering}} for a global constancy test and
#'   nested-CV predictive improvement for independent binary observations.
#' @examples
#' x <- seq(0, 1, length.out = 40)
#' y <- rep(c(0, 1, 0, 1, 1), 8)
#' model <- fit.logistic.trend.filtering(x, y, lambda = c(1, 0.1))
#' predict(model, seq(0, 1, length.out = 10), lambda.index = 2)
#' @export
fit.logistic.trend.filtering <- function(
    x, y, lambda, weights = NULL, covariates = NULL, offset = NULL,
    nuisance.ridge = 0, curve.ridge = 0,
    lambda.selection = c("fixed", "cv"), foldid = NULL,
    control = list(), warn = TRUE, backend = c("auto", "dense", "banded"),
    cv.penalty.scale = c("training_mass", "none")) {
    .fit.univariate.trend.filtering(x, y, lambda, weights, covariates, offset,
        nuisance.ridge, curve.ridge, lambda.selection, foldid, control, warn, backend,
        family = "binomial", call = match.call(), cv.penalty.scale = cv.penalty.scale)
}

.logistic.tf.basis <- function(u, knots) {
    H <- cbind(u, u^2)
    for (j in 3L:(length(knots)-1L))
        H <- cbind(H, 0.5 * (u - knots[j-1L]) * (u - knots[j]) * (u > knots[j]))
    H
}

#' Predict from Logistic Trend Filtering
#' @param object A `logistic.trend.filtering.fit` object.
#' @param newdata Finite numeric predictor locations; default training locations.
#' @param covariates Additional predictor matrix, with the training column order.
#' @param offset Offset vector, required for new locations if used during fitting.
#' @param type `"response"` for probabilities or `"link"` for logits.
#' @param lambda.index Path indices; defaults to the CV selection if available,
#'   otherwise all supplied penalties.
#' @param allow.unconverged Explicitly allow predictions from failed fits.
#' @param ... Reserved for compatibility with generic methods.
#' @return A numeric vector for one path index, otherwise a matrix. Attributes
#'   `outside.range` and `lambda` identify boundary predictions and penalties.
#' @examples
#' x <- seq(0, 1, length.out = 24)
#' fit <- fit.logistic.trend.filtering(x, rep(c(0, 1), 12), lambda = Inf)
#' predict(fit, c(0.25, 0.75))
#' head(fitted(fit))
#' head(residuals(fit))
#' print(fit)
#' @export
predict.logistic.trend.filtering.fit <- function(
    object, newdata = NULL, covariates = NULL, offset = NULL,
    type = c("response", "link"), lambda.index = NULL,
    allow.unconverged = FALSE, ...) {
    .predict.univariate.trend.filtering(object, newdata, covariates, offset,
        type, lambda.index, allow.unconverged, family = "binomial")
}

#' @rdname predict.logistic.trend.filtering.fit
#' @export
fitted.logistic.trend.filtering.fit <- function(object, ...) stats::predict(object, ...)

#' @rdname predict.logistic.trend.filtering.fit
#' @export
residuals.logistic.trend.filtering.fit <- function(object, ...) object$y - stats::fitted(object, ...)

#' @param x A fitted logistic trend-filtering object.
#' @rdname predict.logistic.trend.filtering.fit
#' @method print logistic.trend.filtering.fit
#' @export
print.logistic.trend.filtering.fit <- function(x, ...) {
    cat("Order-2 logistic trend filtering (native C++; ",
        if (is.null(x$backend)) "dense" else x$backend, ")\n", sep = "")
    cat(length(x$x), "observations;", length(x$design$knots), "distinct locations\n")
    print(x$diagnostics[, c("lambda", "converged", "status", "stationarity", "complementarity")], row.names = FALSE)
    if (!is.null(x$selected.index)) cat("Selected lambda:", x$lambda[x$selected.index], "\n")
    invisible(x)
}
