# Pointwise variance and confidence bands for fixed-configuration LPS fits
# (Tier 4 / E4.1). Builds on the E0.2 linear-smoother identity: at a fixed
# configuration (singleton grids, explicit numeric chart dimension) the LPS
# fitted vector is linear in the response, yhat = S %*% y, with S independent
# of y. lps.smoother.matrix() extracts S analytically, row by row, from the
# same local solves the fit performed; lps.pointwise.band() derives
# Var(yhat_i) = sigma^2 * sum_j S_ij^2, df = tr(S),
# sigma.hat^2 = RSS / (n - tr(S)), and the band
# yhat_i +/- z_{(1+level)/2} * sigma * ||S_i.||_2.

.klp.uq.supported.design.basis <- "orthogonal.polynomial.drop"

.klp.uq.validate.fit <- function(object, require.square = FALSE,
                                 caller = "lps.smoother.matrix") {
    if (!inherits(object, "lps")) {
        stop("'", caller, "' requires a fitted \"lps\" object.",
             call. = FALSE)
    }
    if (!identical(object$outcome.family, "gaussian")) {
        stop("'", caller, "' supports outcome.family = \"gaussian\" only; ",
             "this fit used \"", object$outcome.family, "\".", call. = FALSE)
    }
    if (!identical(object$backend.used, "R")) {
        stop("'", caller, "' supports the R backend only; this fit resolved ",
             "to backend \"", object$backend.used, "\".", call. = FALSE)
    }
    if (!identical(object$design.basis, .klp.uq.supported.design.basis)) {
        stop("'", caller, "' supports design.basis = ",
             "\"orthogonal.polynomial.drop\" only; this fit used \"",
             object$design.basis, "\".", call. = FALSE)
    }
    if (!is.data.frame(object$cv.table) || nrow(object$cv.table) != 1L) {
        stop("'", caller, "' requires a fixed configuration: singleton ",
             "support.grid, degree.grid, and kernel.grid (a cv.table with ",
             "exactly one row), so that no data-driven selection occurred. ",
             "The CV-selected pipeline is not a linear smoother.",
             call. = FALSE)
    }
    mode <- .klp.chart.dim.mode(object$requested.chart.dim,
                                object$coordinate.method)
    if (!mode %in% c("ambient", "fixed")) {
        stop("'", caller, "' requires a y-independent, explicitly fixed ",
             "chart dimension: coordinate.method = \"coordinates\", or ",
             "\"local.pca\" with an explicit numeric chart.dim (never ",
             "\"auto\", \"local.auto\", or the NULL default). This fit's ",
             "chart-dimension mode is \"", mode, "\".", call. = FALSE)
    }
    if (!is.null(object$diagnostics$chart.dim$by.anchor) &&
        length(unique(object$diagnostics$chart.dim$by.anchor)) > 1L) {
        stop("'", caller, "' does not support per-anchor chart dimensions.",
             call. = FALSE)
    }
    chart.dim <- object$chart.dim
    if (length(chart.dim) != 1L || !is.numeric(chart.dim) ||
        !is.finite(chart.dim) || chart.dim < 1) {
        stop("'", caller, "' requires a stored numeric chart dimension; ",
             "got ", deparse(chart.dim), ".", call. = FALSE)
    }
    if (any(!is.finite(object$y))) {
        stop("'", caller, "' requires the stored training response to be ",
             "finite.", call. = FALSE)
    }
    if (require.square) {
        same <- identical(dim(object$X.eval), dim(object$X)) &&
            all(object$X.eval == object$X)
        if (!isTRUE(same)) {
            stop("'", caller, "' requires X.eval identical to X: df = tr(S) ",
                 "and RSS are defined on the square training-fit smoother. ",
                 "Refit with the default X.eval.", call. = FALSE)
        }
    }
    invisible(TRUE)
}

# Influence (smoother) row of one local fit, in support coordinates: the
# vector l with prediction = sum_j l_j y_j for the weighted local polynomial
# solve that .klp.fit.intercept.design() performs with
# design.basis = "orthogonal.polynomial.drop". Mirrors that routine branch for
# branch (ok-mask, orthogonal transform, ridge-multiplier loop, condition
# guard, fallbacks); the only new algebra is that, with prediction
# p' (B'WB + ridge)^{-1} B'W y, the influence row is
# l = W B (B'WB + ridge)^{-1} p, computed in the transformed basis actually
# solved. Fallback semantics are reproduced honestly: the weighted-mean
# fallback is the exact linear row w / sum(w); the "na" fallback is an all-NA
# row.
.klp.uq.influence.row <- function(z, weights, degree,
                                  design.drop.tol,
                                  ridge.multiplier.grid,
                                  ridge.condition.max,
                                  unstable.action) {
    support.size <- length(weights)
    fallback <- function() {
        if (identical(unstable.action, "na")) {
            return(rep(NA_real_, support.size))
        }
        total <- sum(weights)
        if (!is.finite(total) || total <= 0) {
            return(rep(NA_real_, support.size))
        }
        weights / total
    }
    ok <- is.finite(weights) & weights > 0
    design <- .klp.get.local.design(
        z = z,
        degree = degree,
        chart.dim = ncol(z),
        design.cache = new.env(parent = emptyenv())
    )
    design.ok <- rowSums(is.finite(design)) == ncol(design)
    ok <- ok & design.ok
    if (sum(ok) < 1L) return(fallback())
    prediction.row <- matrix(c(1, rep(0, ncol(design) - 1L)), nrow = 1L)
    transformed <- .klp.orthogonal.polynomial.transform(
        design = design[ok, , drop = FALSE],
        weights = weights[ok],
        prediction.rows = prediction.row,
        design.drop.tol = design.drop.tol
    )
    if (!isTRUE(transformed$ok)) return(fallback())
    basis <- transformed$design
    pred <- transformed$prediction.rows
    if (nrow(basis) < ncol(basis)) return(fallback())
    xw <- basis * sqrt(weights[ok])
    cross <- crossprod(xw)
    scale <- .klp.local.ridge.scale(cross)
    penalty.base <- diag(1, nrow = ncol(cross))
    for (rho in ridge.multiplier.grid) {
        ridge <- rho * scale
        normal <- cross + ridge * penalty.base
        cond <- .klp.local.design.condition(normal)
        if (is.finite(ridge.condition.max) &&
            (!is.finite(cond) || cond > ridge.condition.max)) {
            next
        }
        solved <- tryCatch(solve(normal, t(pred)),
                           error = function(e) NULL)
        if (is.null(solved) || !all(is.finite(solved))) next
        influence.ok <- weights[ok] * as.numeric(basis %*% solved)
        if (!all(is.finite(influence.ok))) next
        row <- numeric(support.size)
        row[ok] <- influence.ok
        return(row)
    }
    fallback()
}

#' Extract the Linear-Smoother Matrix of a Fixed-Configuration LPS Fit
#'
#' For a fixed configuration (singleton `support.grid`, `degree.grid`, and
#' `kernel.grid`, with an explicit numeric chart dimension) the LPS fitted
#' vector is linear in the response: `fitted = S %*% y` with `S` depending on
#' `X`, the kernel weights, and the configuration, but not on `y`. This
#' function reconstructs `S` analytically, one evaluation row at a time, by
#' rebuilding each local fit's support, kernel weights, local chart, and
#' design through the same internal routines `fit.lps()` used, and reading off
#' the influence row of the local weighted least-squares solve.
#'
#' The extraction refuses configurations where the linearity premise fails or
#' is unsupported: it requires `outcome.family = "gaussian"`, the R backend,
#' `design.basis = "orthogonal.polynomial.drop"`, singleton grids (no
#' data-driven selection: the CV-selected pipeline is *not* a linear
#' smoother), and `coordinate.method = "coordinates"` or `"local.pca"` with an
#' explicit numeric `chart.dim` (never `"auto"` or `"local.auto"`).
#'
#' Self-guard: before returning, the function verifies
#' `max(abs(S %*% y - fitted.values.raw)) <= check.tol` against the fit it was
#' given (and that the `NA` patterns agree), so any divergence between the
#' reconstruction and the estimator is a hard error, never a silent drift.
#'
#' Local fits that fell back are represented honestly: a weighted-mean
#' fallback (`unstable.action = "mean"`) contributes its exact linear row
#' `w / sum(w)`; an `unstable.action = "na"` non-fit contributes an all-`NA`
#' row.
#'
#' @param object A fitted `"lps"` object from [fit.lps()] at a fixed
#'   configuration (see Details).
#' @param check.tol Positive scalar: maximum allowed absolute discrepancy of
#'   the self-guard identity `S %*% y == fitted.values.raw`. Default `1e-10`,
#'   the program's algebraic tolerance.
#' @return A numeric matrix `S` with `nrow(object$X.eval)` rows and
#'   `nrow(object$X)` columns: row `i` holds the weights through which the
#'   training responses enter the prediction at evaluation point `i`.
#' @seealso [lps.pointwise.band()] for pointwise variances and confidence
#'   bands derived from `S`.
#' @examples
#' set.seed(1)
#' n <- 30
#' X <- matrix(runif(2 * n, -1, 1), ncol = 2)
#' y <- sin(pi * X[, 1]) + 0.1 * rnorm(n)
#' fit <- fit.lps(X, y, foldid = rep(1:2, length.out = n),
#'                support.grid = 12L, degree.grid = 1L,
#'                kernel.grid = "tricube", backend = "R",
#'                design.basis = "orthogonal.polynomial.drop",
#'                ridge.multiplier.grid = 0, ridge.condition.max = Inf,
#'                unstable.action = "na")
#' S <- lps.smoother.matrix(fit)
#' max(abs(S %*% y - fit$fitted.values))   # ~1e-15: the linear identity
#' sum(diag(S))                            # effective degrees of freedom
#' @export
lps.smoother.matrix <- function(object, check.tol = 1e-10) {
    .klp.uq.validate.fit(object, require.square = FALSE,
                         caller = "lps.smoother.matrix")
    check.tol <- .klp.validate.positive.scalar(check.tol, "check.tol")
    X <- object$X
    X.eval <- object$X.eval
    n.train <- nrow(X)
    n.eval <- nrow(X.eval)
    support.size <- min(as.integer(object$selected$support.size[[1L]]),
                        n.train)
    degree <- as.integer(object$selected$degree[[1L]])
    kernel <- object$selected$kernel[[1L]]
    S <- matrix(0, nrow = n.eval, ncol = n.train)
    for (i in seq_len(n.eval)) {
        center <- X.eval[i, , drop = TRUE]
        ordered <- .klp.local.order(
            X.train = X,
            center = center,
            support.size = support.size
        )
        weights <- .klp.kernel.weights(ordered$distances, kernel)
        if (!any(weights > 0)) weights[] <- 1
        z <- .klp.local.coordinates(
            X.support = X[ordered$index, , drop = FALSE],
            center = center,
            coordinate.method = object$coordinate.method,
            chart.dim = as.integer(object$chart.dim),
            local.chart.method = object$local.chart.method.effective,
            weights = weights,
            return.chart = FALSE
        )
        row <- .klp.uq.influence.row(
            z = z,
            weights = weights,
            degree = degree,
            design.drop.tol = object$design.drop.tol,
            ridge.multiplier.grid = object$ridge.multiplier.grid,
            ridge.condition.max = object$ridge.condition.max,
            unstable.action = object$unstable.action
        )
        if (anyNA(row)) {
            S[i, ] <- NA_real_
        } else {
            S[i, ordered$index] <- row
        }
    }
    fitted.check <- as.numeric(S %*% object$y)
    fitted.raw <- as.numeric(object$fitted.values.raw)
    if (!identical(is.na(fitted.check), is.na(fitted.raw))) {
        stop("lps.smoother.matrix() self-guard failed: the NA pattern of ",
             "S %*% y does not match the fit's fitted.values.raw. The ",
             "reconstruction does not reproduce this fit; do not use its ",
             "output.", call. = FALSE)
    }
    finite <- !is.na(fitted.raw)
    max.diff <- if (any(finite)) {
        max(abs(fitted.check[finite] - fitted.raw[finite]))
    } else {
        0
    }
    if (!is.finite(max.diff) || max.diff > check.tol) {
        stop("lps.smoother.matrix() self-guard failed: max |S %*% y - ",
             "fitted.values.raw| = ", format(max.diff), " exceeds check.tol ",
             "= ", format(check.tol), ". The reconstruction does not ",
             "reproduce this fit; do not use its output.", call. = FALSE)
    }
    S
}

#' Pointwise Variance and Confidence Band for a Fixed-Configuration LPS Fit
#'
#' Computes, from the analytically extracted linear-smoother matrix `S` of a
#' fixed-configuration LPS fit (see [lps.smoother.matrix()]), the pointwise
#' variance `Var(fitted_i) = sigma^2 * sum_j S_ij^2`, the effective degrees of
#' freedom `df = tr(S)`, the plug-in noise estimate
#' `sigma.hat^2 = RSS / (n - tr(S))`, and the pointwise confidence band
#' `fitted_i +/- z * sigma * ||S_i.||_2` with
#' `z = qnorm(1 - (1 - level) / 2)`.
#'
#' With `sigma` supplied (known noise standard deviation), the band uses it
#' directly (`sigma.source = "known"`); with `sigma = NULL` the plug-in
#' `sigma.hat` is used (`sigma.source = "plug.in"`). `sigma.hat` is reported
#' in both modes. The variance is purely the linear-smoother sampling
#' variance: the band makes no bias correction, so where bias dominates
#' (boundary, high curvature) undercoverage is expected and must be reported,
#' not masked.
#'
#' Requires `X.eval` identical to `X` (the square training-fit smoother, on
#' which `tr(S)` and `RSS` are defined), in addition to all restrictions of
#' [lps.smoother.matrix()]. If any evaluation point's local fit returned `NA`,
#' its variance/band entries are `NA`, and `df`, `rss`, and `sigma.hat` are
#' `NA` as well; supplying a known `sigma` still yields bands at the
#' unaffected points.
#'
#' @param object A fitted `"lps"` object from [fit.lps()] at a fixed
#'   configuration, with `X.eval` identical to `X`.
#' @param sigma Known noise standard deviation (positive scalar), or `NULL`
#'   (default) to use the plug-in `sigma.hat`.
#' @param level Confidence level of the band, a single number strictly
#'   between 0 and 1. Default `0.95`.
#' @param check.tol Passed to [lps.smoother.matrix()]'s self-guard.
#' @return A list of class `"lps.pointwise.band"` with named fields:
#'   `fitted` (the fit's raw fitted values), `se`, `variance`, `lower`,
#'   `upper`, `level`, `z`, `sigma` (the supplied known sigma, or `NA` in
#'   plug-in mode), `sigma.hat`, `sigma.source` (`"known"` or `"plug.in"`),
#'   `df` (`tr(S)`), `rss`, `n.train`, `smoother.row.norm` (`||S_i.||_2`),
#'   and `configuration` (the pinned fit configuration).
#' @seealso [lps.smoother.matrix()]
#' @examples
#' set.seed(1)
#' n <- 30
#' X <- matrix(runif(2 * n, -1, 1), ncol = 2)
#' y <- sin(pi * X[, 1]) + 0.1 * rnorm(n)
#' fit <- fit.lps(X, y, foldid = rep(1:2, length.out = n),
#'                support.grid = 12L, degree.grid = 1L,
#'                kernel.grid = "tricube", backend = "R",
#'                design.basis = "orthogonal.polynomial.drop",
#'                ridge.multiplier.grid = 0, ridge.condition.max = Inf,
#'                unstable.action = "na")
#' band.known <- lps.pointwise.band(fit, sigma = 0.1)
#' band.plugin <- lps.pointwise.band(fit)
#' band.known$df
#' band.plugin$sigma.hat
#' @export
lps.pointwise.band <- function(object, sigma = NULL, level = 0.95,
                               check.tol = 1e-10) {
    .klp.uq.validate.fit(object, require.square = TRUE,
                         caller = "lps.pointwise.band")
    if (!is.null(sigma)) {
        sigma <- .klp.validate.positive.scalar(sigma, "sigma")
    }
    if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
        level <= 0 || level >= 1) {
        stop("'level' must be a single number strictly between 0 and 1.",
             call. = FALSE)
    }
    S <- lps.smoother.matrix(object, check.tol = check.tol)
    n <- nrow(object$X)
    fitted <- as.numeric(object$fitted.values.raw)
    row.sq <- rowSums(S^2)
    df <- sum(diag(S))
    residuals <- object$y - fitted
    rss <- sum(residuals^2)
    sigma.hat <- if (is.finite(df) && is.finite(rss) && (n - df) > 0) {
        sqrt(rss / (n - df))
    } else {
        NA_real_
    }
    sigma.source <- if (is.null(sigma)) "plug.in" else "known"
    sigma.used <- if (is.null(sigma)) sigma.hat else sigma
    if (identical(sigma.source, "plug.in") && !is.finite(sigma.used)) {
        stop("the plug-in sigma.hat is unavailable for this fit (df = ",
             format(df), ", rss = ", format(rss), ", n = ", n, "); supply a ",
             "known 'sigma', or refit with a configuration whose local fits ",
             "all succeed and whose df is below n.", call. = FALSE)
    }
    z <- stats::qnorm(1 - (1 - level) / 2)
    variance <- sigma.used^2 * row.sq
    se <- sqrt(variance)
    out <- list(
        fitted = fitted,
        se = se,
        variance = variance,
        lower = fitted - z * se,
        upper = fitted + z * se,
        level = level,
        z = z,
        sigma = if (identical(sigma.source, "known")) sigma.used else
            NA_real_,
        sigma.hat = sigma.hat,
        sigma.source = sigma.source,
        df = df,
        rss = rss,
        n.train = n,
        smoother.row.norm = sqrt(row.sq),
        configuration = list(
            support.size = as.integer(object$selected$support.size[[1L]]),
            degree = as.integer(object$selected$degree[[1L]]),
            kernel = object$selected$kernel[[1L]],
            coordinate.method = object$coordinate.method,
            chart.dim = as.integer(object$chart.dim),
            design.basis = object$design.basis,
            ridge.multiplier.grid = object$ridge.multiplier.grid,
            ridge.condition.max = object$ridge.condition.max,
            unstable.action = object$unstable.action,
            backend.used = object$backend.used,
            outcome.family = object$outcome.family
        )
    )
    class(out) <- c("lps.pointwise.band", "list")
    out
}
