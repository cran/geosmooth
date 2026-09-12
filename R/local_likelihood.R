#' Fit a Local-Likelihood Density or Bernoulli Smoother
#'
#' Fits a local likelihood model at each evaluation point and reads off one raw
#' fitted value at that evaluation anchor.  Density and Bernoulli workflows
#' should convert the returned fitted field with \code{\link{normalize.density}}
#' when the target is a subject-occupation density.
#'
#' The \code{likelihood.family = "density"} branch uses a local exponential
#' tilt of the chart reference measure.  The \code{"bernoulli"} branch uses a
#' weighted local logistic likelihood and returns fitted probabilities at the
#' evaluation anchors.
#'
#' @param X Numeric matrix with one row per source/support point.
#' @param y Numeric response vector.  For \code{likelihood.family = "density"},
#'   this must be a nonnegative mass/intensity vector with positive total mass.
#'   For \code{"bernoulli"}, values must lie in \code{[0, 1]}.
#' @param X.eval Optional numeric matrix of evaluation points. Defaults to
#'   \code{X}.
#' @param likelihood.family Local likelihood family.
#' @param support.size Number of source points in each local support.
#' @param degree Local chart feature degree.  Supported values are 0, 1, and 2.
#'   The density branch omits an intercept because the intercept is not
#'   identifiable in the normalized local density likelihood.  The Bernoulli
#'   branch includes an intercept in its internal logistic feature map.
#' @param kernel Kernel name. Supported values are \code{"gaussian"},
#'   \code{"tricube"}, \code{"epanechnikov"}, and \code{"triangular"}.
#' @param bandwidth.multiplier Positive multiplier applied to the local support
#'   radius.
#' @param support.grid Optional integer candidate neighborhood sizes for
#'   cross-validation.  CV selection is currently implemented for
#'   \code{likelihood.family = "bernoulli"}.
#' @param degree.grid Optional local polynomial degree candidates for
#'   Bernoulli cross-validation.
#' @param kernel.grid Optional kernel candidates for Bernoulli
#'   cross-validation.
#' @param bandwidth.multiplier.grid Optional bandwidth-multiplier candidates
#'   for Bernoulli cross-validation.
#' @param lambda.ridge.grid Optional ridge-penalty candidates for Bernoulli
#'   cross-validation.
#' @param foldid Optional positive integer vector assigning source rows to
#'   cross-validation folds.
#' @param cv.folds Number of folds used when \code{foldid} is not supplied.
#' @param cv.seed Random seed used to generate folds when \code{foldid} is not
#'   supplied.
#' @param coordinate.method Local coordinate method. \code{"coordinates"} uses
#'   centered ambient coordinates. \code{"local.pca"} projects centered support
#'   points onto a local PCA basis.
#' @param chart.dim Local PCA dimension when \code{coordinate.method =
#'   "local.pca"}. If \code{NULL}, the dimension is
#'   \code{min(ncol(X), support.size - 1)}.  The deployable input-only
#'   policies \code{"auto"} and \code{"local.auto"} use the same local-PCA
#'   dimension diagnostics as \code{\link{fit.lps}}; \code{"auto"} resolves
#'   one global chart dimension, while \code{"local.auto"} resolves one
#'   dimension per evaluation anchor.
#' @param chart.dim.grid Optional candidate chart dimensions for experimental
#'   cross-validation. Numeric grids can be evaluated with
#'   \code{selection.strategy = "sparse_kd"}.
#' @param selection.strategy Candidate-selection strategy. \code{"grid"}
#'   evaluates the requested candidate grid. \code{"sparse_kd"} evaluates a
#'   sparse support-size by chart-dimension skeleton when
#'   \code{chart.dim.grid} is supplied. \code{"plateau_kd"} uses the same
#'   geometry-only support-size and chart-dimension plateau rule as
#'   \code{\link{fit.lps}}.
#' @param chart.dim.max Optional explicit maximum chart dimension for the
#'   sparse coupled candidate family.
#' @param design.margin Nonnegative integer feasibility margin used to screen
#'   local polynomial design size before candidate evaluation.
#' @param auto.chart.support.metric Support system used by \code{chart.dim =
#'   "auto"} or \code{"local.auto"}.  Local-likelihood smoothers currently use
#'   coordinate supports for both coordinate and operator diagnostics.
#' @param auto.chart.selection.metric Which auto chart-dimension diagnostic to
#'   use when both coordinate and operator summaries are requested.
#' @param quadrature.weights Optional positive reference-measure weights.
#'   Defaults to unit weights.
#' @param lambda.ridge Nonnegative ridge penalty on identifiable coefficients.
#' @param min.local.mass Minimum local kernel-weighted mass needed before
#'   attempting a higher-degree density fit.  For Bernoulli fits this is used
#'   only as diagnostic telemetry.
#' @param min.nonzero.mass Minimum number of locally weighted positive-mass
#'   source points needed before attempting a higher-degree local fit.
#' @param fallback Fallback policy for underidentified or failed higher-degree
#'   local fits.  Zero local mass always uses \code{"zero"}.
#' @param optimizer Optimizer for nonzero-degree density fits.  \code{"newton"}
#'   is implemented; \code{"optim"} is accepted and delegates to BFGS.
#' @param max.iter Maximum optimizer iterations.
#' @param tol Convergence tolerance for gradient norm and step norm.
#' @param return.details Logical; if \code{TRUE}, keep per-evaluation
#'   diagnostics.
#'
#' @return A list with class \code{"local_likelihood"} containing
#'   \code{fitted.values}, selected controls, and local solver diagnostics.
#' @examples
#' X <- matrix(seq(0, 1, length.out = 21), ncol = 1)
#' y <- exp(-20 * (X[, 1] - 0.5)^2)
#' y <- y / sum(y)
#' fit.local.likelihood(X, y, likelihood.family = "density",
#'                      support.size = 7L, degree = 0L)
#' @export
fit.local.likelihood <- function(
    X,
    y,
    X.eval = NULL,
    likelihood.family = c("density", "bernoulli"),
    support.size = min(15L, nrow(X)),
    degree = 1L,
    kernel = c("gaussian", "tricube", "epanechnikov", "triangular"),
    bandwidth.multiplier = 1,
    support.grid = NULL,
    degree.grid = NULL,
    kernel.grid = NULL,
    bandwidth.multiplier.grid = NULL,
    lambda.ridge.grid = NULL,
    foldid = NULL,
    cv.folds = 5L,
    cv.seed = 1L,
    coordinate.method = c("coordinates", "local.pca"),
    chart.dim = NULL,
    chart.dim.grid = NULL,
    selection.strategy = c("grid", "sparse_kd", "plateau_kd"),
    chart.dim.max = NULL,
    design.margin = 2L,
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator"),
    quadrature.weights = NULL,
    lambda.ridge = 1e-8,
    min.local.mass = sqrt(.Machine$double.eps),
    min.nonzero.mass = 1L,
    fallback = c("degree0", "zero", "chart_kernel", "na"),
    optimizer = c("newton", "optim"),
    max.iter = 50L,
    tol = 1e-8,
    return.details = TRUE) {

    prepared <- .local.chart.prepare.X.eval(X, X.eval)
    X <- prepared$X
    X.eval <- prepared$X.eval
    n <- nrow(X)
    p <- ncol(X)
    likelihood.family <- match.arg(likelihood.family)
    if (identical(likelihood.family, "density")) {
        y <- .local.chart.validate.response(y, n, nonnegative = TRUE)
        if (sum(y) <= 0) {
            stop("y must have positive total mass for likelihood.family = \"density\".",
                 call. = FALSE)
        }
    } else {
        y <- .local.chart.validate.response(y, n, nonnegative = FALSE)
        if (any(y < 0 | y > 1)) {
            stop("y must have values in [0, 1] for likelihood.family = \"bernoulli\".",
                 call. = FALSE)
        }
    }
    support.size <- .local.chart.validate.support.size(support.size, n)
    degree <- .local.likelihood.validate.degree(degree)
    kernel <- match.arg(kernel)
    cv.requested <- .local.chart.cv.requested(
        foldid = foldid,
        support.grid = support.grid,
        degree.grid = degree.grid,
        kernel.grid = kernel.grid,
        bandwidth.multiplier.grid = bandwidth.multiplier.grid,
        lambda.ridge.grid = lambda.ridge.grid,
        chart.dim.grid = chart.dim.grid
    )
    support.grid <- if (is.null(support.grid)) {
        support.size
    } else {
        .klp.clean.support.grid(support.grid, n)
    }
    degree.grid <- if (is.null(degree.grid)) {
        degree
    } else {
        .klp.clean.degree.grid(degree.grid)
    }
    kernel.grid <- if (is.null(kernel.grid)) {
        kernel
    } else {
        .klp.clean.kernel.grid(kernel.grid)
    }
    bandwidth.multiplier <- .local.chart.validate.positive.scalar(
        bandwidth.multiplier, "bandwidth.multiplier"
    )
    bandwidth.multiplier.grid <- if (is.null(bandwidth.multiplier.grid)) {
        bandwidth.multiplier
    } else {
        .klp.clean.bandwidth.multiplier.grid(bandwidth.multiplier.grid)
    }
    coordinate.method <- match.arg(coordinate.method)
    selection.strategy <- .coupled.kd.selection.strategy(selection.strategy)
    auto.chart.support.metric <- match.arg(auto.chart.support.metric)
    auto.chart.selection.metric <- match.arg(auto.chart.selection.metric)
    design.margin <- as.integer(design.margin)
    if (length(design.margin) != 1L || !is.finite(design.margin) ||
        design.margin < 0L) {
        stop("'design.margin' must be a nonnegative integer scalar.",
             call. = FALSE)
    }
    requested.chart.dim <- chart.dim
    chart.dim.info <- .local.chart.resolve.chart.dim(
        X = X,
        support.size = support.size,
        degree = degree,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric
    )
    chart.dim <- chart.dim.info$chart.dim
    quadrature.weights <- .local.chart.validate.quadrature(
        quadrature.weights, n
    )
    lambda.ridge <- .local.chart.validate.nonnegative.scalar(
        lambda.ridge, "lambda.ridge"
    )
    lambda.ridge.grid <- if (is.null(lambda.ridge.grid)) {
        lambda.ridge
    } else {
        .local.likelihood.clean.lambda.ridge.grid(lambda.ridge.grid)
    }
    min.local.mass <- .local.chart.validate.nonnegative.scalar(
        min.local.mass, "min.local.mass"
    )
    min.nonzero.mass <- .local.chart.validate.positive.integer(
        min.nonzero.mass, "min.nonzero.mass"
    )
    fallback <- match.arg(fallback)
    optimizer <- match.arg(optimizer)
    max.iter <- .local.chart.validate.positive.integer(max.iter, "max.iter")
    tol <- .local.chart.validate.positive.scalar(tol, "tol")

    cv.table <- NULL
    cv.predictions <- NULL
    coupled.plan <- NULL
    coupled.telemetry <- NULL
    if (isTRUE(cv.requested)) {
        foldid <- .klp.prepare.foldid(n, foldid, cv.folds, cv.seed)
        candidate.spec <- .coupled.kd.lps.candidate.spec(
            X = X,
            support.grid = support.grid,
            degree.grid = degree.grid,
            kernel.grid = kernel.grid,
            bandwidth.multiplier.grid = bandwidth.multiplier.grid,
            chart.dim = requested.chart.dim,
            chart.dim.grid = chart.dim.grid,
            coordinate.method = coordinate.method,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric,
            selection.strategy = selection.strategy,
            chart.dim.max = chart.dim.max,
            design.margin = design.margin
        )
        cand <- .local.likelihood.expand.lambda.ridge.candidates(
            candidate.spec$candidates,
            lambda.ridge.grid
        )
        coupled.plan <- candidate.spec$coupled.plan
        coupled.telemetry <- candidate.spec$telemetry
        if (!is.null(coupled.telemetry)) {
            coupled.telemetry$evaluated.candidates <- nrow(cand)
        }
        cv.result <- .local.likelihood.cv.table(
            X = X,
            y = y,
            foldid = foldid,
            cand = cand,
            likelihood.family = likelihood.family,
            coordinate.method = coordinate.method,
            chart.dim = requested.chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric,
            quadrature.weights = quadrature.weights,
            min.local.mass = min.local.mass,
            min.nonzero.mass = min.nonzero.mass,
            fallback = fallback,
            optimizer = optimizer,
            max.iter = max.iter,
            tol = tol
        )
        cv.table <- cv.result$cv.table
        cv.predictions <- cv.result$predictions
        score.column <- if (identical(likelihood.family, "bernoulli")) {
            "cv.brier.observed"
        } else {
            "cv.rmse.observed"
        }
        best.idx <- .local.chart.select.best.idx(cv.table,
                                                 score.column = score.column)
        selected.row <- cv.table[best.idx, , drop = FALSE]
        support.size <- selected.row$support.size[[1L]]
        degree <- selected.row$degree[[1L]]
        kernel <- selected.row$kernel[[1L]]
        bandwidth.multiplier <- selected.row$bandwidth.multiplier[[1L]]
        lambda.ridge <- selected.row$lambda.ridge[[1L]]
        requested.chart.dim <- if ("chart.dim" %in% names(selected.row) &&
                                   !is.na(selected.row$chart.dim[[1L]])) {
            .local.chart.decode.chart.dim(selected.row$chart.dim[[1L]])
        } else {
            requested.chart.dim
        }
        chart.dim.info <- .local.chart.resolve.chart.dim(
            X = X,
            support.size = support.size,
            degree = degree,
            coordinate.method = coordinate.method,
            chart.dim = requested.chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        chart.dim <- chart.dim.info$chart.dim
    }

    ne <- nrow(X.eval)
    fitted <- numeric(ne)
    status <- character(ne)
    local.mass <- numeric(ne)
    n.nonzero.local <- integer(ne)
    effective.support <- integer(ne)
    degree.used <- integer(ne)
    resolved.chart.dim <- integer(ne)
    iterations <- integer(ne)
    objective <- numeric(ne)
    gradient.norm <- numeric(ne)
    normalization.constant <- numeric(ne)
    bandwidth <- numeric(ne)
    fallback.used <- logical(ne)

    for (ii in seq_len(ne)) {
        fit.fun <- if (identical(likelihood.family, "density")) {
            .local.likelihood.density.fit
        } else {
            .local.likelihood.bernoulli.fit
        }
        local.chart.dim <- .local.chart.resolve.eval.chart.dim(
            X = X,
            x0 = X.eval[ii, ],
            support.size = support.size,
            degree = degree,
            coordinate.method = coordinate.method,
            chart.dim = requested.chart.dim,
            summary.dim = chart.dim
        )
        local <- fit.fun(
            X = X,
            y = y,
            x0 = X.eval[ii, ],
            support.size = support.size,
            degree = degree,
            kernel = kernel,
            bandwidth.multiplier = bandwidth.multiplier,
            coordinate.method = coordinate.method,
            chart.dim = local.chart.dim,
            quadrature.weights = quadrature.weights,
            lambda.ridge = lambda.ridge,
            min.local.mass = min.local.mass,
            min.nonzero.mass = min.nonzero.mass,
            fallback = fallback,
            optimizer = optimizer,
            max.iter = max.iter,
            tol = tol
        )
        fitted[[ii]] <- local$value
        status[[ii]] <- local$status
        local.mass[[ii]] <- local$local.mass
        n.nonzero.local[[ii]] <- local$n.nonzero.local
        effective.support[[ii]] <- local$effective.support
        degree.used[[ii]] <- local$degree.used
        resolved.chart.dim[[ii]] <- local$chart.dim
        iterations[[ii]] <- local$iterations
        objective[[ii]] <- local$objective
        gradient.norm[[ii]] <- local$gradient.norm
        normalization.constant[[ii]] <- local$normalization.constant
        bandwidth[[ii]] <- local$bandwidth
        fallback.used[[ii]] <- local$fallback.used
    }

    diagnostics <- list(
        status.counts = as.list(table(factor(
            status,
            levels = c("ok", "zero_mass_fallback", "degree0_fallback",
                       "chart_kernel_fallback", "optimizer_failed",
                       "nonfinite_fit", "constant_response")
        ))),
        fallback.count = sum(fallback.used),
        fallback.fraction = mean(fallback.used),
        min.local.mass = min(local.mass),
        median.local.mass = stats::median(local.mass),
        min.normalization.constant =
            .local.likelihood.finite.min(normalization.constant),
        median.normalization.constant =
            .local.likelihood.finite.median(normalization.constant),
        degree.used.summary = summary(degree.used),
        chart.dim = .local.chart.dimension.telemetry(
            chart.dim.info = chart.dim.info,
            chart.dim.by.anchor = resolved.chart.dim,
            source.path = "fit.local.likelihood.local.chart_resolution"
        )
    )
    if (isTRUE(return.details)) {
        diagnostics$per.eval <- data.frame(
            eval.index = seq_len(ne),
            status = status,
            M.local = local.mass,
            n.nonzero.local = n.nonzero.local,
            support.size = support.size,
            effective.support = effective.support,
            degree.requested = degree,
            degree.used = degree.used,
            lambda.ridge = lambda.ridge,
            iterations = iterations,
            objective = objective,
            gradient.norm = gradient.norm,
            normalization.constant = normalization.constant,
            bandwidth = bandwidth,
            raw.fitted = fitted,
            fallback.used = fallback.used
        )
    }

    out <- list(
        method.id = "local_likelihood",
        likelihood.family = likelihood.family,
        X = X,
        X.eval = X.eval,
        y = y,
        fitted.values = fitted,
        selected = list(
            likelihood.family = likelihood.family,
            support.size = support.size,
            degree = degree,
            kernel = kernel,
            bandwidth.multiplier = bandwidth.multiplier,
            coordinate.method = coordinate.method,
            requested.chart.dim = requested.chart.dim,
            chart.dim = chart.dim,
            auto.chart.dim = chart.dim.info$auto.chart.dim,
            auto.chart.dim.local = chart.dim.info$auto.chart.dim.local,
            chart.dim.mode = chart.dim.info$chart.dim.mode,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric,
            lambda.ridge = lambda.ridge,
            min.local.mass = min.local.mass,
            min.nonzero.mass = min.nonzero.mass,
            fallback = fallback,
            optimizer = optimizer,
            selection.strategy = selection.strategy,
            max.iter = max.iter,
            tol = tol,
            cv.rmse.observed = if (!is.null(cv.table) &&
                                   "cv.rmse.observed" %in%
                                       names(cv.table)) {
                min(cv.table$cv.rmse.observed, na.rm = TRUE)
            } else {
                NA_real_
            },
            cv.brier.observed = if (!is.null(cv.table) &&
                                    "cv.brier.observed" %in%
                                        names(cv.table)) {
                min(cv.table$cv.brier.observed, na.rm = TRUE)
            } else {
                NA_real_
            },
            cv.logloss.observed = if (!is.null(cv.table) &&
                                      "cv.logloss.observed" %in%
                                          names(cv.table)) {
                min(cv.table$cv.logloss.observed, na.rm = TRUE)
            } else {
                NA_real_
            }
        ),
        cv.table = cv.table,
        foldid = if (isTRUE(cv.requested)) foldid else NULL,
        cv.predictions = if (isTRUE(return.details)) cv.predictions else NULL,
        coupled.kd.candidate.plan = if (isTRUE(return.details)) {
            coupled.plan
        } else {
            NULL
        },
        selection.strategy = selection.strategy,
        quadrature.weights = quadrature.weights,
        diagnostics = diagnostics,
        call = match.call()
    )
    out$diagnostics$coupled.kd.selection <- coupled.telemetry
    class(out) <- c("local_likelihood", "list")
    out
}

.local.likelihood.expand.lambda.ridge.candidates <- function(cand,
                                                             lambda.ridge.grid) {
    lambda.ridge.grid <- .local.likelihood.clean.lambda.ridge.grid(
        lambda.ridge.grid
    )
    out <- cand[rep(seq_len(nrow(cand)), each = length(lambda.ridge.grid)),
                ,
                drop = FALSE]
    out$lambda.ridge <- rep(lambda.ridge.grid, times = nrow(cand))
    out$candidate.id <- seq_len(nrow(out))
    rownames(out) <- NULL
    out
}

.local.likelihood.cv.table <- function(X,
                                       y,
                                       foldid,
                                       cand,
                                       likelihood.family,
                                       coordinate.method,
                                       chart.dim,
                                       auto.chart.support.metric,
                                       auto.chart.selection.metric,
                                       quadrature.weights,
                                       min.local.mass,
                                       min.nonzero.mass,
                                       fallback,
                                       optimizer,
                                       max.iter,
                                       tol) {
    pred <- matrix(NA_real_, nrow = length(y), ncol = nrow(cand))
    folds <- sort(unique(foldid))
    candidate.chart.dim <- "chart.dim" %in% names(cand) &&
        any(!is.na(cand$chart.dim))
    for (fold in folds) {
        test <- which(foldid == fold)
        train <- which(foldid != fold)
        X.train <- X[train, , drop = FALSE]
        y.train <- y[train]
        q.train <- quadrature.weights[train]
        for (rr in seq_len(nrow(cand))) {
            chart.dim.fold <- if (identical(coordinate.method,
                                            "coordinates")) {
                NULL
            } else if (isTRUE(candidate.chart.dim)) {
                .local.chart.decode.chart.dim(cand$chart.dim[[rr]])
            } else {
                chart.dim
            }
            fit <- fit.local.likelihood(
                X = X.train,
                y = y.train,
                X.eval = X[test, , drop = FALSE],
                likelihood.family = likelihood.family,
                support.size = cand$support.size[[rr]],
                degree = cand$degree[[rr]],
                kernel = cand$kernel[[rr]],
                bandwidth.multiplier = cand$bandwidth.multiplier[[rr]],
                coordinate.method = coordinate.method,
                chart.dim = chart.dim.fold,
                auto.chart.support.metric = auto.chart.support.metric,
                auto.chart.selection.metric = auto.chart.selection.metric,
                quadrature.weights = q.train,
                lambda.ridge = cand$lambda.ridge[[rr]],
                min.local.mass = min.local.mass,
                min.nonzero.mass = min.nonzero.mass,
                fallback = fallback,
                optimizer = optimizer,
                max.iter = max.iter,
                tol = tol,
                return.details = FALSE
            )
            pred[test, rr] <- fit$fitted.values
        }
    }
    cv.table <- cand
    cv.table$cv.rmse.observed <- vapply(
        seq_len(ncol(pred)),
        function(j) .klp.rmse(pred[, j], y),
        numeric(1L)
    )
    if (identical(likelihood.family, "bernoulli")) {
        cv.table$cv.brier.observed <- cv.table$cv.rmse.observed^2
        cv.table$cv.logloss.observed <- vapply(
            seq_len(ncol(pred)),
            function(j) {
                if (!all(is.finite(pred[, j]))) {
                    return(Inf)
                }
                .klp.logloss(y, pmin(pmax(pred[, j], 1e-15), 1 - 1e-15))
            },
            numeric(1L)
        )
    }
    list(cv.table = cv.table, predictions = pred)
}

.local.likelihood.bernoulli.fit <- function(X,
                                            y,
                                            x0,
                                            support.size,
                                            degree,
                                            kernel,
                                            bandwidth.multiplier,
                                            coordinate.method,
                                            chart.dim,
                                            quadrature.weights,
                                            lambda.ridge,
                                            min.local.mass,
                                            min.nonzero.mass,
                                            fallback,
                                            optimizer,
                                            max.iter,
                                            tol) {
    support <- .local.chart.support(X, x0, support.size)
    idx <- support$idx
    coords <- .local.chart.coordinates(
        centered = support$centered,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim
    )
    distances <- sqrt(rowSums(coords^2))
    kernel.info <- .local.chart.kernel(
        distances = distances,
        kernel = kernel,
        bandwidth.multiplier = bandwidth.multiplier
    )
    r <- kernel.info$weights
    q <- quadrature.weights[idx]
    y.local <- y[idx]
    weights <- q * r
    local.mass <- sum(r * y.local)
    n.nonzero <- sum(y.local > 0 & r > 0)
    effective <- kernel.info$effective.support
    weight.total <- sum(weights)

    if (!is.finite(weight.total) || weight.total <= 0) {
        return(.local.likelihood.density.fallback(
            status = "zero_mass_fallback",
            value = 0,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero,
            effective.support = effective,
            degree.used = 0L,
            chart.dim = ncol(coords),
            iterations = 0L,
            objective = NA_real_,
            gradient.norm = NA_real_,
            normalization.constant = NA_real_,
            bandwidth = kernel.info$bandwidth,
            fallback.used = TRUE
        ))
    }

    raw.features <- .local.chart.feature.matrix(coords, degree)
    features <- cbind(intercept = 1, raw.features)
    if (degree == 0L || ncol(raw.features) == 0L ||
        effective <= ncol(features) + 1L) {
        return(.local.likelihood.bernoulli.degree0.fit(
            weights = weights,
            y.local = y.local,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero,
            effective.support = effective,
            chart.dim = ncol(coords),
            bandwidth = kernel.info$bandwidth,
            status = if (degree == 0L) "ok" else "degree0_fallback",
            fallback.used = degree != 0L
        ))
    }

    solved <- .local.likelihood.bernoulli.solve(
        features = features,
        weights = weights,
        y.local = y.local,
        lambda.ridge = lambda.ridge,
        optimizer = optimizer,
        max.iter = max.iter,
        tol = tol
    )
    if (!isTRUE(solved$converged) || !is.finite(solved$value)) {
        return(.local.likelihood.bernoulli.apply.fallback(
            fallback = fallback,
            weights = weights,
            y.local = y.local,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero,
            effective.support = effective,
            chart.dim = ncol(coords),
            bandwidth = kernel.info$bandwidth,
            iterations = solved$iterations,
            objective = solved$objective,
            gradient.norm = solved$gradient.norm
        ))
    }
    list(
        value = solved$value,
        status = "ok",
        local.mass = local.mass,
        n.nonzero.local = n.nonzero,
        effective.support = effective,
        degree.used = degree,
        chart.dim = ncol(coords),
        iterations = solved$iterations,
        objective = solved$objective,
        gradient.norm = solved$gradient.norm,
        normalization.constant = NA_real_,
        bandwidth = kernel.info$bandwidth,
        fallback.used = FALSE
    )
}

.local.likelihood.validate.degree <- function(degree) {
    degree <- .local.chart.validate.nonnegative.integer(degree, "degree")
    if (!degree %in% 0:2) {
        stop("degree must be 0, 1, or 2 for fit.local.likelihood().",
             call. = FALSE)
    }
    degree
}

.local.likelihood.finite.min <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) {
        return(NA_real_)
    }
    min(x)
}

.local.likelihood.finite.median <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) {
        return(NA_real_)
    }
    stats::median(x)
}

.local.likelihood.density.fit <- function(X,
                                          y,
                                          x0,
                                          support.size,
                                          degree,
                                          kernel,
                                          bandwidth.multiplier,
                                          coordinate.method,
                                          chart.dim,
                                          quadrature.weights,
                                          lambda.ridge,
                                          min.local.mass,
                                          min.nonzero.mass,
                                          fallback,
                                          optimizer,
                                          max.iter,
                                          tol) {
    support <- .local.chart.support(X, x0, support.size)
    idx <- support$idx
    coords <- .local.chart.coordinates(
        centered = support$centered,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim
    )
    distances <- sqrt(rowSums(coords^2))
    kernel.info <- .local.chart.kernel(
        distances = distances,
        kernel = kernel,
        bandwidth.multiplier = bandwidth.multiplier
    )
    r <- kernel.info$weights
    q <- quadrature.weights[idx]
    y.local <- y[idx]
    local.mass <- sum(r * y.local)
    n.nonzero <- sum(y.local > 0 & r > 0)
    base <- q * r

    if (!is.finite(local.mass) || local.mass < min.local.mass) {
        return(.local.likelihood.density.fallback(
            status = "zero_mass_fallback",
            value = 0,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero,
            effective.support = kernel.info$effective.support,
            degree.used = 0L,
            chart.dim = ncol(coords),
            iterations = 0L,
            objective = NA_real_,
            gradient.norm = NA_real_,
            normalization.constant = sum(base),
            bandwidth = kernel.info$bandwidth,
            fallback.used = TRUE
        ))
    }

    features <- .local.chart.feature.matrix(coords, degree)
    if (degree == 0L || ncol(features) == 0L ||
        n.nonzero < min.nonzero.mass ||
        kernel.info$effective.support <= ncol(features) + 1L) {
        return(.local.likelihood.degree0.fit(
            base = base,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero,
            effective.support = kernel.info$effective.support,
            chart.dim = ncol(coords),
            bandwidth = kernel.info$bandwidth,
            status = if (degree == 0L) "ok" else "degree0_fallback",
            fallback.used = degree != 0L
        ))
    }

    solved <- .local.likelihood.density.solve(
        features = features,
        base = base,
        target = r * y.local,
        local.mass = local.mass,
        lambda.ridge = lambda.ridge,
        optimizer = optimizer,
        max.iter = max.iter,
        tol = tol
    )
    if (!isTRUE(solved$converged) || !is.finite(solved$value)) {
        return(.local.likelihood.apply.fallback(
            fallback = fallback,
            base = base,
            r = r,
            q = q,
            y.local = y.local,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero,
            effective.support = kernel.info$effective.support,
            chart.dim = ncol(coords),
            bandwidth = kernel.info$bandwidth,
            iterations = solved$iterations,
            objective = solved$objective,
            gradient.norm = solved$gradient.norm
        ))
    }
    list(
        value = solved$value,
        status = "ok",
        local.mass = local.mass,
        n.nonzero.local = n.nonzero,
        effective.support = kernel.info$effective.support,
        degree.used = degree,
        chart.dim = ncol(coords),
        iterations = solved$iterations,
        objective = solved$objective,
        gradient.norm = solved$gradient.norm,
        normalization.constant = solved$normalization.constant,
        bandwidth = kernel.info$bandwidth,
        fallback.used = FALSE
    )
}

.local.likelihood.degree0.fit <- function(base,
                                          local.mass,
                                          n.nonzero.local,
                                          effective.support,
                                          chart.dim,
                                          bandwidth,
                                          status,
                                          fallback.used) {
    Z <- sum(base)
    value <- if (is.finite(Z) && Z > 0) local.mass / Z else NA_real_
    .local.likelihood.density.fallback(
        status = status,
        value = value,
        local.mass = local.mass,
        n.nonzero.local = n.nonzero.local,
        effective.support = effective.support,
        degree.used = 0L,
        chart.dim = chart.dim,
        iterations = 0L,
        objective = if (is.finite(Z) && Z > 0) local.mass * log(Z) else NA_real_,
        gradient.norm = 0,
        normalization.constant = Z,
        bandwidth = bandwidth,
        fallback.used = fallback.used
    )
}

.local.likelihood.density.fallback <- function(status,
                                               value,
                                               local.mass,
                                               n.nonzero.local,
                                               effective.support,
                                               degree.used,
                                               chart.dim,
                                               iterations,
                                               objective,
                                               gradient.norm,
                                               normalization.constant,
                                               bandwidth,
                                               fallback.used) {
    list(
        value = value,
        status = status,
        local.mass = local.mass,
        n.nonzero.local = n.nonzero.local,
        effective.support = effective.support,
        degree.used = degree.used,
        chart.dim = chart.dim,
        iterations = iterations,
        objective = objective,
        gradient.norm = gradient.norm,
        normalization.constant = normalization.constant,
        bandwidth = bandwidth,
        fallback.used = fallback.used
    )
}

.local.likelihood.apply.fallback <- function(fallback,
                                             base,
                                             r,
                                             q,
                                             y.local,
                                             local.mass,
                                             n.nonzero.local,
                                             effective.support,
                                             chart.dim,
                                             bandwidth,
                                             iterations,
                                             objective,
                                             gradient.norm) {
    if (identical(fallback, "zero")) {
        return(.local.likelihood.density.fallback(
            status = "optimizer_failed",
            value = 0,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero.local,
            effective.support = effective.support,
            degree.used = 0L,
            chart.dim = chart.dim,
            iterations = iterations,
            objective = objective,
            gradient.norm = gradient.norm,
            normalization.constant = sum(base),
            bandwidth = bandwidth,
            fallback.used = TRUE
        ))
    }
    if (identical(fallback, "chart_kernel")) {
        denom <- sum(q * r)
        numer <- sum(r * y.local)
        return(.local.likelihood.density.fallback(
            status = "chart_kernel_fallback",
            value = numer / denom,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero.local,
            effective.support = effective.support,
            degree.used = 0L,
            chart.dim = chart.dim,
            iterations = iterations,
            objective = objective,
            gradient.norm = gradient.norm,
            normalization.constant = denom,
            bandwidth = bandwidth,
            fallback.used = TRUE
        ))
    }
    if (identical(fallback, "na")) {
        return(.local.likelihood.density.fallback(
            status = "nonfinite_fit",
            value = NA_real_,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero.local,
            effective.support = effective.support,
            degree.used = 0L,
            chart.dim = chart.dim,
            iterations = iterations,
            objective = objective,
            gradient.norm = gradient.norm,
            normalization.constant = sum(base),
            bandwidth = bandwidth,
            fallback.used = TRUE
        ))
    }
    .local.likelihood.degree0.fit(
        base = base,
        local.mass = local.mass,
        n.nonzero.local = n.nonzero.local,
        effective.support = effective.support,
        chart.dim = chart.dim,
        bandwidth = bandwidth,
        status = "degree0_fallback",
        fallback.used = TRUE
    )
}

.local.likelihood.bernoulli.degree0.fit <- function(weights,
                                                    y.local,
                                                    local.mass,
                                                    n.nonzero.local,
                                                    effective.support,
                                                    chart.dim,
                                                    bandwidth,
                                                    status,
                                                    fallback.used) {
    denom <- sum(weights)
    value <- if (is.finite(denom) && denom > 0) {
        sum(weights * y.local) / denom
    } else {
        NA_real_
    }
    .local.likelihood.density.fallback(
        status = status,
        value = min(1, max(0, value)),
        local.mass = local.mass,
        n.nonzero.local = n.nonzero.local,
        effective.support = effective.support,
        degree.used = 0L,
        chart.dim = chart.dim,
        iterations = 0L,
        objective = NA_real_,
        gradient.norm = 0,
        normalization.constant = NA_real_,
        bandwidth = bandwidth,
        fallback.used = fallback.used
    )
}

.local.likelihood.bernoulli.apply.fallback <- function(fallback,
                                                       weights,
                                                       y.local,
                                                       local.mass,
                                                       n.nonzero.local,
                                                       effective.support,
                                                       chart.dim,
                                                       bandwidth,
                                                       iterations,
                                                       objective,
                                                       gradient.norm) {
    if (identical(fallback, "zero")) {
        return(.local.likelihood.density.fallback(
            status = "optimizer_failed",
            value = 0,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero.local,
            effective.support = effective.support,
            degree.used = 0L,
            chart.dim = chart.dim,
            iterations = iterations,
            objective = objective,
            gradient.norm = gradient.norm,
            normalization.constant = NA_real_,
            bandwidth = bandwidth,
            fallback.used = TRUE
        ))
    }
    if (identical(fallback, "na")) {
        return(.local.likelihood.density.fallback(
            status = "nonfinite_fit",
            value = NA_real_,
            local.mass = local.mass,
            n.nonzero.local = n.nonzero.local,
            effective.support = effective.support,
            degree.used = 0L,
            chart.dim = chart.dim,
            iterations = iterations,
            objective = objective,
            gradient.norm = gradient.norm,
            normalization.constant = NA_real_,
            bandwidth = bandwidth,
            fallback.used = TRUE
        ))
    }
    .local.likelihood.bernoulli.degree0.fit(
        weights = weights,
        y.local = y.local,
        local.mass = local.mass,
        n.nonzero.local = n.nonzero.local,
        effective.support = effective.support,
        chart.dim = chart.dim,
        bandwidth = bandwidth,
        status = "degree0_fallback",
        fallback.used = TRUE
    )
}

.local.likelihood.bernoulli.solve <- function(features,
                                              weights,
                                              y.local,
                                              lambda.ridge,
                                              optimizer,
                                              max.iter,
                                              tol) {
    m <- ncol(features)
    penalty <- c(0, rep(lambda.ridge, max(0L, m - 1L)))
    if (identical(optimizer, "optim")) {
        return(.local.likelihood.bernoulli.solve.optim(
            features = features,
            weights = weights,
            y.local = y.local,
            penalty = penalty,
            max.iter = max.iter,
            tol = tol
        ))
    }
    beta <- .local.likelihood.bernoulli.initial.beta(features, weights, y.local)
    state <- .local.likelihood.bernoulli.state(
        beta, features, weights, y.local, penalty
    )
    converged <- FALSE
    iter <- 0L
    for (iter in seq_len(max.iter)) {
        if (!is.finite(state$objective) ||
            !all(is.finite(state$gradient)) ||
            !all(is.finite(state$hessian))) {
            break
        }
        grad.norm <- sqrt(sum(state$gradient^2))
        if (grad.norm <= tol) {
            converged <- TRUE
            break
        }
        step <- tryCatch(
            solve(state$hessian, state$gradient),
            error = function(e) rep(NA_real_, length(beta))
        )
        if (any(!is.finite(step))) {
            break
        }
        step.scale <- 1
        accepted <- FALSE
        for (hh in seq_len(30L)) {
            candidate <- beta - step.scale * step
            candidate.state <- .local.likelihood.bernoulli.state(
                candidate, features, weights, y.local, penalty
            )
            if (is.finite(candidate.state$objective) &&
                candidate.state$objective <= state$objective + 1e-10) {
                beta <- candidate
                state <- candidate.state
                accepted <- TRUE
                break
            }
            step.scale <- step.scale / 2
        }
        if (!isTRUE(accepted)) {
            break
        }
        if (sqrt(sum((step.scale * step)^2)) <= tol) {
            converged <- TRUE
            break
        }
    }
    list(
        converged = converged,
        value = .local.likelihood.logistic(beta[[1L]]),
        iterations = as.integer(iter),
        objective = state$objective,
        gradient.norm = sqrt(sum(state$gradient^2)),
        normalization.constant = NA_real_
    )
}

.local.likelihood.bernoulli.initial.beta <- function(features,
                                                     weights,
                                                     y.local) {
    beta <- numeric(ncol(features))
    denom <- sum(weights)
    if (length(beta) < 1L || !is.finite(denom) || denom <= 0) {
        return(beta)
    }
    p0 <- sum(weights * y.local) / denom
    if (!is.finite(p0)) {
        return(beta)
    }
    p0 <- min(1 - 1e-6, max(1e-6, p0))
    beta[[1L]] <- log(p0 / (1 - p0))
    beta
}

.local.likelihood.bernoulli.solve.optim <- function(features,
                                                    weights,
                                                    y.local,
                                                    penalty,
                                                    max.iter,
                                                    tol) {
    fn <- function(beta) {
        .local.likelihood.bernoulli.state(
            beta, features, weights, y.local, penalty
        )$objective
    }
    gr <- function(beta) {
        .local.likelihood.bernoulli.state(
            beta, features, weights, y.local, penalty
        )$gradient
    }
    fit <- tryCatch(
        stats::optim(
            par = numeric(ncol(features)),
            fn = fn,
            gr = gr,
            method = "BFGS",
            control = list(maxit = max.iter, reltol = tol)
        ),
        error = function(e) NULL
    )
    if (is.null(fit)) {
        return(list(
            converged = FALSE,
            value = NA_real_,
            iterations = 0L,
            objective = NA_real_,
            gradient.norm = NA_real_,
            normalization.constant = NA_real_
        ))
    }
    state <- .local.likelihood.bernoulli.state(
        fit$par, features, weights, y.local, penalty
    )
    list(
        converged = fit$convergence == 0L,
        value = .local.likelihood.logistic(fit$par[[1L]]),
        iterations = as.integer(fit$counts[["function"]]),
        objective = state$objective,
        gradient.norm = sqrt(sum(state$gradient^2)),
        normalization.constant = NA_real_
    )
}

.local.likelihood.bernoulli.state <- function(beta,
                                              features,
                                              weights,
                                              y.local,
                                              penalty) {
    eta <- as.numeric(features %*% beta)
    p <- .local.likelihood.logistic(eta)
    objective <- sum(weights * (.local.likelihood.log1pexp(eta) -
                                    y.local * eta)) +
        0.5 * sum(penalty * beta^2)
    gradient <- colSums(features * (weights * (p - y.local))) +
        penalty * beta
    hess.weights <- weights * p * (1 - p)
    hessian <- crossprod(
        features,
        features * hess.weights
    ) + diag(penalty, nrow = length(beta), ncol = length(beta))
    list(
        objective = objective,
        gradient = as.numeric(gradient),
        hessian = hessian
    )
}

.local.likelihood.density.solve <- function(features,
                                            base,
                                            target,
                                            local.mass,
                                            lambda.ridge,
                                            optimizer,
                                            max.iter,
                                            tol) {
    m <- ncol(features)
    if (m == 0L) {
        Z <- sum(base)
        return(list(
            converged = TRUE,
            value = local.mass / Z,
            iterations = 0L,
            objective = local.mass * log(Z),
            gradient.norm = 0,
            normalization.constant = Z
        ))
    }
    if (identical(optimizer, "optim")) {
        return(.local.likelihood.density.solve.optim(
            features = features,
            base = base,
            target = target,
            local.mass = local.mass,
            lambda.ridge = lambda.ridge,
            max.iter = max.iter,
            tol = tol
        ))
    }
    beta <- numeric(m)
    state <- .local.likelihood.density.state(
        beta, features, base, target, local.mass, lambda.ridge
    )
    converged <- FALSE
    iter <- 0L
    for (iter in seq_len(max.iter)) {
        if (!is.finite(state$objective) ||
            !all(is.finite(state$gradient)) ||
            !all(is.finite(state$hessian))) {
            break
        }
        grad.norm <- sqrt(sum(state$gradient^2))
        if (grad.norm <= tol) {
            converged <- TRUE
            break
        }
        step <- tryCatch(
            solve(state$hessian, state$gradient),
            error = function(e) rep(NA_real_, length(beta))
        )
        if (any(!is.finite(step))) {
            break
        }
        step.scale <- 1
        accepted <- FALSE
        for (hh in seq_len(30L)) {
            candidate <- beta - step.scale * step
            candidate.state <- .local.likelihood.density.state(
                candidate, features, base, target, local.mass, lambda.ridge
            )
            if (is.finite(candidate.state$objective) &&
                candidate.state$objective <= state$objective + 1e-10) {
                beta <- candidate
                state <- candidate.state
                accepted <- TRUE
                break
            }
            step.scale <- step.scale / 2
        }
        if (!isTRUE(accepted)) {
            break
        }
        if (sqrt(sum((step.scale * step)^2)) <= tol) {
            converged <- TRUE
            break
        }
    }
    value <- .local.likelihood.anchor.value(local.mass, state$logZ)
    list(
        converged = converged,
        value = value,
        iterations = as.integer(iter),
        objective = state$objective,
        gradient.norm = sqrt(sum(state$gradient^2)),
        normalization.constant = exp(state$logZ)
    )
}

.local.likelihood.density.solve.optim <- function(features,
                                                  base,
                                                  target,
                                                  local.mass,
                                                  lambda.ridge,
                                                  max.iter,
                                                  tol) {
    fn <- function(beta) {
        .local.likelihood.density.state(
            beta, features, base, target, local.mass, lambda.ridge
        )$objective
    }
    gr <- function(beta) {
        .local.likelihood.density.state(
            beta, features, base, target, local.mass, lambda.ridge
        )$gradient
    }
    fit <- tryCatch(
        stats::optim(
            par = numeric(ncol(features)),
            fn = fn,
            gr = gr,
            method = "BFGS",
            control = list(maxit = max.iter, reltol = tol)
        ),
        error = function(e) NULL
    )
    if (is.null(fit)) {
        return(list(
            converged = FALSE,
            value = NA_real_,
            iterations = 0L,
            objective = NA_real_,
            gradient.norm = NA_real_,
            normalization.constant = NA_real_
        ))
    }
    state <- .local.likelihood.density.state(
        fit$par, features, base, target, local.mass, lambda.ridge
    )
    list(
        converged = fit$convergence == 0L,
        value = .local.likelihood.anchor.value(local.mass, state$logZ),
        iterations = as.integer(fit$counts[["function"]]),
        objective = state$objective,
        gradient.norm = sqrt(sum(state$gradient^2)),
        normalization.constant = exp(state$logZ)
    )
}

.local.likelihood.density.state <- function(beta,
                                            features,
                                            base,
                                            target,
                                            local.mass,
                                            lambda.ridge) {
    eta <- as.numeric(features %*% beta)
    log.weights <- rep(-Inf, length(base))
    positive.base <- base > 0
    log.weights[positive.base] <- log(base[positive.base]) +
        eta[positive.base]
    logZ <- .local.likelihood.logsumexp(log.weights)
    if (!is.finite(logZ)) {
        m <- length(beta)
        return(list(
            objective = Inf,
            gradient = rep(Inf, m),
            hessian = matrix(Inf, m, m),
            logZ = logZ
        ))
    }
    prob <- exp(log.weights - logZ)
    mean.features <- colSums(features * prob)
    weighted.features <- features * prob
    second <- crossprod(features, weighted.features)
    cov.features <- second - tcrossprod(mean.features)
    objective <- -sum(target * eta) + local.mass * logZ +
        0.5 * lambda.ridge * sum(beta^2)
    gradient <- -colSums(features * target) + local.mass * mean.features +
        lambda.ridge * beta
    hessian <- local.mass * cov.features +
        diag(lambda.ridge, nrow = length(beta), ncol = length(beta))
    list(
        objective = objective,
        gradient = as.numeric(gradient),
        hessian = hessian,
        logZ = logZ
    )
}

.local.likelihood.anchor.value <- function(local.mass, logZ) {
    if (!is.finite(logZ)) {
        return(NA_real_)
    }
    local.mass * exp(-logZ)
}

.local.likelihood.logsumexp <- function(x) {
    finite <- is.finite(x)
    if (!any(finite)) {
        return(-Inf)
    }
    mx <- max(x[finite])
    mx + log(sum(exp(x[finite] - mx)))
}

.local.likelihood.logistic <- function(x) {
    out <- numeric(length(x))
    pos <- x >= 0
    out[pos] <- 1 / (1 + exp(-x[pos]))
    ex <- exp(x[!pos])
    out[!pos] <- ex / (1 + ex)
    out
}

.local.likelihood.log1pexp <- function(x) {
    out <- numeric(length(x))
    pos <- x > 0
    out[pos] <- x[pos] + log1p(exp(-x[pos]))
    out[!pos] <- log1p(exp(x[!pos]))
    out
}
