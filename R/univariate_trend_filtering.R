# Shared order-2 geometry, loss-family dispatch, CV and prediction.

.fit.univariate.trend.filtering <- function(
    x, y, lambda, weights = NULL, covariates = NULL, offset = NULL,
    nuisance.ridge = 0, curve.ridge = 0,
    lambda.selection = c("fixed", "cv"), foldid = NULL,
    control = list(), warn = TRUE, backend = c("auto", "dense", "banded"),
    family = "binomial", call = NULL,
    cv.penalty.scale = c("training_mass", "none")) {
    gaussian <- identical(family, "gaussian")
    lambda.selection <- match.arg(lambda.selection)
    cv.penalty.scale <- match.arg(cv.penalty.scale)
    backend.requested <- match.arg(backend)
    scalar <- function(v, name, positive = FALSE) {
        if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
            if (positive) v <= 0 else v < 0)
            stop(name, " must be a finite ", if (positive) "positive" else "nonnegative",
                 " scalar.", call. = FALSE)
        v
    }
    defaults <- list(tolerance = 1e-9, gap.tolerance = 1e-10,
                     max.iterations = 1000L, max.workspace.bytes = 1024^3)
    if (!is.list(control) || (length(control) &&
        (is.null(names(control)) || any(!names(control) %in% names(defaults)) ||
         anyDuplicated(names(control))))) stop("Unknown or unnamed control fields.", call. = FALSE)
    control <- utils::modifyList(defaults, control)
    for (nm in names(control)) scalar(control[[nm]], nm, TRUE)
    if (control$max.iterations != floor(control$max.iterations) ||
        control$max.iterations > .Machine$integer.max)
        stop("max.iterations must be a positive integer.", call. = FALSE)
    if (!is.logical(warn) || length(warn) != 1L || is.na(warn))
        stop("warn must be TRUE or FALSE.", call. = FALSE)
    scalar(nuisance.ridge, "nuisance.ridge"); scalar(curve.ridge, "curve.ridge")
    if (!is.numeric(x) || !is.null(dim(x)) || length(x) < 4L || any(!is.finite(x)))
        stop("x must be a finite numeric vector with at least four observations.", call. = FALSE)
    n <- length(x)
    if (!is.numeric(y) || !is.null(dim(y)) || length(y) != n ||
        any(!is.finite(y)) || (!gaussian && any(y < 0 | y > 1)))
        stop(if (gaussian) "y must be a finite numeric vector with length(x) entries." else
            "y must be a finite numeric vector in [0, 1], with length(x) entries.", call. = FALSE)
    if (!is.numeric(lambda) || !length(lambda) || anyNA(lambda) ||
        any(lambda < 0) || any(is.nan(lambda)))
        stop("lambda must contain nonnegative numbers or Inf.", call. = FALSE)
    weights <- if (is.null(weights)) rep(1, n) else weights
    offset.used <- !is.null(offset)
    offset <- if (is.null(offset)) numeric(n) else offset
    for (nm in c("weights", "offset")) {
        v <- get(nm)
        if (!is.numeric(v) || !is.null(dim(v)) || length(v) != n || any(!is.finite(v)))
            stop(nm, " must be a finite numeric vector with length(x) entries.", call. = FALSE)
    }
    if (any(weights <= 0) || !is.finite(sum(weights)))
        stop("weights must be positive with a finite sum.", call. = FALSE)
    if (is.null(covariates)) covariates <- matrix(numeric(), n, 0L)
    if (!is.matrix(covariates) || !is.numeric(covariates) ||
        nrow(covariates) != n || any(!is.finite(covariates)))
        stop("covariates must be a finite numeric matrix with length(x) rows.", call. = FALSE)
    knots <- sort(unique(x))
    if (length(knots) < 4L || !is.finite(diff(range(x))))
        stop("x needs at least four distinct values and a finite range.", call. = FALSE)
    lo <- min(x); span <- max(x) - lo; knots <- (knots - lo)/span
    if (any(diff(knots) <= 0))
        stop("Distinct x values collapsed during normalization.", call. = FALSE)
    p <- length(knots) + ncol(covariates)
    dense.workspace <- 8 * (6 * as.double(n) * p + 12 * as.double(p)^2 +
                               length(lambda) * (3 * as.double(n) + 3 * p))
    backend <- if (backend.requested == "auto") {
        if (length(knots) > 256L || dense.workspace > control$max.workspace.bytes) "banded" else "dense"
    } else backend.requested
    workspace <- if (backend == "dense") {
        dense.workspace
    } else {
        8 * (as.double(n) * (24 + 8 * ncol(covariates)) +
                 length(knots) * (64 + 8 * ncol(covariates)) +
                 16 * as.double(ncol(covariates) + 4)^2 +
                 length(lambda) * (4 * as.double(n) + 4 * p))
    }
    if (!is.finite(workspace) || workspace > control$max.workspace.bytes)
        stop(backend, " solver workspace estimate exceeds max.workspace.bytes.", call. = FALSE)
    u <- (x - lo)/span; w <- weights/sum(weights)
    if (backend == "dense") {
        H <- .logistic.tf.basis(u, knots)
        center <- as.vector(crossprod(w, H))
        H <- sweep(H, 2L, center)
        scales <- sqrt(colSums(H^2 * w))
        if (any(!is.finite(scales)) || any(scales <= 0))
            stop("Degenerate falling-factorial basis.", call. = FALSE)
        H <- sweep(H, 2L, scales, "/")
    } else {
        geometry <- logistic_tf_structured_geometry(as.numeric(u), as.numeric(knots), as.numeric(weights))
        center <- geometry$center; scales <- geometry$scales
    }
    zcenter <- as.vector(crossprod(w, covariates))
    Z <- sweep(covariates, 2L, zcenter)
    zscale <- sqrt(colSums(Z^2 * w))
    if (any(!is.finite(zcenter)) || any(!is.finite(zscale)))
        stop("Nonfinite covariate centering or scaling; rescale covariates first.", call. = FALSE)
    if (any(zscale <= 0)) stop("covariates must not contain constant columns.", call. = FALSE)
    Z <- sweep(Z, 2L, zscale, "/")
    nz <- ncol(Z)
    if (backend == "dense") {
        A <- cbind(1, Z, H)
        # Keep ridge residuals instead of forming their Gram matrix. Large,
        # cancelling falling-factorial coefficients make b'Qb inaccurate even
        # when the fitted curve (and its ridge penalty) is well behaved.
        ridge.factor <- matrix(0, (if (curve.ridge > 0) n else 0L) + nz, p)
        if (nz) ridge.factor[seq_len(nz), 1L + seq_len(nz)] <- diag(sqrt(nuisance.ridge), nz)
        hi <- (nz + 2L):p
        if (curve.ridge > 0) ridge.factor[nz + seq_len(n), hi] <- H * sqrt(curve.ridge * w)
        penalty <- c(numeric(nz + 3L), 1/scales[-c(1L, 2L)])
        if (any(!is.finite(A)) || any(!is.finite(ridge.factor)) || any(!is.finite(penalty)))
            stop("Nonfinite transformed design or penalty.", call. = FALSE)
    }
    if (!gaussian && (sum(weights * y) == 0 || sum(weights * (1-y)) == 0)) {
        stop("No finite unpenalized intercept: y contains only one outcome.", call. = FALSE)
    }
    # Fit stronger penalties first for warm starts, but preserve caller order.
    ord <- order(lambda, decreasing = TRUE)
    fits <- if (backend == "dense") {
        native <- if (gaussian) gaussian_tf_native else logistic_tf_native
        native(A, ridge.factor, as.numeric(y), as.numeric(weights),
            as.numeric(offset), penalty, as.numeric(lambda[ord]),
            control$tolerance, control$gap.tolerance, as.integer(control$max.iterations))
    } else {
        native <- if (gaussian) gaussian_tf_banded_native else logistic_tf_banded_native
        native(as.numeric(u), as.numeric(knots), Z,
            as.numeric(y), as.numeric(weights), as.numeric(offset),
            as.numeric(lambda[ord]), nuisance.ridge, curve.ridge,
            control$tolerance, control$gap.tolerance, as.integer(control$max.iterations))
    }
    fits <- fits[order(ord)]
    diagnostics <- do.call(rbind, lapply(seq_along(fits), function(i) {
        f <- fits[[i]]
        data.frame(lambda = lambda[i], converged = f$converged, status = f$status,
            objective = f$objective, stationarity = f$stationarity,
            complementarity = f$complementarity, dual.feasibility = f$dual.feasibility,
            iterations = f$iterations, barrier.stages = f$barrier.stages,
            barrier.parameter = f$barrier.parameter, backtracks = f$backtracks,
            max.damping = f$max.damping, max.abs.eta = f$max.abs.eta)
    }))
    diagnostics$backend <- backend
    diagnostics$working.knots <- if (backend == "banded")
        vapply(fits, function(f) length(f$active.knots), integer(1)) else NA_integer_
    diagnostics$initialization.iterations <- if (backend == "banded")
        vapply(fits, function(f) as.integer(f$initialization.iterations), integer(1)) else 0L
    out <- list(call = call, backend = backend, x = x, y = y, weights = weights, offset = offset,
        covariates = covariates, offset.used = offset.used, lambda = lambda,
        fits = fits, diagnostics = diagnostics, selected.index = NULL,
        cv.table = NULL, cv.penalty.scale = cv.penalty.scale,
        control = control, nuisance.ridge = nuisance.ridge,
        curve.ridge = curve.ridge, workspace.estimate.bytes = workspace,
        design = list(lo = lo, span = span, knots = knots, center = center,
            scales = scales, zcenter = zcenter, zscale = zscale))
    class(out) <- if (gaussian) "gaussian.trend.filtering.fit" else "logistic.trend.filtering.fit"
    if (gaussian) out$order <- 2L
    if (lambda.selection == "cv") {
        if (!is.numeric(foldid) || length(foldid) != n || anyNA(foldid) ||
            any(!is.finite(foldid)) || any(foldid != floor(foldid)) ||
            length(unique(foldid)) < 2L)
            stop("CV requires finite integer foldid with at least two folds.", call. = FALSE)
        folds <- sort(unique(foldid))
        scores <- matrix(NA_real_, length(folds), length(lambda))
        fold.status <- character(length(folds))
        factors <- numeric(length(folds))
        fold.lambda <- matrix(NA_real_, length(folds), length(lambda))
        for (k in seq_along(folds)) {
            test <- foldid == folds[k]; train <- !test
            factors[k] <- .tf.cv.penalty.factor(sum(weights[train]), sum(weights),
                                                cv.penalty.scale)
            fold.lambda[k, ] <- lambda * factors[k]
            ff <- tryCatch(.fit.univariate.trend.filtering(x[train], y[train], fold.lambda[k, ],
                weights[train], covariates[train, , drop = FALSE], offset[train],
                nuisance.ridge * factors[k], curve.ridge * factors[k],
                control = control, warn = FALSE, backend = backend.requested,
                family = family, cv.penalty.scale = cv.penalty.scale),
                error = function(e) e)
            if (inherits(ff, "error")) { fold.status[k] <- conditionMessage(ff); next }
            eta <- stats::predict(ff, x[test], covariates = covariates[test, , drop = FALSE],
                offset = offset[test], type = "link", lambda.index = seq_along(lambda))
            eta <- matrix(eta, nrow = sum(test))
            loss <- if (gaussian) (y[test] - eta)^2 else
                pmax(eta, 0) - y[test] * eta + log1p(exp(-abs(eta)))
            scores[k, ] <- colSums(loss * weights[test])
            fold.status[k] <- if (all(ff$diagnostics$converged)) "complete" else "nonconverged candidates"
        }
        cv.loss <- colSums(scores)/sum(weights)
        valid <- is.finite(cv.loss) & diagnostics$converged
        out$cv.table <- data.frame(lambda = lambda, log.loss = cv.loss, eligible = valid)
        if (gaussian) names(out$cv.table)[2L] <- "mse"
        out$cv.fold.loss.sums <- scores; out$cv.fold.status <- fold.status; out$foldid <- foldid
        out$cv.fold.lambda <- fold.lambda
        out$cv.fold.penalty.factors <- stats::setNames(factors, folds)
        out$cv.fold.ridge <- data.frame(fold = folds,
            nuisance.ridge = nuisance.ridge * factors, curve.ridge = curve.ridge * factors)
        if (any(valid)) out$selected.index <- which.min(ifelse(valid, cv.loss, Inf))
        else if (warn) warning("No lambda passed every CV fold and the full fit.", call. = FALSE)
    }
    if (warn && any(!diagnostics$converged))
        warning(sum(!diagnostics$converged), " lambda fit(s) failed convergence checks; inspect diagnostics.", call. = FALSE)
    out
}

.predict.univariate.trend.filtering <- function(
    object, newdata = NULL, covariates = NULL, offset = NULL,
    type = c("response", "link"), lambda.index = NULL,
    allow.unconverged = FALSE, family = "binomial") {
    type <- match.arg(type); training <- is.null(newdata)
    if (training) {
        newdata <- object$x
        if (is.null(covariates)) covariates <- object$covariates
        if (is.null(offset)) offset <- object$offset
    }
    if (!is.numeric(newdata) || !is.null(dim(newdata)) || any(!is.finite(newdata)))
        stop("newdata must be a finite numeric vector.", call. = FALSE)
    n <- length(newdata); nz <- length(object$design$zscale)
    if (is.null(covariates)) covariates <- matrix(numeric(), n, 0L)
    if (!is.matrix(covariates) || !is.numeric(covariates) ||
        nrow(covariates) != n || ncol(covariates) != nz || any(!is.finite(covariates)))
        stop("covariates must match training columns and newdata rows.", call. = FALSE)
    if (is.null(offset)) {
        if (object$offset.used && !training) stop("Supply offset for newdata.", call. = FALSE)
        offset <- numeric(n)
    }
    if (!is.numeric(offset) || length(offset) != n || any(!is.finite(offset)))
        stop("offset must be finite with one entry per newdata row.", call. = FALSE)
    if (is.null(lambda.index)) lambda.index <- if (!is.null(object$selected.index))
        object$selected.index else seq_along(object$lambda)
    if (!is.numeric(lambda.index) || !length(lambda.index) || anyNA(lambda.index) ||
        any(!is.finite(lambda.index)) || any(lambda.index != floor(lambda.index)) ||
        any(lambda.index < 1 | lambda.index > length(object$lambda)))
        stop("Invalid lambda.index.", call. = FALSE)
    if (!is.logical(allow.unconverged) || length(allow.unconverged) != 1L || is.na(allow.unconverged))
        stop("allow.unconverged must be TRUE or FALSE.", call. = FALSE)
    d <- object$design; u <- (newdata - d$lo)/d$span
    outside <- u < 0 | u > 1
    Z <- sweep(sweep(covariates, 2L, d$zcenter), 2L, d$zscale, "/")
    banded <- identical(object$backend, "banded")
    if (!banded) {
        H <- sweep(sweep(.logistic.tf.basis(pmin(1, pmax(0, u)), d$knots), 2L, d$center), 2L, d$scales, "/")
        A <- cbind(1, Z, H)
    }
    ans <- vapply(lambda.index, function(j) {
        f <- object$fits[[j]]
        if (!f$converged && !allow.unconverged) return(rep(NA_real_, n))
        eta <- if (banded) {
            logistic_tf_banded_predict(d$knots, f$active.knots, f$compact.coefficients,
                pmin(1, pmax(0, u))) + as.vector(Z %*% f$nuisance.coefficients) + offset
        } else as.vector(A %*% f$coefficients + offset)
        if (type == "response" && family == "binomial") stats::plogis(eta) else eta
    }, numeric(n))
    ans <- matrix(ans, nrow = n, ncol = length(lambda.index))
    colnames(ans) <- paste0("lambda=", object$lambda[lambda.index])
    if (length(lambda.index) == 1L) ans <- as.vector(ans)
    attr(ans, "outside.range") <- outside; attr(ans, "lambda") <- object$lambda[lambda.index]
    ans
}
