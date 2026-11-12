# Simultaneous fixed-direction fitting, with the same joint certificate as backfit.
.ptf.joint.failure <- function(message, kind = "numerical", details = NULL) {
    structure(list(message = message, call = NULL, kind = kind, details = details),
        class = c("projection_joint_failure", "error", "condition"))
}

.ptf.joint.geometry <- function(problem, state) {
    n <- length(problem$y)
    if (length(state$components) != 2L || !identical(dim(state$terms), c(n, 2L)) ||
        !is.matrix(state$directions) || nrow(state$directions) != 2L ||
        ncol(state$directions) != ncol(problem$Z))
        stop("Invalid two-component joint state.", call. = FALSE)
    x <- lapply(1:2, function(j) .ptf.projection(problem, state$directions[j, ]))
    m <- vapply(x, function(v) length(unique(v)), integer(1))
    p <- 1+sum(m-1)
    workspace <- 8*(12*as.double(n)*p+16*p^2+12*n+12*p)
    if (!is.finite(workspace) || workspace > problem$control$joint.max.workspace.bytes)
        stop(.ptf.joint.failure("joint_workspace_refused", "workspace",
            list(workspace.estimate.bytes = workspace, columns = p)))
    blocks <- lapply(1:2, function(j) {
        lo <- min(x[[j]]); span <- diff(range(x[[j]]))
        knots <- (sort(unique(x[[j]]))-lo)/span
        if (any(diff(knots) <= 0))
            stop(.ptf.joint.failure("Distinct projections collapsed during normalization."))
        H <- .logistic.tf.basis((x[[j]]-lo)/span, knots)
        center <- as.vector(crossprod(problem$w, H))
        H <- sweep(H, 2L, center)
        scales <- sqrt(colSums(H^2*problem$w))
        if (any(!is.finite(scales)) || any(scales <= 0))
            stop(.ptf.joint.failure("Degenerate joint falling-factorial basis."))
        H <- sweep(H, 2L, scales, "/")
        previous <- state$components[[j]]$fit$design
        difference <- if (!is.null(previous) && identical(previous$knots, knots))
            list(center = max(abs(center-previous$center)),
                 scales = max(abs(scales-previous$scales))) else NULL
        list(x = x[[j]], H = H, penalty = c(0, 0, 1/scales[-(1:2)]),
            geometry.difference = difference,
            design = list(lo = lo, span = span, knots = knots, center = center,
                          scales = scales, zcenter = numeric(), zscale = numeric()))
    })
    A <- cbind(1, blocks[[1L]]$H, blocks[[2L]]$H)
    indices <- list(1L+seq_len(m[1L]-1L), m[1L]+seq_len(m[2L]-1L))
    R <- matrix(0, if (problem$rho > 0) 2L*n else 0L, p)
    if (problem$rho > 0) for (j in 1:2)
        R[(j-1L)*n+seq_len(n), indices[[j]]] <-
            blocks[[j]]$H*sqrt(sum(problem$weights)*problem$rho*problem$w)
    penalty <- c(0, blocks[[1L]]$penalty, blocks[[2L]]$penalty)
    if (nrow(A) != n || ncol(A) != p || ncol(R) != p || length(penalty) != p ||
        any(A[, 1L] != 1) || any(!is.finite(c(A, R, penalty))) || any(penalty < 0) ||
        length(problem$weights) != n || any(!is.finite(problem$weights)) ||
        any(problem$weights <= 0) || !is.finite(sum(problem$weights)) ||
        length(problem$offset) != n || any(!is.finite(c(problem$y, problem$offset))))
        stop(.ptf.joint.failure("Invalid dense joint geometry or training inputs."))
    list(A = A, R = R, penalty = penalty, blocks = blocks, indices = indices,
         workspace.estimate.bytes = workspace)
}

.ptf.joint.reconstruct <- function(problem, state, geometry, native) {
    p <- ncol(geometry$A); n <- nrow(geometry$A); mass <- sum(problem$weights)
    if (!isTRUE(native$converged) || length(native$coefficients) != p ||
        length(native$dual) != p || length(native$linear.predictors) != n ||
        any(!is.finite(c(native$coefficients, native$dual, native$linear.predictors,
                        native$objective))))
        stop(.ptf.joint.failure(paste0("joint_native_", native$status), details = native))
    if (is.infinite(problem$lambda) &&
        any(c(native$coefficients[geometry$penalty > 0], native$dual[geometry$penalty > 0]) != 0))
        stop(.ptf.joint.failure("Nonzero penalized coefficients at the polynomial limit."))
    raw <- vapply(1:2, function(j)
        as.vector(geometry$blocks[[j]]$H %*% native$coefficients[geometry$indices[[j]]]), numeric(n))
    centers <- colSums(raw*problem$w)
    state$terms <- sweep(raw, 2L, centers)
    state$intercept <- native$coefficients[1L]+sum(centers)
    eta <- problem$offset+state$intercept+rowSums(state$terms)
    control <- utils::modifyList(list(tolerance = 1e-9, gap.tolerance = 1e-10,
        max.iterations = 1000L, max.workspace.bytes = 1024^3), problem$control$inner)
    for (j in 1:2) {
        block <- geometry$blocks[[j]]; index <- geometry$indices[[j]]
        f <- native
        f$coefficients <- c(0, native$coefficients[index])
        f$dual <- c(0, native$dual[index])
        f$separation.direction <- numeric(length(index)+1L)
        f$certificate.scope <- "joint"
        f$linear.predictors <- eta
        f$fitted.values <- if (problem$family == "gaussian") eta else stats::plogis(eta)
        diagnostics <- data.frame(lambda = problem$lambda*mass, converged = TRUE,
            status = native$status, objective = native$objective,
            stationarity = native$stationarity, complementarity = native$complementarity,
            dual.feasibility = native$dual.feasibility, iterations = native$iterations,
            backend = "dense", certificate.scope = "joint")
        fit <- list(call = NULL, backend = "dense", x = block$x, y = problem$y,
            weights = problem$weights,
            offset = eta-raw[, j], covariates = matrix(numeric(), n, 0L),
            offset.used = TRUE, lambda = problem$lambda*mass, fits = list(f),
            diagnostics = diagnostics, selected.index = NULL, cv.table = NULL,
            cv.penalty.scale = "training_mass", control = control, nuisance.ridge = 0,
            curve.ridge = problem$rho*mass,
            workspace.estimate.bytes = geometry$workspace.estimate.bytes, design = block$design)
        class(fit) <- if (problem$family == "gaussian") "gaussian.trend.filtering.fit" else
            "logistic.trend.filtering.fit"
        if (problem$family == "gaussian") fit$order <- 2L
        state$components[[j]] <- list(fit = fit, center = centers[j])
    }
    check <- .ptf.diagnose(problem, state)
    representation.error <- max(abs(eta-native$linear.predictors))
    objective.error <- abs(check$objective-native$objective/mass)
    if (!isTRUE(check$converged) ||
        representation.error > 1e-10*max(1, abs(native$linear.predictors)) ||
        objective.error > 1e-10*max(1, abs(check$objective)))
        stop(.ptf.joint.failure("joint_reconstruction_failed", details = list(
            native = native, checks = check, representation.error = representation.error,
            objective.error = objective.error)))
    state$checks <- check; state$converged <- TRUE; state$status <- "converged"
    state$sweeps <- 0L; state$solver.calls <- 1L
    state$joint.diagnostics <- list(native = native,
        representation.error = representation.error, objective.error = objective.error,
        columns = p, observations = n, workspace.estimate.bytes = geometry$workspace.estimate.bytes,
        geometry.difference = lapply(geometry$blocks, `[[`, "geometry.difference"))
    state
}

.ptf.backfit.dense <- function(problem, state) {
    incoming <- .ptf.diagnose(problem, state)
    if (isTRUE(incoming$converged)) {
        out <- .ptf.backfit.core(problem, state)
        out$joint.work <- list(policy = "dense", direct.entries = 0L,
                              fallback = FALSE, direct.status = "incoming_certified")
        return(out)
    }
    entered <- 0L
    trial <- tryCatch(.ptf.observed(problem$observer, "joint_dense", list(), {
        g <- .ptf.observed(problem$observer, "joint_geometry", list(),
            .ptf.joint.geometry(problem, state))
        ctl <- utils::modifyList(list(tolerance = 1e-9, gap.tolerance = 1e-10,
            max.iterations = 1000L), problem$control$inner)
        native <- .ptf.native.call(problem$observer, problem$before.fit,
            list(backend = "dense", call.type = "joint", components = 2L,
                 observations = nrow(g$A), columns = ncol(g$A)), {
                entered <- 1L
                solve <- if (problem$family == "gaussian") gaussian_tf_native else logistic_tf_native
                solve(g$A, g$R, as.numeric(problem$y), as.numeric(problem$weights),
                    as.numeric(problem$offset), g$penalty, problem$lambda*sum(problem$weights),
                    ctl$tolerance, ctl$gap.tolerance, as.integer(ctl$max.iterations))
            })[[1L]]
        out <- .ptf.observed(problem$observer, "joint_reconstruction", list(),
            .ptf.joint.reconstruct(problem, state, g, native))
        if (out$checks$objective > incoming$objective+
            problem$control$objective.tolerance*max(1, abs(incoming$objective)))
            stop(.ptf.joint.failure("objective_increase", details = out$joint.diagnostics))
        out$objective.history <- c(incoming$objective, out$checks$objective)
        out
    }), error = function(e) {
        .ptf.numerical.error(e)
        if (!inherits(e, "projection_joint_failure")) stop(e)
        e
    })
    if (!inherits(trial, "projection_joint_failure")) {
        trial$joint.work <- list(policy = "dense", direct.entries = entered,
                                fallback = FALSE, direct.status = "converged")
        return(trial)
    }
    fallback <- .ptf.observed(problem$observer, "joint_fallback",
        list(reason = conditionMessage(trial), kind = trial$kind, direct.entries = entered,
             details = trial$details), .ptf.backfit.core(problem, state))
    fallback$solver.calls <- fallback$solver.calls+entered
    fallback$joint.work <- list(policy = "dense", direct.entries = entered,
        fallback = TRUE, direct.status = conditionMessage(trial), failure.kind = trial$kind,
        direct.details = trial$details)
    fallback
}
