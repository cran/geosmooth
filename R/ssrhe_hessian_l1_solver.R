# Numerical ADMM engine shared by Hessian L1 and transported trend filtering.
# The statistical objective is unchanged by automatic rho or proximal damping.

.ssrhe.soft.threshold <- function(x, lambda) {
    sign(x) * pmax(abs(x) - lambda, 0)
}

.ssrhe.hessian.l1.admm <- function(y, weights, D, lambda, solver.args) {
    D <- methods::as(D, "dgCMatrix")
    n <- length(y)
    m <- nrow(D)
    weights <- as.numeric(weights)
    y <- as.numeric(y)
    norm2 <- function(x) sqrt(sum(x * x))
    operator.scale <- if (m) mean(Matrix::rowSums(D^2)) else 0
    if (!is.finite(operator.scale)) {
        stop("The ADMM penalty operator is too large for finite arithmetic.",
             call. = FALSE)
    }
    rho <- solver.args$admm.rho %||% (1 / max(1, operator.scale))
    rho.initial <- rho
    rho.updates <- 0L
    adaptive <- isTRUE(solver.args$admm.adaptive.rho %||% TRUE)
    AtA <- Matrix::crossprod(D)
    factorize <- function(rho) {
        system <- AtA * rho + Matrix::Diagonal(n, x = weights)
        proximal <- max(1e-10, sqrt(.Machine$double.eps) *
                           max(1, mean(Matrix::diag(system))))
        list(
            factor = Matrix::Cholesky(
                system + Matrix::Diagonal(n, x = rep(proximal, n)),
                LDL = FALSE, Imult = 0),
            proximal = proximal
        )
    }
    factored <- factorize(rho)
    wy <- weights * y
    beta <- rep(0, n)
    z <- u <- rep(0, m)
    converged <- FALSE
    primal <- dual <- stationarity <- Inf
    eps.pri <- eps.dual <- NA_real_
    # A finite adaptation window leaves an ordinary fixed-rho iteration
    # afterwards. Relative bounds preserve very small/large operator units.
    rho.lower <- max(.Machine$double.xmin, rho.initial / 1e6)
    rho.upper <- min(.Machine$double.xmax, rho.initial * 1e6)
    for (iter in seq_len(solver.args$admm.maxiter)) {
        previous <- beta
        z.old <- z
        # Center damping at the previous iterate: unlike a fixed ridge
        # penalty, this term vanishes from the limiting optimality equation.
        rhs <- wy + rho * as.vector(Matrix::crossprod(D, z - u)) +
            factored$proximal * previous
        beta <- as.vector(Matrix::solve(factored$factor, rhs))
        Dbeta <- as.vector(D %*% beta)
        z <- .ssrhe.soft.threshold(Dbeta + u, lambda / rho)
        u <- u + Dbeta - z
        primal <- norm2(Dbeta - z)
        dual <- rho * norm2(as.vector(Matrix::crossprod(D, z - z.old)))
        dual.gradient <- as.vector(Matrix::crossprod(D, rho * u))
        stationarity <- norm2(weights * (beta - y) + dual.gradient)
        eps.pri <- sqrt(m) * solver.args$admm.abstol +
            solver.args$admm.reltol * max(norm2(Dbeta), norm2(z))
        eps.dual <- sqrt(n) * solver.args$admm.abstol +
            solver.args$admm.reltol * norm2(dual.gradient)
        if (!all(is.finite(c(beta, primal, dual, stationarity,
                             eps.pri, eps.dual)))) break
        if (primal <= eps.pri && dual <= eps.dual &&
            stationarity <= eps.dual) {
            converged <- TRUE
            break
        }
        if (adaptive && iter %% 25L == 0L && iter <= 1000L &&
            iter < solver.args$admm.maxiter) {
            relative.primal <- primal / max(eps.pri, .Machine$double.eps)
            relative.dual <- dual / max(eps.dual, .Machine$double.eps)
            next.rho <- if (relative.primal > 10 * relative.dual) {
                min(rho.upper, 2 * rho)
            } else if (relative.dual > 10 * relative.primal) {
                max(rho.lower, rho / 2)
            } else rho
            if (next.rho != rho) {
                u <- u * (rho / next.rho)
                rho <- next.rho
                factored <- factorize(rho)
                rho.updates <- rho.updates + 1L
            }
        }
    }
    status <- if (converged) "converged" else if (
        all(is.finite(c(beta, primal, dual, stationarity)))) {
        "iteration_limit"
    } else "nonfinite_iterate"
    attr(beta, "admm") <- list(
        iterations = iter,
        converged = converged,
        status = status,
        primal.residual = primal,
        dual.residual = dual,
        stationarity.residual = stationarity,
        primal.tolerance = eps.pri,
        dual.tolerance = eps.dual,
        rho = rho,
        rho.initial = rho.initial,
        rho.updates = rho.updates,
        adaptive.rho = adaptive,
        proximal.diagonal = factored$proximal
    )
    beta
}

.ssrhe.hessian.l1.admm.grid <- function(y, weights, D, lambda.grid,
                                        solver.args) {
    admm.info <- vector("list", length(lambda.grid))
    beta <- vapply(seq_along(lambda.grid), function(ii) {
        fit <- tryCatch(.ssrhe.hessian.l1.admm(
            y = y, weights = weights, D = D, lambda = lambda.grid[ii],
            solver.args = solver.args
        ), error = function(e) {
            structure(rep(NA_real_, length(y)), admm = list(
                iterations = 0L, converged = FALSE, status = "solver_error",
                message = conditionMessage(e)))
        })
        admm.info[[ii]] <<- attr(fit, "admm")
        as.vector(fit)
    }, numeric(length(y)))
    attr(beta, "admm") <- admm.info
    beta
}

.ssrhe.hessian.l1.candidate.status <- function(beta) {
    finite <- colSums(is.finite(beta)) == nrow(beta)
    out <- data.frame(
        status = ifelse(finite, "path_available", "nonfinite_fit"),
        acceptable = finite,
        converged = rep(NA, ncol(beta)),
        iterations = rep(NA_integer_, ncol(beta)),
        message = rep("", ncol(beta)),
        stringsAsFactors = FALSE
    )
    admm <- attr(beta, "admm")
    if (!is.null(admm)) {
        for (ii in seq_along(admm)) {
            out$converged[ii] <- isTRUE(admm[[ii]]$converged)
            out$status[ii] <- admm[[ii]]$status
            out$iterations[ii] <- admm[[ii]]$iterations
            out$message[ii] <- admm[[ii]]$message %||% ""
        }
        out$acceptable <- finite & out$converged
    }
    out
}
