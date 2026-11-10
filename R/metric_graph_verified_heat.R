# Opt-in, double-precision cross-check of the observed-data heat calculation.
# This is an agreement test, not a forward-error or interval certificate.
.metric.graph.verified.policy <- function(object) {
    identical(object$parameters$numerical.policy, "verified.fallback")
}

.metric.graph.verified.limits <- function() {
    c(candidate.residual = 1e-8, reference.residual = 1e-9,
      operator.error = 1e-7, nullspace.error = 1e-7,
      prediction.error = 1e-7, df.error = 1e-7,
      gcv.scaled.error = 1, gcv.final.scaled.error = 1)
}

.metric.graph.verify.agreement <- function(diagnostics) {
    limits <- .metric.graph.verified.limits()
    x <- unlist(diagnostics[names(limits)], use.names = TRUE)
    if (length(x) != length(limits) || any(!is.finite(x)) ||
        any(x < 0) || any(x > limits)) {
        failed <- names(limits)[!is.finite(x) | x < 0 | x > limits]
        stop(paste0("Verified heat fallback failed numerical agreement: ",
                    paste(failed, collapse = ", "), ". No fit accepted."), call. = FALSE)
    }
    invisible(TRUE)
}

.metric.graph.gcv.tuning <- function(scores, status, grid,
                                     reference.scores = scores,
                                     structural.flat.reason = NULL) {
    use <- which(status == "valid" & is.finite(scores) & grid > 0)
    if (!length(use)) return(list(status = "no_valid_candidate", best.idx = NA_integer_,
                                 lower.endpoint = FALSE, upper.endpoint = FALSE,
                                 bracketed = FALSE, near.minimum.indices = integer(),
                                 selection.allowed = FALSE))
    if (length(reference.scores) != length(scores) ||
        any(!is.finite(reference.scores[use]))) {
        stop("Tuning requires finite reference scores at every valid time.", call. = FALSE)
    }
    # Pointwise agreement and resolving a score difference are separate tests.
    # Retain both per-score budgets, and cover both formulations. These are
    # operational resolution bands, not interval-arithmetic certificates.
    c <- scores[use]
    r <- reference.scores[use]
    ec <- 1e-10 + 1e-7 * abs(c)
    er <- 1e-10 + 1e-7 * abs(r)
    lower <- pmax(0, pmin(c - ec, r - er))
    upper <- pmax(c + ec, r + er)
    resolution <- data.frame(grid.index = use, eta = grid[use],
                             candidate.score = c, reference.score = r,
                             candidate.radius = ec, reference.radius = er,
                             lower = lower, upper = upper)
    # Repeated copies of a time are one alternative; preserve all original
    # indices while enveloping their bands for the time-selection decision.
    times <- unique(grid[use])
    groups <- lapply(times, function(t) which(grid[use] == t))
    lo <- vapply(groups, function(ix) min(lower[ix]), numeric(1))
    hi <- vapply(groups, function(ix) max(upper[ix]), numeric(1))
    possible <- which(lo <= min(hi))
    near <- use[which(grid[use] %in% times[possible])]
    structural <- !is.null(structural.flat.reason)
    flat <- structural || max(lo) <= min(hi)
    single <- length(times) == 1L
    allowed <- single || (!flat && length(possible) == 1L)
    label <- if (single) "single_candidate" else if (flat) {
        "unresolved_flat_gcv"
    } else if (allowed) "unique_grid_minimum" else "multiple_near_minima"
    if (structural) near <- use
    selected <- if (single) 1L else if (allowed) possible else NA_integer_
    best <- if (allowed) use[groups[[selected]][1L]] else NA_integer_
    lower.endpoint <- if (allowed) grid[best] == min(grid[use]) else NA
    upper.endpoint <- if (allowed) grid[best] == max(grid[use]) else NA
    # Strictly positive margin is required to select among distinct times.
    margin <- if (allowed && !single) min(lo[-selected]) - hi[selected] else NA_real_
    list(status = label, best.idx = best, lower.endpoint = lower.endpoint,
         upper.endpoint = upper.endpoint,
         bracketed = allowed && !single && !lower.endpoint && !upper.endpoint,
         near.minimum.indices = near, selection.allowed = allowed,
         structural.flat.reason = structural.flat.reason,
         resolution.rule = "two_form_score_envelope_v1", resolution = resolution,
         separation.margin = margin,
         candidate.minimum.index = use[which.min(c)],
         reference.minimum.index = use[which.min(r)])
}

.metric.graph.require.tuning <- function(tuning) {
    if (is.null(tuning)) return(invisible(NULL))
    if (tuning$status == "unresolved_flat_gcv") {
        stop("GCV is flat or unresolved at the declared numerical resolution: no smoothing time selected. Inspect apply.metric.graph.lowpass.path() or use a prespecified time.",
             call. = FALSE)
    }
    if (tuning$status == "multiple_near_minima") {
        stop("GCV differences are numerically unresolved between competing times: no smoothing time selected. Inspect apply.metric.graph.lowpass.path() or use a prespecified time.",
             call. = FALSE)
    }
}

.metric.graph.verified.reference <- function(basis) {
    if (!isTRUE(basis$spectral$is.complete) ||
        basis$operator$laplacian.type != "unnormalized") {
        stop("verified.fallback requires a complete spectrum and the unnormalized Laplacian.", call. = FALSE)
    }
    # Recompute from the actual operator using base R LAPACK, not the candidate
    # eigenvectors. Scale first to retain behavior under global conductance changes.
    L <- basis$operator$laplacian$matrix
    unit <- max(abs(L))
    if (unit == 0) unit <- 1
    eig <- eigen(as.matrix(L / unit), symmetric = TRUE)
    ix <- rev(seq_along(eig$values))
    ref <- list(eigenvalues = eig$values[ix] * unit,
                eigenvectors = eig$vectors[, ix, drop = FALSE])
    .check.metric.graph.lowpass.spectrum(ref, basis$operator, "full")
}

.metric.graph.verified.heat.path <- function(basis, Y, grid, filter.type) {
    if (filter.type != "heat_kernel") {
        stop("verified.fallback supports only heat_kernel.", call. = FALSE)
    }
    reference <- .metric.graph.verified.reference(basis)
    V <- basis$spectral$eigenvectors
    lambda <- basis$spectral$eigenvalues
    W <- reference$eigenvectors
    mu <- reference$eigenvalues
    n <- nrow(Y)
    p <- ncol(Y)
    J <- length(grid)
    labels <- basis$spectral$components
    nc <- basis$spectral$nullity
    Q <- W[, seq_len(nc), drop = FALSE]
    positive <- which(mu > 0)
    fitted <- array(NA_real_, c(n, J, p))
    gcv <- reference.gcv <- df <- matrix(NA_real_, J, p)
    structural.reasons <- vector("list", p)
    status <- matrix("degenerate_residual_operator", J, p)
    n.obs <- colSums(!is.na(Y))
    heat.error <- .metric.graph.heat.time.diagnostic(basis$spectral, grid, filter.type)
    keys <- apply(!is.na(Y), 2L, function(x) paste(which(x), collapse = ","))
    groups <- split(seq_len(p), factor(keys, levels = unique(keys)))
    masks <- lapply(groups, function(cols) {
        O <- which(!is.na(Y[, cols[1L]]))
        if (length(unique(labels[O])) != nc) {
            stop("Each connected component needs an observed response; an unobserved component is unidentified.", call. = FALSE)
        }
        list(cols = cols, O = O, U = setdiff(seq_len(n), O))
    })
    diagnostics <- list()
    for (j in seq_len(J)) {
        if (grid[j] == 0) {
            fitted[, j, ] <- Y
            df[j, ] <- n.obs
            status[j, ] <- "identity_undefined_gcv"
            next
        }
        a <- -expm1(-grid[j] * lambda)
        b <- -expm1(-grid[j] * mu)
        s <- max(a)
        sr <- max(b)
        if (s == 0 && sr == 0 && nc == n) {
            fitted[, j, ] <- Y
            df[j, ] <- n.obs
            next
        }
        if (!is.finite(s) || !is.finite(sr) || s <= 0 || sr <= 0) {
            stop("Verified heat fallback cannot resolve the heat complement.", call. = FALSE)
        }
        A <- tcrossprod(sweep(V, 2L, a / s, `*`), V)
        A <- (A + t(A)) / 2
        # A = F'F. Retain the square root for an energy-form Schur complement:
        # subtracting A_OO - A_OU A_UU^{-1} A_UO loses relative accuracy when
        # the observed residual is tiny compared with either subtractand.
        F <- if (any(n.obs < n)) sweep(t(V), 1L, sqrt(a / s), `*`) else NULL
        # The alternative route uses (I-H)^+ and an explicit nullspace border.
        P <- NULL
        if (any(n.obs < n)) {
            P <- tcrossprod(sweep(W[, positive, drop = FALSE], 2L,
                                  sr / b[positive], `*`), W[, positive, drop = FALSE])
        }
        for (k in seq_along(masks)) {
            mask <- masks[[k]]
            O <- mask$O
            U <- mask$U
            cols <- mask$cols
            q <- length(O)
            K <- A[O, O, drop = FALSE]
            candidate.residual <- reference.residual <- 0
            condition <- 1
            T <- matrix(0, n, q)
            Tr <- matrix(0, n, q)
            if (length(U)) {
                B <- A[U, U, drop = FALSE]
                C <- A[U, O, drop = FALSE]
                condition <- rcond(B)
                ch <- tryCatch(chol(B), error = function(e) NULL)
                if (is.null(ch)) stop("Verified heat block is singular; no ridge was added.", call. = FALSE)
                X <- backsolve(ch, forwardsolve(t(ch), C))
                candidate.residual <- max(abs(B %*% X - C)) /
                    max(1, max(abs(C)), max(rowSums(abs(B))) * max(abs(X)))
                T[U, ] <- -X
                # For the same harmonic extension X, K = R'R with
                # R = F_O - F_U X. This is the same Schur complement in exact
                # arithmetic, with no ridge, truncation or changed nullspace.
                residual.factor <- F[, O, drop = FALSE] - F[, U, drop = FALSE] %*% X
                K <- crossprod(residual.factor)
                G <- rbind(cbind(P[O, O, drop = FALSE], Q[O, , drop = FALSE]),
                           cbind(t(Q[O, , drop = FALSE]), matrix(0, nc, nc)))
                rhs <- rbind(diag(q), matrix(0, nc, q))
                Z <- tryCatch(solve(G, rhs, tol = 0), error = function(e) NULL)
                if (is.null(Z)) stop("Verified heat reference solve failed.", call. = FALSE)
                reference.residual <- max(abs(G %*% Z - rhs)) /
                    max(1, max(rowSums(abs(G))) * max(abs(Z)))
                Kr <- Z[seq_len(q), , drop = FALSE]
                Tr[U, ] <- P[U, O, drop = FALSE] %*% Kr +
                    Q[U, , drop = FALSE] %*% Z[q + seq_len(nc), , drop = FALSE]
            } else {
                Kr <- tcrossprod(sweep(W, 2L, b / sr, `*`), W)
            }
            K <- (K + t(K)) / 2
            Kr <- (Kr + t(Kr)) / 2
            T[O, ] <- diag(q) - s * K
            Tr[O, ] <- diag(q) - sr * Kr
            tr <- sum(diag(K))
            trr <- sum(diag(Kr))
            d <- q - s * tr
            dr <- q - sr * trr
            values <- Y[O, cols, drop = FALSE]
            # Remove the observed component-constant nullspace before any
            # response multiplication. Relative-to-anchor centering avoids
            # summing a large baseline. Each component needs its own mean;
            # subtracting one global mean is insufficient for disconnected graphs.
            centered <- values
            centers <- matrix(0, nc, length(cols))
            for (component in seq_len(nc)) {
                ix <- which(labels[O] == component)
                anchor <- values[ix[1L], ]
                delta <- sweep(values[ix, , drop = FALSE], 2L, anchor, `-`)
                mean.delta <- colMeans(delta)
                centered[ix, ] <- sweep(delta, 2L, mean.delta, `-`)
                centers[component, ] <- anchor + mean.delta
            }
            magnitude <- apply(abs(centered), 2L, max)
            magnitude[magnitude == 0] <- 1
            scales <- magnitude * sqrt(colMeans(sweep(centered, 2L, magnitude, `/`)^2))
            scales[scales == 0] <- 1
            for (col in seq_along(cols)) {
                reason <- if (q - nc == 1L) "one_dimensional_observed_residual_space" else
                    if (all(centered[, col] == 0)) "componentwise_constant_response" else NULL
                structural.reasons[cols[col]] <- list(reason)
            }
            standardized <- sweep(centered, 2L, scales, `/`)
            pred <- T %*% standardized
            predr <- Tr %*% standardized
            valid <- q > nc && is.finite(tr) && is.finite(trr) &&
                min(tr, trr) > 128 * .Machine$double.eps * q
            score <- scorer <- rep(NA_real_, length(cols))
            score.error <- final.score.error <- 0
            raw.score <- reference.raw.score <- rep(NA_real_, length(cols))
            if (valid) {
                score <- q * colSums((K %*% standardized / tr)^2)
                scorer <- q * colSums((Kr %*% standardized / trr)^2)
                score.error <- max(abs(score - scorer) / (1e-10 + 1e-7 * abs(scorer)))
                # Restore units from the checked quantities, not by multiplying
                # K by uncentered values again. Sequential scaling avoids forming
                # scale^2 when that intermediate would overflow.
                raw.score <- (score * scales) * scales
                reference.raw.score <- (scorer * scales) * scales
                final.score.error <- max(abs(raw.score - reference.raw.score) /
                    (1e-10 + 1e-7 * abs(reference.raw.score)))
            }
            checks <- list(candidate.residual = candidate.residual,
                           reference.residual = reference.residual,
                           operator.error = max(rowSums(abs(T - Tr))),
                           nullspace.error = max(abs(T %*% Q[O, , drop = FALSE] - Q),
                                                 abs(Tr %*% Q[O, , drop = FALSE] - Q)),
                           prediction.error = max(abs(pred - predr)),
                           df.error = abs(d - dr), gcv.scaled.error = score.error,
                           gcv.final.scaled.error = final.score.error)
            .metric.graph.verify.agreement(checks)
            if (q > nc && !valid) {
                stop("Verified heat fallback has unresolved residual degrees of freedom.", call. = FALSE)
            }
            # Restore component means only after applying the checked operator
            # to centered responses. Its nullspace-preservation check remains
            # mandatory; centering does not excuse an inaccurate operator.
            raw <- sweep(pred, 2L, scales, `*`) + centers[labels, , drop = FALSE]
            if (any(!is.finite(raw)) || (valid && any(!is.finite(raw.score)))) {
                stop("Verified heat fallback produced non-finite response results.", call. = FALSE)
            }
            fitted[, j, cols] <- raw
            df[j, cols] <- d
            if (valid) {
                gcv[j, cols] <- raw.score
                reference.gcv[j, cols] <- reference.raw.score
                status[j, cols] <- "valid"
            }
            uncertainty <- heat.error[j] / s + 32 * .Machine$double.eps * n
            old.accept <- is.finite(heat.error[j]) && heat.error[j] <= basis$spectral$verification$tolerance &&
                (!length(U) || (is.finite(condition) && condition > 0 && uncertainty / condition <= 1e-8))
            diagnostics[[length(diagnostics) + 1L]] <- as.data.frame(c(
                list(grid.index = j, eta = grid[j], mask = k, n.observed = q,
                     conservative.guard.passed = old.accept, fallback.used = !old.accept,
                     heat.error = heat.error[j], reciprocal.condition = condition,
                     condition.error = if (length(U)) uncertainty / condition else 0), checks))
        }
    }
    tuning <- lapply(seq_len(p), function(k) .metric.graph.gcv.tuning(
        gcv[, k], status[, k], grid, reference.gcv[, k], structural.reasons[[k]]))
    list(fitted = fitted, gcv = gcv, df = df, status = status, n.observed = n.obs,
         tuning = tuning,
         numerical = list(policy = "verified.fallback", limits = .metric.graph.verified.limits(),
                          reference.gcv.scores = reference.gcv,
                          reference.verification = reference$verification,
                          diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(),
                          response.mask = match(keys, unique(keys)),
                          interpretation = "Double-precision agreement; not a rigorous forward-error certificate."))
}
