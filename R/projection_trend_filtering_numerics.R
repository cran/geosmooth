# Numerical geometry for constrained projection-pursuit proposals.

.ptf.curve.slopes <- function(fit, x = fit$x) {
    d <- fit$design
    if (ncol(fit$covariates)) stop("Curve slopes require a fit without covariates.")
    b <- fit$fits[[1L]]$coefficients[-1L]/d$scales
    v <- (x-d$lo)/d$span
    g <- b[1L]+2*b[2L]*v
    if (length(b) > 2L) for (j in 3:length(b))
        g <- g+b[j]*(v-(d$knots[j-1L]+d$knots[j])/2)*(v > d$knots[j])
    g[v < 0 | v > 1] <- 0
    g/d$span
}

.ptf.curve.diagnostics <- function(component) {
    fit <- component$fit; d <- fit$design
    b <- fit$fits[[1L]]$coefficients[-1L]/d$scales
    j <- if (length(b) > 2L) 3:length(b) else integer()
    threshold <- 1e-10*max(1, abs(b))
    jumps <- b[j]*(d$knots[j]-d$knots[j-1L])/(2*d$span)
    ends <- .ptf.curve.slopes(fit, range(fit$x))
    list(active.knots = sum(abs(b[j]) > threshold), coefficient.threshold = threshold,
         interior.jumps = jumps, boundary.jumps = c(ends[1L], -ends[2L]),
         largest.jump = max(abs(c(jumps, ends))), projection.span = d$span)
}

.ptf.projector.basis <- function(E, dimension, tolerance = 1e-12) {
    n <- nrow(E); V <- matrix(0, n, 0L)
    if (!dimension) return(V)
    for (j in seq_len(n)) {
        v <- E[, j]
        for (pass in 1:2) if (ncol(V)) v <- v-as.vector(V %*% crossprod(V, v))
        size <- sqrt(sum(v^2))
        if (size <= tolerance) next
        v <- v/size
        if (v[which.max(abs(v))] < 0) v <- -v
        V <- cbind(V, v)
        if (ncol(V) == dimension) return(unname(V))
    }
    stop("Could not recover the projector's declared dimension.", call. = FALSE)
}

.ptf.gn.equations <- function(Z, y, eta, g, w, family, curvature.floor = 1e-6) {
    J <- Z*g
    v <- if (family == "gaussian") rep(1, length(y)) else
        stats::plogis(eta)*stats::plogis(-eta)
    vt <- if (family == "gaussian") v else pmax(v, curvature.floor)
    residual <- y-if (family == "gaussian") eta else stats::plogis(eta)
    A <- crossprod(J, J*(w*vt)); b <- as.vector(crossprod(J, w*residual))
    if (any(!is.finite(c(A, b)))) stop("Nonfinite direction surrogate.", call. = FALSE)
    list(A = A, b = b,
         unfloored = crossprod(J, J*(w*v)), J = J, variance = v,
         stabilized.variance = vt, floored = sum(vt > v),
         working.residual = residual/vt, working.weights = w*vt)
}

.ptf.gn.increment <- function(A, b, u, mu = NULL) {
    p <- length(u)
    if (p == 1L) return(list(delta = numeric(p), T = matrix(0, p, 0),
                            mu = NA_real_, status = "rank_one"))
    T <- .ptf.projector.basis(diag(p)-tcrossprod(u), p-1L)
    curvature <- crossprod(T, A %*% T); score <- as.vector(crossprod(T, b))
    if (is.null(mu)) mu <- max(1e-8, 1e-3*sum(diag(curvature))/(p-1L))
    delta <- as.vector(T %*% solve(curvature+diag(mu, p-1L), score))
    threshold <- 1e-12*max(1, sqrt(sum(b^2)))
    list(delta = delta, T = T, mu = mu, tangent.score = sqrt(sum(score^2)),
         score.threshold = threshold,
         status = if (sqrt(sum(score^2)) <= threshold) "zero_tangent_score" else "proposal")
}

.ptf.retract <- function(u, delta, fraction = 1) .ptf.canonical(u+fraction*delta)

.ptf.pca <- function(X, weights = NULL, theta = .9, rank = NULL) {
    X <- .ptf.matrix(X)
    if (!ncol(X) || nrow(X) < 2L) stop("PCA needs predictors and at least two rows.")
    if (is.null(weights)) weights <- rep(1, nrow(X))
    if (length(weights) != nrow(X) || any(!is.finite(weights)) ||
        any(weights <= 0) || !is.finite(sum(weights))) stop("Invalid PCA weights.")
    if (is.null(rank)) {
        .ptf.scalar(theta, "theta", positive = TRUE)
        if (theta > 1) stop("theta must not exceed one.")
    } else .ptf.scalar(rank, "rank", integer = TRUE, positive = TRUE)
    w <- weights/sum(weights); center <- colSums(X*w)
    constant <- vapply(seq_len(ncol(X)), function(j) all(X[, j] == X[1L, j]), logical(1))
    center[constant] <- X[1L, constant]
    Z <- sweep(X, 2, center); C <- crossprod(Z, Z*w)
    eig <- eigen(C, symmetric = TRUE); values <- pmax(0, eig$values)
    blocks <- list(); V <- matrix(0, ncol(X), 0L); a <- 1L
    while (a <= length(values) && values[a] > 0) {
        end <- max(which(seq_along(values) >= a & values[a]-values <= 1e-8*values[a]))
        if (values[end] <= 1e-10*values[1L]) break
        ix <- a:end; E <- tcrossprod(eig$vectors[, ix, drop = FALSE])
        V <- cbind(V, .ptf.projector.basis(E, length(ix)))
        blocks[[length(blocks)+1L]] <- list(indices = ix, projector = E)
        a <- end+1L
    }
    r <- ncol(V); ends <- vapply(blocks, function(b) max(b$indices), integer(1))
    available <- is.null(rank) || rank <= r
    q <- if (!r || !available) 0L else if (is.null(rank))
        if (theta == 1) r else
        ends[which(cumsum(values)[ends] >= theta*sum(values[seq_len(r)]))[1L]] else
        ends[which(ends >= rank)[1L]]
    P <- V[, seq_len(q), drop = FALSE]
    list(scores = Z %*% P, basis = P, center = center, covariance = C,
         score.covariance = crossprod(P, C %*% P), values = values,
         blocks = blocks, numerical.rank = r, requested.rank = rank,
         rank = q, available = available, theta = theta,
         rank.expanded = !is.null(rank) && available && q > rank,
         tolerances = c(block = 1e-8, rank = 1e-10, basis = 1e-12))
}

.ptf.seed.data <- function(problem, state, component = NULL) {
    eta <- problem$offset+state$intercept+rowSums(state$terms)
    if (problem$family == "gaussian") {
        r <- problem$y-eta; w <- problem$w
    } else {
        p <- stats::plogis(eta); v <- pmax(p*stats::plogis(-eta), 1e-6)
        r <- (problem$y-p)/v; w <- problem$w*v
    }
    if (!is.null(component)) r <- r+state$terms[, component]
    zbar <- colSums(problem$Z*w)/sum(w); rbar <- sum(w*r)/sum(w)
    list(Z = sweep(problem$Z, 2, zbar), r = r-rbar, w = w,
         predictor.center = zbar, response.center = rbar)
}

.ptf.sir <- function(Z, r, w, classes = NULL) {
    n <- nrow(Z)
    w <- w/sum(w); S <- 1
    # Treat equality of cumulative masses consistently after weight rescaling.
    # Response values and their exact tie groups are never perturbed.
    mass.tolerance <- 64*.Machine$double.eps
    C <- crossprod(Z, Z*w)/S
    if (is.null(classes)) {
        H <- min(10L, floor(n/15))
        if (H < 2L) return(list(direction = NULL, status = "sir_unavailable", slices = integer(n)))
        values <- sort(unique(r))
        mass <- vapply(values, function(v) sum(w[r == v]), numeric(1))
        bounds <- unique(vapply(seq_len(H-1L), function(h)
            values[which(cumsum(mass)+mass.tolerance >= h*S/H)[1L]], numeric(1)))
        slices <- 1L+rowSums(outer(r, bounds, ">"))
    } else slices <- match(classes, sort(unique(classes)))
    ids <- sort(unique(slices))
    if (length(ids) < 2L) return(list(direction = NULL, status = "sir_unavailable", slices = slices))
    M <- matrix(0, ncol(Z), ncol(Z))
    for (h in ids) {
        take <- slices == h; mass <- sum(w[take])
        m <- colSums(Z[take, , drop = FALSE]*w[take])/mass
        M <- M+(mass/S)*tcrossprod(m)
    }
    ec <- eigen(C, symmetric = TRUE)
    keep <- which(ec$values > 1e-10*max(ec$values))
    if (!length(keep)) return(list(direction = NULL, status = "sir_unavailable", slices = slices))
    W <- sweep(ec$vectors[, keep, drop = FALSE], 2, sqrt(ec$values[keep]), "/")
    em <- eigen(crossprod(W, M %*% W), symmetric = TRUE)
    gamma <- em$values[1L]; direction <- NULL
    if (gamma >= 1e-8) {
        tied <- which(gamma-em$values <= 1e-8*gamma)
        U <- qr.Q(qr(W %*% em$vectors[, tied, drop = FALSE]))
        direction <- .ptf.canonical(.ptf.projector.basis(tcrossprod(U), length(tied))[, 1L])
    }
    list(direction = direction, status = if (is.null(direction)) "sir_unavailable" else "available",
         slices = slices, covariance = C, between = M, eigenvalue = gamma,
         slice.mass.tolerance = mass.tolerance)
}

.ptf.supervised.seeds <- function(problem, state, component = NULL, only.pls = FALSE) {
    s <- .ptf.seed.data(problem, state, component); p <- ncol(s$Z)
    c <- as.vector(crossprod(s$Z, s$w*s$r))
    zero <- sqrt(sum(c^2)) < 1e-10*max(1, sqrt(sum(s$r^2)))
    ridge <- pls <- NULL
    if (!zero) {
        pls <- .ptf.canonical(c)
        if (!only.pls) {
            C <- crossprod(s$Z, s$Z*s$w)
            ridge <- .ptf.canonical(solve(C+diag(1e-3*sum(diag(C))/p, p), c))
        }
    }
    first.binary <- problem$family == "binomial" && all(problem$weights == 1) && !ncol(state$terms) &&
        is.null(component) && all(problem$y %in% c(0, 1))
    sir <- if (only.pls) list(direction = NULL, status = "not_requested") else
        .ptf.sir(s$Z, s$r, s$w, if (first.binary) problem$y else NULL)
    list(directions = list(ridge = ridge, pls = pls, sir = sir$direction),
         zero.seed = zero, sir = sir, centered = s, first.bernoulli = first.binary)
}

.ptf.key.seed <- function(seed, key) {
    as.integer(abs(as.double(digest::digest2int(paste(c(seed, key), collapse = "\034")))) %%
                   .Machine$integer.max)
}

.ptf.gn.starts <- function(problem, state, term, seed, method = "gn", key = character()) {
    p <- ncol(problem$Z)
    seeds <- .ptf.supervised.seeds(problem, state, only.pls = method == "pls")
    stream <- .ptf.key.seed(seed, c(key, term, "random-start"))
    count <- if (method == "pls") as.integer(is.null(seeds$directions$pls)) else 2L*p
    draws <- if (count) .ptf.with.seed(stream,
        matrix(stats::rnorm(p*count), ncol = p, byrow = TRUE)) else matrix(0, 0, p)
    candidates <- list(); origins <- character()
    add <- function(u, origin) {
        if (is.null(u)) return()
        candidates[[length(candidates)+1L]] <<- .ptf.canonical(u)
        origins <<- c(origins, origin)
    }
    if (method == "pls") {
        if (is.null(seeds$directions$pls)) add(draws[1L, ], "random_fallback") else
            add(seeds$directions$pls, "pls")
    } else {
        for (j in seq_len(p)) add(diag(p)[j, ], paste0("axis_", j))
        for (nm in names(seeds$directions)) add(seeds$directions[[nm]], nm)
        for (j in seq_len(nrow(draws))) add(draws[j, ], paste0("random_", j))
    }
    U <- matrix(0, 0, p); source <- list()
    for (j in seq_along(candidates)) {
        u <- candidates[[j]]
        duplicate <- which(vapply(seq_len(nrow(U)), function(k) max(abs(U[k, ]-u)) <= 1e-12, logical(1)))
        if (length(duplicate)) source[[duplicate[1L]]] <- c(source[[duplicate[1L]]], origins[j]) else {
            U <- rbind(U, u); source[[length(source)+1L]] <- origins[j]
        }
    }
    list(directions = U, origins = source, seed = stream, supervised = seeds,
         random.draws = count, rng.kind = RNGkind())
}
