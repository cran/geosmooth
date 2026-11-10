#' Unfold a Local Neighborhood of the L-infinity Simplex
#'
#' Rotate one adjacent face into the query's face before constructing a local
#' chart. This function certifies a conservative one- or two-face chart; it does
#' not force a chart through a corner or through three or more faces.
#'
#' @param X Nonnegative numeric matrix, with support points in rows. Each row
#'   must have a positive maximum; rows are max-normalized internally.
#' @param query Nonnegative query vector with length ncol(X) and positive maximum.
#' @param face.tol Nonnegative tolerance for identifying tied maxima.
#' @param distance.tol Positive absolute tolerance for agreement between unfolded
#'   Euclidean chords and intrinsic distances.
#' @details On reference face i, an adjacent-face j point is mapped to \code{v[i]=1},
#'   \code{v[j]=2-u[i]}, retaining the other coordinates. The query is included when
#'   identifying faces. Tied maxima or more than two faces yield an uncertified
#'   chart. For two faces, the intrinsic support radius must be less than the
#'   query's distance to every other face, and all support/query pair distances
#'   must agree with the unfolded chords. Distances use
#'   \code{dgraphs::linf.simplex.distances()}, including intermediate faces.
#'   Certification costs quadratic work and memory in the local support size.
#'   A failed certificate does not prove that no useful chart exists.
#' @return A list with certified, status, coordinates, query, faces, radius,
#'   clearance and distance.error. Coordinates and query are NULL when the
#'   chart is uncertified. Coordinates retain the constant reference column.
#' @examples
#' X <- rbind(c(1, .9, .2), c(.9, 1, .2), c(1, .8, .3))
#' unfold.simplex.chart(X, c(1, .95, .2))
#' @export
unfold.simplex.chart <- function(X, query, face.tol = 1e-12,
                                 distance.tol = 1e-8) {
    X <- .simplex.chart.normalize(X)
    query <- .simplex.chart.query(query, ncol(X))
    face.tol <- .local.chart.validate.nonnegative.scalar(face.tol, "face.tol")
    distance.tol <- .local.chart.validate.positive.scalar(distance.tol, "distance.tol")
    U <- rbind(query, X)
    active <- abs(U - 1) <= face.tol
    faces <- max.col(U, ties.method = "first")
    ans <- list(certified = FALSE, status = "tied.maxima", coordinates = NULL,
                query = NULL, faces = unique(faces), radius = NA_real_,
                clearance = NA_real_, distance.error = NA_real_)
    if (any(rowSums(active) != 1L)) return(ans)
    if (length(unique(faces)) > 2L) {
        ans$status <- "multiple.faces"
        return(ans)
    }
    D <- dgraphs::linf.simplex.distances(U)
    ans$radius <- max(D[1L, ])
    other <- setdiff(seq_len(ncol(U)), unique(faces))
    ans$clearance <- if (length(other)) min(1 - query[other]) else Inf
    V <- U
    adjacent <- setdiff(unique(faces), faces[1L])
    if (length(adjacent)) {
        hit <- faces == adjacent
        V[hit, adjacent] <- 2 - U[hit, faces[1L]]
        V[hit, faces[1L]] <- 1
        if (ans$radius >= ans$clearance) {
            ans$status <- "corner.clearance"
            return(ans)
        }
    }
    ans$distance.error <- max(abs(as.matrix(stats::dist(V)) - D))
    if (ans$distance.error > distance.tol) {
        ans$status <- "distance.mismatch"
        return(ans)
    }
    ans$certified <- TRUE
    ans$status <- if (length(adjacent)) "two.faces" else "one.face"
    ans$coordinates <- V[-1L, , drop = FALSE]
    ans$query <- V[1L, ]
    ans
}

.simplex.chart.normalize <- function(X) {
    X <- .state.density.validate.X(X)
    if (any(X < 0) || any(apply(X, 1L, max) <= 0))
        stop("Simplex coordinates must be nonnegative with positive row maxima.", call. = FALSE)
    X / apply(X, 1L, max)
}

.simplex.chart.query <- function(query, p) {
    if (!is.numeric(query) || length(query) != p || any(!is.finite(query)) ||
        any(query < 0) || max(query) <= 0)
        stop("query must be a finite nonnegative vector with a positive maximum.", call. = FALSE)
    query / max(query)
}

# Training-only, weighted centering. PLS1 uses R = W(P'W)^-1 to project queries.
.local.binary.rotation <- function(X, y, weights, method, dimension) {
    w <- weights / sum(weights)
    center <- colSums(X * w)
    Z <- sweep(X, 2L, center)
    if (method == "local.pca") {
        sv <- svd(Z * sqrt(w), nu = 0L, nv = min(dimension, ncol(X), nrow(X)))
        d <- min(dimension, sum(sv$d^2 > max(1e-14, sv$d[1L]^2 * 1e-10)))
        rotation <- sv$v[, seq_len(d), drop = FALSE]
    } else {
        E <- Z
        f <- y - sum(w * y)
        W <- P <- matrix(numeric(), nrow = ncol(X), ncol = 0L)
        for (j in seq_len(min(dimension, ncol(X), nrow(X)))) {
            a <- drop(crossprod(E, w * f))
            norm <- sqrt(sum(a^2))
            if (norm < 1e-12) break
            a <- a / norm
            t <- drop(E %*% a)
            tt <- sum(w * t^2)
            if (tt < 1e-14) break
            p <- drop(crossprod(E, w * t)) / tt
            q <- sum(w * f * t) / tt
            W <- cbind(W, a)
            P <- cbind(P, p)
            E <- E - tcrossprod(t, p)
            f <- f - q * t
        }
        rotation <- if (ncol(W)) W %*% solve(crossprod(P, W)) else W
    }
    list(center = center, rotation = rotation)
}

# Ridge slopes and symmetric intercept pseudo-counts at the weighted center.
.local.binary.solve <- function(A, y, w, query, ridge, mass, max.iter, tol) {
    b <- c(log((sum(w * y) + mass) / (sum(w * (1 - y)) + mass)), rep(0, ncol(A) - 1L))
    penalty <- c(0, rep(ridge, ncol(A) - 1L))
    softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
    objective <- function(b) {
        eta <- drop(A %*% b)
        sum(w * (softplus(eta) - y * eta)) + sum(penalty * b^2) / 2 +
            2 * mass * softplus(b[1L]) - mass * b[1L]
    }
    gradient <- function(b) {
        g <- drop(crossprod(A, w * (stats::plogis(A %*% b) - y))) + penalty * b
        g[1L] <- g[1L] + mass * (2 * stats::plogis(b[1L]) - 1)
        g
    }
    status <- "iteration.limit"
    for (iter in seq_len(max.iter)) {
        g <- gradient(b)
        if (max(abs(g)) < tol) { status <- "converged"; break }
        pr <- drop(stats::plogis(A %*% b))
        pi <- stats::plogis(b[1L])
        H <- crossprod(A, A * (w * pr * (1 - pr))) + diag(penalty, ncol(A))
        H[1L, 1L] <- H[1L, 1L] + 2 * mass * pi * (1 - pi)
        step <- tryCatch(solve(H, g), error = function(e) NULL)
        if (is.null(step) || any(!is.finite(step))) { status <- "singular.hessian"; break }
        loss <- objective(b)
        rate <- 1
        accepted <- FALSE
        for (attempt in seq_len(30L)) {
            candidate <- b - rate * step
            trial <- objective(candidate)
            if (is.finite(trial) && (trial <= loss - 1e-4 * rate * sum(g * step) ||
                (abs(trial - loss) <= 32 * .Machine$double.eps * max(1, abs(loss)) &&
                 max(abs(gradient(candidate))) < max(abs(g))))) {
                accepted <- TRUE
                break
            }
            rate <- rate / 2
        }
        if (!accepted) { status <- "line.search.failed"; break }
        b <- candidate
    }
    residual <- max(abs(gradient(b)))
    if (residual < tol) status <- "converged"
    list(value = stats::plogis(sum(query * b)), status = status,
         gradient = residual, coefficients = b)
}
