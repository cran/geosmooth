#' Construct an SSRHE-Style Local Hessian Energy Operator
#'
#' Constructs the discrete local Hessian operator used by the semi-supervised
#' regularization framework of Kim, Steinke, and Hein (2009). The constructor
#' is intended as a package-facing, auditable R/C++ port of the reference
#' behavior: it returns the sparse row operator \eqn{A}, the quadratic energy
#' matrix \eqn{B=A^\top A}, and optionally the supplemental stabilizer
#' \eqn{B_S}.
#'
#' @param X Numeric matrix with one observation per row. Distances for the
#'   default neighborhood search are Euclidean distances in these coordinates.
#' @param k Integer number of nearest neighbors per point, including the point
#'   itself, used when \code{neighborhood.type = "knn"}. This follows the
#'   SSRHE Matlab convention. For adaptive-radius neighborhoods, \code{k} may be
#'   omitted when \code{adaptive.k.scale} is supplied.
#' @param tangent.dim Integer tangent dimension used for the local PCA chart.
#'   Required when \code{tangent.dim.rule = "fixed"}. If \code{NULL} and
#'   \code{tangent.dim.rule = "eigen.cumulative"}, the dimension is selected
#'   locally from the PCA variance ratio.
#' @param nn.index Optional integer matrix of neighbor indices, with
#'   \code{nrow(X)} rows and \code{k} columns. Each row must contain its center
#'   vertex. If \code{NULL}, Euclidean k-NN including self is computed in C++.
#' @param neighborhood.type Local support rule. \code{"knn"} uses the original
#'   rectangular self-including kNN neighborhoods. \code{"adaptive.radius"}
#'   builds variable-size supports from
#'   \code{dgraphs::create.rknn.graph()}.
#'   \code{"supplied"} uses \code{support.index} directly.
#' @param support.index Optional list of integer vectors, one per row of
#'   \code{X}. Each element gives a variable-size local support and must contain
#'   its center vertex.
#' @param adaptive.k.scale Integer local-scale k used by
#'   \code{dgraphs::create.rknn.graph()} when
#'   \code{neighborhood.type = "adaptive.radius"}.
#' @param radius.rule,radius.factor Adaptive-radius graph parameters passed to
#'   \code{dgraphs::create.rknn.graph()}.
#' @param min.support Optional minimum local support size for adaptive-radius
#'   supports. Defaults to the local quadratic design size plus
#'   \code{support.buffer}, clamped to \code{nrow(X)}.
#' @param max.support Optional maximum local support size. Oversized adaptive
#'   supports are trimmed to the closest \code{max.support} vertices while
#'   keeping the center vertex.
#' @param support.buffer Nonnegative integer added to the local quadratic design
#'   size when choosing the default \code{min.support}.
#' @param support.topup How undersized adaptive-radius supports are enlarged.
#'   \code{"nearest"} appends ambient nearest neighbors. \code{"none"} leaves
#'   supports unchanged and lets the C++ operator reject undersized supports.
#' @param tangent.dim.rule Either \code{"fixed"} or \code{"eigen.cumulative"}.
#' @param eigen.tolerance Cumulative local PCA variance threshold used when
#'   \code{tangent.dim.rule = "eigen.cumulative"} and \code{tangent.dim = NULL}.
#' @param derivative.order Integer derivative order of the local SSRHE-style
#'   operator. \code{2L} is the original local Hessian energy. \code{3L} is an
#'   experimental third-derivative energy whose rows estimate unique symmetric
#'   third-derivative tensor components with factorial and tensor-multiplicity
#'   scaling.
#' @param stabilizer Logical. If \code{TRUE}, also construct the supplemental
#'   stabilizer matrix described in the SSRHE supplement. Currently supported
#'   only for \code{derivative.order = 2L}.
#' @param pinv.tol Nonnegative tolerance multiplier for local pseudoinverses.
#' @param local.solver Local least-squares backend used to map local function
#'   values to derivative coefficients. \code{"auto"} is the default: it uses
#'   normal equations when the local design is full rank and has condition
#'   number no larger than \code{normal.equations.max.condition}, and otherwise
#'   falls back to SVD. \code{"normal.equations"} requests the normal-equation
#'   solve for full-rank local designs and falls back only on hard numerical
#'   failures. \code{"svd"} is the most stable reference path. \code{"qr"} uses
#'   pivoted QR and falls back to SVD for rank-deficient local designs.
#' @param normal.equations.max.condition Positive condition-number guard used
#'   by \code{local.solver = "auto"} before accepting the normal-equation
#'   backend. The default is deliberately conservative because normal-equation
#'   solves square the effective condition number and order-3 near-minimum
#'   local supports can be numerically fragile.
#' @param return.A Logical. If \code{TRUE}, return \eqn{A} as a sparse matrix.
#' @param return.B Logical. If \code{TRUE}, return \eqn{B=A^\top A}.
#' @param return.BS Logical. If \code{TRUE} and \code{stabilizer = TRUE}, return
#'   the supplemental stabilizer \eqn{B_S}.
#' @param return.sparse Logical. If \code{TRUE}, attach \pkg{Matrix} sparse
#'   matrices in addition to raw triplets.
#' @param return.local.diagnostics Logical. If \code{TRUE}, compute additional
#'   R-side local chart diagnostics such as chart distortion and boundary
#'   asymmetry. These diagnostics are useful for operator audits, but can be
#'   skipped in fitting and cross-validation paths for speed.
#' @param return.timing Logical. If \code{TRUE}, attach a phase-level elapsed
#'   time table to the returned operator. The timings separate R-side
#'   validation, neighborhood/support construction, native local operator
#'   construction, optional local diagnostics, sparse matrix assembly, and
#'   output finalization. For adaptive-radius neighborhoods, subphase timings
#'   are also attached in \code{neighborhoods$timing}.
#' @param verbose Logical. If \code{TRUE}, print a short native construction
#'   message.
#'
#' @details
#' For each point \eqn{x_i}, the operator constructs a local PCA chart from a
#' self-including local support, projects the support into
#' \eqn{\mathbb R^m}, and centers the local coordinates at \eqn{x_i}. A local
#' quadratic model is then fitted with the intercept fixed at \eqn{f_i}:
#' \deqn{
#'   f(x_j) - f(x_i)
#'   \approx
#'   \sum_{a \le b} h_{ab} z_{ja} z_{jb}
#'   +
#'   \sum_a g_a z_{ja}.
#' }
#' The local linear map from function values to the quadratic coefficients is
#' the Matlab-compatible
#' \deqn{
#'   \mathrm{RegMat} = X^+ - X^+ \mathrm{IndMat},
#' }
#' where \eqn{X} is the reduced quadratic-plus-linear design and
#' \eqn{\mathrm{IndMat}} copies the center value into all local rows. Diagonal
#' Hessian components are scaled by \eqn{\sqrt{2}} and off-diagonal components
#' by \eqn{1}, so that the row energy matches the reference Hessian energy.
#' With \code{derivative.order = 3L}, the reduced local design uses cubic
#' monomials first, followed by quadratic and linear monomials. The returned
#' rows estimate unique symmetric third-derivative components
#' \eqn{\partial_{abc} f} with scale
#' \deqn{
#'   s_{abc}=m_{abc}\sqrt{\mu_{abc}},
#' }
#' where \eqn{m_{abc}} is the factorial monomial-to-derivative multiplier
#' (\eqn{6} for \eqn{aaa}, \eqn{2} for \eqn{aab}, and \eqn{1} for three
#' distinct indices) and \eqn{\mu_{abc}} is the ordered-tensor multiplicity
#' (\eqn{1}, \eqn{3}, or \eqn{6}). Thus \eqn{\|Af\|_2^2} approximates the
#' squared Frobenius norm of the full symmetric third-derivative tensor.
#'
#' The returned \eqn{A} is the stacked sparse matrix of these local derivative
#' coefficient rows. The \eqn{\ell_2} SSRHE penalty is
#' \deqn{
#'   f^\top B f = \|Af\|_2^2,\quad B=A^\top A.
#' }
#' Exposing \eqn{A} is useful for auditability and for future
#' \eqn{\ell_1}-style variants based on \eqn{\|Af\|_1}.
#'
#' Adaptive-radius graph construction uses \pkg{dgraphs}; fixed-k and supplied
#' neighborhoods are package-local geosmooth paths.
#'
#' @return A list of class \code{"ssrhe.hessian.operator"} containing:
#'   \itemize{
#'     \item \code{A}: sparse local Hessian row operator, if requested;
#'     \item \code{B}: sparse quadratic energy matrix \eqn{A^\top A}, if requested;
#'     \item \code{BS}: optional supplemental stabilizer matrix;
#'     \item \code{A.triplet}, \code{BS.triplet}: raw sparse triplets;
#'     \item \code{nn.index}: neighbor index matrix for rectangular kNN
#'       neighborhoods, if used;
#'     \item \code{support.index}: list of actual local supports used by the
#'       backend;
#'     \item \code{neighborhoods}: support-construction diagnostics;
#'     \item \code{row.table}: one row per Hessian component row of \eqn{A};
#'     \item \code{diagnostics}: local design rank, condition, and PCA summaries;
#'     \item \code{local.diagnostics}: if requested, one row per center vertex with support
#'       size, design rank/condition, PCA variance, chart-distortion,
#'       boundary-asymmetry, and curvature-bias proxy diagnostics;
#'     \item \code{timing}: if requested, a phase-level elapsed-time table;
#'     \item \code{native.timing}: if requested, a C++ backend subphase timing
#'       table for native operator construction;
#'     \item \code{parity}: notes documenting reference-behavior conventions.
#'   }
#'
#' @references
#' Kim, K. I., Steinke, F., and Hein, M. (2009). Semi-supervised regression
#' using Hessian energy with an application to semi-supervised dimensionality
#' reduction. \emph{Advances in Neural Information Processing Systems 22}.
#' \url{https://papers.nips.cc/paper_files/paper/2009/hash/f4552671f8909587cf485ea990207f3b-Abstract.html}
#'
#' @examples
#' X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
#' ssrhe.hessian.operator(X, k = 6L, tangent.dim = 1L)
#' @export
ssrhe.hessian.operator <- function(
    X,
    k = NULL,
    tangent.dim,
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    stabilizer = FALSE,
    pinv.tol = sqrt(.Machine$double.eps),
    local.solver = c("auto", "normal.equations", "svd", "qr"),
    normal.equations.max.condition = 1e4,
    return.A = TRUE,
    return.B = TRUE,
    return.BS = stabilizer,
    return.sparse = TRUE,
    return.local.diagnostics = TRUE,
    return.timing = FALSE,
    verbose = FALSE) {

    timing.total.start <- proc.time()[["elapsed"]]
    timing.phase.start <- timing.total.start
    timing.rows <- list()
    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    neighborhood.type <- match.arg(neighborhood.type)
    radius.rule <- match.arg(radius.rule)
    support.topup <- match.arg(support.topup)
    local.solver <- match.arg(local.solver)
    return.timing <- isTRUE(return.timing)

    tangent.dim.rule <- match.arg(tangent.dim.rule)
    if (missing(tangent.dim) || is.null(tangent.dim)) {
        if (identical(tangent.dim.rule, "fixed")) {
            stop("tangent.dim is required when tangent.dim.rule = 'fixed'.",
                 call. = FALSE)
        }
        tangent.dim.cpp <- 0L
    } else {
        tangent.dim.cpp <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        if (tangent.dim.cpp > ncol(X)) {
            stop("tangent.dim must be no larger than ncol(X).",
                 call. = FALSE)
        }
    }

    eigen.tolerance <- .validate.ssrhe.numeric.scalar(eigen.tolerance, "eigen.tolerance")
    if (eigen.tolerance <= 0 || eigen.tolerance > 1) {
        stop("eigen.tolerance must be in (0, 1].", call. = FALSE)
    }
    pinv.tol <- .validate.ssrhe.numeric.scalar(pinv.tol, "pinv.tol")
    if (pinv.tol < 0) {
        stop("pinv.tol must be nonnegative.", call. = FALSE)
    }
    normal.equations.max.condition <- .validate.ssrhe.numeric.scalar(
        normal.equations.max.condition, "normal.equations.max.condition")
    if (normal.equations.max.condition <= 0) {
        stop("normal.equations.max.condition must be positive.",
             call. = FALSE)
    }
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    stabilizer <- isTRUE(stabilizer)
    if (derivative.order == 3L && stabilizer) {
        stop("stabilizer is currently only supported for derivative.order = 2.",
             call. = FALSE)
    }
    return.A <- isTRUE(return.A)
    return.B <- isTRUE(return.B)
    return.BS <- isTRUE(return.BS)
    return.sparse <- isTRUE(return.sparse)
    return.local.diagnostics <- isTRUE(return.local.diagnostics)
    verbose <- isTRUE(verbose)

    if ((return.A || return.B || return.BS || return.sparse) &&
        !requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required to attach sparse SSRHE operators.",
             call. = FALSE)
    }
    if (!return.sparse && (return.A || return.B || return.BS)) {
        stop("return.sparse must be TRUE when returning A, B, or BS matrices.",
             call. = FALSE)
    }

    if (return.timing) {
        timing.rows[["validation"]] <- .ssrhe.operator.timing.row(
            phase = "validation",
            elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start
        )
        timing.phase.start <- proc.time()[["elapsed"]]
    }

    neighborhood <- .prepare.ssrhe.neighborhood(
        X = X,
        k = k,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim = if (missing(tangent.dim)) NULL else tangent.dim,
        tangent.dim.rule = tangent.dim.rule,
        derivative.order = derivative.order,
        return.timing = return.timing
    )

    if (return.timing) {
        timing.rows[["neighborhood"]] <- .ssrhe.operator.timing.row(
            phase = "neighborhood",
            elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start
        )
        timing.phase.start <- proc.time()[["elapsed"]]
    }

    if (!missing(tangent.dim) && !is.null(tangent.dim) &&
        tangent.dim.cpp > min(min(neighborhood$support.size), ncol(X))) {
        stop("tangent.dim must be no larger than the smallest local support size and ncol(X).",
             call. = FALSE)
    }

    raw <- rcpp_ssrhe_hessian_operator(
        X,
        as.integer(neighborhood$k),
        as.integer(tangent.dim.cpp),
        neighborhood$nn.index,
        neighborhood$support.index,
        tangent.dim.rule,
        as.double(eigen.tolerance),
        as.integer(derivative.order),
        as.logical(stabilizer),
        as.double(pinv.tol),
        local.solver,
        as.double(normal.equations.max.condition),
        as.logical(verbose)
    )

    if (return.timing) {
        timing.rows[["native.operator"]] <- .ssrhe.operator.timing.row(
            phase = "native.operator",
            elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start
        )
        timing.phase.start <- proc.time()[["elapsed"]]
    }

    raw$row.table <- as.data.frame(raw$row.table, stringsAsFactors = FALSE)
    raw$diagnostics <- as.data.frame(raw$diagnostics, stringsAsFactors = FALSE)
    raw$native.timing <- as.data.frame(raw$native.timing, stringsAsFactors = FALSE)
    raw$local.diagnostics <- if (return.local.diagnostics) {
        .ssrhe.local.geometry.diagnostics(
            X = X,
            support.index = raw$support.index,
            diagnostics = raw$diagnostics
        )
    } else {
        NULL
    }

    if (return.timing) {
        timing.rows[["local.diagnostics"]] <- .ssrhe.operator.timing.row(
            phase = "local.diagnostics",
            elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start
        )
        timing.phase.start <- proc.time()[["elapsed"]]
    }

    A <- .ssrhe.triplet.to.sparse(raw$A.triplet)
    B <- Matrix::crossprod(A)
    BS <- NULL
    if (stabilizer && !is.null(raw$BS.triplet)) {
        BS <- .ssrhe.triplet.to.sparse(raw$BS.triplet)
        BS <- (BS + Matrix::t(BS)) / 2
    }

    if (return.timing) {
        timing.rows[["sparse.assembly"]] <- .ssrhe.operator.timing.row(
            phase = "sparse.assembly",
            elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start
        )
        timing.phase.start <- proc.time()[["elapsed"]]
    }

    out <- list(
        A = if (return.A) A else NULL,
        B = if (return.B) B else NULL,
        BS = if (return.BS && stabilizer) BS else NULL,
        A.triplet = raw$A.triplet,
        BS.triplet = raw$BS.triplet,
        nn.index = raw$nn.index,
        support.index = raw$support.index,
        neighborhoods = neighborhood$metadata,
        row.table = raw$row.table,
        diagnostics = raw$diagnostics,
        local.diagnostics = raw$local.diagnostics,
        native.timing = if (return.timing) raw$native.timing else NULL,
        timing = NULL,
        parity = raw$parity,
        parameters = raw$parameters
    )
    out$parameters$neighborhood.type <- neighborhood$type
    out$parameters$k <- neighborhood$k
    out$parameters$adaptive.k.scale <- neighborhood$adaptive.k.scale
    out$parameters$radius.rule <- neighborhood$radius.rule
    out$parameters$radius.factor <- neighborhood$radius.factor
    out$parameters$derivative.order <- derivative.order
    out$parameters$local.solver <- local.solver
    out$parameters$normal.equations.max.condition <- normal.equations.max.condition

    if (return.timing) {
        timing.rows[["output.finalization"]] <- .ssrhe.operator.timing.row(
            phase = "output.finalization",
            elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start
        )
        timing <- do.call(rbind, timing.rows)
        timing$total.elapsed.sec <- proc.time()[["elapsed"]] - timing.total.start
        timing$fraction.of.total <- timing$elapsed.sec / timing$total.elapsed.sec
        rownames(timing) <- NULL
        out$timing <- timing
    }

    class(out) <- c("ssrhe.hessian.operator", "list")
    out
}

.ssrhe.operator.timing.row <- function(phase, elapsed.sec) {
    data.frame(
        phase = phase,
        elapsed.sec = as.numeric(elapsed.sec),
        stringsAsFactors = FALSE
    )
}

.validate.ssrhe.X <- function(X) {
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    storage.mode(X) <- "double"
    if (!is.numeric(X) || nrow(X) < 2L || ncol(X) < 1L) {
        stop("X must be a numeric matrix with at least two rows and one column.",
             call. = FALSE)
    }
    if (any(!is.finite(X))) {
        stop("X must contain only finite values.", call. = FALSE)
    }
    X
}

.validate.ssrhe.positive.integer <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 ||
        abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be a positive integer scalar.", name), call. = FALSE)
    }
    as.integer(round(x))
}

.validate.ssrhe.derivative.order <- function(x) {
    x <- .validate.ssrhe.positive.integer(x, "derivative.order")
    if (!x %in% c(2L, 3L)) {
        stop("derivative.order must be either 2L or 3L.", call. = FALSE)
    }
    x
}

.validate.ssrhe.numeric.scalar <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x)) {
        stop(sprintf("%s must be a finite numeric scalar.", name), call. = FALSE)
    }
    as.double(x)
}

.validate.ssrhe.nn.index <- function(nn.index, n, k) {
    if (is.null(nn.index)) {
        return(NULL)
    }
    if (!is.matrix(nn.index)) {
        nn.index <- as.matrix(nn.index)
    }
    storage.mode(nn.index) <- "integer"
    if (!identical(dim(nn.index), c(n, k))) {
        stop("nn.index must be an nrow(X) by k integer matrix.", call. = FALSE)
    }
    if (any(is.na(nn.index)) || any(nn.index < 1L) || any(nn.index > n)) {
        stop("nn.index must contain integer indices in 1:nrow(X).", call. = FALSE)
    }
    has.center <- vapply(seq_len(n), function(i) any(nn.index[i, ] == i), logical(1))
    if (!all(has.center)) {
        stop("Each row of nn.index must contain its center vertex.", call. = FALSE)
    }
    nn.index
}

.validate.ssrhe.support.index <- function(support.index, n) {
    if (is.null(support.index)) {
        return(NULL)
    }
    if (is.matrix(support.index)) {
        support.index <- lapply(seq_len(nrow(support.index)), function(i) {
            support.index[i, !is.na(support.index[i, ])]
        })
    }
    if (!is.list(support.index) || length(support.index) != n) {
        stop("support.index must be a list with length nrow(X).", call. = FALSE)
    }
    out <- vector("list", n)
    for (i in seq_len(n)) {
        ids <- as.integer(support.index[[i]])
        if (length(ids) < 2L || any(is.na(ids)) || any(ids < 1L) || any(ids > n)) {
            stop("Each support.index element must contain at least two valid vertex indices.",
                 call. = FALSE)
        }
        ids <- unique(ids)
        if (!i %in% ids) {
            stop("Each support.index element must contain its center vertex.",
                 call. = FALSE)
        }
        out[[i]] <- ids
    }
    out
}

.ssrhe.design.ncol <- function(m, derivative.order = 2L) {
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    q2 <- m * (m + 1L) / 2L
    if (derivative.order == 2L) {
        q2 + m
    } else {
        q3 <- m * (m + 1L) * (m + 2L) / 6L
        q3 + q2 + m
    }
}

.ssrhe.default.min.support <- function(tangent.dim,
                                       tangent.dim.rule,
                                       ambient.dim,
                                       support.buffer,
                                       n,
                                       derivative.order = 2L) {
    support.buffer <- .validate.ssrhe.nonnegative.integer(support.buffer,
                                                         "support.buffer")
    if (!is.null(tangent.dim)) {
        m <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
    } else if (identical(tangent.dim.rule, "eigen.cumulative")) {
        m <- ambient.dim
    } else {
        stop("tangent.dim is required to choose the default min.support.",
             call. = FALSE)
    }
    min(n, .ssrhe.design.ncol(m, derivative.order = derivative.order) + support.buffer)
}

#' Build Candidate Support Profiles for SSRHE Adaptive-Radius Tuning
#'
#' Constructs a compact grid of adaptive-radius support profiles for use with
#' \code{support.selection = "cv"} in SSRHE fitting functions. Each row gives an
#' \code{adaptive.k.scale} value and a requested \code{min.support}. The grid is
#' intentionally small by default because support selection nests operator
#' construction inside response cross-validation.
#'
#' @param n Number of observations.
#' @param tangent.dim Tangent dimension used by the SSRHE local polynomial
#'   design.
#' @param derivative.order SSRHE derivative order, currently \code{2L} or
#'   \code{3L}.
#' @param support.buffer Nonnegative integer added to the local design size when
#'   constructing the smallest candidate support.
#' @param max.candidates Maximum number of candidate support profiles to return.
#'
#' @return A data frame with columns \code{adaptive.k.scale},
#'   \code{min.support}, and \code{max.support}. \code{max.support} is
#'   \code{NA} by default, meaning no truncation.
#' @examples
#' ssrhe.support.grid(n = 50L, tangent.dim = 2L, max.candidates = 4L)
#' @export
ssrhe.support.grid <- function(n,
                               tangent.dim,
                               derivative.order = 2L,
                               support.buffer = 2L,
                               max.candidates = 8L) {
    n <- .validate.ssrhe.positive.integer(n, "n")
    tangent.dim <- .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    support.buffer <- .validate.ssrhe.nonnegative.integer(support.buffer,
                                                         "support.buffer")
    max.candidates <- .validate.ssrhe.positive.integer(max.candidates,
                                                      "max.candidates")
    if (n < 3L) {
        stop("n must be at least 3.", call. = FALSE)
    }
    base <- min(n - 1L,
                .ssrhe.design.ncol(tangent.dim, derivative.order) +
                    support.buffer)
    candidates <- unique(as.integer(round(c(
        base,
        base + c(3L, 5L, 10L),
        1.5 * base,
        2 * base,
        3 * base,
        4 * base,
        floor(0.15 * n),
        floor(0.25 * n)
    ))))
    candidates <- sort(unique(pmax(2L, pmin(candidates, n - 1L))))
    if (length(candidates) > max.candidates) {
        keep <- unique(round(seq(1, length(candidates),
                                 length.out = max.candidates)))
        candidates <- candidates[keep]
    }
    adaptive.k.scale <- pmax(2L, floor(candidates / 3L))
    adaptive.k.scale <- pmin(adaptive.k.scale, n - 1L)
    unique(data.frame(
        adaptive.k.scale = as.integer(adaptive.k.scale),
        min.support = as.integer(candidates),
        max.support = NA_integer_
    ))
}

.validate.ssrhe.support.grid <- function(support.grid,
                                         n,
                                         tangent.dim,
                                         derivative.order,
                                         support.buffer,
                                         max.candidates) {
    if (is.null(support.grid)) {
        return(ssrhe.support.grid(
            n = n,
            tangent.dim = tangent.dim,
            derivative.order = derivative.order,
            support.buffer = support.buffer,
            max.candidates = max.candidates
        ))
    }
    if (!is.data.frame(support.grid)) {
        stop("support.grid must be a data frame.", call. = FALSE)
    }
    required <- c("adaptive.k.scale", "min.support")
    if (!all(required %in% names(support.grid))) {
        stop("support.grid must contain adaptive.k.scale and min.support columns.",
             call. = FALSE)
    }
    grid <- support.grid[, unique(c(required, intersect("max.support", names(support.grid)))),
                         drop = FALSE]
    if (!"max.support" %in% names(grid)) {
        grid$max.support <- NA_integer_
    }
    grid$adaptive.k.scale <- as.integer(round(grid$adaptive.k.scale))
    grid$min.support <- as.integer(round(grid$min.support))
    grid$max.support <- ifelse(is.na(grid$max.support), NA_integer_,
                               as.integer(round(grid$max.support)))
    bad <- is.na(grid$adaptive.k.scale) | grid$adaptive.k.scale < 1L |
        grid$adaptive.k.scale >= n |
        is.na(grid$min.support) | grid$min.support < 2L |
        grid$min.support > n |
        (!is.na(grid$max.support) &
             (grid$max.support < grid$min.support | grid$max.support > n))
    if (any(bad)) {
        stop("support.grid contains invalid adaptive.k.scale/min.support/max.support values.",
             call. = FALSE)
    }
    grid <- unique(grid)
    if (nrow(grid) > max.candidates) {
        grid <- grid[seq_len(max.candidates), , drop = FALSE]
    }
    rownames(grid) <- NULL
    grid
}

.ssrhe.support.diagnostics <- function(operator) {
    support.size <- operator$neighborhoods$support.size
    topup <- operator$neighborhoods$n.topup
    truncated <- operator$neighborhoods$n.truncated
    data.frame(
        support.size.min = min(support.size),
        support.size.median = stats::median(support.size),
        support.size.max = max(support.size),
        support.topup.median = stats::median(topup, na.rm = TRUE),
        support.topup.max = max(topup, na.rm = TRUE),
        support.truncated.max = max(truncated, na.rm = TRUE),
        operator.rows = operator$A.triplet$dim[[1L]],
        operator.cols = operator$A.triplet$dim[[2L]],
        stringsAsFactors = FALSE
    )
}

.ssrhe.local.geometry.diagnostics <- function(X, support.index, diagnostics) {
    n <- nrow(X)
    if (!is.list(support.index) || length(support.index) != n) {
        stop("Internal error: invalid SSRHE support.index.", call. = FALSE)
    }
    out <- diagnostics
    support.size <- out$k
    if (is.null(support.size)) {
        support.size <- lengths(support.index)
    }
    out$support.size <- support.size
    out$design.rank.deficiency <- pmax(0L, out$design.ncol - out$design.rank)
    out$pca.variance.explained <- out$local.variance.ratio
    out$pca.discarded.variance.ratio <- pmax(0, 1 - out$local.variance.ratio)

    chart.distortion <- rep(NA_real_, n)
    boundary.asymmetry <- rep(NA_real_, n)
    curvature.bias.proxy <- out$pca.discarded.variance.ratio

    for (center in seq_len(n)) {
        ids <- as.integer(support.index[[center]])
        if (!length(ids) || !center %in% ids) next
        local <- X[ids, , drop = FALSE]
        centered <- sweep(local, 2L, colMeans(local), "-")
        svd.local <- tryCatch(svd(centered), error = function(e) NULL)
        if (is.null(svd.local) || !length(svd.local$d)) next
        m <- as.integer(out$tangent.dim[center])
        m <- max(1L, min(m, ncol(svd.local$v), length(svd.local$d)))
        coords <- centered %*% svd.local$v[, seq_len(m), drop = FALSE]
        base.local <- match(center, ids)
        coords <- sweep(coords, 2L, coords[base.local, ], "-")

        radius <- sqrt(max(rowSums(coords^2), na.rm = TRUE))
        centroid.norm <- sqrt(sum(colMeans(coords)^2))
        boundary.asymmetry[center] <- if (is.finite(radius) && radius > 0) {
            centroid.norm / radius
        } else {
            0
        }

        if (length(ids) >= 3L) {
            d.ambient <- stats::dist(local)
            d.chart <- stats::dist(coords)
            d.ambient <- as.numeric(d.ambient)
            d.chart <- as.numeric(d.chart)
            ok <- is.finite(d.ambient) & is.finite(d.chart) & d.ambient > 0
            denom <- mean(d.ambient[ok])
            chart.distortion[center] <- if (any(ok) && is.finite(denom) && denom > 0) {
                sqrt(mean((d.chart[ok] - d.ambient[ok])^2)) / denom
            } else {
                0
            }
        } else {
            chart.distortion[center] <- 0
        }
    }

    out$chart.distortion <- chart.distortion
    out$boundary.asymmetry <- boundary.asymmetry
    out$curvature.bias.proxy <- curvature.bias.proxy
    out
}

.validate.ssrhe.nonnegative.integer <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 ||
        abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be a nonnegative integer scalar.", name),
             call. = FALSE)
    }
    as.integer(round(x))
}

.prepare.ssrhe.neighborhood <- function(X,
                                        k,
                                        nn.index,
                                        neighborhood.type,
                                        support.index,
                                        adaptive.k.scale,
                                        radius.rule,
                                        radius.factor,
                                        min.support,
                                        max.support,
                                        support.buffer,
                                        support.topup,
                                        tangent.dim,
                                        tangent.dim.rule,
                                        derivative.order = 2L,
                                        return.timing = FALSE) {
    return.timing <- isTRUE(return.timing)
    timing.start <- proc.time()[["elapsed"]]
    timing.phase.start <- timing.start
    timing.rows <- list()
    add.timing <- function(subphase) {
        if (return.timing) {
            timing.rows[[length(timing.rows) + 1L]] <<-
                .ssrhe.operator.timing.row(
                    phase = paste0("neighborhood.", subphase),
                    elapsed.sec = proc.time()[["elapsed"]] - timing.phase.start
                )
            timing.phase.start <<- proc.time()[["elapsed"]]
        }
        invisible(NULL)
    }
    finalize.timing <- function() {
        if (!return.timing) {
            return(NULL)
        }
        timing <- do.call(rbind, timing.rows)
        if (is.null(timing)) {
            timing <- data.frame(phase = character(),
                                 elapsed.sec = numeric(),
                                 stringsAsFactors = FALSE)
        }
        total.elapsed <- proc.time()[["elapsed"]] - timing.start
        timing$total.elapsed.sec <- total.elapsed
        timing$fraction.of.total <- if (total.elapsed > 0) {
            timing$elapsed.sec / total.elapsed
        } else {
            0
        }
        rownames(timing) <- NULL
        timing
    }

    n <- nrow(X)
    if (!is.null(support.index) &&
        identical(neighborhood.type, "knn") &&
        is.null(nn.index) && is.null(k)) {
        neighborhood.type <- "supplied"
    }

    if (identical(neighborhood.type, "knn")) {
        if (is.null(k)) {
            if (!is.null(nn.index)) {
                k <- ncol(as.matrix(nn.index))
            } else {
                stop("k is required when neighborhood.type = 'knn'.",
                     call. = FALSE)
            }
        }
        k <- .validate.ssrhe.positive.integer(k, "k")
        if (k < 2L || k > n) {
            stop("k must be between 2 and nrow(X).", call. = FALSE)
        }
        nn.index <- .validate.ssrhe.nn.index(nn.index, n, k)
        support.size <- rep(k, n)
        add.timing("validation")
        return(list(
            type = "knn",
            k = k,
            nn.index = nn.index,
            support.index = NULL,
            support.size = support.size,
            adaptive.k.scale = NA_integer_,
            radius.rule = NA_character_,
            radius.factor = NA_real_,
            metadata = list(
                type = "knn",
                support.size = support.size,
                n.topup = rep(0L, n),
                n.truncated = rep(0L, n),
                timing = finalize.timing()
            )
        ))
    }

    if (identical(neighborhood.type, "supplied")) {
        support.index <- .validate.ssrhe.support.index(support.index, n)
        support.size <- lengths(support.index)
        add.timing("validation")
        return(list(
            type = "supplied",
            k = 0L,
            nn.index = NULL,
            support.index = support.index,
            support.size = support.size,
            adaptive.k.scale = NA_integer_,
            radius.rule = NA_character_,
            radius.factor = NA_real_,
            metadata = list(
                type = "supplied",
                support.size = support.size,
                n.topup = rep(0L, n),
                n.truncated = rep(0L, n),
                timing = finalize.timing()
            )
        ))
    }

    if (!identical(neighborhood.type, "adaptive.radius")) {
        stop("Unknown SSRHE neighborhood type.", call. = FALSE)
    }
    if (is.null(adaptive.k.scale)) {
        if (is.null(k)) {
            stop("adaptive.k.scale is required when neighborhood.type = 'adaptive.radius'.",
                 call. = FALSE)
        }
        adaptive.k.scale <- k
    }
    adaptive.k.scale <- .validate.ssrhe.positive.integer(adaptive.k.scale,
                                                        "adaptive.k.scale")
    if (adaptive.k.scale >= n) {
        stop("adaptive.k.scale must be smaller than nrow(X).", call. = FALSE)
    }
    radius.factor <- .validate.ssrhe.numeric.scalar(radius.factor, "radius.factor")
    if (radius.factor <= 0) {
        stop("radius.factor must be positive.", call. = FALSE)
    }
    if (is.null(min.support)) {
        min.support <- .ssrhe.default.min.support(
            tangent.dim = tangent.dim,
            tangent.dim.rule = tangent.dim.rule,
            ambient.dim = ncol(X),
            support.buffer = support.buffer,
            n = n,
            derivative.order = derivative.order
        )
    } else {
        min.support <- .validate.ssrhe.positive.integer(min.support,
                                                       "min.support")
    }
    if (min.support < 2L || min.support > n) {
        stop("min.support must be between 2 and nrow(X).", call. = FALSE)
    }
    if (!is.null(max.support)) {
        max.support <- .validate.ssrhe.positive.integer(max.support,
                                                       "max.support")
        if (max.support < min.support || max.support > n) {
            stop("max.support must be between min.support and nrow(X).",
                 call. = FALSE)
        }
    }
    add.timing("validation")

    graph <- dgraphs::create.rknn.graph(
        X = X,
        type = "adaptive.radius",
        k.scale = adaptive.k.scale,
        radius.factor = radius.factor,
        radius.rule = radius.rule,
        prune.method = "none",
        connect.components = FALSE,
        return.timing = return.timing,
        graph.detail = "minimal"
    )
    add.timing("create.graph")

    graph.adj <- .geosmooth.graph.payload(graph)$adj.list
    support.index <- lapply(seq_len(n), function(i) {
        unique(c(i, graph.adj[[i]]))
    })
    add.timing("initial.supports")

    n.topup <- integer(n)
    if (identical(support.topup, "nearest")) {
        undersized <- which(lengths(support.index) < min.support)
        for (i in undersized) {
            ids <- support.index[[i]]
            need <- min.support - length(ids)
            add <- .ssrhe.nearest.outside.support(
                X = X,
                center = i,
                exclude = ids,
                n.add = need
            )
            ids <- unique(c(ids, add))
            n.topup[[i]] <- max(0L, length(ids) - length(support.index[[i]]))
            support.index[[i]] <- ids
        }
    }
    add.timing("topup")

    n.truncated <- integer(n)
    for (i in seq_len(n)) {
        ids <- support.index[[i]]
        if (!is.null(max.support) && length(ids) > max.support) {
            non.center <- setdiff(ids, i)
            d <- rowSums((t(t(X[non.center, , drop = FALSE]) - X[i, ]))^2)
            keep <- non.center[order(d, non.center)][seq_len(max.support - 1L)]
            n.truncated[[i]] <- length(ids) - max.support
            ids <- c(i, keep)
        } else if (ids[[1L]] != i) {
            ids <- c(i, setdiff(ids, i))
        }
        support.index[[i]] <- as.integer(ids)
    }
    add.timing("truncate.reorder")

    support.index <- .validate.ssrhe.support.index(support.index, n)
    support.size <- lengths(support.index)
    add.timing("final.validation")
    list(
        type = "adaptive.radius",
        k = 0L,
        nn.index = NULL,
        support.index = support.index,
        support.size = support.size,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        metadata = list(
            type = "adaptive.radius",
            graph = graph,
            adaptive.k.scale = adaptive.k.scale,
            radius.rule = radius.rule,
            radius.factor = radius.factor,
            min.support = min.support,
            max.support = max.support,
            support.topup = support.topup,
            support.size = support.size,
            n.topup = n.topup,
            n.truncated = n.truncated,
            sigma = graph$sigma %||% graph$metadata$sigma,
            graph.timing = graph$timing %||% graph$metadata$timing,
            timing = finalize.timing()
        )
    )
}

.ssrhe.nearest.outside.support <- function(X, center, exclude, n.add) {
    if (n.add <= 0L) {
        return(integer())
    }
    n <- nrow(X)
    d <- rowSums((t(t(X) - X[center, ]))^2)
    d[unique(c(center, exclude))] <- Inf
    candidates <- order(d, seq_len(n))
    candidates <- candidates[is.finite(d[candidates])]
    as.integer(candidates[seq_len(min(n.add, length(candidates)))])
}

.ssrhe.triplet.to.sparse <- function(triplet) {
    Matrix::sparseMatrix(
        i = triplet$i,
        j = triplet$j,
        x = triplet$x,
        dims = triplet$dim
    )
}

#' @rdname geosmooth-print-methods
#' @method print ssrhe.hessian.operator
#' @export
print.ssrhe.hessian.operator <- function(x, ...) {
    n <- if (!is.null(x$B)) ncol(x$B) else x$A.triplet$dim[[2]]
    nrow.A <- x$A.triplet$dim[[1]]
    cat("SSRHE Hessian operator\n")
    cat("  vertices:", n, "\n")
    cat("  operator rows:", nrow.A, "\n")
    cat("  neighborhood.type:", x$parameters$neighborhood.type %||% "knn", "\n")
    if (identical(x$parameters$neighborhood.type, "knn")) {
        cat("  k:", x$parameters$k, "\n")
    } else if (!is.null(x$neighborhoods$support.size)) {
        cat("  support size:",
            paste(range(x$neighborhoods$support.size), collapse = "-"),
            "\n")
    }
    cat("  tangent.dim.rule:", x$parameters$tangent.dim.rule, "\n")
    cat("  derivative.order:", x$parameters$derivative.order %||% 2L, "\n")
    cat("  stabilizer:", isTRUE(x$parameters$stabilizer), "\n")
    invisible(x)
}

#' Fit SSRHE-Style Hessian-Energy Regression
#'
#' Fits the \eqn{\ell_2} Hessian-energy regularized estimator associated with
#' \code{\link{ssrhe.hessian.operator}}. This is a direct SSRHE-style
#' comparator for Hessian smoothing on point clouds; it is not a replacement
#' for graph trend-filtering regression and is distinct from
#' \eqn{\ell_1}-adaptive graph trend filtering.
#'
#' @inheritParams ssrhe.hessian.operator
#' @param y Numeric response vector of length \code{nrow(X)} or numeric matrix
#'   with \code{nrow(X)} rows. Matrix columns are fit as separate responses.
#'   \code{NA} values are allowed and are treated as unobserved by setting the
#'   corresponding observation weight to zero.
#' @param lambda1 Nonnegative Hessian-energy penalty multiplier for
#'   \eqn{f^\top B f = \|Af\|_2^2}.
#' @param lambda2 Nonnegative supplemental-stabilizer multiplier for
#'   \eqn{f^\top B_S f}. If positive, \code{stabilizer} must be \code{TRUE}.
#'   Currently supported only for \code{derivative.order = 2L}.
#' @param weights Optional nonnegative observation weights. May be \code{NULL},
#'   a vector of length \code{nrow(X)}, or a matrix with the same dimensions as
#'   \code{y}.
#' @param ridge Nonnegative diagonal ridge added to the linear system for
#'   numerical stabilization.
#' @param lambda.selection Penalty selection: \code{"fixed"} (the default)
#'   uses \code{lambda1} and \code{lambda2}; \code{"cv"} holds out observed
#'   labels; \code{"gcv"} uses generalized cross-validation. Grid inputs
#'   require an explicit \code{"cv"} or \code{"gcv"} selection mode.
#' @param lambda1.grid,lambda2.grid Nonnegative penalty candidates for CV or
#'   GCV. \code{lambda1.grid} is required in either selection mode;
#'   \code{lambda2.grid = NULL} means \code{0}. Do not also supply scalar
#'   \code{lambda1} or \code{lambda2}. Unless explicitly set, \code{stabilizer}
#'   is enabled when any supplemental penalty candidate is positive.
#' @param cv.control Named list used only with \code{lambda.selection = "cv"}.
#'   Entries: \code{foldid = NULL} (one fold label per row of \code{X};
#'   nonpositive or missing labels are ignored), \code{cv.folds = 5L},
#'   \code{loss = "mse"} (or \code{"mae"}), and
#'   \code{selection = "min"} (or \code{"one.se"}, choosing the largest
#'   total penalty within one standard error of the minimum). Support-search
#'   entries are \code{support.selection = "rule"} (or \code{"cv"}),
#'   \code{support.grid = NULL}, and \code{support.max.candidates = 8L}.
#'   Fold count is ignored when \code{foldid} is supplied.
#' @param gcv.control Named list used only with \code{lambda.selection = "gcv"}.
#'   Entries: \code{trace.method = "exact"} (or \code{"hutchinson"}),
#'   \code{trace.n.probes = 50L}, \code{trace.seed = NULL},
#'   \code{support.selection = "rule"} (or \code{"gcv"}),
#'   \code{support.grid = NULL}, and \code{support.max.candidates = 8L}.
#'   Unknown, unnamed, or duplicated control entries are errors. Controls
#'   for a different selection mode are rejected rather than ignored.
#'
#' @details
#' For each response column, this function solves
#' \deqn{
#'   (W + \lambda_1 B + \lambda_2 B_S + \epsilon I)\hat f = Wy,
#' }
#' where \eqn{W} is the diagonal matrix of observation weights and
#' \eqn{\epsilon} is \code{ridge}. Fully observed data with
#' \code{lambda1 = lambda2 = ridge = 0} therefore reproduce the observed
#' response exactly. Missing responses are excluded from the data-fit term by
#' setting their weights to zero. The Matlab SSRHE semi-supervised convention is
#' represented by a \code{0}/\code{1} labeled-indicator \code{weights} vector, or
#' equivalently by setting unlabeled responses to \code{NA}: labeled vertices
#' contribute unit diagonal data-fit terms and unlabeled vertices contribute no
#' data-fit term.
#'
#' The lower-level operator is matched to the Kim--Steinke--Hein SSRHE Matlab
#' construction: self-including kNN neighborhoods, local PCA charts, a
#' fixed-intercept local quadratic fit, doubled diagonal Hessian components,
#' and optional supplemental stabilizer. The package exposes both \eqn{A} and
#' \eqn{B=A^\top A}; the fitted \eqn{\ell_2} estimator uses \eqn{B}.
#'
#' @section Penalty selection:
#' Fixed fitting supports response vectors or matrices. CV and GCV currently
#' select penalties for one response vector at a time. Label CV removes each
#' validation fold from the data-fit term and scores its held-out observed
#' positive-weight labels. It supports missing responses. GCV requires a fully
#' observed response and strictly positive observation weights, and minimizes
#' \deqn{\frac{n^{-1}\sum_i w_i(y_i-\hat f_i)^2}
#'   {(1-\mathrm{tr}(S_\lambda)/n)^2},}
#' where \eqn{S_\lambda=(W+\lambda_1 B+\lambda_2 B_S+\epsilon I)^{-1}W}.
#' The exact trace is deterministic. The Hutchinson option uses Rademacher
#' probes, reports trace standard errors, and is reproducible with
#' \code{gcv.control = list(trace.method = "hutchinson", trace.seed = ...)}.
#'
#' Both selection modes can also search adaptive-radius support profiles:
#' set \code{support.selection} to the same mode in its control list. Each row
#' of \code{support.grid} contains \code{adaptive.k.scale}, \code{min.support},
#' and optionally \code{max.support}. A \code{NULL} grid is constructed by
#' \code{\link{ssrhe.support.grid}}. Each profile requires a new operator and
#' a full penalty search, so this can substantially increase runtime.
#'
#' @section Migration:
#' The former \code{fit.ssrhe.hessian.regression.cv()} and
#' \code{fit.ssrhe.hessian.regression.gcv()} entry points have been retired.
#' Use this function with \code{lambda.selection = "cv"} or \code{"gcv"}.
#' Move \code{fold.id} to \code{cv.control$foldid}, \code{nfolds} to
#' \code{cv.control$cv.folds}, and \code{loss}/\code{selection} to
#' \code{cv.control}. Move \code{gcv.trace.method}, \code{gcv.trace.n.probes},
#' and \code{gcv.trace.seed} to \code{trace.method}, \code{trace.n.probes},
#' and \code{trace.seed} in \code{gcv.control}. Support search settings also
#' belong in the active control list; both former support candidate limits
#' are now named \code{support.max.candidates}.
#'
#' @return A list of class \code{"ssrhe.hessian.fit"} containing fitted values,
#'   residuals, input response, weights, lambda parameters, objective/energy
#'   diagnostics, the reused \code{operator}, solver metadata, the call, and
#'   \code{lambda.selection}. CV adds class \code{"ssrhe.hessian.cv.fit"},
#'   \code{cv.table}, \code{fold.id}, and \code{selection}; GCV adds class
#'   \code{"ssrhe.hessian.gcv.fit"}, \code{gcv.table}, \code{gcv.trace}, and
#'   \code{selection}. Support searches retain \code{selected.support} and
#'   the corresponding \code{support.cv.table} or \code{support.gcv.table}.
#'
#' @references
#' Kim, K. I., Steinke, F., and Hein, M. (2009). Semi-supervised regression
#' using Hessian energy with an application to semi-supervised dimensionality
#' reduction. \emph{Advances in Neural Information Processing Systems 22}.
#' \url{https://papers.nips.cc/paper_files/paper/2009/hash/f4552671f8909587cf485ea990207f3b-Abstract.html}
#'
#' @examples
#' X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
#' fit.ssrhe.hessian.regression(
#'   X, sin(2 * pi * X[, 1]), k = 6L, tangent.dim = 1L, lambda1 = 0.1
#' )
#' fit.ssrhe.hessian.regression(
#'   X, sin(2 * pi * X[, 1]), k = 6L, tangent.dim = 1L,
#'   lambda.selection = "cv", lambda1.grid = c(0.01, 0.1),
#'   cv.control = list(cv.folds = 2L)
#' )
#' fit.ssrhe.hessian.regression(
#'   X, sin(2 * pi * X[, 1]), k = 6L, tangent.dim = 1L,
#'   lambda.selection = "gcv", lambda1.grid = c(0.01, 0.1)
#' )
#' @aliases fit.ssrhe.hessian.regression.cv fit.ssrhe.hessian.regression.gcv
#' @section Choosing controls:
#' **First fit:** supply `X`, `y`, neighborhood size `k`, `tangent.dim`, and
#' `lambda1`; use `lambda.selection = "fixed"` for a specified penalty.
#'
#' **Selection:** `lambda.selection`, `lambda1.grid`, and `lambda2.grid` select
#' penalties. Fold settings belong in `cv.control`; trace settings belong in
#' `gcv.control`. Label CV and GCV are different criteria. Support selection
#' rebuilds the operator for each candidate and multiplies this work.
#'
#' **Geometry:** neighborhood, adaptive-radius, tangent and derivative controls
#' determine the operator. Reuse it when fitting multiple responses.
#'
#' **Numerics:** ridge and local conditioning controls stabilize solves.
#' Exact GCV traces require more solves than stochastic trace estimates; the
#' latter report sampling uncertainty and need a reproducible probe seed.
#'
#' **Diagnostics:** inspect the selection table, energies, observed-label mask,
#' and solver metadata. [fitted()] and [residuals()] preserve response matrices.
#'
#' @export
fit.ssrhe.hessian.regression <- function(
    X,
    y,
    k = NULL,
    tangent.dim,
    lambda1,
    lambda2 = 0,
    weights = NULL,
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    stabilizer = lambda2 > 0,
    pinv.tol = sqrt(.Machine$double.eps),
    local.solver = c("auto", "normal.equations", "svd", "qr"),
    normal.equations.max.condition = 1e4,
    ridge = 0,
    return.A = TRUE,
    return.local.diagnostics = FALSE,
    return.timing = FALSE,
    verbose = FALSE,
    lambda.selection = c("fixed", "cv", "gcv"),
    lambda1.grid = NULL,
    lambda2.grid = NULL,
    cv.control = list(),
    gcv.control = list()) {

    lambda.selection <- match.arg(lambda.selection)
    cv <- .ssrhe.selection.control(cv.control, list(
        foldid = NULL, cv.folds = 5L, loss = "mse", selection = "min",
        support.selection = "rule", support.grid = NULL,
        support.max.candidates = 8L
    ), "cv.control")
    gcv <- .ssrhe.selection.control(gcv.control, list(
        trace.method = "exact", trace.n.probes = 50L, trace.seed = NULL,
        support.selection = "rule", support.grid = NULL,
        support.max.candidates = 8L
    ), "gcv.control")
    if (length(cv.control) && lambda.selection != "cv") {
        stop("cv.control requires lambda.selection = 'cv'.", call. = FALSE)
    }
    if (length(gcv.control) && lambda.selection != "gcv") {
        stop("gcv.control requires lambda.selection = 'gcv'.", call. = FALSE)
    }
    if (lambda.selection != "fixed") {
        if (!missing(lambda1) || !missing(lambda2)) {
            stop("Supply penalty grids, not scalar lambda1/lambda2, when lambda.selection is 'cv' or 'gcv'.",
                 call. = FALSE)
        }
        if (is.null(lambda1.grid)) {
            stop("lambda1.grid is required when lambda.selection is 'cv' or 'gcv'.",
                 call. = FALSE)
        }
        lambda1.grid <- .validate.ssrhe.lambda.grid(lambda1.grid, "lambda1.grid")
        if (is.null(lambda2.grid)) lambda2.grid <- 0
        lambda2.grid <- .validate.ssrhe.lambda.grid(lambda2.grid, "lambda2.grid")
        if (missing(stabilizer)) stabilizer <- any(lambda2.grid > 0)
        common <- list(
            X = X, y = y, k = k,
            lambda1.grid = lambda1.grid, lambda2.grid = lambda2.grid,
            weights = weights, nn.index = nn.index,
            neighborhood.type = match.arg(neighborhood.type),
            support.index = support.index, adaptive.k.scale = adaptive.k.scale,
            radius.rule = match.arg(radius.rule), radius.factor = radius.factor,
            min.support = min.support, max.support = max.support,
            support.buffer = support.buffer, support.topup = match.arg(support.topup),
            tangent.dim.rule = match.arg(tangent.dim.rule),
            eigen.tolerance = eigen.tolerance, derivative.order = derivative.order,
            stabilizer = stabilizer, pinv.tol = pinv.tol,
            local.solver = match.arg(local.solver),
            normal.equations.max.condition = normal.equations.max.condition,
            ridge = ridge, return.A = return.A,
            return.local.diagnostics = return.local.diagnostics,
            return.timing = return.timing, verbose = verbose
        )
        # Preserve the operator's existing missing-dimension behavior.
        if (!missing(tangent.dim)) common["tangent.dim"] <- list(tangent.dim)
        if (lambda.selection == "cv") {
            out <- do.call(.fit.ssrhe.hessian.regression.cv, c(common, list(
                fold.id = cv$foldid, nfolds = cv$cv.folds,
                loss = cv$loss, selection = cv$selection,
                support.selection = cv$support.selection,
                support.grid = cv$support.grid,
                support.cv.max.candidates = cv$support.max.candidates
            )))
        } else {
            out <- do.call(.fit.ssrhe.hessian.regression.gcv, c(common, list(
                gcv.trace.method = gcv$trace.method,
                gcv.trace.n.probes = gcv$trace.n.probes,
                gcv.trace.seed = gcv$trace.seed,
                support.selection = gcv$support.selection,
                support.grid = gcv$support.grid,
                support.gcv.max.candidates = gcv$support.max.candidates
            )))
        }
        attr(out, "call") <- match.call()
        out$lambda.selection <- lambda.selection
        return(out)
    }
    if (!is.null(lambda1.grid) || !is.null(lambda2.grid)) {
        stop("Penalty grids require lambda.selection = 'cv' or 'gcv'.", call. = FALSE)
    }
    if (missing(lambda1)) {
        stop("lambda1 is required when lambda.selection = 'fixed'.", call. = FALSE)
    }
    X <- .validate.ssrhe.X(X)
    lambda1 <- .validate.ssrhe.nonnegative.scalar(lambda1, "lambda1")
    lambda2 <- .validate.ssrhe.nonnegative.scalar(lambda2, "lambda2")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    local.solver <- match.arg(local.solver)
    normal.equations.max.condition <- .validate.ssrhe.numeric.scalar(
        normal.equations.max.condition, "normal.equations.max.condition")
    if (normal.equations.max.condition <= 0) {
        stop("normal.equations.max.condition must be positive.",
             call. = FALSE)
    }
    stabilizer <- isTRUE(stabilizer)
    if (lambda2 > 0 && !stabilizer) {
        stop("lambda2 > 0 requires stabilizer = TRUE.", call. = FALSE)
    }
    if (derivative.order == 3L && (lambda2 > 0 || stabilizer)) {
        stop("lambda2/stabilizer is currently only supported for derivative.order = 2.",
             call. = FALSE)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        derivative.order = derivative.order,
        stabilizer = stabilizer,
        pinv.tol = pinv.tol,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        return.A = return.A,
        return.B = TRUE,
        return.BS = stabilizer,
        return.sparse = TRUE,
        return.local.diagnostics = return.local.diagnostics,
        return.timing = return.timing,
        verbose = verbose
    )

    out <- .fit.ssrhe.hessian.from.operator(
        operator = operator,
        y = y,
        lambda1 = lambda1,
        lambda2 = lambda2,
        weights = weights,
        ridge = ridge,
        verbose = verbose
    )
    out$X <- X
    out$lambda.selection <- lambda.selection
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.fit", "list")
    out
}

# Validate control names before any geometry is built; retain explicit NULLs.
.ssrhe.selection.control <- function(control, defaults, name) {
    if (!is.list(control) || is.data.frame(control) ||
        (length(control) && (is.null(names(control)) ||
         anyNA(names(control)) || any(!nzchar(names(control))) ||
         anyDuplicated(names(control))))) {
        stop(name, " must be a named list with unique, nonempty names.", call. = FALSE)
    }
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) {
        stop("Unknown ", name, " entries: ", paste(unknown, collapse = ", "),
             call. = FALSE)
    }
    defaults[names(control)] <- control
    defaults
}

#' Refit SSRHE-Style Hessian-Energy Regression
#'
#' Reuses the operator from \code{\link{fit.ssrhe.hessian.regression}} to fit
#' new responses or new fixed penalty weights without rebuilding local PCA
#' neighborhoods or Hessian-energy matrices.
#'
#' @param object A \code{"ssrhe.hessian.fit"} object.
#' @param y New numeric response vector or matrix with one row per vertex.
#'   If \code{NULL}, the original response is reused.
#' @inheritParams fit.ssrhe.hessian.regression
#'
#' @param ... Additional arguments are not currently supported.
#' @return A list of class \code{"ssrhe.hessian.refit"}.
#' @examples
#' X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
#' fit <- fit.ssrhe.hessian.regression(
#'   X, X[, 1]^2, k = 6L, tangent.dim = 1L, lambda1 = 0.1
#' )
#' refit(fit, y = X[, 1]^3)
#' @seealso \code{\link{refit}} for supported fitted objects and migration.
#' @method refit ssrhe.hessian.fit
#' @aliases refit.ssrhe.hessian.regression
#' @export
refit.ssrhe.hessian.fit <- function(object,
                                    y = NULL,
                                    lambda1 = object$lambda$lambda1,
                                    lambda2 = object$lambda$lambda2,
                                    weights = NULL,
                                    ridge = object$lambda$ridge,
                                    verbose = FALSE, ...) {
    .geosmooth.check.dots(...)
    if (!inherits(object, "ssrhe.hessian.fit")) {
        stop("object must be a 'ssrhe.hessian.fit' object.", call. = FALSE)
    }
    if (is.null(object$operator) ||
        !inherits(object$operator, "ssrhe.hessian.operator")) {
        stop("object does not contain a reusable SSRHE operator.",
             call. = FALSE)
    }
    if (is.null(y)) {
        y <- object$y
    }
    if (is.null(weights)) {
        weights <- object$weights
    }
    out <- .fit.ssrhe.hessian.from.operator(
        operator = object$operator,
        y = y,
        lambda1 = lambda1,
        lambda2 = lambda2,
        weights = weights,
        ridge = ridge,
        verbose = verbose
    )
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.refit", "list")
    out
}

# Internal CV selection engine; public controls are validated and mapped above.
.fit.ssrhe.hessian.regression.cv <- function(
    X,
    y,
    k = NULL,
    tangent.dim,
    lambda1.grid,
    lambda2.grid = 0,
    weights = NULL,
    nfolds = 5L,
    fold.id = NULL,
    loss = c("mse", "mae"),
    selection = c("min", "one.se"),
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    stabilizer = any(lambda2.grid > 0),
    pinv.tol = sqrt(.Machine$double.eps),
    local.solver = c("auto", "normal.equations", "svd", "qr"),
    normal.equations.max.condition = 1e4,
    ridge = 0,
    return.A = TRUE,
    return.local.diagnostics = FALSE,
    return.timing = FALSE,
    support.selection = c("rule", "cv"),
    support.grid = NULL,
    support.cv.max.candidates = 8L,
    verbose = FALSE) {

    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for SSRHE CV regression.",
             call. = FALSE)
    }
    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("fit.ssrhe.hessian.regression(lambda.selection = 'cv') currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    lambda1.grid <- .validate.ssrhe.lambda.grid(lambda1.grid, "lambda1.grid")
    lambda2.grid <- .validate.ssrhe.lambda.grid(lambda2.grid, "lambda2.grid")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    loss <- match.arg(loss)
    selection <- match.arg(selection)
    support.selection <- match.arg(support.selection)
    neighborhood.type <- match.arg(neighborhood.type)
    local.solver <- match.arg(local.solver)
    normal.equations.max.condition <- .validate.ssrhe.numeric.scalar(
        normal.equations.max.condition, "normal.equations.max.condition")
    if (normal.equations.max.condition <= 0) {
        stop("normal.equations.max.condition must be positive.",
             call. = FALSE)
    }
    support.cv.max.candidates <- .validate.ssrhe.positive.integer(
        support.cv.max.candidates, "support.cv.max.candidates")
    stabilizer <- isTRUE(stabilizer)
    if (any(lambda2.grid > 0) && !stabilizer) {
        stop("positive lambda2.grid values require stabilizer = TRUE.",
             call. = FALSE)
    }
    if (derivative.order == 3L && (any(lambda2.grid > 0) || stabilizer)) {
        stop("lambda2.grid/stabilizer is currently only supported for derivative.order = 2.",
             call. = FALSE)
    }
    if (identical(support.selection, "cv")) {
        if (!identical(neighborhood.type, "adaptive.radius")) {
            stop("support.selection = 'cv' is currently supported only for neighborhood.type = 'adaptive.radius'.",
                 call. = FALSE)
        }
        tangent.dim.grid <- if (missing(tangent.dim) || is.null(tangent.dim)) {
            ncol(X)
        } else {
            .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        }
        support.grid <- .validate.ssrhe.support.grid(
            support.grid = support.grid,
            n = n,
            tangent.dim = tangent.dim.grid,
            derivative.order = derivative.order,
            support.buffer = support.buffer,
            max.candidates = support.cv.max.candidates
        )
        fold.id.cv <- .prepare.ssrhe.cv.folds(
            observed = as.vector(y.info$observed & W > 0),
            nfolds = nfolds,
            fold.id = fold.id
        )
        candidate.fits <- vector("list", nrow(support.grid))
        candidate.rows <- vector("list", nrow(support.grid))
        for (ii in seq_len(nrow(support.grid))) {
            max.support.ii <- support.grid$max.support[ii]
            if (is.na(max.support.ii)) max.support.ii <- NULL
            elapsed <- system.time({
                cand <- tryCatch(
                    .fit.ssrhe.hessian.regression.cv(
                        X = X,
                        y = y,
                        k = k,
                        tangent.dim = if (missing(tangent.dim)) NULL else tangent.dim,
                        lambda1.grid = lambda1.grid,
                        lambda2.grid = lambda2.grid,
                        weights = weights,
                        nfolds = nfolds,
                        fold.id = fold.id.cv,
                        loss = loss,
                        selection = selection,
                        support.selection = "rule",
                        nn.index = nn.index,
                        neighborhood.type = "adaptive.radius",
                        support.index = support.index,
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        radius.rule = radius.rule,
                        radius.factor = radius.factor,
                        min.support = support.grid$min.support[ii],
                        max.support = max.support.ii,
                        support.buffer = support.buffer,
                        support.topup = support.topup,
                        tangent.dim.rule = tangent.dim.rule,
                        eigen.tolerance = eigen.tolerance,
                        derivative.order = derivative.order,
                        stabilizer = stabilizer,
                        pinv.tol = pinv.tol,
                        local.solver = local.solver,
                        normal.equations.max.condition = normal.equations.max.condition,
                        ridge = ridge,
                        return.A = return.A,
                        return.local.diagnostics = return.local.diagnostics,
                        return.timing = return.timing,
                        verbose = verbose
                    ),
                    error = function(e) e
                )
            })[["elapsed"]]
            if (inherits(cand, "error")) {
                candidate.rows[[ii]] <- data.frame(
                    candidate = ii,
                    status = "error",
                    adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                    min.support = support.grid$min.support[ii],
                    max.support = support.grid$max.support[ii],
                    cv.mean = Inf,
                    lambda1 = NA_real_,
                    lambda2 = NA_real_,
                    runtime.sec = unname(elapsed),
                    message = conditionMessage(cand),
                    stringsAsFactors = FALSE
                )
            } else {
                candidate.fits[[ii]] <- cand
                diag <- .ssrhe.support.diagnostics(cand$operator)
                candidate.rows[[ii]] <- cbind(
                    data.frame(
                        candidate = ii,
                        status = "ok",
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        min.support = support.grid$min.support[ii],
                        max.support = support.grid$max.support[ii],
                        cv.mean = cand$selection$cv.mean,
                        cv.se = cand$selection$cv.se,
                        lambda1 = cand$selection$lambda1,
                        lambda2 = cand$selection$lambda2,
                        runtime.sec = unname(elapsed),
                        message = "",
                        stringsAsFactors = FALSE
                    ),
                    diag
                )
            }
        }
        support.cv.table <- do.call(rbind, candidate.rows)
        if (!any(is.finite(support.cv.table$cv.mean))) {
            stop("All SSRHE support-CV candidates failed.", call. = FALSE)
        }
        selected.support.idx <- which.min(support.cv.table$cv.mean)
        final <- candidate.fits[[selected.support.idx]]
        final$support.selection <- "cv"
        final$support.cv.table <- support.cv.table
        final$selected.support <- support.grid[selected.support.idx, , drop = FALSE]
        final$selection$support.index <- selected.support.idx
        final$selection$support.cv.mean <- support.cv.table$cv.mean[selected.support.idx]
        attr(final, "call") <- match.call()
        class(final) <- c("ssrhe.hessian.cv.fit", "ssrhe.hessian.fit", "list")
        return(final)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        derivative.order = derivative.order,
        stabilizer = stabilizer,
        pinv.tol = pinv.tol,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        return.A = return.A,
        return.B = TRUE,
        return.BS = stabilizer,
        return.sparse = TRUE,
        return.local.diagnostics = return.local.diagnostics,
        return.timing = return.timing,
        verbose = verbose
    )

    fold.id <- .prepare.ssrhe.cv.folds(
        observed = as.vector(y.info$observed & W > 0),
        nfolds = nfolds,
        fold.id = fold.id
    )
    folds <- sort(unique(fold.id[fold.id > 0]))
    grid <- expand.grid(lambda1 = lambda1.grid,
                        lambda2 = lambda2.grid,
                        KEEP.OUT.ATTRS = FALSE)
    fold.loss <- matrix(NA_real_, nrow = nrow(grid), ncol = length(folds))
    colnames(fold.loss) <- paste0("fold", folds)

    Y <- as.vector(y.info$Y)
    for (g in seq_len(nrow(grid))) {
        for (ff in seq_along(folds)) {
            validation <- fold.id == folds[ff]
            train.weights <- as.vector(W)
            train.weights[validation] <- 0
            fold.fit <- tryCatch(
                .fit.ssrhe.hessian.from.operator(
                    operator = operator,
                    y = Y,
                    lambda1 = grid$lambda1[g],
                    lambda2 = grid$lambda2[g],
                    weights = train.weights,
                    ridge = ridge,
                    verbose = verbose
                ),
                error = function(e) NULL
            )
            if (!is.null(fold.fit)) {
                err <- fold.fit$fitted.values[validation] - Y[validation]
                wv <- as.vector(W)[validation]
                if (identical(loss, "mse")) {
                    fold.loss[g, ff] <- sum(wv * err^2) / sum(wv)
                } else {
                    fold.loss[g, ff] <- sum(wv * abs(err)) / sum(wv)
                }
            }
        }
    }

    cv.mean <- rowMeans(fold.loss, na.rm = TRUE)
    cv.se <- apply(fold.loss, 1L, stats::sd, na.rm = TRUE) /
        sqrt(rowSums(is.finite(fold.loss)))
    cv.mean[!is.finite(cv.mean)] <- Inf
    cv.se[!is.finite(cv.se)] <- Inf
    cv.table <- data.frame(
        lambda1 = grid$lambda1,
        lambda2 = grid$lambda2,
        cv.mean = cv.mean,
        cv.se = cv.se,
        n.successful.folds = rowSums(is.finite(fold.loss)),
        total.penalty = grid$lambda1 + grid$lambda2,
        stringsAsFactors = FALSE
    )
    cv.table <- cbind(cv.table, as.data.frame(fold.loss, check.names = FALSE))
    if (!any(is.finite(cv.table$cv.mean))) {
        stop("All SSRHE cross-validation fits failed.", call. = FALSE)
    }
    min.idx <- which.min(cv.table$cv.mean)
    if (identical(selection, "one.se")) {
        threshold <- cv.table$cv.mean[min.idx] + cv.table$cv.se[min.idx]
        eligible <- which(cv.table$cv.mean <= threshold)
        selected.idx <- eligible[which.max(cv.table$total.penalty[eligible])]
    } else {
        selected.idx <- min.idx
    }

    final <- .fit.ssrhe.hessian.from.operator(
        operator = operator,
        y = Y,
        lambda1 = cv.table$lambda1[selected.idx],
        lambda2 = cv.table$lambda2[selected.idx],
        weights = as.vector(W),
        ridge = ridge,
        verbose = verbose
    )
    final$X <- X
    final$support.selection <- "rule"
    final$selected.support <- data.frame(
        adaptive.k.scale = operator$parameters$adaptive.k.scale,
        min.support = operator$neighborhoods$min.support %||% NA_integer_,
        max.support = operator$neighborhoods$max.support %||% NA_integer_
    )
    final$cv.table <- cv.table
    final$fold.id <- fold.id
    final$selection <- list(
        rule = selection,
        loss = loss,
        selected.index = selected.idx,
        min.index = min.idx,
        lambda1 = cv.table$lambda1[selected.idx],
        lambda2 = cv.table$lambda2[selected.idx],
        cv.mean = cv.table$cv.mean[selected.idx],
        cv.se = cv.table$cv.se[selected.idx],
        nfolds = length(folds)
    )
    attr(final, "call") <- match.call()
    class(final) <- c("ssrhe.hessian.cv.fit", "ssrhe.hessian.fit", "list")
    final
}

# Internal GCV selection engine; public controls are validated and mapped above.
.fit.ssrhe.hessian.regression.gcv <- function(
    X,
    y,
    k = NULL,
    tangent.dim,
    lambda1.grid,
    lambda2.grid = 0,
    weights = NULL,
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    stabilizer = any(lambda2.grid > 0),
    pinv.tol = sqrt(.Machine$double.eps),
    local.solver = c("auto", "normal.equations", "svd", "qr"),
    normal.equations.max.condition = 1e4,
    ridge = 0,
    return.A = TRUE,
    return.local.diagnostics = FALSE,
    return.timing = FALSE,
    support.selection = c("rule", "gcv"),
    support.grid = NULL,
    support.gcv.max.candidates = 8L,
    gcv.trace.method = c("exact", "hutchinson"),
    gcv.trace.n.probes = 50L,
    gcv.trace.seed = NULL,
    verbose = FALSE) {

    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for SSRHE GCV regression.",
             call. = FALSE)
    }
    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("fit.ssrhe.hessian.regression(lambda.selection = 'gcv') currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    if (!all(y.info$observed)) {
        stop("fit.ssrhe.hessian.regression(lambda.selection = 'gcv') currently requires a fully observed response.",
             call. = FALSE)
    }
    if (any(W <= 0)) {
        stop("fit.ssrhe.hessian.regression(lambda.selection = 'gcv') currently requires strictly positive weights.",
             call. = FALSE)
    }
    lambda1.grid <- .validate.ssrhe.lambda.grid(lambda1.grid, "lambda1.grid")
    lambda2.grid <- .validate.ssrhe.lambda.grid(lambda2.grid, "lambda2.grid")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    support.selection <- match.arg(support.selection)
    neighborhood.type <- match.arg(neighborhood.type)
    local.solver <- match.arg(local.solver)
    normal.equations.max.condition <- .validate.ssrhe.numeric.scalar(
        normal.equations.max.condition, "normal.equations.max.condition")
    if (normal.equations.max.condition <= 0) {
        stop("normal.equations.max.condition must be positive.",
             call. = FALSE)
    }
    gcv.trace.method <- match.arg(gcv.trace.method)
    gcv.trace.n.probes <- .validate.ssrhe.positive.integer(
        gcv.trace.n.probes, "gcv.trace.n.probes")
    gcv.trace.seed <- .validate.ssrhe.optional.seed(
        gcv.trace.seed, "gcv.trace.seed")
    support.gcv.max.candidates <- .validate.ssrhe.positive.integer(
        support.gcv.max.candidates, "support.gcv.max.candidates")
    stabilizer <- isTRUE(stabilizer)
    if (any(lambda2.grid > 0) && !stabilizer) {
        stop("positive lambda2.grid values require stabilizer = TRUE.",
             call. = FALSE)
    }
    if (derivative.order == 3L && (any(lambda2.grid > 0) || stabilizer)) {
        stop("lambda2.grid/stabilizer is currently only supported for derivative.order = 2.",
             call. = FALSE)
    }

    if (identical(support.selection, "gcv")) {
        if (!identical(neighborhood.type, "adaptive.radius")) {
            stop("support.selection = 'gcv' is currently supported only for neighborhood.type = 'adaptive.radius'.",
                 call. = FALSE)
        }
        tangent.dim.grid <- if (missing(tangent.dim) || is.null(tangent.dim)) {
            ncol(X)
        } else {
            .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        }
        support.grid <- .validate.ssrhe.support.grid(
            support.grid = support.grid,
            n = n,
            tangent.dim = tangent.dim.grid,
            derivative.order = derivative.order,
            support.buffer = support.buffer,
            max.candidates = support.gcv.max.candidates
        )
        candidate.fits <- vector("list", nrow(support.grid))
        candidate.rows <- vector("list", nrow(support.grid))
        for (ii in seq_len(nrow(support.grid))) {
            max.support.ii <- support.grid$max.support[ii]
            if (is.na(max.support.ii)) max.support.ii <- NULL
            elapsed <- system.time({
                cand <- tryCatch(
                    .fit.ssrhe.hessian.regression.gcv(
                        X = X,
                        y = y,
                        k = k,
                        tangent.dim = if (missing(tangent.dim)) NULL else tangent.dim,
                        lambda1.grid = lambda1.grid,
                        lambda2.grid = lambda2.grid,
                        weights = weights,
                        support.selection = "rule",
                        nn.index = nn.index,
                        neighborhood.type = "adaptive.radius",
                        support.index = support.index,
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        radius.rule = radius.rule,
                        radius.factor = radius.factor,
                        min.support = support.grid$min.support[ii],
                        max.support = max.support.ii,
                        support.buffer = support.buffer,
                        support.topup = support.topup,
                        tangent.dim.rule = tangent.dim.rule,
                        eigen.tolerance = eigen.tolerance,
                        derivative.order = derivative.order,
                        stabilizer = stabilizer,
                        pinv.tol = pinv.tol,
                        local.solver = local.solver,
                        normal.equations.max.condition = normal.equations.max.condition,
                        ridge = ridge,
                        return.A = return.A,
                        return.local.diagnostics = return.local.diagnostics,
                        return.timing = return.timing,
                        gcv.trace.method = gcv.trace.method,
                        gcv.trace.n.probes = gcv.trace.n.probes,
                        gcv.trace.seed = gcv.trace.seed,
                        verbose = verbose
                    ),
                    error = function(e) e
                )
            })[["elapsed"]]
            if (inherits(cand, "error")) {
                candidate.rows[[ii]] <- data.frame(
                    candidate = ii,
                    status = "error",
                    adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                    min.support = support.grid$min.support[ii],
                    max.support = support.grid$max.support[ii],
                    gcv = Inf,
                    rss = NA_real_,
                    trace.S = NA_real_,
                    trace.se = NA_real_,
                    trace.method = gcv.trace.method,
                    trace.n.probes = if (identical(gcv.trace.method, "hutchinson")) {
                        gcv.trace.n.probes
                    } else {
                        NA_integer_
                    },
                    lambda1 = NA_real_,
                    lambda2 = NA_real_,
                    runtime.sec = unname(elapsed),
                    message = conditionMessage(cand),
                    stringsAsFactors = FALSE
                )
            } else {
                candidate.fits[[ii]] <- cand
                diag <- .ssrhe.support.diagnostics(cand$operator)
                candidate.rows[[ii]] <- cbind(
                    data.frame(
                        candidate = ii,
                        status = "ok",
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        min.support = support.grid$min.support[ii],
                        max.support = support.grid$max.support[ii],
                        gcv = cand$selection$gcv,
                        rss = cand$selection$rss,
                        trace.S = cand$selection$trace.S,
                        trace.se = cand$selection$trace.se,
                        trace.method = cand$selection$trace.method,
                        trace.n.probes = cand$selection$trace.n.probes,
                        lambda1 = cand$selection$lambda1,
                        lambda2 = cand$selection$lambda2,
                        runtime.sec = unname(elapsed),
                        message = "",
                        stringsAsFactors = FALSE
                    ),
                    diag
                )
            }
        }
        support.gcv.table <- do.call(rbind, candidate.rows)
        if (!any(is.finite(support.gcv.table$gcv))) {
            stop("All SSRHE support-GCV candidates failed.", call. = FALSE)
        }
        selected.support.idx <- which.min(support.gcv.table$gcv)
        final <- candidate.fits[[selected.support.idx]]
        final$support.selection <- "gcv"
        final$support.gcv.table <- support.gcv.table
        final$selected.support <- support.grid[selected.support.idx, , drop = FALSE]
        final$selection$support.index <- selected.support.idx
        final$selection$support.gcv <- support.gcv.table$gcv[selected.support.idx]
        attr(final, "call") <- match.call()
        class(final) <- c("ssrhe.hessian.gcv.fit", "ssrhe.hessian.fit", "list")
        return(final)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        derivative.order = derivative.order,
        stabilizer = stabilizer,
        pinv.tol = pinv.tol,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        return.A = return.A,
        return.B = TRUE,
        return.BS = stabilizer,
        return.sparse = TRUE,
        return.local.diagnostics = return.local.diagnostics,
        return.timing = return.timing,
        verbose = verbose
    )

    final <- .fit.ssrhe.hessian.gcv.from.operator(
        operator = operator,
        y = as.vector(y.info$Y),
        lambda1.grid = lambda1.grid,
        lambda2.grid = lambda2.grid,
        weights = as.vector(W),
        ridge = ridge,
        trace.method = gcv.trace.method,
        trace.n.probes = gcv.trace.n.probes,
        trace.seed = gcv.trace.seed,
        verbose = verbose
    )
    final$gcv.trace <- list(
        method = gcv.trace.method,
        n.probes = if (identical(gcv.trace.method, "hutchinson")) {
            gcv.trace.n.probes
        } else {
            NA_integer_
        },
        seed = gcv.trace.seed
    )
    final$X <- X
    final$support.selection <- "rule"
    final$selected.support <- data.frame(
        adaptive.k.scale = operator$parameters$adaptive.k.scale,
        min.support = operator$neighborhoods$min.support %||% NA_integer_,
        max.support = operator$neighborhoods$max.support %||% NA_integer_
    )
    attr(final, "call") <- match.call()
    class(final) <- c("ssrhe.hessian.gcv.fit", "ssrhe.hessian.fit", "list")
    final
}

.fit.ssrhe.hessian.gcv.from.operator <- function(operator,
                                                 y,
                                          lambda1.grid,
                                          lambda2.grid,
                                                 weights,
                                                 ridge,
                                                 trace.method = c("exact", "hutchinson"),
                                                 trace.n.probes = 50L,
                                                 trace.seed = NULL,
                                                 verbose = FALSE) {
    n <- ncol(operator$B)
    trace.method <- match.arg(trace.method)
    trace.n.probes <- .validate.ssrhe.positive.integer(
        trace.n.probes, "trace.n.probes")
    trace.seed <- .validate.ssrhe.optional.seed(trace.seed, "trace.seed")
    if (length(y) != n) stop("Internal error: y length mismatch.", call. = FALSE)
    if (length(weights) != n) {
        stop("Internal error: weights length mismatch.", call. = FALSE)
    }
    if (any(!is.finite(y))) {
        stop("GCV requires a finite response vector.", call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights <= 0)) {
        stop("GCV requires finite strictly positive weights.", call. = FALSE)
    }
    grid <- expand.grid(lambda1 = lambda1.grid,
                        lambda2 = lambda2.grid,
                        KEEP.OUT.ATTRS = FALSE)
    rows <- vector("list", nrow(grid))
    fits <- vector("list", nrow(grid))
    for (g in seq_len(nrow(grid))) {
        fit <- tryCatch(
            .fit.ssrhe.hessian.from.operator(
                operator = operator,
                y = y,
                lambda1 = grid$lambda1[g],
                lambda2 = grid$lambda2[g],
                weights = weights,
                ridge = ridge,
                verbose = verbose
            ),
            error = function(e) e
        )
        if (inherits(fit, "error")) {
            rows[[g]] <- data.frame(
                lambda1 = grid$lambda1[g],
                lambda2 = grid$lambda2[g],
                rss = NA_real_,
                trace.S = NA_real_,
                trace.se = NA_real_,
                trace.method = trace.method,
                trace.n.probes = if (identical(trace.method, "hutchinson")) {
                    trace.n.probes
                } else {
                    NA_integer_
                },
                edf = NA_real_,
                gcv = Inf,
                status = "error",
                message = conditionMessage(fit),
                stringsAsFactors = FALSE
            )
            next
        }
        system <- .ssrhe.hessian.system.matrix(
            B = operator$B,
            BS = operator$BS,
            weights = weights,
            lambda1 = grid$lambda1[g],
            lambda2 = grid$lambda2[g],
            ridge = ridge
        )
        grid.seed <- .ssrhe.trace.grid.seed(trace.seed, g)
        trace.info <- tryCatch(
            .ssrhe.hessian.smoother.trace(
                system = system,
                weights = weights,
                method = trace.method,
                n.probes = trace.n.probes,
                seed = grid.seed,
                verbose = verbose
            ),
            error = function(e) e
        )
        if (inherits(trace.info, "error")) {
            rows[[g]] <- data.frame(
                lambda1 = grid$lambda1[g],
                lambda2 = grid$lambda2[g],
                rss = NA_real_,
                trace.S = NA_real_,
                trace.se = NA_real_,
                trace.method = trace.method,
                trace.n.probes = if (identical(trace.method, "hutchinson")) {
                    trace.n.probes
                } else {
                    NA_integer_
                },
                edf = NA_real_,
                gcv = Inf,
                status = "error",
                message = conditionMessage(trace.info),
                stringsAsFactors = FALSE
            )
            next
        }
        trace.S <- trace.info$trace
        rss <- sum(weights * (y - fit$fitted.values)^2)
        denom <- (1 - trace.S / n)^2
        gcv <- if (denom > 0) (rss / n) / denom else Inf
        rows[[g]] <- data.frame(
            lambda1 = grid$lambda1[g],
            lambda2 = grid$lambda2[g],
            rss = rss,
            trace.S = trace.S,
            trace.se = trace.info$trace.se,
            trace.method = trace.info$method,
            trace.n.probes = trace.info$n.probes,
            edf = trace.S,
            gcv = gcv,
            status = "ok",
            message = "",
            stringsAsFactors = FALSE
        )
        fits[[g]] <- fit
    }
    gcv.table <- do.call(rbind, rows)
    if (!any(is.finite(gcv.table$gcv))) {
        stop("All SSRHE GCV fits failed.", call. = FALSE)
    }
    selected.idx <- which.min(gcv.table$gcv)
    final <- fits[[selected.idx]]
    final$gcv.table <- gcv.table
    final$selection <- list(
        rule = "gcv",
        selected.index = selected.idx,
        lambda1 = gcv.table$lambda1[selected.idx],
        lambda2 = gcv.table$lambda2[selected.idx],
        gcv = gcv.table$gcv[selected.idx],
        rss = gcv.table$rss[selected.idx],
        trace.S = gcv.table$trace.S[selected.idx],
        trace.se = gcv.table$trace.se[selected.idx],
        trace.method = gcv.table$trace.method[selected.idx],
        trace.n.probes = gcv.table$trace.n.probes[selected.idx],
        n = n
    )
    final
}

.fit.ssrhe.hessian.from.operator <- function(operator,
                                             y,
                                      lambda1,
                                      lambda2,
                                             weights,
                                             ridge,
                                             verbose = FALSE) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for SSRHE regression solves.",
             call. = FALSE)
    }
    if (is.null(operator$B)) {
        stop("operator must contain B. Rebuild with return.B = TRUE.",
             call. = FALSE)
    }
    n <- ncol(operator$B)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    Y.raw <- y.info$Y
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    Y.clean <- Y.raw
    Y.clean[!is.finite(Y.clean)] <- 0
    lambda1 <- .validate.ssrhe.nonnegative.scalar(lambda1, "lambda1")
    lambda2 <- .validate.ssrhe.nonnegative.scalar(lambda2, "lambda2")
    ridge <- .validate.ssrhe.nonnegative.scalar(ridge, "ridge")
    if (lambda2 > 0 && is.null(operator$BS)) {
        stop("lambda2 > 0 requires an operator with BS. Rebuild with stabilizer = TRUE.",
             call. = FALSE)
    }

    n.responses <- ncol(Y.clean)
    Y.hat <- matrix(NA_real_, nrow = n, ncol = n.responses)
    residuals <- matrix(NA_real_, nrow = n, ncol = n.responses)
    data.loss <- hessian.energy <- stabilizer.energy <- objective <- numeric(n.responses)
    solve.method <- character(n.responses)

    common.weights <- n.responses == 1L || all(vapply(seq_len(n.responses), function(j) {
        identical(W[, j], W[, 1L])
    }, logical(1)))

    if (common.weights) {
        solved <- .solve.ssrhe.hessian.system(
            B = operator$B,
            BS = operator$BS,
            weights = W[, 1L],
            Y = Y.clean,
            lambda1 = lambda1,
            lambda2 = lambda2,
            ridge = ridge,
            verbose = verbose
        )
        Y.hat[,] <- solved$fitted
        solve.method[] <- solved$method
    } else {
        for (j in seq_len(n.responses)) {
            solved <- .solve.ssrhe.hessian.system(
                B = operator$B,
                BS = operator$BS,
                weights = W[, j],
                Y = Y.clean[, j, drop = FALSE],
                lambda1 = lambda1,
                lambda2 = lambda2,
                ridge = ridge,
                verbose = verbose
            )
            Y.hat[, j] <- solved$fitted[, 1L]
            solve.method[j] <- solved$method
        }
    }

    for (j in seq_len(n.responses)) {
        observed <- y.info$observed[, j] & W[, j] > 0
        residuals[observed, j] <- Y.raw[observed, j] - Y.hat[observed, j]
        residuals[!observed, j] <- NA_real_
        data.loss[j] <- 0.5 * sum(W[, j] * (Y.clean[, j] - Y.hat[, j])^2)
        hessian.energy[j] <- as.numeric(Matrix::crossprod(Y.hat[, j], operator$B %*% Y.hat[, j]))
        if (!is.null(operator$BS)) {
            stabilizer.energy[j] <- as.numeric(Matrix::crossprod(Y.hat[, j], operator$BS %*% Y.hat[, j]))
        } else {
            stabilizer.energy[j] <- 0
        }
        objective[j] <- data.loss[j] +
            0.5 * lambda1 * hessian.energy[j] +
            0.5 * lambda2 * stabilizer.energy[j]
    }

    if (!is.null(y.info$col.names)) {
        colnames(Y.hat) <- y.info$col.names
        colnames(residuals) <- y.info$col.names
        names(data.loss) <- y.info$col.names
        names(hessian.energy) <- y.info$col.names
        names(stabilizer.energy) <- y.info$col.names
        names(objective) <- y.info$col.names
    }

    single <- n.responses == 1L
    list(
        fitted.values = if (single) as.vector(Y.hat) else Y.hat,
        residuals = if (single) as.vector(residuals) else residuals,
        y = if (single) as.vector(Y.raw) else Y.raw,
        weights = if (single) as.vector(W) else W,
        lambda = list(lambda1 = lambda1, lambda2 = lambda2, ridge = ridge),
        objective = if (single) objective[1L] else objective,
        energies = list(
            data.loss = if (single) data.loss[1L] else data.loss,
            hessian = if (single) hessian.energy[1L] else hessian.energy,
            stabilizer = if (single) stabilizer.energy[1L] else stabilizer.energy
        ),
        operator = operator,
        n.responses = n.responses,
        solver = list(method = if (length(unique(solve.method)) == 1L) solve.method[1L] else solve.method),
        observed = if (single) as.vector(y.info$observed & W > 0) else y.info$observed & W > 0
    )
}

.solve.ssrhe.hessian.system <- function(B, BS, weights, Y,
                                        lambda1, lambda2, ridge,
                                        verbose = FALSE) {
    n <- nrow(B)
    if (length(weights) != n) stop("Internal error: weights length mismatch.", call. = FALSE)
    system <- .ssrhe.hessian.system.matrix(
        B = B,
        BS = BS,
        weights = weights,
        lambda1 = lambda1,
        lambda2 = lambda2,
        ridge = ridge
    )
    rhs <- weights * Y

    fit <- tryCatch({
        as.matrix(Matrix::solve(system, rhs))
    }, error = function(e) {
        if (isTRUE(verbose)) {
            message("Sparse solve failed; retrying with dense solve().")
        }
        tryCatch(
            solve(as.matrix(system), as.matrix(rhs)),
            error = function(e2) {
                stop("SSRHE linear solve failed: ", conditionMessage(e2),
                     call. = FALSE)
            }
        )
    })
    list(fitted = fit, method = "linear.solve")
}

.ssrhe.hessian.system.matrix <- function(B, BS, weights,
                                         lambda1, lambda2, ridge) {
    n <- nrow(B)
    system <- lambda1 * B
    if (lambda2 > 0) {
        if (is.null(BS)) {
            stop("lambda2 > 0 requires BS.", call. = FALSE)
        }
        system <- system + lambda2 * BS
    }
    if (ridge > 0) {
        system <- system + Matrix::Diagonal(n, ridge)
    }
    if (any(weights > 0)) {
        system <- system + Matrix::Diagonal(n, weights)
    }
    system
}

.ssrhe.hessian.smoother.trace <- function(system,
                                          weights,
                                          method = c("exact", "hutchinson"),
                                          n.probes = 50L,
                                          seed = NULL,
                                          verbose = FALSE) {
    method <- match.arg(method)
    if (identical(method, "exact")) {
        trace <- .ssrhe.hessian.smoother.trace.exact(system, weights,
                                                     verbose = verbose)
        return(list(
            trace = trace,
            trace.se = NA_real_,
            method = "exact",
            n.probes = NA_integer_
        ))
    }
    .ssrhe.hessian.smoother.trace.hutchinson(
        system = system,
        weights = weights,
        n.probes = n.probes,
        seed = seed,
        verbose = verbose
    )
}

.ssrhe.hessian.smoother.trace.exact <- function(system,
                                                weights,
                                                verbose = FALSE) {
    rhs <- Matrix::Diagonal(length(weights), x = weights)
    S <- tryCatch({
        Matrix::solve(system, rhs)
    }, error = function(e) {
        if (isTRUE(verbose)) {
            message("Sparse smoother-trace solve failed; retrying with dense solve().")
        }
        tryCatch(
            solve(as.matrix(system), as.matrix(rhs)),
            error = function(e2) {
                stop("SSRHE smoother-trace solve failed: ", conditionMessage(e2),
                     call. = FALSE)
            }
        )
    })
    if (inherits(S, "Matrix")) {
        sum(as.numeric(Matrix::diag(S)))
    } else {
        sum(diag(S))
    }
}

.ssrhe.hessian.smoother.trace.hutchinson <- function(system,
                                                     weights,
                                                     n.probes,
                                                     seed = NULL,
                                                     verbose = FALSE) {
    n.probes <- .validate.ssrhe.positive.integer(n.probes, "n.probes")
    seed <- .validate.ssrhe.optional.seed(seed, "seed")
    n <- length(weights)
    probes <- .ssrhe.with.optional.seed(seed, {
        matrix(
            sample(c(-1, 1), n * n.probes, replace = TRUE),
            nrow = n,
            ncol = n.probes
        )
    })
    rhs <- weights * probes
    solved <- tryCatch({
        as.matrix(Matrix::solve(system, rhs))
    }, error = function(e) {
        if (isTRUE(verbose)) {
            message("Sparse Hutchinson trace solve failed; retrying with dense solve().")
        }
        tryCatch(
            solve(as.matrix(system), as.matrix(rhs)),
            error = function(e2) {
                stop("SSRHE Hutchinson trace solve failed: ",
                     conditionMessage(e2), call. = FALSE)
            }
        )
    })
    trace.samples <- colSums(probes * solved)
    list(
        trace = mean(trace.samples),
        trace.se = if (n.probes > 1L) {
            stats::sd(trace.samples) / sqrt(n.probes)
        } else {
            NA_real_
        },
        method = "hutchinson",
        n.probes = n.probes
    )
}

.ssrhe.with.optional.seed <- function(seed, expr) {
    if (is.null(seed)) {
        return(force(expr))
    }
    had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had.seed) {
        old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
        if (had.seed) {
            assign(".Random.seed", old.seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    force(expr)
}

.ssrhe.trace.grid.seed <- function(seed, grid.index) {
    if (is.null(seed)) return(NULL)
    as.integer(((seed + grid.index - 1L) %% 2147483646L) + 1L)
}

.prepare.ssrhe.response.matrix <- function(y, n, name) {
    is.matrix.input <- is.matrix(y) || inherits(y, "Matrix") ||
        (is.data.frame(y) && ncol(y) > 1L)
    if (is.matrix.input) {
        Y <- if (inherits(y, "Matrix")) as.matrix(y) else as.matrix(y)
        if (nrow(Y) != n) stop(sprintf("nrow(%s) must be %d.", name, n), call. = FALSE)
        col.names <- colnames(Y)
    } else {
        if (length(y) != n) stop(sprintf("%s must have length %d.", name, n), call. = FALSE)
        Y <- matrix(y, ncol = 1L)
        col.names <- NULL
    }
    storage.mode(Y) <- "double"
    if (any(is.infinite(Y))) {
        stop(sprintf("%s cannot contain infinite values.", name), call. = FALSE)
    }
    observed <- is.finite(Y)
    list(Y = Y, observed = observed, col.names = col.names)
}

.prepare.ssrhe.weight.matrix <- function(weights, y.info, n) {
    p <- ncol(y.info$Y)
    if (is.null(weights)) {
        W <- matrix(1, nrow = n, ncol = p)
    } else if (is.matrix(weights) || inherits(weights, "Matrix") || is.data.frame(weights)) {
        W <- if (inherits(weights, "Matrix")) as.matrix(weights) else as.matrix(weights)
        if (!identical(dim(W), dim(y.info$Y))) {
            stop("weights matrix must have the same dimensions as y.", call. = FALSE)
        }
    } else {
        if (length(weights) != n) stop("weights vector must have length nrow(X).", call. = FALSE)
        W <- matrix(as.double(weights), nrow = n, ncol = p)
    }
    storage.mode(W) <- "double"
    if (any(!is.finite(W)) || any(W < 0)) {
        stop("weights must be finite and nonnegative.", call. = FALSE)
    }
    W[!y.info$observed] <- 0
    W
}

.validate.ssrhe.nonnegative.scalar <- function(x, name) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
        stop(sprintf("%s must be a finite nonnegative numeric scalar.", name),
             call. = FALSE)
    }
    as.double(x)
}

.validate.ssrhe.optional.seed <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 ||
        abs(x - round(x)) > .Machine$double.eps^0.5) {
        stop(sprintf("%s must be NULL or a finite nonnegative integer scalar.", name),
             call. = FALSE)
    }
    as.integer(round(x))
}

.validate.ssrhe.lambda.grid <- function(x, name) {
    if (!is.numeric(x) || !length(x) ||
        any(is.na(x)) || any(!is.finite(x)) || any(x < 0)) {
        stop(sprintf("%s must be a nonempty finite nonnegative numeric vector.",
                     name), call. = FALSE)
    }
    sort(unique(as.double(x)))
}

.prepare.ssrhe.cv.folds <- function(observed, nfolds, fold.id) {
    n <- length(observed)
    if (is.null(fold.id)) {
        nfolds <- .validate.ssrhe.positive.integer(nfolds, "nfolds")
        observed.idx <- which(observed)
        if (length(observed.idx) < 2L) {
            stop("At least two observed positive-weight labels are required for CV.",
                 call. = FALSE)
        }
        nfolds <- min(nfolds, length(observed.idx))
        out <- integer(n)
        out[observed.idx] <- rep(seq_len(nfolds), length.out = length(observed.idx))
    } else {
        if (length(fold.id) != n) {
            stop("fold.id must have length nrow(X).", call. = FALSE)
        }
        out <- as.integer(fold.id)
        out[is.na(out) | out < 1L | !observed] <- 0L
    }
    folds <- sort(unique(out[out > 0L]))
    if (length(folds) < 2L) {
        stop("At least two nonempty validation folds are required.", call. = FALSE)
    }
    for (ff in folds) {
        n.validation <- sum(out == ff)
        n.training <- sum(observed & out != ff)
        if (n.validation < 1L || n.training < 1L) {
            stop("Each validation fold must leave at least one training label.",
                 call. = FALSE)
        }
    }
    out
}

#' Fit SSRHE-Style Hessian L1 Regression
#'
#' Fits an \eqn{\ell_1}-adaptive SSRHE Hessian comparator using the row
#' operator \eqn{A} from \code{\link{ssrhe.hessian.operator}}:
#' \deqn{
#'   \widehat f_\lambda =
#'   \arg\min_f
#'   \left\{
#'     \frac{1}{2}\sum_i w_i(y_i-f_i)^2
#'     + \lambda \|Af\|_1
#'   \right\}.
#' }
#'
#' @inheritParams ssrhe.hessian.operator
#' @param y Numeric response vector of length \code{nrow(X)}. Matrix responses
#'   are not yet supported for the \eqn{\ell_1} path.
#' @param lambda.grid Optional nonnegative lambda grid. For
#'   \code{lambda.selection = "fixed"}, this must contain exactly one value. For
#'   \code{lambda.selection = "cv"}, \code{NULL} builds a default grid from the
#'   fitted full-data generalized-lasso path.
#' @param lambda.selection \code{"cv"} for observed-label K-fold
#'   cross-validation or \code{"fixed"} for a supplied fixed lambda.
#' @param weights Optional nonnegative observation weights. Missing responses
#'   automatically receive zero weight.
#' @param n.lambda Number of default lambda candidates when
#'   \code{lambda.grid = NULL} and \code{lambda.selection = "cv"}.
#' @param nfolds Number of validation folds over observed positive-weight
#'   labels. Ignored when \code{fold.id} is supplied.
#' @param fold.id Optional integer fold assignments of length \code{nrow(X)}.
#' @param loss Validation loss, currently \code{"mse"} or \code{"mae"}.
#' @param selection Selection rule. \code{"min"} chooses the smallest mean
#'   validation loss. \code{"one.se"} chooses the largest lambda within one
#'   standard error of the minimum.
#' @param support.selection Support-profile selection rule.
#'   \code{"rule"} uses the supplied \code{adaptive.k.scale},
#'   \code{min.support}, and \code{max.support}. \code{"cv"} is currently
#'   supported for \code{neighborhood.type = "adaptive.radius"} with
#'   \code{lambda.selection = "cv"} and chooses among rows of
#'   \code{support.grid} by an outer response cross-validation loop.
#' @param support.grid Optional data frame of support profiles with columns
#'   \code{adaptive.k.scale}, \code{min.support}, and optional
#'   \code{max.support}. If \code{NULL}, \code{\link{ssrhe.support.grid}} builds
#'   a compact default grid.
#' @param support.cv.max.candidates Maximum number of support profiles to try
#'   when \code{support.selection = "cv"}.
#' @param solver Solver backend. \code{"genlasso"} uses the generalized-lasso
#'   path to propose a grid when needed, then refines every candidate with ADMM
#'   on the original objective. \code{"admm"} skips the path and requires an
#'   explicit grid. \code{"auto"} also proceeds when path construction fails,
#'   provided an explicit grid is supplied.
#' @param row.scaling Optional row scaling for the penalty matrix before
#'   solving. \code{"l2"} scales nonzero rows of \eqn{A} to unit Euclidean norm.
#' @param admm.rho,admm.maxiter,admm.abstol,admm.reltol ADMM controls used by
#'   all backends, including refinement of generalized-lasso candidates.
#'   \code{admm.rho = NULL} chooses an initial algorithm penalty of
#'   \code{1 / max(1, mean(rowSums(D^2)))}, where \code{D} is the penalty
#'   operator after any requested row scaling. A positive numeric value sets
#'   the initial penalty explicitly. This is separate from statistical
#'   regularization \code{lambda.grid}. The iteration limit is a hard cap.
#' @param admm.adaptive.rho Logical; adapt the algorithm penalty during the
#'   first 1000 iterations by balancing tolerance-normalized primal and dual
#'   residuals. Set \code{FALSE} to hold the initial penalty fixed.
#' @param maxsteps,minlam,approx,rtol,btol,eps,verbose Controls passed to
#'   \pkg{genlasso}.
#'
#' @details
#' This function exposes the SSRHE local polynomial derivative estimator as an
#' \eqn{\ell_1} penalty, rather than as the original SSRHE \eqn{\ell_2}
#' Hessian-energy penalty used by
#' \code{\link{fit.ssrhe.hessian.regression}}. With
#' \code{derivative.order = 2L}, \eqn{A} contains local quadratic/Hessian rows.
#' With \code{derivative.order = 3L}, \eqn{A} contains local cubic
#' third-derivative tensor rows with the same symmetric component scaling used
#' by \code{\link{ssrhe.hessian.operator}}. The existing \eqn{\ell_2} fit
#' solves a linear system involving \eqn{B=A^\top A}. This function instead
#' keeps the rows of \eqn{A} and penalizes their absolute values. It is
#' therefore closer in spirit to graph trend filtering, while retaining SSRHE's
#' local-PCA derivative construction.
#'
#' The default solver backend is \pkg{genlasso}. For larger or numerically
#' fragile third-order operators, \code{solver = "admm"} fits the requested
#' lambda values directly without computing a generalized-lasso path, and
#' all backends use ADMM for the returned coefficients. \code{solver = "auto"}
#' can additionally recover from path construction errors with an explicit grid.
#' Cross-validation removes held-out labels
#' from the data-fit term by setting their weights to zero, then scores
#' predictions on those held-out labels.
#'
#' ADMM uses proximal numerical stabilization, which does not add a ridge
#' penalty to the requested objective. Convergence requires primal feasibility,
#' the dual stopping residual, and stationarity of the original weighted L1
#' objective to meet their reported tolerances. A finite fixed-lambda iterate
#' that reaches its cap is returned with a warning and
#' \code{solver$converged = FALSE}; it is not a completed fit. Inspect
#' \code{solver$admm} before using it, and consider a larger iteration limit
#' or an adaptive penalty. Non-finite selected fits are errors.
#'
#' CV candidates are eligible only when every fold has an acceptable solve.
#' ADMM iterates that do not converge are not scored. The CV result retains
#' \code{fold.status}, \code{fold.converged}, \code{fold.diagnostics},
#' \code{n.successful.folds}, and \code{eligible}. If no candidate qualifies,
#' an error of class \code{ssrhe_l1_cv_error} retains fold diagnostics for
#' inspection with \code{tryCatch}. A selected full-data ADMM fit that fails
#' to converge is also an error in CV mode.
#'
#' A generalized-lasso path alone can be inaccurate for nearly rank-deficient
#' operators. It is retained in \code{path} for diagnostics, but its coefficients
#' are not the returned fit. \code{solver$backend = "genlasso_admm_refined"}
#' identifies this route; \code{solver$path.objective} records the unrefined
#' objectives at the grid values. Backend warnings, including its numerical
#' ridge for underdetermined designs, are retained in \code{solver$warnings}.
#' CV uses the original objective without that additional ridge. The default
#' ADMM tolerances are absolute \code{1e-6} and relative \code{1e-5}; convergence
#' describes numerical stopping at those tolerances, not predictive accuracy.
#' With \code{row.scaling = "l2"}, the reported objective uses the scaled
#' penalty; \code{energies$penalty.l1} is that penalty norm, while
#' \code{energies$hessian.l1} retains the unscaled operator norm.
#'
#' With \code{support.selection = "cv"}, each adaptive-radius support candidate
#' builds a fresh SSRHE operator and runs the usual lambda CV with shared fold
#' assignments. The selected fit is the candidate with the smallest selected CV
#' error. This is useful when the default support-size rule is too rigid, but it
#' can be much slower than tuning lambda for one fixed operator.
#'
#' @return A list of class \code{"ssrhe.hessian.l1.fit"} containing fitted
#'   values, residuals, selected lambda, the SSRHE operator, generalized-lasso
#'   path metadata, lambda-grid fitted values, and CV diagnostics when
#'   requested.
#'
#' @references
#' Kim, K. I., Steinke, F., and Hein, M. (2009). Semi-supervised regression
#' using Hessian energy with an application to semi-supervised dimensionality
#' reduction. \emph{Advances in Neural Information Processing Systems 22}.
#' \url{https://papers.nips.cc/paper_files/paper/2009/hash/f4552671f8909587cf485ea990207f3b-Abstract.html}
#'
#' @examples
#' X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
#' fit.ssrhe.hessian.l1.regression(
#'   X, sin(2 * pi * X[, 1]), k = 6L, tangent.dim = 1L,
#'   lambda.grid = 0.05, lambda.selection = "fixed", solver = "admm"
#' )
#' @section Choosing controls:
#' **First fit:** specify `X`, `y`, `k`, `tangent.dim`, and `lambda.grid`.
#' For a single supplied penalty set `lambda.selection = "fixed"` and
#' `solver = "admm"` to avoid constructing a path.
#'
#' **Selection:** `fold.id`, `nfolds`, `loss`, and `selection` control label CV.
#' A grid requires one solve per penalty per fold; support CV also rebuilds
#' geometry. Incomplete candidates are excluded, not silently scored.
#'
#' **Geometry:** neighborhood, tangent and derivative controls define `A`.
#' `row.scaling` changes the penalty itself and thus the statistical objective.
#'
#' **Numerics:** `admm.*` controls stopping and algorithmic penalties for all
#' backends. Path controls (`maxsteps`, `rtol`, `btol`, `eps`) affect the optional
#' proposal path. They do not certify the final original-objective solve.
#'
#' **Diagnostics:** `solver$admm`, candidate/fold status and objective components
#' describe numerical acceptance. `path` and `beta.grid` can use substantial
#' memory; use a compact explicit grid with ADMM when no full path is needed.
#'
#' @export
fit.ssrhe.hessian.l1.regression <- function(
    X,
    y,
    k = NULL,
    tangent.dim,
    lambda.grid = NULL,
    lambda.selection = c("cv", "fixed"),
    weights = NULL,
    n.lambda = 40L,
    nfolds = 5L,
    fold.id = NULL,
    loss = c("mse", "mae"),
    selection = c("min", "one.se"),
    nn.index = NULL,
    neighborhood.type = c("knn", "adaptive.radius", "supplied"),
    support.index = NULL,
    adaptive.k.scale = NULL,
    radius.rule = c("geomean", "max", "min"),
    radius.factor = 1.25,
    min.support = NULL,
    max.support = NULL,
    support.buffer = 2L,
    support.topup = c("nearest", "none"),
    tangent.dim.rule = c("fixed", "eigen.cumulative"),
    eigen.tolerance = 0.95,
    derivative.order = 2L,
    pinv.tol = sqrt(.Machine$double.eps),
    local.solver = c("auto", "normal.equations", "svd", "qr"),
    normal.equations.max.condition = 1e4,
    solver = c("genlasso", "admm", "auto"),
    row.scaling = c("none", "l2"),
    admm.rho = NULL,
    admm.maxiter = 2000L,
    admm.abstol = 1e-6,
    admm.reltol = 1e-5,
    maxsteps = 2000L,
    minlam = 0,
    approx = FALSE,
    rtol = 1e-7,
    btol = 1e-7,
    eps = 1e-4,
    support.selection = c("rule", "cv"),
    support.grid = NULL,
    support.cv.max.candidates = 8L,
    return.local.diagnostics = FALSE,
    return.timing = FALSE,
    verbose = FALSE,
    admm.adaptive.rho = TRUE) {

    solver <- match.arg(solver)
    if (!identical(solver, "admm") &&
        !requireNamespace("genlasso", quietly = TRUE)) {
        stop("Package 'genlasso' is required for fit.ssrhe.hessian.l1.regression().",
             call. = FALSE)
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for fit.ssrhe.hessian.l1.regression().",
             call. = FALSE)
    }
    X <- .validate.ssrhe.X(X)
    n <- nrow(X)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("fit.ssrhe.hessian.l1.regression currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    lambda.selection <- match.arg(lambda.selection)
    loss <- match.arg(loss)
    selection <- match.arg(selection)
    support.selection <- match.arg(support.selection)
    neighborhood.type <- match.arg(neighborhood.type)
    local.solver <- match.arg(local.solver)
    normal.equations.max.condition <- .validate.ssrhe.numeric.scalar(
        normal.equations.max.condition, "normal.equations.max.condition")
    if (normal.equations.max.condition <= 0) {
        stop("normal.equations.max.condition must be positive.",
             call. = FALSE)
    }
    derivative.order <- .validate.ssrhe.derivative.order(derivative.order)
    support.cv.max.candidates <- .validate.ssrhe.positive.integer(
        support.cv.max.candidates, "support.cv.max.candidates")
    solver.args <- .validate.ssrhe.hessian.l1.solver.args(
        n.lambda = n.lambda,
        nfolds = nfolds,
        fold.id = fold.id,
        observed = as.vector(y.info$observed & W > 0),
        maxsteps = maxsteps,
        minlam = minlam,
        approx = approx,
        rtol = rtol,
        btol = btol,
        eps = eps,
        solver = solver,
        row.scaling = row.scaling,
        admm.rho = admm.rho,
        admm.maxiter = admm.maxiter,
        admm.abstol = admm.abstol,
        admm.reltol = admm.reltol,
        admm.adaptive.rho = admm.adaptive.rho,
        verbose = verbose
    )
    lambda.grid <- .validate.ssrhe.hessian.l1.lambda.grid(
        lambda.grid,
        lambda.selection
    )
    if (identical(support.selection, "cv")) {
        if (!identical(neighborhood.type, "adaptive.radius")) {
            stop("support.selection = 'cv' is currently supported only for neighborhood.type = 'adaptive.radius'.",
                 call. = FALSE)
        }
        if (!identical(lambda.selection, "cv")) {
            stop("support.selection = 'cv' requires lambda.selection = 'cv'.",
                 call. = FALSE)
        }
        tangent.dim.grid <- if (missing(tangent.dim) || is.null(tangent.dim)) {
            ncol(X)
        } else {
            .validate.ssrhe.positive.integer(tangent.dim, "tangent.dim")
        }
        support.grid <- .validate.ssrhe.support.grid(
            support.grid = support.grid,
            n = n,
            tangent.dim = tangent.dim.grid,
            derivative.order = derivative.order,
            support.buffer = support.buffer,
            max.candidates = support.cv.max.candidates
        )
        fold.id.cv <- solver.args$fold.id %||% .prepare.ssrhe.cv.folds(
            observed = as.vector(y.info$observed & W > 0),
            nfolds = nfolds,
            fold.id = NULL
        )
        candidate.fits <- vector("list", nrow(support.grid))
        candidate.rows <- vector("list", nrow(support.grid))
        for (ii in seq_len(nrow(support.grid))) {
            max.support.ii <- support.grid$max.support[ii]
            if (is.na(max.support.ii)) max.support.ii <- NULL
            elapsed <- system.time({
                cand <- tryCatch(
                    fit.ssrhe.hessian.l1.regression(
                        X = X,
                        y = y,
                        k = k,
                        tangent.dim = if (missing(tangent.dim)) NULL else tangent.dim,
                        lambda.grid = lambda.grid,
                        lambda.selection = "cv",
                        weights = weights,
                        n.lambda = n.lambda,
                        nfolds = nfolds,
                        fold.id = fold.id.cv,
                        loss = loss,
                        selection = selection,
                        support.selection = "rule",
                        nn.index = nn.index,
                        neighborhood.type = "adaptive.radius",
                        support.index = support.index,
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        radius.rule = radius.rule,
                        radius.factor = radius.factor,
                        min.support = support.grid$min.support[ii],
                        max.support = max.support.ii,
                        support.buffer = support.buffer,
                        support.topup = support.topup,
                        tangent.dim.rule = tangent.dim.rule,
                        eigen.tolerance = eigen.tolerance,
                        derivative.order = derivative.order,
                        pinv.tol = pinv.tol,
                        local.solver = local.solver,
                        normal.equations.max.condition = normal.equations.max.condition,
                        solver = solver,
                        row.scaling = row.scaling,
                        admm.rho = admm.rho,
                        admm.maxiter = admm.maxiter,
                        admm.abstol = admm.abstol,
                        admm.reltol = admm.reltol,
                        admm.adaptive.rho = admm.adaptive.rho,
                        maxsteps = maxsteps,
                        minlam = minlam,
                        approx = approx,
                        rtol = rtol,
                        btol = btol,
                        eps = eps,
                        return.local.diagnostics = return.local.diagnostics,
                        return.timing = return.timing,
                        verbose = verbose
                    ),
                    error = function(e) e
                )
            })[["elapsed"]]
            if (inherits(cand, "error")) {
                candidate.rows[[ii]] <- data.frame(
                    candidate = ii,
                    status = "error",
                    adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                    min.support = support.grid$min.support[ii],
                    max.support = support.grid$max.support[ii],
                    cv.error = Inf,
                    lambda = NA_real_,
                    runtime.sec = unname(elapsed),
                    message = conditionMessage(cand),
                    stringsAsFactors = FALSE
                )
            } else {
                candidate.fits[[ii]] <- cand
                diag <- .ssrhe.support.diagnostics(cand$operator)
                candidate.rows[[ii]] <- cbind(
                    data.frame(
                        candidate = ii,
                        status = "ok",
                        adaptive.k.scale = support.grid$adaptive.k.scale[ii],
                        min.support = support.grid$min.support[ii],
                        max.support = support.grid$max.support[ii],
                        cv.error = cand$cv$error.selected,
                        cv.se = cand$cv$se[cand$cv$selected.idx],
                        lambda = cand$lambda,
                        runtime.sec = unname(elapsed),
                        message = "",
                        stringsAsFactors = FALSE
                    ),
                    diag
                )
            }
        }
        support.cv.table <- do.call(rbind, candidate.rows)
        if (!any(is.finite(support.cv.table$cv.error))) {
            stop("All SSRHE Hessian L1 support-CV candidates failed.",
                 call. = FALSE)
        }
        selected.support.idx <- which.min(support.cv.table$cv.error)
        out <- candidate.fits[[selected.support.idx]]
        out$support.selection <- "cv"
        out$support.cv.table <- support.cv.table
        out$selected.support <- support.grid[selected.support.idx, , drop = FALSE]
        out$selection$support.index <- selected.support.idx
        out$selection$support.cv.error <- support.cv.table$cv.error[selected.support.idx]
        attr(out, "call") <- match.call()
        class(out) <- c("ssrhe.hessian.l1.fit", "list")
        return(out)
    }

    operator <- ssrhe.hessian.operator(
        X = X,
        k = k,
        tangent.dim = tangent.dim,
        nn.index = nn.index,
        neighborhood.type = neighborhood.type,
        support.index = support.index,
        adaptive.k.scale = adaptive.k.scale,
        radius.rule = radius.rule,
        radius.factor = radius.factor,
        min.support = min.support,
        max.support = max.support,
        support.buffer = support.buffer,
        support.topup = support.topup,
        tangent.dim.rule = tangent.dim.rule,
        eigen.tolerance = eigen.tolerance,
        derivative.order = derivative.order,
        stabilizer = FALSE,
        pinv.tol = pinv.tol,
        local.solver = local.solver,
        normal.equations.max.condition = normal.equations.max.condition,
        return.A = TRUE,
        return.B = FALSE,
        return.BS = FALSE,
        return.sparse = TRUE,
        return.local.diagnostics = return.local.diagnostics,
        return.timing = return.timing,
        verbose = verbose
    )

    out <- .fit.ssrhe.hessian.l1.from.operator(
        operator = operator,
        y = as.vector(y.info$Y),
        weights = as.vector(W),
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        n.lambda = n.lambda,
        fold.id = solver.args$fold.id,
        loss = loss,
        selection = selection,
        solver.args = solver.args,
        verbose = verbose
    )
    out$X <- X
    out$support.selection <- "rule"
    out$selected.support <- data.frame(
        adaptive.k.scale = operator$parameters$adaptive.k.scale,
        min.support = operator$neighborhoods$min.support %||% NA_integer_,
        max.support = operator$neighborhoods$max.support %||% NA_integer_
    )
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.l1.fit", "list")
    out
}

#' Refit SSRHE-Style Hessian L1 Regression
#'
#' Reuses the \eqn{A} operator from
#' \code{\link{fit.ssrhe.hessian.l1.regression}} to fit a new response or a new
#' lambda grid without rebuilding local neighborhoods or local PCA Hessian
#' rows.
#'
#' @param object A \code{"ssrhe.hessian.l1.fit"} object.
#' @param y Optional new numeric response vector. If \code{NULL}, the
#'   original response is reused.
#' @inheritParams fit.ssrhe.hessian.l1.regression
#'
#' @param ... Additional arguments are not currently supported.
#' @return A list of class \code{"ssrhe.hessian.l1.refit"}.
#' @examples
#' X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
#' fit <- fit.ssrhe.hessian.l1.regression(
#'   X, X[, 1]^2, k = 6L, tangent.dim = 1L,
#'   lambda.grid = 0.05, lambda.selection = "fixed", solver = "admm"
#' )
#' refit(
#'   fit, y = X[, 1]^3, solver = "admm"
#' )
#' @seealso \code{\link{refit}} for supported fitted objects and migration.
#' @method refit ssrhe.hessian.l1.fit
#' @aliases refit.ssrhe.hessian.l1.regression
#' @export
refit.ssrhe.hessian.l1.fit <- function(
    object,
    y = NULL,
    lambda.grid = object$lambda,
    lambda.selection = c("fixed", "cv"),
    weights = NULL,
    n.lambda = 40L,
    nfolds = 5L,
    fold.id = NULL,
    loss = c("mse", "mae"),
    selection = c("min", "one.se"),
    solver = c("genlasso", "admm", "auto"),
    row.scaling = c("none", "l2"),
    admm.rho = NULL,
    admm.maxiter = 2000L,
    admm.abstol = 1e-6,
    admm.reltol = 1e-5,
    maxsteps = 2000L,
    minlam = 0,
    approx = FALSE,
    rtol = 1e-7,
    btol = 1e-7,
    eps = 1e-4,
    verbose = FALSE, admm.adaptive.rho = TRUE, ...) {
    .geosmooth.check.dots(...)

    if (!inherits(object, "ssrhe.hessian.l1.fit")) {
        stop("object must be a 'ssrhe.hessian.l1.fit' object.",
             call. = FALSE)
    }
    if (is.null(object$operator) ||
        !inherits(object$operator, "ssrhe.hessian.operator")) {
        stop("object does not contain a reusable SSRHE operator.",
             call. = FALSE)
    }
    if (is.null(object$operator$A)) {
        stop("object operator does not contain A.", call. = FALSE)
    }
    if (is.null(y)) {
        y <- object$y
    }
    n <- length(object$fitted.values)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("refit() for Hessian L1 currently supports one response vector.",
             call. = FALSE)
    }
    if (is.null(weights)) {
        weights <- object$weights
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    lambda.selection <- match.arg(lambda.selection)
    loss <- match.arg(loss)
    selection <- match.arg(selection)
    solver <- match.arg(solver)
    if (!identical(solver, "admm") &&
        !requireNamespace("genlasso", quietly = TRUE)) {
        stop("Package 'genlasso' is required for the selected SSRHE Hessian L1 solver.",
             call. = FALSE)
    }
    solver.args <- .validate.ssrhe.hessian.l1.solver.args(
        n.lambda = n.lambda,
        nfolds = nfolds,
        fold.id = fold.id,
        observed = as.vector(y.info$observed & W > 0),
        maxsteps = maxsteps,
        minlam = minlam,
        approx = approx,
        rtol = rtol,
        btol = btol,
        eps = eps,
        solver = solver,
        row.scaling = row.scaling,
        admm.rho = admm.rho,
        admm.maxiter = admm.maxiter,
        admm.abstol = admm.abstol,
        admm.reltol = admm.reltol,
        admm.adaptive.rho = admm.adaptive.rho,
        verbose = verbose
    )
    lambda.grid <- .validate.ssrhe.hessian.l1.lambda.grid(
        lambda.grid,
        lambda.selection
    )

    out <- .fit.ssrhe.hessian.l1.from.operator(
        operator = object$operator,
        y = as.vector(y.info$Y),
        weights = as.vector(W),
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        n.lambda = n.lambda,
        fold.id = solver.args$fold.id,
        loss = loss,
        selection = selection,
        solver.args = solver.args,
        verbose = verbose
    )
    attr(out, "call") <- match.call()
    class(out) <- c("ssrhe.hessian.l1.refit", "ssrhe.hessian.l1.fit", "list")
    out
}

.fit.ssrhe.hessian.l1.from.operator <- function(operator,
                                                y,
                                                weights,
                                         lambda.grid,
                                         lambda.selection,
                                                n.lambda,
                                                fold.id,
                                                loss,
                                                selection,
                                                solver.args,
                                                verbose = FALSE) {
    if (is.null(operator$A)) {
        stop("operator must contain A. Rebuild with return.A = TRUE.",
             call. = FALSE)
    }
    n <- ncol(operator$A)
    y.info <- .prepare.ssrhe.response.matrix(y, n, "y")
    if (ncol(y.info$Y) != 1L) {
        stop("SSRHE Hessian L1 fitting currently supports one response vector.",
             call. = FALSE)
    }
    W <- .prepare.ssrhe.weight.matrix(weights, y.info, n)
    y.raw <- as.vector(y.info$Y)
    weights <- as.vector(W)
    y.clean <- y.raw
    y.clean[!is.finite(y.clean)] <- 0

    D.info <- .ssrhe.hessian.l1.penalty.matrix(
        operator$A,
        row.scaling = solver.args$row.scaling
    )
    diagnostics <- .ssrhe.hessian.l1.diagnostics(
        A = operator$A,
        D.info = D.info
    )
    backend <- solver.args$solver
    path.info <- NULL
    path <- NULL
    if (!identical(backend, "admm")) {
        path.info <- tryCatch(
            .fit.ssrhe.hessian.l1.path(
                y = y.clean,
                weights = weights,
                D = D.info$D,
                solver.args = solver.args,
                use.svd = D.info$svd
            ),
            error = function(e) {
                if (identical(backend, "auto")) {
                    structure(list(error = e), class = "ssrhe.l1.path.error")
                } else {
                    stop(e)
                }
            }
        )
        if (!inherits(path.info, "ssrhe.l1.path.error")) {
            path <- path.info$path
        }
        if (is.null(lambda.grid)) {
            if (inherits(path.info, "ssrhe.l1.path.error")) {
                stop("lambda.grid is required when solver = 'auto' falls back before a genlasso path is available.",
                     call. = FALSE)
            }
            lambda.grid <- .ssrhe.hessian.l1.default.lambda.grid(path, n.lambda)
        }
    } else if (is.null(lambda.grid)) {
        stop("lambda.grid is required when solver = 'admm'.", call. = FALSE)
    }
    if (identical(lambda.selection, "fixed") && length(lambda.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda value.",
             call. = FALSE)
    }

    cv <- NULL
    if (identical(lambda.selection, "cv")) {
        if (is.null(fold.id)) {
            fold.id <- .prepare.ssrhe.cv.folds(
                observed = as.vector(y.info$observed & weights > 0),
                nfolds = solver.args$nfolds,
                fold.id = NULL
            )
        }
        cv <- .ssrhe.hessian.l1.cv(
            y = y.clean,
            weights = weights,
            D = D.info$D,
            lambda.grid = lambda.grid,
            fold.id = fold.id,
            loss = loss,
            selection = selection,
            solver.args = solver.args,
            use.svd = D.info$svd
        )
        selected.idx <- cv$selected.idx
    } else {
        selected.idx <- 1L
    }

    # A finite generalized-lasso path is not an optimality certificate. Its
    # rank decisions (and its ridge for missing labels) can change the answer.
    # Solve the original weighted objective at every proposed lambda, including
    # when the caller explicitly requests the path backend.
    path.beta <- if (!is.null(path)) {
        .ssrhe.hessian.l1.coef.matrix(path, lambda.grid, n)
    } else NULL
    beta.grid <- .ssrhe.hessian.l1.admm.grid(
        y = y.clean, weights = weights, D = D.info$D.sparse,
        lambda.grid = lambda.grid, solver.args = solver.args
    )
    fit.source <- if (identical(backend, "admm")) "admm" else if (
        inherits(path.info, "ssrhe.l1.path.error")) "admm_fallback" else
        "genlasso_admm_refined"
    path.objective <- NULL
    if (!is.null(path.beta)) {
        path.objective <- vapply(seq_along(lambda.grid), function(i) {
            b <- path.beta[, i]
            if (any(!is.finite(b))) return(NA_real_)
            .5 * sum(weights * (y.clean - b)^2) +
                lambda.grid[i] * sum(abs(D.info$D.sparse %*% b))
        }, numeric(1))
    }
    candidate.status <- .ssrhe.hessian.l1.candidate.status(beta.grid)
    fitted <- beta.grid[, selected.idx]
    if (!all(is.finite(fitted))) {
        stop("The selected Hessian L1 fit failed: ",
             candidate.status$status[selected.idx], ". ",
             candidate.status$message[selected.idx], call. = FALSE)
    }
    if (!candidate.status$acceptable[selected.idx]) {
        message <- paste0(
            "The selected Hessian L1 ADMM fit did not converge (",
            candidate.status$status[selected.idx], ", ",
            candidate.status$iterations[selected.idx], " iterations). ",
            "Inspect solver$admm; consider adaptive rho or a larger admm.maxiter.")
        if (identical(lambda.selection, "cv")) stop(message, call. = FALSE)
        warning(message, call. = FALSE)
    }
    observed <- as.vector(y.info$observed & weights > 0)
    residuals <- rep(NA_real_, n)
    residuals[observed] <- y.raw[observed] - fitted[observed]
    data.loss <- 0.5 * sum(weights * (y.clean - fitted)^2)
    hessian.l1 <- sum(abs(as.vector(operator$A %*% fitted)))
    penalty.l1 <- sum(abs(as.vector(D.info$D.sparse %*% fitted)))
    lambda <- lambda.grid[selected.idx]

    list(
        fitted.values = as.vector(fitted),
        residuals = residuals,
        y = y.raw,
        weights = weights,
        lambda = lambda,
        lambda.grid = lambda.grid,
        lambda.selection = lambda.selection,
        objective = data.loss + lambda * penalty.l1,
        energies = list(data.loss = data.loss, hessian.l1 = hessian.l1,
                        penalty.l1 = penalty.l1),
        operator = operator,
        beta.grid = beta.grid,
        path = path,
        cv = cv,
        fold.id = fold.id,
        observed = observed,
        n.responses = 1L,
        diagnostics = diagnostics,
        solver = list(
            backend = fit.source,
            status = candidate.status$status[selected.idx],
            converged = candidate.status$converged[selected.idx],
            candidate.status = candidate.status,
            warnings = if (is.null(path.info)) character() else
                path.info$warnings %||% character(),
            path.error = if (inherits(path.info, "ssrhe.l1.path.error"))
                conditionMessage(path.info$error) else NULL,
            requested = backend,
            path.objective = path.objective,
            representation = D.info$representation,
            svd = D.info$svd,
            row.scaling = D.info$row.scaling,
            n.observed = if (is.null(path.info) ||
                inherits(path.info, "ssrhe.l1.path.error")) {
                sum(observed)
            } else {
                path.info$n.observed
            },
            admm = if (grepl("admm", fit.source, fixed = TRUE)) {
                attr(beta.grid, "admm")
            } else {
                NULL
            }
        ),
        selection = list(
            rule = selection,
            loss = loss,
            selected.index = selected.idx,
            lambda = lambda
        )
    )
}

.ssrhe.hessian.l1.penalty.matrix <- function(A,
                                             row.scaling = c("none", "l2")) {
    row.scaling <- match.arg(row.scaling)
    A.sparse <- methods::as(A, "dgCMatrix")
    row.scale <- rep(1, nrow(A.sparse))
    if (identical(row.scaling, "l2")) {
        row.norm <- sqrt(Matrix::rowSums(A.sparse^2))
        nz <- is.finite(row.norm) & row.norm > 0
        row.scale[nz] <- 1 / row.norm[nz]
        A.sparse <- Matrix::Diagonal(x = row.scale) %*% A.sparse
        A.sparse <- methods::as(A.sparse, "dgCMatrix")
    }
    use.svd <- nrow(A) >= ncol(A)
    list(
        D = if (use.svd) as.matrix(A.sparse) else A.sparse,
        D.sparse = A.sparse,
        svd = use.svd,
        representation = if (use.svd) "dense" else "sparse",
        row.scaling = row.scaling,
        row.scale = row.scale
    )
}

.ssrhe.hessian.l1.diagnostics <- function(A, D.info) {
    A.sparse <- methods::as(A, "dgCMatrix")
    scaled <- D.info$D.sparse
    row.norm <- sqrt(Matrix::rowSums(A.sparse^2))
    scaled.row.norm <- sqrt(Matrix::rowSums(scaled^2))
    list(
        nrow = nrow(A.sparse),
        ncol = ncol(A.sparse),
        nnzero = Matrix::nnzero(A.sparse),
        density = Matrix::nnzero(A.sparse) / (nrow(A.sparse) * ncol(A.sparse)),
        row.norm = list(
            min = suppressWarnings(min(row.norm, na.rm = TRUE)),
            median = stats::median(row.norm),
            max = suppressWarnings(max(row.norm, na.rm = TRUE)),
            zero = sum(!is.finite(row.norm) | row.norm == 0)
        ),
        scaled.row.norm = list(
            min = suppressWarnings(min(scaled.row.norm, na.rm = TRUE)),
            median = stats::median(scaled.row.norm),
            max = suppressWarnings(max(scaled.row.norm, na.rm = TRUE)),
            zero = sum(!is.finite(scaled.row.norm) | scaled.row.norm == 0)
        ),
        row.scaling = D.info$row.scaling,
        representation = D.info$representation,
        svd = D.info$svd
    )
}

.fit.ssrhe.hessian.l1.path <- function(y, weights, D, solver.args,
                                       use.svd = FALSE) {
    path.warnings <- character()
    capture.path <- function(expr) {
        withCallingHandlers(expr, warning = function(w) {
            path.warnings <<- c(path.warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    }
    observed <- is.finite(y) & is.finite(weights) & weights > 0
    if (!any(observed)) {
        stop("At least one observed positive-weight response is required.",
             call. = FALSE)
    }
    if (all(observed) && all(abs(weights - 1) < sqrt(.Machine$double.eps))) {
        path <- capture.path(genlasso::genlasso(
            y = y,
            D = D,
            approx = solver.args$approx,
            maxsteps = solver.args$maxsteps,
            minlam = solver.args$minlam,
            rtol = solver.args$rtol,
            btol = solver.args$btol,
            eps = solver.args$eps,
            verbose = solver.args$verbose,
            svd = use.svd
        ))
    } else {
        obs <- which(observed)
        sw <- sqrt(weights[obs])
        X.design <- as.matrix(Matrix::sparseMatrix(
            i = seq_along(obs),
            j = obs,
            x = sw,
            dims = c(length(obs), length(y)),
            giveCsparse = TRUE
        ))
        path <- capture.path(genlasso::genlasso(
            y = sw * y[obs],
            X = X.design,
            D = D,
            approx = solver.args$approx,
            maxsteps = solver.args$maxsteps,
            minlam = solver.args$minlam,
            rtol = solver.args$rtol,
            btol = solver.args$btol,
            eps = solver.args$eps,
            verbose = solver.args$verbose,
            svd = use.svd
        ))
    }
    list(path = path, n.observed = sum(observed), warnings = unique(path.warnings))
}

.ssrhe.hessian.l1.default.lambda.grid <- function(path, n.lambda) {
    lambda.path <- path$lambda
    lambda.path <- lambda.path[is.finite(lambda.path) & lambda.path >= 0]
    lambda.max <- suppressWarnings(max(lambda.path, na.rm = TRUE))
    if (!is.finite(lambda.max) || lambda.max <= 0) {
        lambda.max <- 1
    }
    n.lambda <- max(2L, as.integer(n.lambda))
    lambda.positive <- lambda.path[lambda.path > 0]
    complete.path <- isTRUE(path$completepath)
    lambda.min <- if (complete.path || !length(lambda.positive)) {
        lambda.max * 1e-4
    } else {
        max(min(lambda.positive) * (1 + 1e-8), lambda.max * 1e-4)
    }
    grid <- exp(seq(log(lambda.max), log(lambda.min), length.out = n.lambda - 1L))
    if (complete.path) grid <- c(grid, 0)
    unique(grid)
}

.ssrhe.hessian.l1.coef.matrix <- function(path, lambda.grid, n) {
    beta <- tryCatch(
        as.matrix(stats::coef(path, lambda = lambda.grid)$beta),
        error = function(e) NULL
    )
    if (is.null(beta) ||
        !identical(dim(beta), c(as.integer(n), length(lambda.grid)))) {
        beta <- vapply(lambda.grid, function(lambda) {
            tryCatch(
                as.vector(stats::coef(path, lambda = lambda)$beta),
                error = function(e) rep(NA_real_, n)
            )
        }, numeric(n))
    }
    beta
}

.ssrhe.hessian.l1.cv <- function(y, weights, D, lambda.grid, fold.id, loss,
                                 selection, solver.args, use.svd = FALSE) {
    n <- length(y)
    folds <- sort(unique(fold.id[fold.id > 0L]))
    cv.errors <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda.grid))
    rownames(cv.errors) <- paste0("Fold", folds)
    colnames(cv.errors) <- format(signif(lambda.grid, 6), scientific = TRUE)

    fold.status <- matrix("solver_error", nrow(cv.errors), ncol(cv.errors),
                          dimnames = dimnames(cv.errors))
    fold.converged <- matrix(NA, nrow(cv.errors), ncol(cv.errors),
                             dimnames = dimnames(cv.errors))
    fold.diagnostics <- vector("list", length(folds))
    names(fold.diagnostics) <- rownames(cv.errors)

    for (ii in seq_along(folds)) {
        fold <- folds[ii]
        test <- which(fold.id == fold)
        train.weights <- weights
        train.weights[test] <- 0
        failure <- character()
        path.warnings <- character()
        failed <- function(e) {
            failure <<- c(failure, conditionMessage(e))
            NULL
        }
        # The full-data path already proposed the grid. Holding out labels
        # introduces a rank-deficient design; score only original-objective
        # ADMM solutions, never genlasso's ridge-augmented fold coefficients.
        beta <- tryCatch(.ssrhe.hessian.l1.admm.grid(
            y = y, weights = train.weights, D = D,
            lambda.grid = lambda.grid, solver.args = solver.args
        ), error = failed)
        if (is.null(beta)) {
            fold.diagnostics[[ii]] <- list(errors = failure, warnings = path.warnings)
            next
        }
        status <- .ssrhe.hessian.l1.candidate.status(beta)
        fold.status[ii, ] <- status$status
        fold.converged[ii, ] <- status$converged
        fold.diagnostics[[ii]] <- list(
            candidate.status = status, admm = attr(beta, "admm"),
            errors = failure, warnings = path.warnings
        )
        beta.test <- beta[test, , drop = FALSE]
        finite.pred <- colSums(is.finite(beta.test)) == nrow(beta.test)
        target <- matrix(y[test], nrow = length(test), ncol = length(lambda.grid))
        w.test <- matrix(weights[test], nrow = length(test), ncol = length(lambda.grid))
        if (identical(loss, "mse")) {
            fold.error <- colSums(w.test * (target - beta.test)^2) / colSums(w.test)
        } else {
            fold.error <- colSums(w.test * abs(target - beta.test)) / colSums(w.test)
        }
        fold.error[!finite.pred | !status$acceptable] <- NA_real_
        cv.errors[ii, ] <- fold.error
    }

    n.successful.folds <- colSums(is.finite(cv.errors))
    eligible.candidates <- n.successful.folds == length(folds)
    mean.error <- colMeans(cv.errors)
    mean.error[!eligible.candidates | !is.finite(mean.error)] <- Inf
    se <- apply(cv.errors, 2L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(Inf)
        stats::sd(x) / sqrt(length(x))
    })
    best.idx <- which.min(mean.error)
    if (!is.finite(mean.error[best.idx])) {
        counts <- table(fold.status)
        reason <- paste(paste0(names(counts), "=", as.integer(counts)), collapse = ", ")
        condition <- structure(list(
            message = paste0("No Hessian L1 CV candidate succeeded on every fold (",
                             reason, "). Inspect fold diagnostics or increase admm.maxiter."),
            call = NULL, fold.status = fold.status,
            fold.diagnostics = fold.diagnostics, fold.errors = cv.errors,
            n.successful.folds = n.successful.folds, lambda.grid = lambda.grid
        ), class = c("ssrhe_l1_cv_error", "error", "condition"))
        stop(condition)
    }
    if (identical(selection, "one.se")) {
        threshold <- mean.error[best.idx] + se[best.idx]
        eligible <- which(eligible.candidates & is.finite(mean.error) &
                              mean.error <= threshold)
        selected.idx <- eligible[which.max(lambda.grid[eligible])]
    } else {
        selected.idx <- best.idx
    }
    list(
        fold.id = fold.id,
        lambda.grid = lambda.grid,
        fold.errors = cv.errors,
        fold.status = fold.status,
        fold.converged = fold.converged,
        fold.diagnostics = fold.diagnostics,
        n.successful.folds = n.successful.folds,
        eligible = eligible.candidates,
        mean.error = mean.error,
        se = se,
        best.idx = best.idx,
        selected.idx = selected.idx,
        lambda.min = lambda.grid[best.idx],
        lambda = lambda.grid[selected.idx],
        error.min = mean.error[best.idx],
        error.selected = mean.error[selected.idx],
        selection = selection,
        loss = loss
    )
}

.validate.ssrhe.hessian.l1.lambda.grid <- function(lambda.grid,
                                            lambda.selection) {
    if (is.null(lambda.grid)) {
        if (identical(lambda.selection, "fixed")) {
            stop("lambda.grid is required when lambda.selection = 'fixed'.",
                 call. = FALSE)
        }
        return(NULL)
    }
    lambda.grid <- .validate.ssrhe.lambda.grid(lambda.grid, "lambda.grid")
    if (identical(lambda.selection, "fixed") && length(lambda.grid) != 1L) {
        stop("lambda.selection = 'fixed' requires exactly one lambda value.",
             call. = FALSE)
    }
    lambda.grid
}

.validate.ssrhe.hessian.l1.solver.args <- function(n.lambda,
                                                   nfolds,
                                                   fold.id,
                                                   observed,
                                                   maxsteps,
                                                   minlam,
                                                   approx,
                                                   rtol,
                                                   btol,
                                                   eps,
                                                   solver = c("genlasso", "admm", "auto"),
                                                   row.scaling = c("none", "l2"),
                                                   admm.rho = NULL,
                                                   admm.maxiter = 2000L,
                                                   admm.abstol = 1e-6,
                                                   admm.reltol = 1e-5,
                                                   verbose,
                                                   admm.adaptive.rho = TRUE) {
    n.lambda <- .validate.ssrhe.positive.integer(n.lambda, "n.lambda")
    nfolds <- .validate.ssrhe.positive.integer(nfolds, "nfolds")
    maxsteps <- .validate.ssrhe.positive.integer(maxsteps, "maxsteps")
    solver <- match.arg(solver)
    row.scaling <- match.arg(row.scaling)
    minlam <- .validate.ssrhe.nonnegative.scalar(minlam, "minlam")
    rtol <- .validate.ssrhe.nonnegative.scalar(rtol, "rtol")
    btol <- .validate.ssrhe.nonnegative.scalar(btol, "btol")
    eps <- .validate.ssrhe.nonnegative.scalar(eps, "eps")
    if (!is.null(admm.rho)) {
        admm.rho <- .validate.ssrhe.numeric.scalar(admm.rho, "admm.rho")
    }
    if (!is.null(admm.rho) && (!is.finite(admm.rho) || admm.rho <= 0)) {
        stop("admm.rho must be a finite positive numeric scalar.",
             call. = FALSE)
    }
    if (!is.logical(admm.adaptive.rho) || length(admm.adaptive.rho) != 1L ||
        is.na(admm.adaptive.rho)) {
        stop("admm.adaptive.rho must be TRUE or FALSE.", call. = FALSE)
    }
    admm.maxiter <- .validate.ssrhe.positive.integer(admm.maxiter, "admm.maxiter")
    admm.abstol <- .validate.ssrhe.nonnegative.scalar(admm.abstol, "admm.abstol")
    admm.reltol <- .validate.ssrhe.nonnegative.scalar(admm.reltol, "admm.reltol")
    if (!is.null(fold.id)) {
        fold.id <- .prepare.ssrhe.cv.folds(observed, nfolds, fold.id)
    }
    list(
        n.lambda = n.lambda,
        nfolds = nfolds,
        fold.id = fold.id,
        maxsteps = maxsteps,
        minlam = minlam,
        approx = isTRUE(approx),
        rtol = rtol,
        btol = btol,
        eps = eps,
        solver = solver,
        row.scaling = row.scaling,
        admm.rho = admm.rho,
        admm.maxiter = admm.maxiter,
        admm.abstol = admm.abstol,
        admm.reltol = admm.reltol,
        admm.adaptive.rho = admm.adaptive.rho,
        verbose = isTRUE(verbose)
    )
}

#' @rdname geosmooth-print-methods
#' @method print ssrhe.hessian.l1.fit
#' @export
print.ssrhe.hessian.l1.fit <- function(x, ...) {
    cat("SSRHE Hessian L1 regression fit\n")
    if (!is.null(x$solver$status)) {
        cat("  solver:", x$solver$backend, "\n")
        cat("  status:", x$solver$status,
            if (identical(x$solver$converged, FALSE)) "(NOT CONVERGED)" else "",
            "\n")
        if (length(x$solver$warnings)) {
            cat("  backend warnings:", length(x$solver$warnings),
                "(see solver$warnings)\n")
        }
    }
    cat("  responses:", x$n.responses, "\n")
    cat("  lambda.selection:", x$lambda.selection, "\n")
    cat("  lambda:", format(x$lambda, digits = 4), "\n")
    cat("  objective:", format(x$objective, digits = 4), "\n")
    if (!is.null(x$cv)) {
        cat("  CV", x$cv$loss, ":",
            format(x$cv$error.selected, digits = 4), "\n")
    }
    invisible(x)
}

#' @rdname geosmooth-print-methods
#' @method print ssrhe.hessian.fit
#' @export
print.ssrhe.hessian.fit <- function(x, ...) {
    cat("SSRHE Hessian regression fit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  lambda1:", format(x$lambda$lambda1, digits = 4), "\n")
    cat("  lambda2:", format(x$lambda$lambda2, digits = 4), "\n")
    cat("  objective:", paste(format(x$objective, digits = 4), collapse = ", "), "\n")
    invisible(x)
}

#' @rdname geosmooth-print-methods
#' @method print ssrhe.hessian.refit
#' @export
print.ssrhe.hessian.refit <- function(x, ...) {
    cat("SSRHE Hessian regression refit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  lambda1:", format(x$lambda$lambda1, digits = 4), "\n")
    cat("  lambda2:", format(x$lambda$lambda2, digits = 4), "\n")
    invisible(x)
}

#' @rdname geosmooth-print-methods
#' @method print ssrhe.hessian.cv.fit
#' @export
print.ssrhe.hessian.cv.fit <- function(x, ...) {
    cat("SSRHE Hessian regression CV fit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  selected lambda1:", format(x$selection$lambda1, digits = 4), "\n")
    cat("  selected lambda2:", format(x$selection$lambda2, digits = 4), "\n")
    cat("  CV", x$selection$loss, ":",
        format(x$selection$cv.mean, digits = 4), "\n")
    invisible(x)
}

#' @rdname geosmooth-print-methods
#' @method print ssrhe.hessian.gcv.fit
#' @export
print.ssrhe.hessian.gcv.fit <- function(x, ...) {
    cat("SSRHE Hessian regression GCV fit\n")
    cat("  responses:", x$n.responses, "\n")
    cat("  selected lambda1:", format(x$selection$lambda1, digits = 4), "\n")
    cat("  selected lambda2:", format(x$selection$lambda2, digits = 4), "\n")
    cat("  GCV:", format(x$selection$gcv, digits = 4), "\n")
    cat("  trace(S):", format(x$selection$trace.S, digits = 4), "\n")
    invisible(x)
}
