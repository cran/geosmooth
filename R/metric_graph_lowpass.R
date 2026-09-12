#' Construct a Metric-Conductance Graph Low-Pass Operator
#'
#' Builds a weighted graph Laplacian by transforming metric edge lengths into
#' conductances. This operator is a direct metric-conductance comparator to
#' the legacy rdgraph regression smoother, not a replacement for the
#' Riemannian-complex/overlap-density smoother.
#'
#' @param adj.list List of integer neighbor vectors using 1-based vertex indices.
#' @param weight.list List of positive metric edge lengths parallel to
#'   \code{adj.list}. These are interpreted as edge lengths, not conductances.
#' @param conductance.rule Character scalar. One of
#'   \code{"inverse.length.power"}, \code{"exp.length"},
#'   \code{"exp.length.squared"}, or \code{"self.tuned.gaussian"}.
#' @param conductance.epsilon Positive numeric regularizer used by inverse-power
#'   conductances and as a local-scale floor.
#' @param conductance.alpha Positive numeric exponent for
#'   \code{"inverse.length.power"}.
#' @param conductance.sigma Optional positive global scale for exponential rules.
#'   If \code{NULL}, it is selected by \code{conductance.sigma.rule}.
#' @param conductance.sigma.rule Rule for selecting a global scale when needed.
#' @param conductance.sigma.quantile Quantile used when
#'   \code{conductance.sigma.rule = "edge.quantile"}.
#' @param conductance.local.k Positive integer local incident-edge order
#'   statistic for \code{"self.tuned.gaussian"}.
#' @param laplacian.type Laplacian operator. \code{"unnormalized"} uses the
#'   weighted graph Laplacian \eqn{L = D - C}. \code{"symmetric.normalized"}
#'   uses \eqn{L_{\mathrm{sym}} = I - D^{-1/2} C D^{-1/2}}.
#' @param return.sparse Logical. If \code{TRUE}, attach a \pkg{Matrix} sparse
#'   Laplacian when \pkg{Matrix} is available.
#' @param verbose Logical. Reserved for future diagnostic messages.
#'
#' @details
#' The legacy rdgraph regression precomputed-graph path uses
#' supplied \code{weight.list} values as edge lengths for neighborhood ordering.
#' The spectral conductance in that smoother is overlap-density based:
#' \deqn{c_e^\rho = 1 / \max(\rho_1(e), 10^{-10}),}
#' where \eqn{\rho_1(e)} is computed by the Riemannian-complex density machinery.
#'
#' This function instead constructs conductances directly from metric lengths:
#' \deqn{c_{ij} = \phi(\ell_{ij}).}
#' Supported phase-1 transforms are
#' \deqn{c_{ij}=(\ell_{ij}+\epsilon)^{-\alpha},}
#' \deqn{c_{ij}=\exp(-\ell_{ij}/\sigma),}
#' \deqn{c_{ij}=\exp(-\ell_{ij}^{2}/\sigma^{2}),}
#' and
#' \deqn{c_{ij}=\exp(-\ell_{ij}^{2}/(\sigma_i\sigma_j)).}
#'
#' For \code{laplacian.type = "symmetric.normalized"}, smoothing is performed
#' in the Euclidean eigenbasis of \eqn{L_{\mathrm{sym}}}. The null vector is
#' proportional to \eqn{\sqrt{d}}, not the constant vector, so this mode does
#' not preserve constant responses in the same way as the unnormalized
#' Laplacian.
#'
#' @return A list of class \code{"metric.graph.lowpass.operator"} containing
#'   the edge table, conductances, degree vector, Laplacian triplets, summaries,
#'   and optionally a sparse Laplacian matrix.
#' @examples
#' adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
#' lengths <- list(1, c(1, 1), c(1, 1), 1)
#' metric.graph.lowpass.operator(adj, lengths)
#' @export
metric.graph.lowpass.operator <- function(
    adj.list,
    weight.list,
    conductance.rule = c("inverse.length.power", "exp.length",
                         "exp.length.squared", "self.tuned.gaussian"),
    conductance.epsilon = 1e-8,
    conductance.alpha = 1,
    conductance.sigma = NULL,
    conductance.sigma.rule = c("edge.quantile", "median", "local.k"),
    conductance.sigma.quantile = 0.75,
    conductance.local.k = 5L,
    laplacian.type = c("unnormalized", "symmetric.normalized"),
    return.sparse = TRUE,
    verbose = FALSE) {

    graph <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
    args <- .validate.metric.graph.lowpass.operator.args(
        conductance.rule = conductance.rule,
        conductance.epsilon = conductance.epsilon,
        conductance.alpha = conductance.alpha,
        conductance.sigma = conductance.sigma,
        conductance.sigma.rule = conductance.sigma.rule,
        conductance.sigma.quantile = conductance.sigma.quantile,
        conductance.local.k = conductance.local.k,
        laplacian.type = laplacian.type,
        verbose = verbose
    )

    out <- rcpp_metric_graph_lowpass_operator(
        graph$adj.list.0based,
        graph$weight.list.cpp,
        args$conductance.rule,
        args$conductance.epsilon,
        args$conductance.alpha,
        args$conductance.sigma,
        args$conductance.sigma.rule,
        args$conductance.sigma.quantile,
        args$conductance.local.k,
        args$laplacian.type
    )

    out$graph$adj.list <- graph$adj.list
    out$graph$weight.list <- graph$weight.list
    out <- .attach.metric.graph.lowpass.laplacian(out, return.sparse)
    class(out) <- c("metric.graph.lowpass.operator", "list")
    out
}

#' Construct a Reusable Metric Graph Low-Pass Spectral Basis
#'
#' Builds a metric-conductance graph Laplacian and computes the low-frequency
#' eigensystem used by graph low-pass filters. The resulting object is
#' response-independent and can be reused for many responses and filter
#' parameters.
#'
#' @inheritParams metric.graph.lowpass.operator
#' @param n.eigenpairs Positive integer number of eigenpairs to compute.
#' @param eigen.solver \code{"auto"}, \code{"sparse"}, or \code{"dense"}.
#' @param dense.eigen.threshold Exact dense threshold for auto mode.
#' @param dense.fallback.threshold Maximum graph size for emergency dense
#'   fallback when sparse decomposition fails and fallback is allowed.
#' @param dense.fallback \code{"auto"}, \code{"never"}, or \code{"always"}.
#'
#' @details
#' The basis is complete only when it contains one eigenvector per graph
#' vertex. A truncated basis represents only the retained low-frequency
#' subspace. The largest retained eigenvalue supplies a conservative proxy for
#' bounding the contribution of omitted modes because all omitted eigenvalues
#' are at least as large.
#'
#' @return A list of class \code{"metric.graph.lowpass.basis"} containing the
#'   graph operator, eigenvalues, eigenvectors, solver metadata, and spectral
#'   completeness diagnostics.
#' @examples
#' adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
#' lengths <- list(1, c(1, 1), c(1, 1), 1)
#' metric.graph.lowpass.basis(adj, lengths, n.eigenpairs = 4L,
#'                            eigen.solver = "dense")
#' @export
metric.graph.lowpass.basis <- function(
    adj.list,
    weight.list,
    conductance.rule = c("inverse.length.power", "exp.length",
                         "exp.length.squared", "self.tuned.gaussian"),
    conductance.epsilon = 1e-8,
    conductance.alpha = 1,
    conductance.sigma = NULL,
    conductance.sigma.rule = c("edge.quantile", "median", "local.k"),
    conductance.sigma.quantile = 0.75,
    conductance.local.k = 5L,
    laplacian.type = c("unnormalized", "symmetric.normalized"),
    n.eigenpairs = 50L,
    eigen.solver = c("auto", "sparse", "dense"),
    dense.eigen.threshold = 200L,
    dense.fallback.threshold = 5000L,
    dense.fallback = c("auto", "never", "always"),
    verbose = FALSE) {

    graph <- .validate.metric.graph.lowpass.graph(adj.list, weight.list)
    n <- length(graph$adj.list)
    args <- .validate.metric.graph.lowpass.operator.args(
        conductance.rule = conductance.rule,
        conductance.epsilon = conductance.epsilon,
        conductance.alpha = conductance.alpha,
        conductance.sigma = conductance.sigma,
        conductance.sigma.rule = conductance.sigma.rule,
        conductance.sigma.quantile = conductance.sigma.quantile,
        conductance.local.k = conductance.local.k,
        laplacian.type = laplacian.type,
        verbose = verbose
    )
    solver <- .validate.metric.graph.lowpass.solver.args(
        n.eigenpairs = n.eigenpairs,
        n = n,
        eigen.solver = eigen.solver,
        dense.eigen.threshold = dense.eigen.threshold,
        dense.fallback.threshold = dense.fallback.threshold,
        dense.fallback = dense.fallback
    )

    raw <- rcpp_metric_graph_lowpass_spectrum(
        graph$adj.list.0based,
        graph$weight.list.cpp,
        args$conductance.rule,
        args$conductance.epsilon,
        args$conductance.alpha,
        args$conductance.sigma,
        args$conductance.sigma.rule,
        args$conductance.sigma.quantile,
        args$conductance.local.k,
        args$laplacian.type,
        solver$n.eigenpairs,
        solver$eigen.solver,
        solver$dense.eigen.threshold,
        solver$dense.fallback.threshold,
        solver$dense.fallback,
        as.logical(verbose)
    )

    operator <- raw$operator
    operator$graph$adj.list <- graph$adj.list
    operator$graph$weight.list <- graph$weight.list
    operator <- .attach.metric.graph.lowpass.laplacian(operator, TRUE)
    class(operator) <- c("metric.graph.lowpass.operator", "list")

    spectral <- raw$spectral
    spectral$eigenvalues <- as.double(spectral$eigenvalues)
    spectral$eigenvectors <- as.matrix(spectral$eigenvectors)
    .require.metric.graph.lowpass.finite(
        spectral$eigenvalues, "spectral$eigenvalues"
    )
    .require.metric.graph.lowpass.finite(
        spectral$eigenvectors, "spectral$eigenvectors"
    )
    spectral$n.vertices <- n
    spectral$n.eigenpairs <- length(spectral$eigenvalues)
    spectral$is.complete <- spectral$n.eigenpairs == n
    spectral$largest.retained.eigenvalue <- max(spectral$eigenvalues)
    positive <- spectral$eigenvalues[
        spectral$eigenvalues > .Machine$double.eps
    ]
    spectral$smallest.positive.eigenvalue <- if (length(positive)) {
        min(positive)
    } else {
        NA_real_
    }
    spectral$omitted.attenuation.bound.type <- if (spectral$is.complete) {
        "complete_spectrum"
    } else {
        "largest_retained_eigenvalue_proxy"
    }

    out <- list(
        graph = operator$graph,
        operator = operator,
        conductance = operator$conductance,
        laplacian = operator$laplacian,
        laplacian.type = operator$laplacian.type,
        spectral = spectral,
        parameters = c(args, solver)
    )
    class(out) <- c("metric.graph.lowpass.basis", "list")
    out
}

#' Construct a Graph Heat-Time Grid
#'
#' Creates a positive heat-time grid from a
#' \code{"metric.graph.lowpass.basis"} object.
#'
#' @param basis A \code{"metric.graph.lowpass.basis"} object.
#' @param rule Grid rule. \code{"w1_inverse_spectrum"} requires a complete
#'   basis and reproduces the W1 inverse-spectrum grid.
#'   \code{"spectral_guarded"} raises the lower endpoint for a truncated basis
#'   until the conservative omitted-mode attenuation bound is no larger than
#'   \code{truncation.tol}, and extends the upper endpoint until the slowest
#'   retained positive mode is attenuated to \code{equilibrium.tol}.
#' @param n.initial Number of positive grid values.
#' @param include.zero Logical. Include the exact no-smoothing endpoint.
#' @param truncation.tol Positive tolerance smaller than one for the
#'   conservative omitted-mode attenuation bound.
#' @param equilibrium.tol Positive tolerance smaller than one for the slowest
#'   positive retained mode at the upper endpoint.
#'
#' @return A numeric vector with grid-construction metadata stored as
#'   attributes.
#' @examples
#' adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
#' lengths <- list(1, c(1, 1), c(1, 1), 1)
#' basis <- metric.graph.lowpass.basis(adj, lengths, n.eigenpairs = 4L,
#'                                     eigen.solver = "dense")
#' metric.graph.heat.eta.grid(basis, n.initial = 6L)
#' @export
metric.graph.heat.eta.grid <- function(
    basis,
    rule = c("spectral_guarded", "w1_inverse_spectrum"),
    n.initial = 40L,
    include.zero = FALSE,
    truncation.tol = 1e-4,
    equilibrium.tol = 1e-4) {

    .validate.metric.graph.lowpass.basis(basis)
    rule <- match.arg(rule)
    if (rule == "w1_inverse_spectrum" &&
        !isTRUE(basis$spectral$is.complete)) {
        stop(
            "rule = 'w1_inverse_spectrum' requires a complete spectral basis.",
            call. = FALSE
        )
    }
    n.initial <- .validate.positive.integer.scalar(n.initial, "n.initial")
    if (n.initial < 2L) stop("n.initial must be at least 2.", call. = FALSE)
    include.zero <- .validate.logical.scalar(include.zero, "include.zero")
    truncation.tol <- .validate.unit.interval.open(
        truncation.tol, "truncation.tol"
    )
    equilibrium.tol <- .validate.unit.interval.open(
        equilibrium.tol, "equilibrium.tol"
    )

    positive <- sort(basis$spectral$eigenvalues[
        basis$spectral$eigenvalues > .Machine$double.eps
    ])
    if (!length(positive)) {
        stop("The basis has no positive retained eigenvalue.", call. = FALSE)
    }
    lambda.slow <- min(positive)
    lambda.cut <- max(positive)

    if (rule == "w1_inverse_spectrum") {
        eta.min <- 1 / lambda.cut
        eta.max <- 1 / lambda.slow
    } else {
        eta.min <- if (isTRUE(basis$spectral$is.complete)) {
            1 / lambda.cut
        } else {
            log(1 / truncation.tol) / lambda.cut
        }
        eta.max <- log(1 / equilibrium.tol) / lambda.slow
    }
    if (!is.finite(eta.min) || !is.finite(eta.max) ||
        eta.min <= 0 || eta.max <= eta.min) {
        stop(
            "The requested heat-time rule did not produce ordered positive endpoints.",
            call. = FALSE
        )
    }

    grid <- exp(seq(log(eta.min), log(eta.max), length.out = n.initial))
    if (include.zero) grid <- c(0, grid)
    attr(grid, "rule") <- rule
    attr(grid, "eta.min") <- eta.min
    attr(grid, "eta.max") <- eta.max
    attr(grid, "truncation.tol") <- truncation.tol
    attr(grid, "equilibrium.tol") <- equilibrium.tol
    attr(grid, "basis.complete") <- isTRUE(basis$spectral$is.complete)
    attr(grid, "lambda.slow") <- lambda.slow
    attr(grid, "lambda.cut") <- lambda.cut
    grid
}

#' Propose a Guarded Lower Graph Heat Time
#'
#' Proposes one lower positive heat time after an external search controller
#' has classified the current lower endpoint. The helper enforces a geometric
#' expansion step, an identity-departure floor, a maximum number of expansion
#' rounds, and the truncated-basis resolution certificate. It does not decide
#' whether an endpoint is competitive; fold-level or cohort-level search logic
#' remains the caller's responsibility.
#'
#' @param basis A \code{"metric.graph.lowpass.basis"} object.
#' @param eta.grid Current finite positive heat-time grid. The exact identity
#'   time \code{eta = 0} is intentionally excluded from competitive search.
#' @param endpoint.status External lower-endpoint classification:
#'   \code{"active"} proposes one extension, while \code{"inactive"} and
#'   \code{"unresolved"} return the grid unchanged.
#' @param expansion.factor Finite numeric factor greater than one. An active
#'   proposal is \code{min(eta.grid) / expansion.factor} before applying the
#'   identity floor.
#' @param max.expansions Nonnegative integer maximum number of lower-extension
#'   rounds.
#' @param expansions.completed Nonnegative integer number of extension rounds
#'   already admitted under the same search contract.
#' @param identity.departure Positive number smaller than one. The proposed
#'   time cannot be smaller than
#'   \code{-log(1 - identity.departure) / lambda.max}, where
#'   \code{lambda.max} is the largest retained eigenvalue. For a complete
#'   basis this is the exact largest graph-Laplacian eigenvalue.
#' @param truncation.tol Positive tolerance smaller than one for the
#'   conservative omitted-mode attenuation bound.
#' @param unresolved.action Action when an active proposal is not certified by
#'   a truncated basis: \code{"error"} or \code{"mark"}.
#'
#' @return A list of class \code{"metric.graph.heat.lower.extension"} with the
#'   augmented or unchanged grid and proposal telemetry. A proposal is added
#'   only when \code{admitted} is \code{TRUE}.
#' @examples
#' adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
#' lengths <- list(1, c(1, 1), c(1, 1), 1)
#' basis <- metric.graph.lowpass.basis(adj, lengths, n.eigenpairs = 4L,
#'                                     eigen.solver = "dense")
#' eta <- metric.graph.heat.eta.grid(basis, n.initial = 6L)
#' metric.graph.heat.extend.lower(basis, eta, endpoint.status = "active")
#' @export
metric.graph.heat.extend.lower <- function(
    basis,
    eta.grid,
    endpoint.status = c("active", "inactive", "unresolved"),
    expansion.factor = 3,
    max.expansions = 3L,
    expansions.completed = 0L,
    identity.departure = 0.01,
    truncation.tol = 1e-4,
    unresolved.action = c("error", "mark")) {

    .validate.metric.graph.lowpass.basis(basis)
    endpoint.status <- match.arg(endpoint.status)
    unresolved.action <- match.arg(unresolved.action)
    if (!is.numeric(eta.grid) || length(eta.grid) < 1L ||
        any(!is.finite(eta.grid)) || any(eta.grid <= 0)) {
        stop("eta.grid must be a finite positive numeric vector.", call. = FALSE)
    }
    eta.grid <- sort(unique(as.double(eta.grid)))
    if (!is.numeric(expansion.factor) || length(expansion.factor) != 1L ||
        is.na(expansion.factor) || !is.finite(expansion.factor) ||
        expansion.factor <= 1) {
        stop(
            "expansion.factor must be a finite numeric scalar greater than one.",
            call. = FALSE
        )
    }
    max.expansions <- .validate.nonnegative.integer.scalar(
        max.expansions, "max.expansions"
    )
    expansions.completed <- .validate.nonnegative.integer.scalar(
        expansions.completed, "expansions.completed"
    )
    identity.departure <- .validate.unit.interval.open(
        identity.departure, "identity.departure"
    )
    truncation.tol <- .validate.unit.interval.open(
        truncation.tol, "truncation.tol"
    )

    positive <- basis$spectral$eigenvalues[
        basis$spectral$eigenvalues > .Machine$double.eps
    ]
    if (!length(positive)) {
        stop("The basis has no positive retained eigenvalue.", call. = FALSE)
    }
    lambda.max <- max(positive)
    identity.floor <- -log1p(-identity.departure) / lambda.max
    complete <- isTRUE(basis$spectral$is.complete)

    result <- list(
        eta.grid = eta.grid,
        previous.eta.min = min(eta.grid),
        proposed.eta.raw = NA_real_,
        proposed.eta = NA_real_,
        identity.floor = identity.floor,
        identity.departure = identity.departure,
        lambda.max = lambda.max,
        lambda.max.type = if (complete) {
            "complete_spectrum"
        } else {
            "largest_retained_eigenvalue_proxy"
        },
        omitted.attenuation.bound = NA_real_,
        truncation.tol = truncation.tol,
        basis.complete = complete,
        endpoint.status = endpoint.status,
        expansion.factor = as.double(expansion.factor),
        max.expansions = max.expansions,
        expansions.completed = expansions.completed,
        expansion.round = NA_integer_,
        admitted = FALSE,
        status = switch(
            endpoint.status,
            inactive = "inactive_endpoint",
            unresolved = "unresolved_endpoint",
            active = "pending"
        )
    )

    if (endpoint.status != "active") {
        class(result) <- c("metric.graph.heat.lower.extension", "list")
        return(result)
    }
    if (expansions.completed >= max.expansions) {
        result$status <- "expansion_cap_reached"
        class(result) <- c("metric.graph.heat.lower.extension", "list")
        return(result)
    }

    proposed.raw <- min(eta.grid) / expansion.factor
    proposed <- max(proposed.raw, identity.floor)
    attenuation <- exp(-proposed * lambda.max)
    resolved <- complete ||
        attenuation <= truncation.tol * (1 + 1e-10)
    result$proposed.eta.raw <- proposed.raw
    result$proposed.eta <- proposed
    result$omitted.attenuation.bound <- attenuation
    result$expansion.round <- expansions.completed + 1L

    lower.endpoint.tol <- 16 * .Machine$double.eps *
        max(.Machine$double.xmin, abs(min(eta.grid)), abs(proposed))
    if (proposed >= min(eta.grid) - lower.endpoint.tol) {
        result$status <- "identity_floor_reached"
        class(result) <- c("metric.graph.heat.lower.extension", "list")
        return(result)
    }

    if (!resolved) {
        result$status <- "unresolved_truncated_basis"
        if (unresolved.action == "error") {
            stop(
                paste0(
                    "The lower heat-time proposal is not resolved by the ",
                    "truncated basis; compute more eigenpairs or use a complete ",
                    "basis before expanding."
                ),
                call. = FALSE
            )
        }
        class(result) <- c("metric.graph.heat.lower.extension", "list")
        return(result)
    }

    duplicate <- any(
        abs(eta.grid - proposed) <=
            16 * .Machine$double.eps *
                pmax(.Machine$double.xmin, abs(eta.grid), abs(proposed))
    )
    if (duplicate) {
        result$status <- "identity_floor_reached"
        class(result) <- c("metric.graph.heat.lower.extension", "list")
        return(result)
    }

    result$eta.grid <- sort(c(eta.grid, proposed))
    result$admitted <- TRUE
    result$status <- if (proposed > proposed.raw) {
        "admitted_at_identity_floor"
    } else {
        "admitted"
    }
    class(result) <- c("metric.graph.heat.lower.extension", "list")
    result
}

#' Apply a Metric Graph Low-Pass Filter Path
#'
#' Applies every requested low-pass parameter to one response or a matrix of
#' responses using a reusable spectral basis.
#'
#' @param basis A \code{"metric.graph.lowpass.basis"} object.
#' @param y Numeric response vector or matrix with one row per graph vertex.
#' @param eta.grid Numeric filter-parameter grid.
#' @param filter.type Spectral low-pass filter family.
#' @param block.size Optional number of response columns processed together.
#' @param truncation.tol Positive tolerance smaller than one used to classify
#'   truncated-basis candidates as spectrally resolved.
#' @param unresolved.action Action when a truncated-basis candidate exceeds
#'   \code{truncation.tol}: \code{"warn"}, \code{"error"}, or \code{"allow"}.
#' @param exact.zero Logical. For heat filtering, return the input exactly at
#'   \code{eta = 0}. This avoids treating a truncated spectral projection as
#'   the identity.
#'
#' @return A list of class \code{"metric.graph.lowpass.path"}. For one response,
#'   \code{fitted.values} is an \eqn{N} by \eqn{J} matrix. For multiple
#'   responses it is an \eqn{N} by \eqn{J} by \eqn{S} array.
#' @examples
#' adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
#' lengths <- list(1, c(1, 1), c(1, 1), 1)
#' basis <- metric.graph.lowpass.basis(adj, lengths, n.eigenpairs = 4L,
#'                                     eigen.solver = "dense")
#' apply.metric.graph.lowpass.path(basis, 1:4, eta.grid = c(0, 0.1, 1))
#' @export
apply.metric.graph.lowpass.path <- function(
    basis,
    y,
    eta.grid,
    filter.type = c("heat_kernel", "tikhonov", "cubic_spline",
                    "gaussian", "exponential", "butterworth"),
    block.size = NULL,
    truncation.tol = 1e-4,
    unresolved.action = c("warn", "error", "allow"),
    exact.zero = TRUE) {

    .validate.metric.graph.lowpass.basis(basis)
    filter.type <- match.arg(filter.type)
    unresolved.action <- match.arg(unresolved.action)
    exact.zero <- .validate.logical.scalar(exact.zero, "exact.zero")
    truncation.tol <- .validate.unit.interval.open(
        truncation.tol, "truncation.tol"
    )
    block.size <- .validate.optional.block.size(block.size)
    y.info <- .prepare.metric.graph.lowpass.response.matrix(
        y, basis$spectral$n.vertices
    )
    Y <- y.info$Y
    eta.grid <- .prepare.metric.graph.lowpass.eta.grid(
        eta.grid, basis$spectral$eigenvalues, filter.type, length(eta.grid)
    )

    V <- basis$spectral$eigenvectors
    eigenvalues <- basis$spectral$eigenvalues
    weights <- compute.filter.weights.matrix(
        eigenvalues, eta.grid, filter.type
    )
    Vt.Y <- crossprod(V, Y)
    .require.metric.graph.lowpass.finite(Vt.Y, "spectral response coefficients")

    lambda.cut <- max(eigenvalues)
    cutoff.weights <- as.numeric(compute.filter.weights.matrix(
        lambda.cut, eta.grid, filter.type
    ))
    complete <- isTRUE(basis$spectral$is.complete)
    exact.identity <- filter.type == "heat_kernel" & eta.grid == 0 & exact.zero
    resolved <- complete |
        cutoff.weights <= truncation.tol * (1 + 1e-10) |
        exact.identity
    unresolved <- which(!resolved)
    if (length(unresolved) && unresolved.action != "allow") {
        msg <- paste0(
            length(unresolved), " candidate(s) exceed the truncated-basis ",
            "attenuation tolerance; increase n.eigenpairs or choose a filter ",
            "parameter that more strongly attenuates omitted modes."
        )
        if (unresolved.action == "error") stop(msg, call. = FALSE)
        warning(msg, call. = FALSE)
    }

    n <- nrow(Y)
    n.eta <- length(eta.grid)
    n.responses <- ncol(Y)
    fitted <- array(NA_real_, dim = c(n, n.eta, n.responses))
    response.blocks <- .make.metric.graph.lowpass.block.index(
        n.responses, block.size
    )
    for (j in seq_len(n.eta)) {
        if (exact.identity[[j]]) {
            fitted[, j, ] <- Y
            next
        }
        for (cols in response.blocks) {
            fitted[, j, cols] <- V %*% (
                weights[, j] * Vt.Y[, cols, drop = FALSE]
            )
        }
    }
    .require.metric.graph.lowpass.finite(fitted, "path fitted values")

    response.names <- y.info$col.names
    eta.names <- format(eta.grid, digits = 10, trim = TRUE)
    if (n.responses == 1L) {
        fitted.out <- matrix(fitted[, , 1L], nrow = n, ncol = n.eta)
        colnames(fitted.out) <- eta.names
    } else {
        dimnames(fitted) <- list(
            NULL, eta = eta.names,
            response = response.names %||% paste0("response", seq_len(n.responses))
        )
        fitted.out <- fitted
    }
    effective.df <- colSums(weights)
    effective.df[exact.identity] <- n
    resolution <- data.frame(
        eta = eta.grid,
        retained.cutoff.weight = cutoff.weights,
        spectrally.resolved = resolved,
        exact.identity = exact.identity,
        stringsAsFactors = FALSE
    )
    out <- list(
        fitted.values = fitted.out,
        eta.grid = eta.grid,
        filter.type = filter.type,
        effective.df = effective.df,
        resolution = resolution,
        n.responses = n.responses,
        response.names = response.names,
        basis = basis,
        parameters = list(
            block.size = block.size,
            truncation.tol = truncation.tol,
            unresolved.action = unresolved.action,
            exact.zero = exact.zero
        )
    )
    class(out) <- c("metric.graph.lowpass.path", "list")
    out
}

.fit.metric.graph.lowpass.eta.search <- function(
    basis,
    y,
    y.spectral,
    V,
    filter.type,
    eta.grid,
    eta.search,
    eta.expansion.factor,
    eta.max.expansions,
    eta.identity.departure,
    eta.truncation.tol) {

    if (!is.numeric(eta.expansion.factor) ||
        length(eta.expansion.factor) != 1L ||
        is.na(eta.expansion.factor) ||
        !is.finite(eta.expansion.factor) ||
        eta.expansion.factor <= 1) {
        stop(
            "eta.expansion.factor must be a finite numeric scalar greater than one.",
            call. = FALSE
        )
    }
    eta.max.expansions <- .validate.nonnegative.integer.scalar(
        eta.max.expansions, "eta.max.expansions"
    )
    eta.identity.departure <- .validate.unit.interval.open(
        eta.identity.departure, "eta.identity.departure"
    )
    eta.truncation.tol <- .validate.unit.interval.open(
        eta.truncation.tol, "eta.truncation.tol"
    )
    settings <- list(
        expansion.factor = as.double(eta.expansion.factor),
        max.expansions = eta.max.expansions,
        identity.departure = eta.identity.departure,
        truncation.tol = eta.truncation.tol
    )
    evaluate <- function(grid) {
        weights <- compute.filter.weights.matrix(
            basis$spectral$eigenvalues, grid, filter.type
        )
        .require.metric.graph.lowpass.finite(weights, "filter weights")
        gcv <- .select.eta.gcv.single(
            y, y.spectral, V, weights, grid
        )
        .validate.metric.graph.lowpass.gcv.result(
            gcv.result = gcv,
            n = length(y),
            n.eta = length(grid),
            context = "fit.metric.graph.lowpass() eta search"
        )
        list(filter.weights.matrix = weights, gcv.result = gcv)
    }

    if (eta.search == "fixed") {
        selected <- evaluate(eta.grid)
        lower.idx <- which.min(eta.grid)
        rounds <- data.frame(
            round = 0L,
            eta.min = min(eta.grid),
            eta.max = max(eta.grid),
            n.candidates = length(eta.grid),
            best.idx = selected$gcv.result$best.idx,
            eta.optimal = selected$gcv.result$eta.optimal,
            gcv.optimal = selected$gcv.result$gcv.min,
            lower.endpoint.selected =
                selected$gcv.result$best.idx == lower.idx,
            extension.status = "not_requested",
            stringsAsFactors = FALSE
        )
        return(list(
            eta.grid = eta.grid,
            filter.weights.matrix = selected$filter.weights.matrix,
            gcv.result = selected$gcv.result,
            search = list(
                mode = eta.search,
                status = "fixed_grid",
                expansions.completed = 0L,
                settings = settings,
                rounds = rounds
            )
        ))
    }

    if (filter.type != "heat_kernel") {
        stop(
            "eta.search = 'guarded.gcv' requires filter.type = 'heat_kernel'.",
            call. = FALSE
        )
    }
    if (any(eta.grid <= 0)) {
        stop(
            paste0(
                "eta.search = 'guarded.gcv' requires positive competitive ",
                "times; eta = 0 may be evaluated separately as an identity ",
                "diagnostic."
            ),
            call. = FALSE
        )
    }

    eta.grid <- sort(unique(as.double(eta.grid)))
    expansions.completed <- 0L
    round.rows <- list()
    repeat {
        selected <- evaluate(eta.grid)
        lower.selected <- selected$gcv.result$best.idx == 1L
        extension.status <- if (lower.selected) "pending" else "not_requested"
        if (lower.selected) {
            extension <- metric.graph.heat.extend.lower(
                basis = basis,
                eta.grid = eta.grid,
                endpoint.status = "active",
                expansion.factor = eta.expansion.factor,
                max.expansions = eta.max.expansions,
                expansions.completed = expansions.completed,
                identity.departure = eta.identity.departure,
                truncation.tol = eta.truncation.tol,
                unresolved.action = "error"
            )
            extension.status <- extension$status
        }
        round.rows[[length(round.rows) + 1L]] <- data.frame(
            round = expansions.completed,
            eta.min = min(eta.grid),
            eta.max = max(eta.grid),
            n.candidates = length(eta.grid),
            best.idx = selected$gcv.result$best.idx,
            eta.optimal = selected$gcv.result$eta.optimal,
            gcv.optimal = selected$gcv.result$gcv.min,
            lower.endpoint.selected = lower.selected,
            extension.status = extension.status,
            stringsAsFactors = FALSE
        )
        if (!lower.selected) {
            search.status <- "interior_optimum"
            break
        }
        if (!extension$admitted) {
            search.status <- extension$status
            break
        }
        eta.grid <- extension$eta.grid
        expansions.completed <- expansions.completed + 1L
    }

    list(
        eta.grid = eta.grid,
        filter.weights.matrix = selected$filter.weights.matrix,
        gcv.result = selected$gcv.result,
        search = list(
            mode = eta.search,
            status = search.status,
            expansions.completed = expansions.completed,
            settings = settings,
            rounds = do.call(rbind, round.rows)
        )
    )
}

#' Fit Metric-Conductance Graph Low-Pass Regression
#'
#' Fits graph-spectral low-pass regression on a supplied graph by transforming
#' metric edge lengths into conductances and smoothing the response in the
#' eigenbasis of the weighted graph Laplacian.
#'
#' @inheritParams metric.graph.lowpass.operator
#' @param y Numeric response vector of length \code{length(adj.list)}.
#' @param n.eigenpairs Positive integer number of eigenpairs to compute.
#' @param filter.type Spectral low-pass filter family.
#' @param eta.grid Optional numeric eta grid. Values must be positive except
#'   for \code{filter.type = "heat_kernel"}, where \code{eta = 0} is allowed
#'   and gives the identity/no-smoothing limit. If \code{NULL}, the existing
#'   package helper \code{generate.eta.grid()} is used.
#' @param n.candidates Number of eta candidates when \code{eta.grid = NULL}.
#' @param eta.search Heat-time search controller. \code{"fixed"} evaluates the
#'   supplied or generated grid once. \code{"guarded.gcv"} is available only
#'   for \code{filter.type = "heat_kernel"} and proposes one lower positive
#'   time whenever GCV selects the current lower endpoint.
#' @param eta.expansion.factor Geometric divisor used by guarded GCV lower-time
#'   expansion.
#' @param eta.max.expansions Maximum number of guarded GCV lower-time
#'   extensions.
#' @param eta.identity.departure Smallest allowed departure from the identity
#'   at the largest retained graph-Laplacian eigenvalue.
#' @param eta.truncation.tol Omitted-mode attenuation tolerance used to certify
#'   guarded proposals for a truncated basis.
#' @param eigen.solver \code{"auto"}, \code{"sparse"}, or \code{"dense"}.
#'   \code{"auto"} uses dense decomposition only for
#'   \code{n <= dense.eigen.threshold}, then sparse-first.
#' @param dense.eigen.threshold Exact dense threshold for auto mode. Default
#'   \code{200L}, intended for small reference/testing problems.
#' @param dense.fallback.threshold Maximum graph size for emergency dense
#'   fallback when sparse decomposition fails and fallback is allowed.
#' @param dense.fallback \code{"auto"}, \code{"never"}, or \code{"always"}.
#'
#' @return A list of class \code{"metric.graph.lowpass.fit"}.
#'
#' @references
#' Gajer, P. and Ravel, J. (2025). Adaptive geometric regression for
#' high-dimensional structured data. arXiv:2511.03817.
#' \doi{10.48550/arXiv.2511.03817}
#'
#' @examples
#' adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
#' lengths <- list(1, c(1, 1), c(1, 1), 1)
#' fit.metric.graph.lowpass(
#'   adj, lengths, y = c(0, 0.2, 1.8, 2), n.eigenpairs = 4L,
#'   eta.grid = c(0.1, 1), eigen.solver = "dense"
#' )
#' @export
fit.metric.graph.lowpass <- function(
    adj.list,
    weight.list,
    y,
    conductance.rule = c("inverse.length.power", "exp.length",
                         "exp.length.squared", "self.tuned.gaussian"),
    conductance.epsilon = 1e-8,
    conductance.alpha = 1,
    conductance.sigma = NULL,
    conductance.sigma.rule = c("edge.quantile", "median", "local.k"),
    conductance.sigma.quantile = 0.75,
    conductance.local.k = 5L,
    laplacian.type = c("unnormalized", "symmetric.normalized"),
    n.eigenpairs = 50L,
    filter.type = c("heat_kernel", "tikhonov", "cubic_spline",
                    "gaussian", "exponential", "butterworth"),
    eta.grid = NULL,
    n.candidates = 40L,
    eigen.solver = c("auto", "sparse", "dense"),
    dense.eigen.threshold = 200L,
    dense.fallback.threshold = 5000L,
    dense.fallback = c("auto", "never", "always"),
    verbose = FALSE,
    eta.search = c("fixed", "guarded.gcv"),
    eta.expansion.factor = 3,
    eta.max.expansions = 3L,
    eta.identity.departure = 0.01,
    eta.truncation.tol = 1e-4) {

    basis <- metric.graph.lowpass.basis(
        adj.list = adj.list,
        weight.list = weight.list,
        conductance.rule = conductance.rule,
        conductance.epsilon = conductance.epsilon,
        conductance.alpha = conductance.alpha,
        conductance.sigma = conductance.sigma,
        conductance.sigma.rule = conductance.sigma.rule,
        conductance.sigma.quantile = conductance.sigma.quantile,
        conductance.local.k = conductance.local.k,
        laplacian.type = laplacian.type,
        n.eigenpairs = n.eigenpairs,
        eigen.solver = eigen.solver,
        dense.eigen.threshold = dense.eigen.threshold,
        dense.fallback.threshold = dense.fallback.threshold,
        dense.fallback = dense.fallback,
        verbose = verbose
    )
    n <- basis$spectral$n.vertices
    y <- .validate.metric.graph.lowpass.response(y, n, "y")
    filter.type <- match.arg(filter.type)
    n.candidates <- .validate.positive.integer.scalar(n.candidates, "n.candidates")
    eta.search <- match.arg(eta.search)

    operator <- basis$operator
    spectral <- basis$spectral
    eigenvalues <- spectral$eigenvalues
    V <- spectral$eigenvectors

    eta.grid <- .prepare.metric.graph.lowpass.eta.grid(
        eta.grid = eta.grid,
        eigenvalues = eigenvalues,
        filter.type = filter.type,
        n.candidates = n.candidates
    )
    y.spectral <- as.vector(crossprod(V, y))
    .require.metric.graph.lowpass.finite(y.spectral, "spectral response coefficients")
    eta.selection <- .fit.metric.graph.lowpass.eta.search(
        basis = basis,
        y = y,
        y.spectral = y.spectral,
        V = V,
        filter.type = filter.type,
        eta.grid = eta.grid,
        eta.search = eta.search,
        eta.expansion.factor = eta.expansion.factor,
        eta.max.expansions = eta.max.expansions,
        eta.identity.departure = eta.identity.departure,
        eta.truncation.tol = eta.truncation.tol
    )
    eta.grid <- eta.selection$eta.grid
    filter.weights.matrix <- eta.selection$filter.weights.matrix
    gcv.result <- eta.selection$gcv.result
    .validate.metric.graph.lowpass.gcv.result(
        gcv.result = gcv.result,
        n = length(y),
        n.eta = length(eta.grid),
        context = "fit.metric.graph.lowpass()"
    )
    best.idx <- gcv.result$best.idx

    spectral$filtered.eigenvalues <- filter.weights.matrix[, best.idx]
    .require.metric.graph.lowpass.finite(spectral$filtered.eigenvalues,
                                         "spectral$filtered.eigenvalues")
    spectral$eta.optimal <- gcv.result$eta.optimal
    spectral$filter.type <- filter.type
    spectral$n.eigenpairs <- length(eigenvalues)

    fitted.grid <- V %*% (y.spectral * filter.weights.matrix)
    .require.metric.graph.lowpass.finite(fitted.grid, "candidate fitted values")
    rss.grid <- colSums((y - fitted.grid)^2)
    df.grid <- colSums(filter.weights.matrix)
    gcv.scores <- rss.grid / pmax(length(y) - df.grid, 1e-10)^2
    .require.metric.graph.lowpass.finite(gcv.scores, "GCV scores")

    residuals <- y - gcv.result$y.hat
    .require.metric.graph.lowpass.finite(gcv.result$y.hat, "fitted values")
    .require.metric.graph.lowpass.finite(residuals, "residuals")
    result <- list(
        fitted.values = as.vector(gcv.result$y.hat),
        residuals = as.vector(residuals),
        y = y,
        graph = operator$graph,
        operator = operator,
        conductance = operator$conductance,
        laplacian = operator$laplacian,
        laplacian.type = operator$laplacian.type,
        spectral = spectral,
        gcv = list(
            eta.grid = eta.grid,
            gcv.scores = as.vector(gcv.scores),
            eta.optimal = gcv.result$eta.optimal,
            gcv.optimal = gcv.result$gcv.min,
            effective.df = gcv.result$effective.df,
            best.idx = best.idx,
            search = eta.selection$search
        ),
        parameters = c(
            basis$parameters,
            list(
                filter.type = filter.type,
                eta.search = eta.search,
                eta.expansion.factor =
                    eta.selection$search$settings$expansion.factor,
                eta.max.expansions =
                    eta.selection$search$settings$max.expansions,
                eta.identity.departure =
                    eta.selection$search$settings$identity.departure,
                eta.truncation.tol =
                    eta.selection$search$settings$truncation.tol
            )
        ),
        timing = NULL
    )
    attr(result, "call") <- match.call()
    class(result) <- c("metric.graph.lowpass.fit", "list")
    result
}

#' Refit Metric-Conductance Graph Low-Pass Regression
#'
#' Reuses a fitted metric graph low-pass eigensystem to smooth new responses.
#'
#' @param fitted.model A \code{"metric.graph.lowpass.fit"} object.
#' @param y.new Numeric vector or matrix with one row per graph vertex.
#' @param per.column.gcv Logical. If \code{TRUE}, select eta independently for
#'   each response column using the cached eigenbasis.
#' @param eta.grid Optional positive numeric eta grid for per-column GCV.
#' @param n.candidates Number of eta candidates when \code{eta.grid = NULL}.
#' @param n.cores Number of cores for per-column GCV. Phase 1 uses sequential
#'   processing if optional parallel packages are unavailable.
#' @param block.size Optional block size for fixed-eta multi-column refits.
#' @param verbose Logical progress flag.
#'
#' @return A list of class \code{"metric.graph.lowpass.refit"}.
#' @examples
#' adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
#' lengths <- list(1, c(1, 1), c(1, 1), 1)
#' fit <- fit.metric.graph.lowpass(
#'   adj, lengths, y = 1:4, n.eigenpairs = 4L,
#'   eta.grid = c(0.1, 1), eigen.solver = "dense"
#' )
#' refit.metric.graph.lowpass(fit, y.new = 4:1)
#' @export
refit.metric.graph.lowpass <- function(fitted.model,
                                       y.new,
                                       per.column.gcv = FALSE,
                                       eta.grid = NULL,
                                       n.candidates = 40L,
                                       n.cores = 1L,
                                       block.size = NULL,
                                       verbose = FALSE) {
    if (!inherits(fitted.model, "metric.graph.lowpass.fit")) {
        stop("fitted.model must be a 'metric.graph.lowpass.fit' object.")
    }
    spectral <- fitted.model$spectral
    V <- spectral$eigenvectors
    eigenvalues <- spectral$eigenvalues
    if (is.null(V) || !is.matrix(V)) stop("fitted.model$spectral$eigenvectors must be a matrix.")
    if (is.null(eigenvalues) || !is.numeric(eigenvalues)) {
        stop("fitted.model$spectral$eigenvalues must be numeric.")
    }
    .require.metric.graph.lowpass.finite(V, "fitted.model$spectral$eigenvectors")
    .require.metric.graph.lowpass.finite(eigenvalues, "fitted.model$spectral$eigenvalues")

    n <- nrow(V)
    y.info <- .prepare.metric.graph.lowpass.response.matrix(y.new, n)
    Y <- y.info$Y
    n.responses <- ncol(Y)
    col.names <- y.info$col.names

    n.candidates <- .validate.positive.integer.scalar(n.candidates, "n.candidates")
    n.cores <- .validate.positive.integer.scalar(n.cores, "n.cores")
    block.size <- .validate.optional.block.size(block.size)
    filter.type <- spectral$filter.type
    if (is.null(filter.type)) filter.type <- "heat_kernel"

    if (isTRUE(per.column.gcv)) {
        eta.grid <- .prepare.metric.graph.lowpass.eta.grid(
            eta.grid = eta.grid,
            eigenvalues = eigenvalues,
            filter.type = filter.type,
            n.candidates = n.candidates
        )
        filter.weights.matrix <- compute.filter.weights.matrix(eigenvalues, eta.grid, filter.type)
        .require.metric.graph.lowpass.finite(filter.weights.matrix, "refit filter weights")
        trace.S.all <- colSums(filter.weights.matrix)
        .require.metric.graph.lowpass.finite(trace.S.all, "refit effective degrees of freedom grid")
        Vt.Y <- crossprod(V, Y)
        .require.metric.graph.lowpass.finite(Vt.Y, "refit spectral response coefficients")

        Y.hat <- matrix(0, nrow = n, ncol = n.responses)
        eta.optimal <- numeric(n.responses)
        gcv.scores <- numeric(n.responses)
        effective.df <- numeric(n.responses)
        best.idx <- integer(n.responses)

        if (verbose && n.responses > 1L) {
            message(sprintf("Selecting eta via GCV for %d response(s).", n.responses))
        }
        for (j in seq_len(n.responses)) {
            gcv.result <- .select.eta.gcv.single(
                Y[, j], Vt.Y[, j], V, filter.weights.matrix, eta.grid
            )
            .validate.metric.graph.lowpass.gcv.result(
                gcv.result = gcv.result,
                n = n,
                n.eta = length(eta.grid),
                context = "refit.metric.graph.lowpass(per.column.gcv = TRUE)"
            )
            Y.hat[, j] <- gcv.result$y.hat
            eta.optimal[j] <- gcv.result$eta.optimal
            gcv.scores[j] <- gcv.result$gcv.min
            effective.df[j] <- gcv.result$effective.df
            best.idx[j] <- gcv.result$best.idx
        }
        .require.metric.graph.lowpass.finite(Y.hat, "refit fitted values")
        .require.metric.graph.lowpass.finite(eta.optimal, "refit eta.optimal")
        .require.metric.graph.lowpass.finite(gcv.scores, "refit GCV scores")
        .require.metric.graph.lowpass.finite(effective.df, "refit effective degrees of freedom")
        if (!is.null(col.names)) {
            colnames(Y.hat) <- col.names
            names(eta.optimal) <- col.names
            names(gcv.scores) <- col.names
            names(effective.df) <- col.names
        }
        residuals <- Y - Y.hat
        .require.metric.graph.lowpass.finite(residuals, "refit residuals")
        out <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = TRUE,
            filter.type = filter.type,
            eta.optimal = if (n.responses == 1L) eta.optimal[1] else eta.optimal,
            gcv.scores = if (n.responses == 1L) gcv.scores[1] else gcv.scores,
            effective.df = if (n.responses == 1L) effective.df[1] else effective.df,
            eta.grid = eta.grid,
            best.idx = if (n.responses == 1L) best.idx[1] else best.idx,
            n.cores.used = 1L
        )
        class(out) <- c("metric.graph.lowpass.refit", "list")
        return(out)
    }

    f.lambda <- spectral$filtered.eigenvalues
    eta.used <- spectral$eta.optimal
    if (is.null(f.lambda) || length(f.lambda) != ncol(V)) {
        stop("fitted.model$spectral$filtered.eigenvalues is missing or has the wrong length.")
    }
    .require.metric.graph.lowpass.finite(f.lambda,
                                         "fitted.model$spectral$filtered.eigenvalues")
    .require.metric.graph.lowpass.finite(eta.used, "fitted.model$spectral$eta.optimal")
    block.index <- .make.metric.graph.lowpass.block.index(n.responses, block.size)
    Y.hat <- matrix(0, nrow = n, ncol = n.responses)
    residuals <- matrix(0, nrow = n, ncol = n.responses)

    for (cols in block.index) {
        y.block <- Y[, cols, drop = FALSE]
        Vt.Y.block <- crossprod(V, y.block)
        Y.hat.block <- V %*% (f.lambda * Vt.Y.block)
        Y.hat[, cols] <- Y.hat.block
        residuals[, cols] <- y.block - Y.hat.block
    }
    .require.metric.graph.lowpass.finite(Y.hat, "refit fitted values")
    .require.metric.graph.lowpass.finite(residuals, "refit residuals")
    if (!is.null(col.names)) {
        colnames(Y.hat) <- col.names
        colnames(residuals) <- col.names
    }
    out <- list(
        fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
        residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
        n.responses = n.responses,
        per.column.gcv = FALSE,
        eta.used = eta.used,
        block.size.used = if (length(block.index) > 1L) block.size else NA_integer_
    )
    class(out) <- c("metric.graph.lowpass.refit", "list")
    out
}

.validate.metric.graph.lowpass.graph <- function(adj.list, weight.list) {
    if (!is.list(adj.list)) stop("adj.list must be a list.")
    if (!is.list(weight.list)) stop("weight.list must be a list.")
    n <- length(adj.list)
    if (n < 2L) stop("adj.list must contain at least two vertices.")
    if (length(weight.list) != n) stop("adj.list and weight.list must have the same length.")

    adj.norm <- vector("list", n)
    weight.norm <- vector("list", n)
    tol <- 1e-10
    for (i in seq_len(n)) {
        nbrs <- adj.list[[i]]
        wts <- weight.list[[i]]
        if (!is.numeric(nbrs) && !is.integer(nbrs)) stop(sprintf("adj.list[[%d]] must be numeric/integer.", i))
        if (!is.numeric(wts)) stop(sprintf("weight.list[[%d]] must be numeric.", i))
        nbrs <- as.integer(nbrs)
        wts <- as.double(wts)
        if (length(nbrs) != length(wts)) {
            stop(sprintf("Length mismatch at vertex %d: adj.list=%d, weight.list=%d.",
                         i, length(nbrs), length(wts)))
        }
        if (anyNA(nbrs)) stop(sprintf("adj.list[[%d]] contains NA.", i))
        if (any(nbrs < 1L | nbrs > n)) stop(sprintf("adj.list[[%d]] has indices outside 1..n.", i))
        if (any(nbrs == i)) stop(sprintf("adj.list[[%d]] contains self-loops.", i))
        if (anyDuplicated(nbrs)) stop(sprintf("adj.list[[%d]] contains duplicate neighbors.", i))
        if (any(!is.finite(wts))) stop(sprintf("weight.list[[%d]] contains non-finite values.", i))
        if (any(wts <= 0)) stop(sprintf("weight.list[[%d]] contains non-positive values.", i))
        adj.norm[[i]] <- nbrs
        weight.norm[[i]] <- wts
    }
    for (i in seq_len(n)) {
        for (idx in seq_along(adj.norm[[i]])) {
            j <- adj.norm[[i]][idx]
            rev.idx <- match(i, adj.norm[[j]])
            if (is.na(rev.idx)) {
                stop(sprintf("Graph must be undirected: edge %d -> %d has no reciprocal entry.", i, j))
            }
            w.ij <- weight.norm[[i]][idx]
            w.ji <- weight.norm[[j]][rev.idx]
            if (abs(w.ij - w.ji) > tol * max(1, abs(w.ij), abs(w.ji))) {
                stop(sprintf("Reciprocal edge weights mismatch for (%d, %d): %.12g vs %.12g.",
                             i, j, w.ij, w.ji))
            }
        }
    }
    list(
        adj.list = adj.norm,
        weight.list = weight.norm,
        adj.list.0based = lapply(adj.norm, function(v) as.integer(v - 1L)),
        weight.list.cpp = lapply(weight.norm, as.double)
    )
}

.validate.metric.graph.lowpass.operator.args <- function(conductance.rule,
                                                         conductance.epsilon,
                                                         conductance.alpha,
                                                         conductance.sigma,
                                                         conductance.sigma.rule,
                                                         conductance.sigma.quantile,
                                                         conductance.local.k,
                                                         laplacian.type,
                                                         verbose) {
    conductance.rule <- match.arg(
        conductance.rule,
        choices = c("inverse.length.power", "exp.length",
                    "exp.length.squared", "self.tuned.gaussian")
    )
    conductance.sigma.rule <- match.arg(
        conductance.sigma.rule,
        choices = c("edge.quantile", "median", "local.k")
    )
    laplacian.type <- match.arg(laplacian.type, choices = c("unnormalized", "symmetric.normalized"))
    if (!is.numeric(conductance.epsilon) || length(conductance.epsilon) != 1L ||
        !is.finite(conductance.epsilon) || conductance.epsilon <= 0) {
        stop("conductance.epsilon must be a finite positive numeric scalar.")
    }
    if (!is.numeric(conductance.alpha) || length(conductance.alpha) != 1L ||
        !is.finite(conductance.alpha) || conductance.alpha <= 0) {
        stop("conductance.alpha must be a finite positive numeric scalar.")
    }
    if (is.null(conductance.sigma)) {
        conductance.sigma <- NaN
    } else if (!is.numeric(conductance.sigma) || length(conductance.sigma) != 1L ||
               !is.finite(conductance.sigma) || conductance.sigma <= 0) {
        stop("conductance.sigma must be NULL or a finite positive numeric scalar.")
    }
    if (!is.numeric(conductance.sigma.quantile) || length(conductance.sigma.quantile) != 1L ||
        !is.finite(conductance.sigma.quantile) ||
        conductance.sigma.quantile < 0 || conductance.sigma.quantile > 1) {
        stop("conductance.sigma.quantile must be in [0, 1].")
    }
    conductance.local.k <- .validate.positive.integer.scalar(conductance.local.k, "conductance.local.k")
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.")
    list(
        conductance.rule = conductance.rule,
        conductance.epsilon = as.double(conductance.epsilon),
        conductance.alpha = as.double(conductance.alpha),
        conductance.sigma = as.double(conductance.sigma),
        conductance.sigma.rule = conductance.sigma.rule,
        conductance.sigma.quantile = as.double(conductance.sigma.quantile),
        conductance.local.k = as.integer(conductance.local.k),
        laplacian.type = laplacian.type,
        verbose = as.logical(verbose)
    )
}

.validate.metric.graph.lowpass.solver.args <- function(n.eigenpairs,
                                                       n,
                                                       eigen.solver,
                                                       dense.eigen.threshold,
                                                       dense.fallback.threshold,
                                                       dense.fallback) {
    n.eigenpairs <- .validate.positive.integer.scalar(n.eigenpairs, "n.eigenpairs")
    if (n.eigenpairs > n) {
        warning("n.eigenpairs exceeds number of vertices; using n.", call. = FALSE)
        n.eigenpairs <- n
    }
    eigen.solver <- match.arg(eigen.solver, choices = c("auto", "sparse", "dense"))
    dense.fallback <- match.arg(dense.fallback, choices = c("auto", "never", "always"))
    dense.eigen.threshold <- .validate.positive.integer.scalar(dense.eigen.threshold, "dense.eigen.threshold")
    dense.fallback.threshold <- .validate.positive.integer.scalar(dense.fallback.threshold, "dense.fallback.threshold")
    list(
        n.eigenpairs = as.integer(n.eigenpairs),
        eigen.solver = eigen.solver,
        dense.eigen.threshold = as.integer(dense.eigen.threshold),
        dense.fallback.threshold = as.integer(dense.fallback.threshold),
        dense.fallback = dense.fallback
    )
}

.validate.positive.integer.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 1 || x != floor(x)) {
        stop(sprintf("%s must be a positive integer scalar.", name))
    }
    as.integer(x)
}

.validate.nonnegative.integer.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 0 || x != floor(x)) {
        stop(sprintf("%s must be a nonnegative integer scalar.", name))
    }
    as.integer(x)
}

.validate.logical.scalar <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
    }
    x
}

.validate.unit.interval.open <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x <= 0 || x >= 1) {
        stop(
            sprintf("%s must be a finite numeric scalar strictly between 0 and 1.", name),
            call. = FALSE
        )
    }
    as.double(x)
}

.validate.metric.graph.lowpass.basis <- function(basis) {
    if (!inherits(basis, "metric.graph.lowpass.basis")) {
        stop("basis must be a 'metric.graph.lowpass.basis' object.", call. = FALSE)
    }
    spectral <- basis$spectral
    if (is.null(spectral$eigenvalues) || !is.numeric(spectral$eigenvalues) ||
        is.null(spectral$eigenvectors) || !is.matrix(spectral$eigenvectors) ||
        ncol(spectral$eigenvectors) != length(spectral$eigenvalues) ||
        nrow(spectral$eigenvectors) != spectral$n.vertices) {
        stop("basis contains an invalid eigensystem.", call. = FALSE)
    }
    .require.metric.graph.lowpass.finite(
        spectral$eigenvalues, "basis$spectral$eigenvalues"
    )
    .require.metric.graph.lowpass.finite(
        spectral$eigenvectors, "basis$spectral$eigenvectors"
    )
    invisible(TRUE)
}

.validate.metric.graph.lowpass.response <- function(y, n, name) {
    if (!is.numeric(y) || length(y) != n) {
        stop(sprintf("%s must be a numeric vector of length %d.", name, n))
    }
    y <- as.double(y)
    if (any(!is.finite(y))) stop(sprintf("%s cannot contain NA/NaN/Inf.", name))
    y
}

.prepare.metric.graph.lowpass.response.matrix <- function(y.new, n) {
    is.matrix.input <- is.matrix(y.new) || inherits(y.new, "Matrix")
    if (is.matrix.input) {
        if (nrow(y.new) != n) stop(sprintf("nrow(y.new) must be %d.", n))
        col.names <- colnames(y.new)
        Y <- if (inherits(y.new, "Matrix")) as.matrix(y.new) else y.new
        if (!is.numeric(Y)) stop("y.new must be numeric.")
    } else {
        Y <- matrix(.validate.metric.graph.lowpass.response(y.new, n, "y.new"), ncol = 1)
        col.names <- NULL
    }
    storage.mode(Y) <- "double"
    if (any(!is.finite(Y))) stop("y.new cannot contain NA/NaN/Inf.")
    list(Y = Y, col.names = col.names)
}

.prepare.metric.graph.lowpass.eta.grid <- function(eta.grid, eigenvalues, filter.type, n.candidates) {
    if (is.null(eta.grid)) {
        eta.grid <- generate.eta.grid(eigenvalues, filter.type, n.candidates)
    } else {
        allow.zero <- identical(filter.type, "heat_kernel")
        bad <- !is.numeric(eta.grid) || length(eta.grid) < 1L ||
            any(!is.finite(eta.grid)) ||
            if (allow.zero) any(eta.grid < 0) else any(eta.grid <= 0)
        if (bad) {
            if (allow.zero) {
                stop("eta.grid must be NULL or a finite nonnegative numeric vector for heat_kernel.")
            }
            stop("eta.grid must be NULL or a positive finite numeric vector.")
        }
        eta.grid <- as.double(eta.grid)
    }
    .require.metric.graph.lowpass.finite(eta.grid, "eta.grid")
    eta.grid
}

.require.metric.graph.lowpass.finite <- function(x, name) {
    if (is.null(x) || length(x) < 1L || any(!is.finite(as.numeric(x)))) {
        stop(sprintf("%s contains non-finite values.", name), call. = FALSE)
    }
    invisible(TRUE)
}

.validate.metric.graph.lowpass.gcv.result <- function(gcv.result, n, n.eta, context) {
    if (!is.list(gcv.result)) {
        stop(sprintf("%s returned an invalid GCV result.", context), call. = FALSE)
    }
    if (is.null(gcv.result$y.hat) || length(gcv.result$y.hat) != n) {
        stop(sprintf("%s returned fitted values with the wrong length.", context),
             call. = FALSE)
    }
    .require.metric.graph.lowpass.finite(gcv.result$y.hat,
                                         sprintf("%s fitted values", context))
    .require.metric.graph.lowpass.finite(gcv.result$eta.optimal,
                                         sprintf("%s selected eta", context))
    .require.metric.graph.lowpass.finite(gcv.result$gcv.min,
                                         sprintf("%s selected GCV score", context))
    .require.metric.graph.lowpass.finite(gcv.result$effective.df,
                                         sprintf("%s effective degrees of freedom", context))
    best.idx <- gcv.result$best.idx
    if (length(best.idx) != 1L || is.na(best.idx) || !is.finite(best.idx) ||
        best.idx != floor(best.idx) || best.idx < 1L || best.idx > n.eta) {
        stop(sprintf("%s selected an invalid eta index.", context), call. = FALSE)
    }
    invisible(TRUE)
}

.attach.metric.graph.lowpass.laplacian <- function(out, return.sparse) {
    if (isTRUE(return.sparse)) {
        if (requireNamespace("Matrix", quietly = TRUE)) {
            trip <- out$laplacian
            out$laplacian$matrix <- Matrix::sparseMatrix(
                i = trip$i + 1L,
                j = trip$j + 1L,
                x = trip$x,
                dims = trip$dim,
                giveCsparse = TRUE
            )
        } else {
            warning("Matrix is unavailable; returning Laplacian triplets only.", call. = FALSE)
        }
    }
    out
}

.validate.optional.block.size <- function(block.size) {
    if (is.null(block.size)) return(NA_integer_)
    .validate.positive.integer.scalar(block.size, "block.size")
}

.make.metric.graph.lowpass.block.index <- function(p, block.size) {
    if (is.na(block.size) || p <= 1L) return(list(seq_len(p)))
    split(seq_len(p), ceiling(seq_len(p) / block.size))
}

generate.eta.grid <- function(eigenvalues, filter.type, n.candidates = 40L) {
    n.candidates <- max(5L, as.integer(n.candidates))
    ev <- sort(as.numeric(eigenvalues))
    ev <- ev[is.finite(ev) & ev > .Machine$double.eps]
    if (length(ev) == 0L) {
        return(seq(1e-3, 1, length.out = n.candidates))
    }

    lo <- 1 / max(ev)
    hi <- 1 / min(ev)
    if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi <= lo) {
        return(seq(1e-3, 1, length.out = n.candidates))
    }

    exp(seq(log(max(lo, 1e-6)), log(hi), length.out = n.candidates))
}

compute.filter.weights.matrix <- function(eigenvalues, eta.grid, filter.type) {
    m <- length(eigenvalues)
    n.eta <- length(eta.grid)
    weights <- matrix(0, nrow = m, ncol = n.eta)

    for (k in seq_len(n.eta)) {
        eta <- eta.grid[k]
        weights[, k] <- switch(filter.type,
            "heat_kernel" = exp(-eta * eigenvalues),
            "tikhonov" = 1.0 / (1.0 + eta * eigenvalues),
            "cubic_spline" = 1.0 / (1.0 + eta * eigenvalues^2),
            "gaussian" = exp(-eta * eigenvalues^2),
            "exponential" = exp(-eta * sqrt(pmax(eigenvalues, 0))),
            "butterworth" = {
                x <- eigenvalues / eta
                1.0 / (1.0 + x^4)
            },
            exp(-eta * eigenvalues)
        )
    }

    weights
}

.select.eta.gcv.single <- function(y.obs, y.spectral, V,
                                  filter.weights.matrix, eta.grid) {
    n <- length(y.obs)
    filtered.spectral <- y.spectral * filter.weights.matrix
    Y.hat.all <- V %*% filtered.spectral
    residuals.all <- y.obs - Y.hat.all
    rss <- colSums(residuals.all^2)
    trace.S <- colSums(filter.weights.matrix)
    denom <- pmax(n - trace.S, 1e-10)
    gcv.scores <- rss / (denom^2)
    best.idx <- which.min(gcv.scores)

    list(
        y.hat = Y.hat.all[, best.idx],
        eta.optimal = eta.grid[best.idx],
        gcv.min = gcv.scores[best.idx],
        effective.df = trace.S[best.idx],
        best.idx = best.idx
    )
}

#' @method print metric.graph.lowpass.fit
#' @export
print.metric.graph.lowpass.fit <- function(x, ...) {
    cat("\nMetric-Conductance Graph Low-Pass Fit\n")
    cat("====================================\n\n")
    cat(sprintf("Vertices: %d; edges: %d\n", x$graph$n.vertices, x$graph$n.edges))
    cat(sprintf("Conductance rule: %s\n", x$conductance$rule))
    cat(sprintf("Laplacian: %s\n", x$laplacian.type))
    cat(sprintf("Filter: %s; eta: %.4e; GCV: %.4e; effective df: %.2f\n",
                x$spectral$filter.type, x$spectral$eta.optimal,
                x$gcv$gcv.optimal, x$gcv$effective.df))
    invisible(x)
}

#' @method print metric.graph.lowpass.refit
#' @export
print.metric.graph.lowpass.refit <- function(x, ...) {
    cat("\nMetric-Conductance Graph Low-Pass Refit\n")
    cat("======================================\n\n")
    cat(sprintf("Responses: %d\n", x$n.responses))
    if (isTRUE(x$per.column.gcv)) {
        cat(sprintf("Per-column GCV: yes; filter: %s\n", x$filter.type))
    } else {
        cat(sprintf("Per-column GCV: no; eta used: %.4e\n", x$eta.used))
    }
    invisible(x)
}
