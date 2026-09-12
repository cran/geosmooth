.local.chart.prepare.X.eval <- function(X, X.eval = NULL) {
    X <- .state.density.validate.X(X)
    if (is.null(X.eval)) {
        X.eval <- X
    } else {
        X.eval <- .state.density.validate.X(X.eval)
        if (ncol(X.eval) != ncol(X)) {
            stop("X.eval must have the same number of columns as X.",
                 call. = FALSE)
        }
    }
    list(X = X, X.eval = X.eval)
}

.local.chart.validate.response <- function(y, n, name = "y",
                                           nonnegative = FALSE) {
    if (!is.numeric(y) || length(y) != n || any(!is.finite(y))) {
        stop(name, " must be a finite numeric vector with length nrow(X).",
             call. = FALSE)
    }
    y <- as.numeric(y)
    if (isTRUE(nonnegative) && any(y < 0)) {
        stop(name, " must be nonnegative.", call. = FALSE)
    }
    y
}

.local.chart.validate.quadrature <- function(quadrature.weights, n) {
    if (is.null(quadrature.weights)) {
        return(rep(1, n))
    }
    if (!is.numeric(quadrature.weights) ||
        length(quadrature.weights) != n ||
        any(!is.finite(quadrature.weights)) ||
        any(quadrature.weights <= 0)) {
        stop("quadrature.weights must be NULL or a positive finite numeric vector with length nrow(X).",
             call. = FALSE)
    }
    as.numeric(quadrature.weights)
}

.local.chart.validate.support.size <- function(support.size, n) {
    if (!is.numeric(support.size) || length(support.size) != 1L ||
        !is.finite(support.size)) {
        stop("support.size must be one finite integer.", call. = FALSE)
    }
    support.size <- as.integer(round(support.size))
    if (support.size < 1L) {
        stop("support.size must be at least 1.", call. = FALSE)
    }
    min(support.size, n)
}

.local.chart.validate.positive.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
        stop(name, " must be one positive finite numeric value.",
             call. = FALSE)
    }
    as.numeric(x)
}

.local.chart.validate.nonnegative.scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
        stop(name, " must be one finite nonnegative numeric value.",
             call. = FALSE)
    }
    as.numeric(x)
}

.local.chart.validate.nonnegative.integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
        stop(name, " must be one finite nonnegative integer.", call. = FALSE)
    }
    x <- as.integer(round(x))
    if (x < 0L) {
        stop(name, " must be nonnegative.", call. = FALSE)
    }
    x
}

.local.chart.validate.positive.integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
        stop(name, " must be one finite positive integer.", call. = FALSE)
    }
    x <- as.integer(round(x))
    if (x < 1L) {
        stop(name, " must be positive.", call. = FALSE)
    }
    x
}

.local.chart.validate.chart.dim <- function(chart.dim,
                                            coordinate.method,
                                            p,
                                            support.size) {
    if (identical(coordinate.method, "coordinates")) {
        if (!is.null(chart.dim)) {
            warning("chart.dim is ignored when coordinate.method = \"coordinates\".",
                    call. = FALSE)
        }
        return(p)
    }
    if (is.null(chart.dim)) {
        return(max(1L, min(p, support.size - 1L)))
    }
    if (!is.numeric(chart.dim) || length(chart.dim) != 1L ||
        !is.finite(chart.dim)) {
        stop("chart.dim must be NULL or one finite integer for local.pca.",
             call. = FALSE)
    }
    chart.dim <- as.integer(round(chart.dim))
    if (chart.dim < 1L) {
        stop("chart.dim must be at least 1.", call. = FALSE)
    }
    min(chart.dim, p, support.size)
}

.local.chart.requested.chart.dim.label <- function(chart.dim) {
    if (is.null(chart.dim)) {
        return(NA_character_)
    }
    if (is.character(chart.dim)) {
        return(as.character(chart.dim[[1L]]))
    }
    as.character(as.integer(round(chart.dim[[1L]])))
}

.local.chart.dimension.telemetry <- function(chart.dim.info,
                                             chart.dim.by.anchor = NULL,
                                             n.anchor = NULL,
                                             source.path = NULL) {
    if (is.null(chart.dim.by.anchor)) {
        if (is.null(n.anchor)) {
            n.anchor <- 1L
        }
        chart.dim.by.anchor <- rep(as.integer(chart.dim.info$chart.dim),
                                   as.integer(n.anchor))
    }
    chart.dim.by.anchor <- as.integer(chart.dim.by.anchor)
    if (!length(chart.dim.by.anchor) ||
        anyNA(chart.dim.by.anchor) ||
        any(chart.dim.by.anchor < 1L)) {
        stop("chart.dim.by.anchor must contain positive integer dimensions.",
             call. = FALSE)
    }
    list(
        requested = .local.chart.requested.chart.dim.label(
            chart.dim.info$requested.chart.dim
        ),
        mode = chart.dim.info$chart.dim.mode,
        resolved = chart.dim.info$chart.dim,
        by.anchor = chart.dim.by.anchor,
        summary = list(
            n.anchor = length(chart.dim.by.anchor),
            min = min(chart.dim.by.anchor),
            max = max(chart.dim.by.anchor),
            median = stats::median(chart.dim.by.anchor),
            n.unique = length(unique(chart.dim.by.anchor))
        ),
        auto = isTRUE(chart.dim.info$auto.chart.dim),
        local.auto = isTRUE(chart.dim.info$auto.chart.dim.local),
        auto.diagnostics = chart.dim.info$auto.chart.dim.diagnostics,
        support.metric = chart.dim.info$auto.chart.support.metric,
        selection.metric = chart.dim.info$auto.chart.selection.metric,
        source.path = source.path
    )
}

.local.chart.dim.mode <- function(chart.dim, coordinate.method) {
    .klp.chart.dim.mode(
        chart.dim = chart.dim,
        coordinate.method = coordinate.method
    )
}

.local.chart.resolve.chart.dim <- function(
    X,
    support.size,
    degree,
    coordinate.method,
    chart.dim,
    feasibility.contract = c("local_polynomial", "chart_kernel"),
    geometry.margin = 0L,
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator")) {

    feasibility.contract <- match.arg(feasibility.contract)
    auto.chart.support.metric <- match.arg(auto.chart.support.metric)
    auto.chart.selection.metric <- match.arg(auto.chart.selection.metric)
    p <- ncol(X)
    support.size <- .local.chart.validate.support.size(support.size, nrow(X))

    if (identical(coordinate.method, "coordinates")) {
        if (!is.null(chart.dim)) {
            warning("chart.dim is ignored when coordinate.method = \"coordinates\".",
                    call. = FALSE)
        }
        return(list(
            chart.dim = p,
            requested.chart.dim = chart.dim,
            chart.dim.mode = "ambient",
            auto.chart.dim = FALSE,
            auto.chart.dim.local = FALSE,
            auto.chart.dim.diagnostics = NULL,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        ))
    }

    if (is.null(chart.dim)) {
        dim <- max(1L, min(p, support.size - 1L))
        return(list(
            chart.dim = as.integer(dim),
            requested.chart.dim = chart.dim,
            chart.dim.mode = "ambient.default",
            auto.chart.dim = FALSE,
            auto.chart.dim.local = FALSE,
            auto.chart.dim.diagnostics = NULL,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        ))
    }

    if (identical(chart.dim, "auto") || identical(chart.dim, "local.auto")) {
        diagnostics <- .local.pca.auto.chart.dim.with.metric(
            X = X,
            support.size = support.size,
            degree = degree,
            feasibility.contract = feasibility.contract,
            geometry.margin = geometry.margin,
            max.anchors = if (identical(chart.dim, "local.auto")) {
                nrow(X)
            } else {
                60L
            },
            operator.support.metric = "coordinates",
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        dim <- max(1L, min(
            p,
            .local.pca.chart.dim.support.cap(
                n.support = support.size,
                degree = degree,
                ambient.dim = p,
                feasibility.contract = feasibility.contract,
                geometry.margin = geometry.margin
            ),
            diagnostics$chart.dim
        ))
        return(list(
            chart.dim = as.integer(dim),
            requested.chart.dim = chart.dim,
            chart.dim.mode = .local.chart.dim.mode(chart.dim,
                                                   coordinate.method),
            auto.chart.dim = TRUE,
            auto.chart.dim.local = identical(chart.dim, "local.auto"),
            auto.chart.dim.diagnostics = diagnostics,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        ))
    }

    if (!is.numeric(chart.dim) || length(chart.dim) != 1L ||
        !is.finite(chart.dim)) {
        stop("chart.dim must be NULL, one finite integer, \"auto\", ",
             "or \"local.auto\" for local.pca.",
             call. = FALSE)
    }
    dim <- as.integer(round(chart.dim))
    if (dim < 1L) {
        stop("chart.dim must be at least 1.", call. = FALSE)
    }
    dim <- min(
        dim,
        .local.pca.chart.dim.support.cap(
            n.support = support.size,
            degree = degree,
            ambient.dim = p,
            feasibility.contract = feasibility.contract,
            geometry.margin = geometry.margin
        )
    )
    list(
        chart.dim = as.integer(dim),
        requested.chart.dim = chart.dim,
        chart.dim.mode = "fixed",
        auto.chart.dim = FALSE,
        auto.chart.dim.local = FALSE,
        auto.chart.dim.diagnostics = NULL,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric
    )
}

.local.chart.resolve.eval.chart.dim <- function(
    X,
    x0,
    support.size,
    degree,
    coordinate.method,
    chart.dim,
    summary.dim,
    feasibility.contract = c("local_polynomial", "chart_kernel"),
    geometry.margin = 0L) {

    feasibility.contract <- match.arg(feasibility.contract)
    if (!identical(coordinate.method, "local.pca") ||
        !identical(chart.dim, "local.auto")) {
        return(as.integer(summary.dim))
    }
    ordered <- .klp.local.order(
        X.train = X,
        center = x0,
        support.size = support.size
    )
    dim <- .klp.local.auto.chart.dim.from.order(
        X.train = X,
        center = x0,
        ordered = ordered,
        support.size = support.size,
        degree = degree
    )
    if (identical(feasibility.contract, "chart_kernel")) {
        support <- .local.chart.support(X, x0, support.size)
        sv <- tryCatch(
            svd(support$centered, nu = 0L, nv = 0L)$d,
            error = function(e) numeric(0)
        )
        row <- .local.pca.auto.chart.dim.from.singular.values(
            sv = sv,
            n.support = nrow(support$centered),
            degree = 0L,
            ambient.dim = ncol(X),
            feasibility.contract = "chart_kernel",
            geometry.margin = geometry.margin
        )
        dim <- row$selected.local.dim[[1L]]
    }
    if (!is.finite(dim) || dim < 1L) {
        dim <- as.integer(summary.dim)
    }
    as.integer(max(1L, min(
        .local.pca.chart.dim.support.cap(
            n.support = support.size,
            degree = degree,
            ambient.dim = ncol(X),
            feasibility.contract = feasibility.contract,
            geometry.margin = geometry.margin
        ),
        dim
    )))
}

.chart.kernel.resolve.chart.dim <- function(
    X,
    support.size,
    coordinate.method,
    chart.dim,
    geometry.margin = 0L,
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator")) {

    .local.chart.resolve.chart.dim(
        X = X,
        support.size = support.size,
        degree = 0L,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        feasibility.contract = "chart_kernel",
        geometry.margin = geometry.margin,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric
    )
}

.chart.kernel.resolve.eval.chart.dim <- function(
    X,
    x0,
    support.size,
    coordinate.method,
    chart.dim,
    summary.dim,
    geometry.margin = 0L) {

    .local.chart.resolve.eval.chart.dim(
        X = X,
        x0 = x0,
        support.size = support.size,
        degree = 0L,
        coordinate.method = coordinate.method,
        chart.dim = chart.dim,
        summary.dim = summary.dim,
        feasibility.contract = "chart_kernel",
        geometry.margin = geometry.margin
    )
}

.local.chart.support <- function(X, x0, support.size) {
    d2 <- rowSums((X - matrix(x0, nrow = nrow(X), ncol = ncol(X),
                              byrow = TRUE))^2)
    ord <- order(d2, seq_along(d2))
    idx <- ord[seq_len(support.size)]
    list(
        idx = idx,
        centered = sweep(X[idx, , drop = FALSE], 2L, x0, "-"),
        ambient.distances = sqrt(d2[idx])
    )
}

.local.chart.coordinates <- function(centered, coordinate.method, chart.dim) {
    if (identical(coordinate.method, "coordinates")) {
        return(centered)
    }
    d <- min(chart.dim, ncol(centered), nrow(centered))
    if (d < 1L) {
        return(matrix(0, nrow = nrow(centered), ncol = 1L))
    }
    sv <- svd(centered, nu = 0L, nv = d)
    centered %*% sv$v[, seq_len(d), drop = FALSE]
}

.local.chart.kernel <- function(distances, kernel, bandwidth.multiplier = 1) {
    weights <- .klp.kernel.weights(
        distances = distances,
        kernel = kernel,
        bandwidth.multiplier = bandwidth.multiplier
    )
    h <- max(distances[is.finite(distances)], 0)
    if (!is.finite(h) || h <= 0) {
        h <- 1
    }
    list(
        weights = weights,
        bandwidth = bandwidth.multiplier * h,
        effective.support = sum(weights > 0)
    )
}

.local.chart.feature.matrix <- function(coords, degree) {
    degree <- .local.chart.validate.nonnegative.integer(degree, "degree")
    n <- nrow(coords)
    d <- ncol(coords)
    if (degree == 0L) {
        return(matrix(numeric(0), nrow = n, ncol = 0L))
    }
    out <- coords
    if (degree >= 2L) {
        quad <- vector("list", d * (d + 1L) / 2L)
        kk <- 0L
        for (aa in seq_len(d)) {
            for (bb in aa:d) {
                kk <- kk + 1L
                quad[[kk]] <- coords[, aa] * coords[, bb]
            }
        }
        out <- cbind(out, do.call(cbind, quad))
    }
    out
}
