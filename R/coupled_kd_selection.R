.coupled.kd.design.ncol <- function(chart.dim, degree) {
    as.integer(choose(as.integer(chart.dim) + as.integer(degree),
                      as.integer(degree)))
}

.coupled.kd.max.feasible.chart.dim <- function(support.size,
                                               degree,
                                               ambient.dim,
                                               design.margin = 2L) {
    support.size <- as.integer(support.size)
    degree <- as.integer(degree)
    ambient.dim <- as.integer(ambient.dim)
    design.margin <- as.integer(design.margin)
    if (!is.finite(support.size) || support.size < 1L ||
        !is.finite(degree) || degree < 0L ||
        !is.finite(ambient.dim) || ambient.dim < 1L ||
        !is.finite(design.margin) || design.margin < 0L) {
        stop("Invalid coupled k-d feasibility inputs.", call. = FALSE)
    }
    dims <- seq_len(ambient.dim)
    ok <- choose(dims + degree, degree) + design.margin <= support.size
    if (!any(ok)) {
        return(NA_integer_)
    }
    as.integer(max(dims[ok]))
}

.coupled.kd.resolve.auto.dim <- function(auto.chart.dim,
                                         support.size,
                                         degree) {
    if (is.null(auto.chart.dim)) {
        stop("'auto' chart-dimension seeds require 'auto.chart.dim'.",
             call. = FALSE)
    }
    if (is.function(auto.chart.dim)) {
        value <- auto.chart.dim(support.size = support.size, degree = degree)
    } else if (is.data.frame(auto.chart.dim)) {
        required <- c("support.size", "degree", "chart.dim")
        if (!all(required %in% names(auto.chart.dim))) {
            stop("'auto.chart.dim' data frames must contain support.size, ",
                 "degree, and chart.dim columns.", call. = FALSE)
        }
        hit <- auto.chart.dim[
            as.integer(auto.chart.dim$support.size) == as.integer(support.size) &
                as.integer(auto.chart.dim$degree) == as.integer(degree),
            ,
            drop = FALSE
        ]
        if (nrow(hit) != 1L) {
            stop("'auto.chart.dim' must provide exactly one seed for ",
                 "each support.size/degree pair.", call. = FALSE)
        }
        value <- hit$chart.dim[[1L]]
    } else if (!is.null(names(auto.chart.dim))) {
        key <- paste(as.integer(support.size), as.integer(degree), sep = ":")
        if (!key %in% names(auto.chart.dim)) {
            key <- as.character(as.integer(support.size))
        }
        if (!key %in% names(auto.chart.dim)) {
            stop("'auto.chart.dim' named vectors must contain support or ",
                 "support:degree keys.", call. = FALSE)
        }
        value <- auto.chart.dim[[key]]
    } else if (length(auto.chart.dim) == 1L) {
        value <- auto.chart.dim[[1L]]
    } else {
        stop("'auto.chart.dim' must be a scalar, named vector, data frame, ",
             "or function.", call. = FALSE)
    }
    value <- as.integer(value)
    if (length(value) != 1L || !is.finite(value) || value < 1L) {
        stop("'auto.chart.dim' resolved to an invalid dimension.",
             call. = FALSE)
    }
    value
}

.coupled.kd.source.rank <- function(source) {
    match(source, c("numeric", "auto_seed", "plateau_geometry", "guard",
                    "manual", "local_auto_policy"))
}

.coupled.kd.stage.rank <- function(stage) {
    match(stage, c("skeleton", "local_refine", "guard", "full_reference"))
}

.coupled.kd.reuse.key <- function(support.size,
                                  kernel,
                                  max.chart.dim,
                                  reuse.type) {
    if (identical(reuse.type, "weighted")) {
        paste(as.integer(support.size), as.character(kernel),
              as.integer(max.chart.dim), sep = "\r")
    } else {
        paste(as.integer(support.size), as.integer(max.chart.dim),
              sep = "\r")
    }
}

.coupled.kd.candidate.table <- function(support.grid,
                                        degree.grid,
                                        chart.dim.grid,
                                        ambient.dim,
                                        kernel.grid = "gaussian",
                                        bandwidth.multiplier.grid = 1,
                                        stage = "full_reference",
                                        chart.dim.max = NULL,
                                        design.margin = 2L,
                                        auto.chart.dim = NULL,
                                        reuse.type = c("weighted", "chart")) {
    reuse.type <- match.arg(reuse.type)
    stage <- match.arg(stage, c("skeleton", "local_refine", "guard",
                               "full_reference"))
    ambient.dim <- as.integer(ambient.dim)
    if (length(ambient.dim) != 1L || !is.finite(ambient.dim) ||
        ambient.dim < 1L) {
        stop("'ambient.dim' must be a positive integer scalar.", call. = FALSE)
    }
    design.margin <- as.integer(design.margin)
    if (length(design.margin) != 1L || !is.finite(design.margin) ||
        design.margin < 0L) {
        stop("'design.margin' must be a nonnegative integer scalar.",
             call. = FALSE)
    }
    support.grid <- .klp.clean.support.grid(support.grid, n = Inf)
    degree.grid <- .klp.clean.degree.grid(degree.grid)
    kernel.grid <- .klp.clean.kernel.grid(kernel.grid)
    bandwidth.multiplier.grid <- .klp.clean.bandwidth.multiplier.grid(
        bandwidth.multiplier.grid
    )
    chart.grid <- .local.chart.clean.chart.dim.grid(chart.dim.grid)
    if (is.null(chart.dim.max)) {
        decoded <- lapply(chart.grid$chart.dim, function(x) {
            tryCatch(.local.chart.decode.chart.dim(x), error = function(e) NA)
        })
        numeric.dim <- unlist(lapply(decoded, function(x) {
            if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
                as.integer(x)
            } else {
                NA_integer_
            }
        }), use.names = FALSE)
        chart.dim.max <- if (any(is.finite(numeric.dim))) {
            max(numeric.dim, na.rm = TRUE)
        } else {
            ambient.dim
        }
    }
    chart.dim.max <- as.integer(chart.dim.max)
    if (length(chart.dim.max) != 1L || !is.finite(chart.dim.max) ||
        chart.dim.max < 1L) {
        stop("'chart.dim.max' must be a positive integer scalar.",
             call. = FALSE)
    }

    rows <- vector("list", length(support.grid) * length(degree.grid) *
                       length(kernel.grid) *
                       length(bandwidth.multiplier.grid) *
                       nrow(chart.grid))
    rr <- 0L
    for (support.size in support.grid) {
        for (degree in degree.grid) {
            d.feasible.max <- .coupled.kd.max.feasible.chart.dim(
                support.size = support.size,
                degree = degree,
                ambient.dim = ambient.dim,
                design.margin = design.margin
            )
            d.hi <- suppressWarnings(min(chart.dim.max, ambient.dim,
                                         d.feasible.max, na.rm = TRUE))
            if (!is.finite(d.hi)) {
                d.hi <- NA_integer_
            }
            for (kernel in kernel.grid) {
                for (bandwidth.multiplier in bandwidth.multiplier.grid) {
                    for (ii in seq_len(nrow(chart.grid))) {
                        label <- chart.grid$chart.dim[[ii]]
                        decoded <- .local.chart.decode.chart.dim(label)
                        source <- "numeric"
                        raw <- NA_integer_
                        clipped <- NA_integer_
                        seed.clipped <- FALSE
                        skip.reason <- NA_character_
                        feasible <- TRUE
                        if (identical(decoded, "local.auto")) {
                            source <- "local_auto_policy"
                            feasible <- FALSE
                            skip.reason <- "local_auto_separate_policy"
                        } else if (identical(decoded, "auto")) {
                            source <- "auto_seed"
                            raw <- .coupled.kd.resolve.auto.dim(
                                auto.chart.dim = auto.chart.dim,
                                support.size = support.size,
                                degree = degree
                            )
                            clipped <- if (is.finite(d.hi)) {
                                as.integer(min(raw, d.hi))
                            } else {
                                NA_integer_
                            }
                            seed.clipped <- is.finite(clipped) &&
                                !identical(as.integer(raw),
                                           as.integer(clipped))
                        } else if (is.null(decoded)) {
                            source <- "manual"
                            raw <- ambient.dim
                            clipped <- raw
                        } else {
                            raw <- as.integer(decoded)
                            clipped <- raw
                        }
                        design.ncol <- if (is.finite(clipped)) {
                            .coupled.kd.design.ncol(clipped, degree)
                        } else {
                            NA_integer_
                        }
                        if (isTRUE(feasible)) {
                            if (!is.finite(clipped) || clipped < 1L) {
                                feasible <- FALSE
                                skip.reason <- "chart_dim_unavailable"
                            } else if (clipped > ambient.dim) {
                                feasible <- FALSE
                                skip.reason <- "chart_dim_exceeds_ambient"
                            } else if (clipped > chart.dim.max) {
                                feasible <- FALSE
                                skip.reason <- "chart_dim_exceeds_max"
                            } else if (!is.finite(d.feasible.max) ||
                                       clipped > d.feasible.max ||
                                       design.ncol + design.margin >
                                           support.size) {
                                feasible <- FALSE
                                skip.reason <- "design_underdetermined"
                            }
                        }
                        rr <- rr + 1L
                        rows[[rr]] <- data.frame(
                            candidate.id = NA_integer_,
                            stage = stage,
                            support.size = as.integer(support.size),
                            chart.dim = if (is.finite(clipped)) {
                                as.character(as.integer(clipped))
                            } else {
                                label
                            },
                            chart.dim.source = source,
                            chart.dim.raw = raw,
                            chart.dim.clipped = clipped,
                            chart.dim.seed.clipped = seed.clipped,
                            chart.dim.max = as.integer(chart.dim.max),
                            kernel = as.character(kernel),
                            degree = as.integer(degree),
                            bandwidth.multiplier =
                                as.numeric(bandwidth.multiplier),
                            design.ncol = design.ncol,
                            design.margin = as.integer(design.margin),
                            feasible = feasible,
                            skip.reason = skip.reason,
                            reuse.key = NA_character_,
                            score = NA_real_,
                            elapsed.sec = NA_real_,
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
    }
    out <- do.call(rbind, rows[seq_len(rr)])
    out$stage.rank <- .coupled.kd.stage.rank(out$stage)
    out$source.rank <- .coupled.kd.source.rank(out$chart.dim.source)
    out$eval.key <- paste(out$support.size, out$degree, out$kernel,
                          format(out$bandwidth.multiplier, digits = 17L),
                          out$chart.dim.clipped, out$feasible,
                          out$skip.reason, sep = "\r")
    out <- out[order(out$support.size, out$degree, out$kernel,
                     out$bandwidth.multiplier, out$chart.dim.clipped,
                     out$stage.rank, out$source.rank,
                     out$chart.dim.raw, na.last = TRUE), , drop = FALSE]
    out <- out[!duplicated(out$eval.key), , drop = FALSE]

    feasible <- out$feasible & is.finite(out$chart.dim.clipped)
    group.cols <- if (identical(reuse.type, "weighted")) {
        c("support.size", "kernel")
    } else {
        "support.size"
    }
    out$reuse.chart.dim.max <- NA_integer_
    if (any(feasible)) {
        split.key <- do.call(
            paste,
            c(out[feasible, group.cols, drop = FALSE], sep = "\r")
        )
        max.by.group <- tapply(out$chart.dim.clipped[feasible], split.key,
                               max, na.rm = TRUE)
        all.key <- do.call(
            paste,
            c(out[, group.cols, drop = FALSE], sep = "\r")
        )
        out$reuse.chart.dim.max[feasible] <- as.integer(max.by.group[all.key[feasible]])
        out$reuse.key[feasible] <- mapply(
            .coupled.kd.reuse.key,
            support.size = out$support.size[feasible],
            kernel = out$kernel[feasible],
            max.chart.dim = out$reuse.chart.dim.max[feasible],
            MoreArgs = list(reuse.type = reuse.type),
            USE.NAMES = FALSE
        )
    }
    out <- out[, c(
        "candidate.id", "stage", "support.size", "chart.dim",
        "chart.dim.source", "chart.dim.raw", "chart.dim.clipped",
        "chart.dim.seed.clipped", "chart.dim.max", "kernel", "degree",
        "bandwidth.multiplier", "design.ncol", "design.margin",
        "feasible", "skip.reason", "reuse.key", "reuse.chart.dim.max",
        "score", "elapsed.sec"
    )]
    rownames(out) <- NULL
    out$candidate.id <- seq_len(nrow(out))
    out
}

.coupled.kd.numeric.chart.dim.vector <- function(chart.dim) {
    unlist(
        lapply(
            chart.dim,
            function(x) {
                decoded <- tryCatch(.local.chart.decode.chart.dim(x),
                                    error = function(e) NULL)
                if (is.numeric(decoded) && length(decoded) == 1L &&
                    is.finite(decoded) && decoded >= 1L) {
                    as.integer(decoded)
                } else {
                    NA_integer_
                }
            }
        ),
        use.names = FALSE
    )
}

.coupled.kd.reuse.plan <- function(candidates,
                                   reuse.type = c("weighted", "chart")) {
    reuse.type <- match.arg(reuse.type)
    if (!is.data.frame(candidates) || !nrow(candidates) ||
        !"chart.dim" %in% names(candidates) ||
        !"support.size" %in% names(candidates)) {
        return(data.frame())
    }
    if (identical(reuse.type, "weighted") &&
        !"kernel" %in% names(candidates)) {
        return(data.frame())
    }
    dim <- .coupled.kd.numeric.chart.dim.vector(candidates$chart.dim)
    ok <- is.finite(dim)
    if ("feasible" %in% names(candidates)) {
        ok <- ok & !is.na(candidates$feasible) &
            as.logical(candidates$feasible)
    }
    if (!any(ok)) {
        return(data.frame())
    }
    group.cols <- if (identical(reuse.type, "weighted")) {
        c("support.size", "kernel")
    } else {
        "support.size"
    }
    tab <- unique(candidates[ok, group.cols, drop = FALSE])
    tab$max.chart.dim <- NA_integer_
    tab$n.candidates <- NA_integer_
    tab$candidate.ids <- NA_character_
    for (ii in seq_len(nrow(tab))) {
        same <- ok & candidates$support.size == tab$support.size[[ii]]
        if (identical(reuse.type, "weighted")) {
            same <- same & candidates$kernel == tab$kernel[[ii]]
        }
        tab$max.chart.dim[[ii]] <- max(dim[same], na.rm = TRUE)
        tab$n.candidates[[ii]] <- sum(same)
        tab$candidate.ids[[ii]] <- paste(candidates$candidate.id[same],
                                         collapse = ",")
    }
    tab$reuse.chart.dim.max <- tab$max.chart.dim
    tab$reuse.key <- mapply(
        .coupled.kd.reuse.key,
        support.size = tab$support.size,
        kernel = if (identical(reuse.type, "weighted")) tab$kernel else "",
        max.chart.dim = tab$max.chart.dim,
        MoreArgs = list(reuse.type = reuse.type),
        USE.NAMES = FALSE
    )
    rownames(tab) <- NULL
    tab
}

.coupled.kd.lookup.reuse.plan <- function(reuse.plan,
                                          support.size,
                                          chart.dim,
                                          kernel = "gaussian",
                                          reuse.type =
                                              c("weighted", "chart")) {
    reuse.type <- match.arg(reuse.type)
    if (is.null(reuse.plan) || !is.data.frame(reuse.plan) ||
        !nrow(reuse.plan)) {
        return(NULL)
    }
    chart.dim <- as.integer(chart.dim)
    if (length(chart.dim) != 1L || !is.finite(chart.dim) || chart.dim < 1L) {
        return(NULL)
    }
    row <- reuse.plan[
        as.integer(reuse.plan$support.size) == as.integer(support.size),
        ,
        drop = FALSE
    ]
    if (identical(reuse.type, "weighted")) {
        if (!"kernel" %in% names(row)) {
            return(NULL)
        }
        row <- row[as.character(row$kernel) == as.character(kernel), ,
                   drop = FALSE]
    }
    if (nrow(row) != 1L || row$max.chart.dim[[1L]] < chart.dim) {
        return(NULL)
    }
    row
}

.coupled.kd.shared.local.pca.supports <- function(
        X,
        support.size,
        chart.dim,
        kernel = "gaussian",
        coordinate.method = "local.pca",
        reuse.plan,
        cache.env,
        reuse.type = c("weighted", "chart")) {
    reuse.type <- match.arg(reuse.type)
    if (!identical(coordinate.method, "local.pca")) {
        return(NULL)
    }
    if (is.null(cache.env)) {
        cache.env <- new.env(parent = emptyenv())
    }
    row <- .coupled.kd.lookup.reuse.plan(
        reuse.plan = reuse.plan,
        support.size = support.size,
        chart.dim = chart.dim,
        kernel = kernel,
        reuse.type = reuse.type
    )
    if (is.null(row)) {
        return(NULL)
    }
    key <- row$reuse.key[[1L]]
    if (exists(key, envir = cache.env, inherits = FALSE)) {
        return(get(key, envir = cache.env, inherits = FALSE))
    }
    native.kernel <- if (identical(reuse.type, "weighted")) {
        as.character(kernel)
    } else {
        "gaussian"
    }
    supports <- tryCatch(
        rcpp_ps_lps_local_pca_supports(
            X = as.matrix(X),
            support_size = as.integer(support.size),
            chart_dim_by_anchor = rep(as.integer(row$max.chart.dim[[1L]]),
                                      nrow(X)),
            kernel = native.kernel
        ),
        error = function(e) NULL
    )
    if (!is.null(supports)) {
        assign(key, supports, envir = cache.env)
    }
    supports
}

.coupled.kd.local.pca.support.cache <- function(
        X,
        candidates,
        reuse.type = c("weighted", "chart"),
        cache.env = new.env(parent = emptyenv())) {
    reuse.type <- match.arg(reuse.type)
    reuse.plan <- .coupled.kd.reuse.plan(candidates, reuse.type = reuse.type)
    list(
        reuse.type = reuse.type,
        reuse.plan = reuse.plan,
        cache.env = cache.env,
        get = function(support.size,
                       chart.dim,
                       kernel = "gaussian",
                       coordinate.method = "local.pca") {
            .coupled.kd.shared.local.pca.supports(
                X = X,
                support.size = support.size,
                chart.dim = chart.dim,
                kernel = kernel,
                coordinate.method = coordinate.method,
                reuse.plan = reuse.plan,
                cache.env = cache.env,
                reuse.type = reuse.type
            )
        }
    )
}

.coupled.kd.selection.strategy <- function(selection.strategy = "grid") {
    selection.strategy <- selection.strategy %||% "grid"
    match.arg(selection.strategy, c("grid", "sparse_kd", "plateau_kd"))
}

.coupled.kd.numeric.chart.dim.range <- function(chart.dim.grid,
                                                chart.dim.max = NULL,
                                                ambient.dim) {
    cleaned <- .local.chart.clean.chart.dim.grid(chart.dim.grid)
    decoded <- lapply(cleaned$chart.dim, function(x) {
        tryCatch(.local.chart.decode.chart.dim(x), error = function(e) NULL)
    })
    numeric.dim <- sort(unique(unlist(lapply(decoded, function(x) {
        if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
            as.integer(x)
        } else {
            NA_integer_
        }
    }), use.names = FALSE)))
    numeric.dim <- numeric.dim[is.finite(numeric.dim) & numeric.dim >= 1L]
    if (!length(numeric.dim)) {
        stop("'plateau_kd' requires at least one numeric chart dimension in ",
             "'chart.dim.grid'.", call. = FALSE)
    }
    d.max <- min(max(numeric.dim), as.integer(ambient.dim))
    if (!is.null(chart.dim.max)) {
        d.max <- min(d.max, as.integer(chart.dim.max))
    }
    as.integer(max(1L, d.max))
}

.coupled.kd.plateau.dimension <- function(X,
                                          anchor.index,
                                          ordered.index,
                                          support.size,
                                          degree,
                                          chart.dim.max,
                                          variance.threshold,
                                          design.margin) {
    support.size <- as.integer(min(support.size, length(ordered.index)))
    if (!is.finite(support.size) || support.size < 2L) {
        return(1L)
    }
    idx <- ordered.index[seq_len(support.size)]
    centered <- sweep(X[idx, , drop = FALSE], 2L,
                      X[anchor.index, , drop = TRUE], "-")
    sv <- tryCatch(svd(centered, nu = 0L, nv = 0L)$d,
                   error = function(e) numeric(0))
    sv <- sv[is.finite(sv) & sv > 0]
    d.cap <- .coupled.kd.max.feasible.chart.dim(
        support.size = support.size,
        degree = degree,
        ambient.dim = min(ncol(X), chart.dim.max),
        design.margin = design.margin
    )
    if (!is.finite(d.cap) || d.cap < 1L || !length(sv)) {
        return(1L)
    }
    max.dim <- min(as.integer(d.cap), as.integer(chart.dim.max), length(sv),
                   ncol(X), max(1L, support.size - 1L))
    if (!is.finite(max.dim) || max.dim < 1L) {
        return(1L)
    }
    energy <- sv^2
    total <- sum(energy)
    if (!is.finite(total) || total <= 0) {
        return(1L)
    }
    d <- which(cumsum(energy) / total >= variance.threshold)[[1L]]
    as.integer(max(1L, min(max.dim, d)))
}

.coupled.kd.plateau.plan <- function(X,
                                     support.grid,
                                     degree.grid,
                                     kernel.grid,
                                     bandwidth.multiplier.grid,
                                     chart.dim.grid,
                                     chart.dim.max = NULL,
                                     design.margin = 2L,
                                     variance.threshold = 0.95,
                                     max.anchors = 60L,
                                     reuse.type = c("weighted", "chart")) {
    reuse.type <- match.arg(reuse.type)
    X <- as.matrix(X)
    support.grid <- .klp.clean.support.grid(support.grid, nrow(X))
    degree.grid <- .klp.clean.degree.grid(degree.grid)
    kernel.grid <- .klp.clean.kernel.grid(kernel.grid)
    bandwidth.multiplier.grid <- .klp.clean.bandwidth.multiplier.grid(
        bandwidth.multiplier.grid
    )
    d.grid.max <- .coupled.kd.numeric.chart.dim.range(
        chart.dim.grid = chart.dim.grid,
        chart.dim.max = chart.dim.max,
        ambient.dim = ncol(X)
    )
    max.anchors <- as.integer(max.anchors)
    if (!is.finite(max.anchors) || max.anchors < 1L) max.anchors <- 60L
    anchor.index <- if (nrow(X) <= max.anchors) {
        seq_len(nrow(X))
    } else {
        unique(as.integer(round(seq(1, nrow(X), length.out = max.anchors))))
    }
    D <- as.matrix(stats::dist(X))
    selected.rows <- list()
    diagnostic.rows <- list()
    rr <- 0L
    dd <- 0L
    for (degree in degree.grid) {
        anchor.k <- integer(length(anchor.index))
        anchor.d <- integer(length(anchor.index))
        for (aa in seq_along(anchor.index)) {
            ii <- anchor.index[[aa]]
            ord <- order(D[ii, ], na.last = NA)
            dims <- vapply(
                support.grid,
                function(k) .coupled.kd.plateau.dimension(
                    X = X,
                    anchor.index = ii,
                    ordered.index = ord,
                    support.size = k,
                    degree = degree,
                    chart.dim.max = d.grid.max,
                    variance.threshold = variance.threshold,
                    design.margin = design.margin
                ),
                integer(1L)
            )
            d0 <- dims[[1L]]
            stable <- dims == d0
            break.at <- which(!stable)
            last.idx <- if (length(break.at)) break.at[[1L]] - 1L else length(stable)
            last.idx <- max(1L, last.idx)
            anchor.k[[aa]] <- as.integer(support.grid[[last.idx]])
            anchor.d[[aa]] <- as.integer(d0)
            dd <- dd + 1L
            diagnostic.rows[[dd]] <- data.frame(
                degree = as.integer(degree),
                anchor.index = as.integer(ii),
                plateau.support.size = as.integer(anchor.k[[aa]]),
                plateau.chart.dim = as.integer(anchor.d[[aa]]),
                min.support.size = as.integer(support.grid[[1L]]),
                max.support.size = as.integer(support.grid[[length(support.grid)]]),
                n.support.grid = length(support.grid),
                dim.path = paste(as.integer(dims), collapse = ","),
                stringsAsFactors = FALSE
            )
        }
        k.med <- stats::median(anchor.k, na.rm = TRUE)
        k.selected <- support.grid[which.min(abs(support.grid - k.med))]
        d.selected <- as.integer(round(stats::median(anchor.d, na.rm = TRUE)))
        d.cap <- .coupled.kd.max.feasible.chart.dim(
            support.size = k.selected,
            degree = degree,
            ambient.dim = min(ncol(X), d.grid.max),
            design.margin = design.margin
        )
        d.selected <- as.integer(max(1L, min(d.selected, d.grid.max, d.cap)))
        for (kernel in kernel.grid) {
            for (bandwidth.multiplier in bandwidth.multiplier.grid) {
                rr <- rr + 1L
                selected.rows[[rr]] <- data.frame(
                    candidate.id = rr,
                    stage = "skeleton",
                    support.size = as.integer(k.selected),
                    chart.dim = as.character(d.selected),
                    chart.dim.source = "plateau_geometry",
                    chart.dim.raw = as.integer(d.selected),
                    chart.dim.clipped = as.integer(d.selected),
                    chart.dim.seed.clipped = FALSE,
                    chart.dim.max = as.integer(d.grid.max),
                    kernel = as.character(kernel),
                    degree = as.integer(degree),
                    bandwidth.multiplier = as.numeric(bandwidth.multiplier),
                    design.ncol = .coupled.kd.design.ncol(d.selected, degree),
                    design.margin = as.integer(design.margin),
                    feasible = TRUE,
                    skip.reason = NA_character_,
                    reuse.key = NA_character_,
                    reuse.chart.dim.max = as.integer(d.selected),
                    score = NA_real_,
                    elapsed.sec = NA_real_,
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    cand <- do.call(rbind, selected.rows)
    cand$stage.rank <- .coupled.kd.stage.rank(cand$stage)
    cand$source.rank <- .coupled.kd.source.rank(cand$chart.dim.source)
    cand <- cand[order(cand$support.size, cand$degree, cand$kernel,
                       cand$bandwidth.multiplier), , drop = FALSE]
    cand$candidate.id <- seq_len(nrow(cand))
    cand$reuse.key <- mapply(
        .coupled.kd.reuse.key,
        support.size = cand$support.size,
        kernel = cand$kernel,
        max.chart.dim = cand$reuse.chart.dim.max,
        MoreArgs = list(reuse.type = reuse.type),
        USE.NAMES = FALSE
    )
    cand$chart.dim.rank <- match(
        cand$chart.dim,
        unique(cand$chart.dim[order(cand$chart.dim.clipped)])
    )
    diagnostics <- do.call(rbind, diagnostic.rows)
    rownames(diagnostics) <- NULL
    plan <- cand[, c(
        "candidate.id", "stage", "support.size", "chart.dim",
        "chart.dim.source", "chart.dim.raw", "chart.dim.clipped",
        "chart.dim.seed.clipped", "chart.dim.max", "kernel", "degree",
        "bandwidth.multiplier", "design.ncol", "design.margin",
        "feasible", "skip.reason", "reuse.key", "reuse.chart.dim.max",
        "score", "elapsed.sec"
    ), drop = FALSE]
    reuse.plan <- .coupled.kd.reuse.plan(cand, reuse.type = reuse.type)
    list(
        candidates = cand[, c("candidate.id", "support.size", "degree",
                              "kernel", "bandwidth.multiplier", "chart.dim",
                              "chart.dim.source", "chart.dim.raw",
                              "chart.dim.clipped", "chart.dim.seed.clipped",
                              "chart.dim.max", "design.ncol",
                              "design.margin", "reuse.key",
                              "reuse.chart.dim.max", "chart.dim.rank"),
                          drop = FALSE],
        coupled.plan = plan,
        diagnostics = diagnostics,
        telemetry = list(
            selection.strategy = "plateau_kd",
            coupled.chart.dim.search = TRUE,
            geometry.only = TRUE,
            planned.candidates = length(support.grid) * length(degree.grid) *
                length(kernel.grid) * length(bandwidth.multiplier.grid),
            evaluated.candidates = nrow(cand),
            skipped.candidates = 0L,
            reuse.groups = nrow(reuse.plan),
            reuse.type = reuse.type,
            support.grid.planned = support.grid,
            support.grid.evaluated = sort(unique(cand$support.size)),
            chart.dim.grid.planned = chart.dim.grid,
            chart.dim.evaluated = sort(unique(cand$chart.dim.clipped)),
            chart.dim.max = chart.dim.max,
            plateau.variance.threshold = variance.threshold,
            plateau.max.anchors = max.anchors,
            plateau.anchor.diagnostics = diagnostics,
            design.margin = design.margin
        )
    )
}

.coupled.kd.sparse.support.grid <- function(support.grid) {
    support.grid <- sort(unique(as.integer(support.grid)))
    if (length(support.grid) <= 3L) {
        return(support.grid)
    }
    mid <- support.grid[[ceiling(length(support.grid) / 2)]]
    sort(unique(c(support.grid[[1L]], mid, support.grid[[length(support.grid)]])))
}

.coupled.kd.sparse.chart.dim.grid <- function(chart.dim.grid,
                                              chart.dim.max = NULL) {
    cleaned <- .local.chart.clean.chart.dim.grid(chart.dim.grid)
    decoded <- lapply(cleaned$chart.dim, function(x) {
        tryCatch(.local.chart.decode.chart.dim(x), error = function(e) NULL)
    })
    numeric.dim <- sort(unique(unlist(lapply(decoded, function(x) {
        if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
            as.integer(x)
        } else {
            NA_integer_
        }
    }), use.names = FALSE)))
    numeric.dim <- numeric.dim[is.finite(numeric.dim)]
    special <- unique(cleaned$chart.dim[!vapply(decoded, is.numeric,
                                               logical(1L))])
    out <- integer()
    if (length(numeric.dim)) {
        d.hi <- max(numeric.dim)
        if (!is.null(chart.dim.max)) {
            d.hi <- min(d.hi, as.integer(chart.dim.max))
        }
        out <- sort(unique(c(numeric.dim[numeric.dim %in% c(1L, 2L)], d.hi)))
    }
    unique(c(as.character(out), special))
}

.coupled.kd.has.special.chart.dim <- function(chart.dim.grid,
                                              chart.dim = NULL) {
    cleaned <- .local.chart.clean.chart.dim.grid(
        chart.dim.grid,
        chart.dim = chart.dim
    )
    decoded <- lapply(cleaned$chart.dim, function(x) {
        tryCatch(.local.chart.decode.chart.dim(x), error = function(e) NA)
    })
    any(vapply(
        decoded,
        function(x) !(is.numeric(x) || is.null(x)),
        logical(1L)
    ))
}

.coupled.kd.auto.chart.dim.function <- function(X,
                                                coordinate.method,
                                                auto.chart.support.metric,
                                                auto.chart.selection.metric) {
    force(X)
    force(coordinate.method)
    force(auto.chart.support.metric)
    force(auto.chart.selection.metric)
    function(support.size, degree) {
        info <- .klp.resolve.chart.dim(
            X = X,
            support.size = support.size,
            degree = degree,
            coordinate.method = coordinate.method,
            chart.dim = "auto",
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        as.integer(info$chart.dim)
    }
}

.coupled.kd.lps.candidate.spec <- function(X,
                                           support.grid,
                                           degree.grid,
                                           kernel.grid,
                                           bandwidth.multiplier.grid,
                                           chart.dim = NULL,
                                           chart.dim.grid = NULL,
                                           coordinate.method,
                                           auto.chart.support.metric,
                                           auto.chart.selection.metric,
                                           selection.strategy = "grid",
                                           chart.dim.max = NULL,
                                           design.margin = 2L,
                                           reuse.type = c("weighted",
                                                          "chart")) {
    reuse.type <- match.arg(reuse.type)
    selection.strategy <- .coupled.kd.selection.strategy(selection.strategy)
    support.grid <- .klp.clean.support.grid(support.grid, nrow(X))
    degree.grid <- .klp.clean.degree.grid(degree.grid)
    kernel.grid <- .klp.clean.kernel.grid(kernel.grid)
    bandwidth.multiplier.grid <- .klp.clean.bandwidth.multiplier.grid(
        bandwidth.multiplier.grid
    )
    if (is.null(chart.dim.grid) && identical(selection.strategy, "plateau_kd")) {
        stop("'plateau_kd' requires 'chart.dim.grid' so the geometry-only ",
             "selector is bounded by an explicit numeric dimension universe.",
             call. = FALSE)
    }
    if (!identical(coordinate.method, "local.pca") &&
        identical(selection.strategy, "plateau_kd")) {
        stop("'plateau_kd' requires coordinate.method = 'local.pca'.",
             call. = FALSE)
    }
    if (identical(selection.strategy, "plateau_kd")) {
        return(.coupled.kd.plateau.plan(
            X = X,
            support.grid = support.grid,
            degree.grid = degree.grid,
            kernel.grid = kernel.grid,
            bandwidth.multiplier.grid = bandwidth.multiplier.grid,
            chart.dim.grid = chart.dim.grid,
            chart.dim.max = chart.dim.max,
            design.margin = design.margin,
            reuse.type = reuse.type
        ))
    }
    chart.dim.grid <- if (is.null(chart.dim.grid)) {
        NULL
    } else if (identical(selection.strategy, "sparse_kd")) {
        .coupled.kd.sparse.chart.dim.grid(
            chart.dim.grid,
            chart.dim.max = chart.dim.max
        )
    } else {
        chart.dim.grid
    }
    support.eval.grid <- if (identical(selection.strategy, "sparse_kd")) {
        .coupled.kd.sparse.support.grid(support.grid)
    } else {
        support.grid
    }
    if (is.null(chart.dim.grid)) {
        cand <- expand.grid(
            support.size = support.eval.grid,
            degree = degree.grid,
            kernel = kernel.grid,
            bandwidth.multiplier = bandwidth.multiplier.grid,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        cand$candidate.id <- seq_len(nrow(cand))
        return(list(
            candidates = cand,
            coupled.plan = NULL,
            telemetry = list(
                selection.strategy = "grid",
                coupled.chart.dim.search = FALSE,
                planned.candidates = nrow(cand),
                evaluated.candidates = nrow(cand),
                skipped.candidates = 0L,
                reuse.groups = 0L
            )
        ))
    }
    if (!identical(coordinate.method, "local.pca")) {
        stop("'chart.dim.grid' requires coordinate.method = 'local.pca'.",
             call. = FALSE)
    }
    auto.fun <- .coupled.kd.auto.chart.dim.function(
        X = X,
        coordinate.method = coordinate.method,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric
    )
    plan <- .coupled.kd.candidate.table(
        support.grid = support.eval.grid,
        degree.grid = degree.grid,
        chart.dim.grid = chart.dim.grid,
        ambient.dim = ncol(X),
        kernel.grid = kernel.grid,
        bandwidth.multiplier.grid = bandwidth.multiplier.grid,
        stage = if (identical(selection.strategy, "sparse_kd")) {
            "skeleton"
        } else {
            "full_reference"
        },
        chart.dim.max = chart.dim.max,
        design.margin = design.margin,
        auto.chart.dim = auto.fun,
        reuse.type = reuse.type
    )
    feasible <- plan$feasible & is.finite(plan$chart.dim.clipped)
    cand <- plan[feasible, , drop = FALSE]
    if (!nrow(cand)) {
        stop("Coupled support/chart-dimension selection produced no ",
             "feasible candidates.", call. = FALSE)
    }
    cand <- cand[, c("candidate.id", "support.size", "degree", "kernel",
                     "bandwidth.multiplier", "chart.dim",
                     "chart.dim.source", "chart.dim.raw",
                     "chart.dim.clipped", "chart.dim.seed.clipped",
                     "chart.dim.max", "design.ncol", "design.margin",
                     "reuse.key", "reuse.chart.dim.max"),
                 drop = FALSE]
    cand$chart.dim <- as.character(cand$chart.dim.clipped)
    cand$chart.dim.rank <- match(
        cand$chart.dim,
        unique(cand$chart.dim[order(cand$chart.dim.clipped)])
    )
    cand$candidate.id <- seq_len(nrow(cand))
    reuse.plan <- .coupled.kd.reuse.plan(cand, reuse.type = reuse.type)
    list(
        candidates = cand,
        coupled.plan = plan,
        telemetry = list(
            selection.strategy = selection.strategy,
            coupled.chart.dim.search = TRUE,
            planned.candidates = nrow(plan),
            evaluated.candidates = nrow(cand),
            skipped.candidates = sum(!feasible),
            reuse.groups = nrow(reuse.plan),
            reuse.type = reuse.type,
            support.grid.planned = support.grid,
            support.grid.evaluated = sort(unique(cand$support.size)),
            chart.dim.grid.planned = chart.dim.grid,
            chart.dim.evaluated = sort(unique(cand$chart.dim.clipped)),
            chart.dim.max = chart.dim.max,
            design.margin = design.margin
        )
    )
}

.coupled.kd.chart.candidate.spec <- function(X,
                                             support.grid,
                                             kernel.grid,
                                             bandwidth.multiplier.grid,
                                             chart.dim = NULL,
                                             chart.dim.grid = NULL,
                                             coordinate.method,
                                             auto.chart.support.metric,
                                             auto.chart.selection.metric,
                                             selection.strategy = "grid",
                                             chart.dim.max = NULL,
                                             geometry.margin = 0L) {
    geometry.margin <- as.integer(geometry.margin)
    if (length(geometry.margin) != 1L || !is.finite(geometry.margin) ||
        geometry.margin < 0L) {
        stop("'geometry.margin' must be a nonnegative integer scalar.",
             call. = FALSE)
    }
    out <- .coupled.kd.lps.candidate.spec(
        X = X,
        support.grid = support.grid,
        degree.grid = 0L,
        kernel.grid = kernel.grid,
        bandwidth.multiplier.grid = bandwidth.multiplier.grid,
        chart.dim = chart.dim,
        chart.dim.grid = chart.dim.grid,
        coordinate.method = coordinate.method,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
        selection.strategy = selection.strategy,
        chart.dim.max = chart.dim.max,
        design.margin = 0L,
        reuse.type = "chart"
    )
    apply.contract <- function(tab) {
        if (is.null(tab)) return(NULL)
        tab$geometry.rank.cap <- pmax(
            0L,
            pmin(ncol(X), as.integer(tab$support.size) - 1L -
                     geometry.margin)
        )
        tab$geometry.margin <- geometry.margin
        tab$feasibility.contract <- "chart_kernel"
        if ("chart.dim.clipped" %in% names(tab)) {
            finite <- is.finite(tab$chart.dim.clipped)
            tab$chart.dim.clipped[finite] <- pmin(
                tab$chart.dim.clipped[finite],
                tab$geometry.rank.cap[finite]
            )
            tab$chart.dim[finite] <- as.character(
                as.integer(tab$chart.dim.clipped[finite])
            )
            bad <- !is.finite(tab$geometry.rank.cap) |
                tab$geometry.rank.cap < 1L |
                !is.finite(tab$chart.dim.clipped) |
                tab$chart.dim.clipped < 1L
            if ("feasible" %in% names(tab)) {
                tab$feasible[bad] <- FALSE
            }
            if ("skip.reason" %in% names(tab)) {
                tab$skip.reason[bad] <- "chart_geometry_rank_infeasible"
            }
            if ("reuse.chart.dim.max" %in% names(tab)) {
                tab$reuse.chart.dim.max[finite] <- pmin(
                    tab$reuse.chart.dim.max[finite],
                    tab$geometry.rank.cap[finite]
                )
            }
            if ("reuse.key" %in% names(tab)) {
                usable <- finite & tab$chart.dim.clipped >= 1L
                tab$reuse.key[usable] <- mapply(
                    .coupled.kd.reuse.key,
                    support.size = tab$support.size[usable],
                    kernel = tab$kernel[usable],
                    max.chart.dim = tab$reuse.chart.dim.max[usable],
                    MoreArgs = list(reuse.type = "chart"),
                    USE.NAMES = FALSE
                )
            }
        }
        tab$degree <- NULL
        tab$design.ncol <- NULL
        tab$design.margin <- NULL
        tab
    }
    out$candidates <- apply.contract(out$candidates)
    out$coupled.plan <- apply.contract(out$coupled.plan)
    out$telemetry$feasibility.contract <- "chart_kernel"
    out$telemetry$geometry.margin <- geometry.margin
    out$telemetry$design.margin <- NULL
    out
}
