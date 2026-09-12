# Shared helpers used by geosmooth split-package code.

`%||%` <- function(a, b) if (is.null(a)) b else a

.geosmooth.ge1.missing.native <- function(feature, phase = "GE2") {
    stop(
        feature, " is not available in geosmooth GE1. ",
        "It is scheduled for ", phase, ".",
        call. = FALSE
    )
}

.geosmooth.graph.geodesic.stage <- function(stage) {
    match.arg(
        stage,
        c("final", "raw", "raw.repaired", "pruned", "pruned.repaired",
          "repaired.pruned")
    )
}

.graph.geodesic.fields <- function(graph, stage = "final") {
    if (!is.list(graph)) {
        stop("'graph' must be a list-like graph object.", call. = FALSE)
    }
    stage <- .geosmooth.graph.geodesic.stage(stage)
    candidates <- switch(
        stage,
        final = list(c("adj_list", "weight_list"),
                     c("adj.list", "weight.list")),
        raw = list(c("raw_adj_list", "raw_weight_list"),
                   c("raw.adj.list", "raw.weight.list")),
        raw.repaired = list(c("raw_repaired_adj_list", "raw_repaired_weight_list"),
                            c("raw.repaired.adj.list", "raw.repaired.weight.list")),
        pruned = list(c("pruned_adj_list", "pruned_weight_list"),
                      c("pruned.adj.list", "pruned.weight.list")),
        pruned.repaired = list(c("pruned_repaired_adj_list",
                                 "pruned_repaired_weight_list"),
                               c("pruned.repaired.adj.list",
                                 "pruned.repaired.weight.list")),
        repaired.pruned = list(c("repaired_pruned_adj_list",
                                 "repaired_pruned_weight_list"),
                               c("repaired.pruned.adj.list",
                                 "repaired.pruned.weight.list"))
    )
    for (pair in candidates) {
        if (!is.null(graph[[pair[[1L]]]]) &&
            !is.null(graph[[pair[[2L]]]])) {
            return(list(adj = pair[[1L]], weight = pair[[2L]], stage = stage))
        }
    }
    stop(
        "Could not extract graph adjacency and weights for graph.stage = '",
        stage, "'.",
        call. = FALSE
    )
}

.validate.graph.geodesic.payload <- function(adj.list, weight.list, fields) {
    if (!is.list(fields) ||
        !all(c("adj", "weight", "stage") %in% names(fields))) {
        stop("Internal graph field metadata is malformed.", call. = FALSE)
    }
    tryCatch(
        .validate.metric.graph.lowpass.graph(adj.list, weight.list),
        error = function(e) {
            stop(
                "Invalid graph-geodesic payload for graph.stage = '",
                fields$stage, "' using fields '", fields$adj, "' and '",
                fields$weight, "': ", conditionMessage(e),
                call. = FALSE
            )
        }
    )
    invisible(TRUE)
}

.geosmooth.shortest.path <- function(graph, edge.lengths, vertices) {
    dgraphs::shortest.path(graph, edge.lengths, vertices)
}

.exact.knn.index <- function(X, k) {
    n <- nrow(X)
    out <- matrix(NA_integer_, nrow = n, ncol = k)
    for (i in seq_len(n)) {
        d <- rowSums((t(t(X) - X[i, ]))^2)
        d[[i]] <- Inf
        out[i, ] <- order(d, seq_len(n))[seq_len(k)]
    }
    out
}
