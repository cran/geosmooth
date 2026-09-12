.local.chart.cv.requested <- function(foldid,
                                      support.grid = NULL,
                                      degree.grid = NULL,
                                      kernel.grid = NULL,
                                      bandwidth.multiplier.grid = NULL,
                                      lambda.ridge.grid = NULL,
                                      chart.dim.grid = NULL) {
    !is.null(foldid) ||
        !is.null(support.grid) ||
        !is.null(degree.grid) ||
        !is.null(kernel.grid) ||
        !is.null(bandwidth.multiplier.grid) ||
        !is.null(lambda.ridge.grid) ||
        !is.null(chart.dim.grid)
}

.local.likelihood.clean.lambda.ridge.grid <- function(lambda.ridge.grid) {
    out <- sort(unique(as.numeric(lambda.ridge.grid)))
    out <- out[is.finite(out) & out >= 0]
    if (!length(out)) {
        stop("'lambda.ridge.grid' must contain at least one finite ",
             "nonnegative value.", call. = FALSE)
    }
    out
}

.local.chart.dim.label <- function(chart.dim) {
    if (is.null(chart.dim)) {
        return("NULL")
    }
    if (is.character(chart.dim) && length(chart.dim) == 1L) {
        if (chart.dim %in% c("auto", "local.auto")) {
            return(chart.dim)
        }
        suppressWarnings(value <- as.integer(chart.dim))
        if (!is.na(value) && value >= 1L &&
            identical(as.character(value), chart.dim)) {
            return(as.character(value))
        }
    }
    if (is.numeric(chart.dim) && length(chart.dim) == 1L &&
        is.finite(chart.dim) && chart.dim >= 1L &&
        chart.dim == as.integer(chart.dim)) {
        return(as.character(as.integer(chart.dim)))
    }
    stop("'chart.dim' candidates must be NULL, positive integers, ",
         "'auto', or 'local.auto'.", call. = FALSE)
}

.local.chart.clean.chart.dim.grid <- function(chart.dim.grid,
                                              chart.dim = NULL) {
    if (is.null(chart.dim.grid)) {
        chart.dim.grid <- list(chart.dim)
    }
    if (!is.list(chart.dim.grid)) {
        chart.dim.grid <- as.list(chart.dim.grid)
    }
    labels <- vapply(
        chart.dim.grid,
        .local.chart.dim.label,
        character(1L)
    )
    labels <- unique(labels)
    if (!length(labels)) {
        stop("'chart.dim.grid' must contain at least one chart-dimension ",
             "candidate.", call. = FALSE)
    }
    data.frame(
        chart.dim = labels,
        chart.dim.rank = seq_along(labels),
        stringsAsFactors = FALSE
    )
}

.local.chart.decode.chart.dim <- function(label) {
    label <- as.character(label[[1L]])
    if (identical(label, "NULL")) {
        return(NULL)
    }
    if (label %in% c("auto", "local.auto")) {
        return(label)
    }
    suppressWarnings(value <- as.integer(label))
    if (!is.na(value) && value >= 1L &&
        identical(as.character(value), label)) {
        return(value)
    }
    stop("Invalid chart-dimension candidate label '", label, "'.",
         call. = FALSE)
}

.local.chart.select.best.idx <- function(cv.table,
                                         score.column = "cv.rmse.observed",
                                         tolerance = 1e-12) {
    if (!score.column %in% names(cv.table)) {
        stop("CV table does not contain selection score column '",
             score.column, "'.", call. = FALSE)
    }
    score <- cv.table[[score.column]]
    finite <- is.finite(score)
    if (!any(finite)) {
        stop("No candidate has a finite selection score in '",
             score.column, "'.", call. = FALSE)
    }
    best <- min(score[finite])
    eligible <- which(
        finite &
            score <= best + max(tolerance, tolerance * abs(best))
    )
    column.or.default <- function(name, default) {
        if (name %in% names(cv.table)) cv.table[[name]][eligible] else default
    }
    eligible[order(
        column.or.default("support.size", rep(0, length(eligible))),
        column.or.default("degree", rep(0, length(eligible))),
        column.or.default("chart.dim.rank", rep(0, length(eligible))),
        column.or.default("chart.dim", rep("", length(eligible))),
        column.or.default("kernel", rep("", length(eligible))),
        column.or.default("bandwidth.multiplier", rep(1, length(eligible))),
        column.or.default("lambda.sync", rep(0, length(eligible))),
        column.or.default("lambda.ridge", rep(0, length(eligible))),
        column.or.default("walk.step", rep(0, length(eligible))),
        column.or.default("affinity.method", rep("", length(eligible))),
        column.or.default("affinity.scale", rep(Inf, length(eligible))),
        column.or.default("affinity.epsilon", rep(1e-12, length(eligible))),
        score[eligible]
    )[[1L]]]
}
