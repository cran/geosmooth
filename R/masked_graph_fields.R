#' Fit a Graph Field from Partially Observed Responses
#'
#' Estimate a conditional mean by weighted least squares with an unnormalized
#' graph Laplacian penalty. Repeated observations at a vertex remain separate.
#' Unobserved components return NA rather than an arbitrary constant.
#' @param operator An unnormalized `metric.graph.lowpass.operator`.
#' @param vertex Integer vertex index for each observed response.
#' @param y Finite numeric observed responses (no missing-value placeholders).
#' @param lambda Positive smoothing parameter.
#' @param weights Positive observation weights; normalized to sum to one.
#' @return A list with vertex fitted values, support flags, observation weights,
#'   objective, and relative linear-system residual.
#' @examples
#' op <- metric.graph.lowpass.operator(list(2L, c(1L, 3L), 2L),
#'                                    list(1, c(1, 1), 1))
#' fit <- fit.masked.graph.field(op, vertex = c(1, 3), y = c(2, 8), lambda = 0.5)
#' fit$fitted
#' @export
fit.masked.graph.field <- function(operator, vertex, y, lambda, weights = NULL) {
    L <- .field.laplacian(operator)
    n <- nrow(L)
    .field.observations(vertex, y, n)
    if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
        stop("lambda must be positive and finite.")
    if (is.null(weights)) weights <- rep(1, length(y))
    if (length(weights) != length(y) || any(!is.finite(weights) | weights <= 0))
        stop("weights must be positive finite observation weights.")
    weights <- weights / sum(weights)
    A <- Matrix::sparseMatrix(i = seq_along(y), j = vertex, x = 1,
                              dims = c(length(y), n))
    q <- as.numeric(Matrix::crossprod(A, weights))
    b <- as.numeric(Matrix::crossprod(A, weights * y))
    component <- .field.components(operator$graph$adj.list)
    supported <- component %in% unique(component[vertex])
    ids <- which(supported)
    H <- Matrix::Diagonal(n, q) + lambda * L
    f <- rep(NA_real_, n)
    f[ids] <- as.numeric(Matrix::solve(Matrix::forceSymmetric(H[ids, ids, drop = FALSE]), b[ids]))
    residual <- as.numeric(H[ids, ids, drop = FALSE] %*% f[ids]) - b[ids]
    rel <- sqrt(sum(residual^2)) / max(sqrt(sum(b[ids]^2)), .Machine$double.eps)
    if (any(!is.finite(f[ids])) || rel > 1e-7)
        stop("Masked field solve failed its numerical residual check.")
    penalty <- as.numeric(Matrix::crossprod(f[ids], L[ids, ids, drop = FALSE] %*% f[ids]))
    list(fitted = f, supported = supported, component = component,
         lambda = lambda, weights = weights,
         objective = sum(weights * (y - f[vertex])^2) + lambda * penalty,
         relative.residual = rel)
}

.field.laplacian <- function(operator) {
    if (!inherits(operator, "metric.graph.lowpass.operator"))
        stop("operator must be a metric.graph.lowpass.operator.")
    if (!identical(operator$laplacian.type, "unnormalized"))
        stop("Use an unnormalized Laplacian.")
    L <- operator$laplacian$matrix
    if (is.null(L) || !inherits(L, "sparseMatrix") || nrow(L) < 1L)
        stop("operator must contain a sparse Laplacian.")
    if (max(abs(Matrix::rowSums(L))) > 1e-8 * max(1, max(abs(L))))
        stop("Use an unnormalized, constant-preserving Laplacian.")
    L
}
.field.observations <- function(vertex, y, n) {
    if (!is.numeric(y) || !length(y) || any(!is.finite(y)) ||
        length(vertex) != length(y) || any(!is.finite(vertex)) ||
        any(vertex != as.integer(vertex) | vertex < 1 | vertex > n))
        stop("Supply finite observed y and matching valid integer vertex indices.")
}
.field.components <- function(adj) {
    n <- length(adj); component <- integer(n); label <- 0L
    queue <- integer(n)
    for (root in seq_len(n)) {
        if (component[root]) next
        label <- label + 1L; component[root] <- label
        head <- 1L; tail <- 1L; queue[1L] <- root
        while (head <= tail) {
            neighbours <- adj[[queue[head]]]; head <- head + 1L
            new <- unique(neighbours[component[neighbours] == 0L])
            if (length(new)) {
                component[new] <- label
                queue[tail + seq_along(new)] <- new; tail <- tail + length(new)
            }
        }
    }
    component
}
.field.weights <- function(group) {
    group <- as.character(group)
    1 / (length(unique(group)) * as.numeric(table(group)[group]))
}
.field.folds <- function(group, k, seed) {
    groups <- sort(unique(as.character(group)))
    if (length(k) != 1L || !is.finite(k) || k != as.integer(k) || k < 2 || k > length(groups))
        stop("Fold count must be between two and the number of groups.")
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv) else
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv))
    set.seed(seed)
    assignment <- sample(rep(seq_len(k), length.out = length(groups)))
    assignment[match(as.character(group), groups)]
}
.field.score <- function(y, pred, group) {
    ok <- is.finite(pred)
    if (!any(ok)) return(NA_real_)
    sum(.field.weights(group[ok]) * (y[ok] - pred[ok])^2)
}
.field.tune <- function(operator, vertex, y, group, lambdas, folds) {
    predictions <- matrix(NA_real_, length(y), length(lambdas))
    for (fold in sort(unique(folds))) {
        train <- folds != fold; test <- !train
        for (j in seq_along(lambdas)) {
            fit <- fit.masked.graph.field(operator, vertex[train], y[train],
                lambdas[j], .field.weights(group[train]))
            predictions[test, j] <- fit$fitted[vertex[test]]
        }
    }
    scores <- vapply(seq_along(lambdas), function(j)
        .field.score(y, predictions[, j], group), numeric(1))
    if (!any(is.finite(scores))) stop("No supported inner-fold predictions.")
    best <- which.min(replace(scores, !is.finite(scores), Inf))
    list(lambda = lambdas[best], scores = scores,
         coverage = colMeans(is.finite(predictions)),
         boundary = best %in% c(1L, length(lambdas)))
}

#' Participant-Grouped Nested Cross-Validation for a Masked Graph Field
#'
#' Select smoothing within outer training folds, predict held-out groups, and
#' finally tune and refit on all observations. Graph construction is label-free
#' and fixed (transductive prediction). Unsupported predictions remain NA.
#' @inheritParams fit.masked.graph.field
#' @param group Participant or leakage-block identifier per observation.
#' @param lambdas Positive finite candidate penalties.
#' @param outer.folds,inner.folds Numbers of grouped folds.
#' @param seed Fold seed, independent of worker scheduling.
#' @param progress Optional function called after each outer fold with fold and
#'   selected penalty; intended for experiment progress records.
#' @param fold.id Optional preassigned outer fold per observation; a group must
#'   occupy exactly one fold.
#' @return Outer predictions, baseline predictions, scores on supported cases,
#'   coverage, tuning records, fold assignments and final full-data field.
#'   Scores do not account for selecting a graph after comparing graphs.
#' @examples
#' op <- metric.graph.lowpass.operator(list(2L, c(1L, 3L), 2L),
#'                                    list(1, c(1, 1), 1))
#' cv <- cv.masked.graph.field(op, rep(c(1, 3), 6), 1:12, letters[1:12],
#'     lambdas = c(0.01, 1), outer.folds = 3, inner.folds = 2)
#' cv$coverage
#' cv$final$fitted
#' @export
cv.masked.graph.field <- function(operator, vertex, y, group,
    lambdas = 10^seq(-6, 2, length.out = 17), outer.folds = 5L,
    inner.folds = 4L, seed = 1L, fold.id = NULL, progress = NULL) {
    .field.observations(vertex, y, nrow(.field.laplacian(operator)))
    if (length(group) != length(y) || anyNA(group) || any(!nzchar(as.character(group))))
        stop("Supply one nonmissing group identifier per observation.")
    if (!is.null(progress) && !is.function(progress)) stop("progress must be a function or NULL.")
    group <- as.character(group)
    if (!is.numeric(lambdas) || !length(lambdas) || any(!is.finite(lambdas) | lambdas <= 0))
        stop("lambdas must be positive finite values.")
    lambdas <- sort(unique(lambdas))
    if (is.null(fold.id)) fold.id <- .field.folds(group, outer.folds, seed)
    if (length(fold.id) != length(y) || anyNA(fold.id) || length(unique(fold.id)) < 2L ||
        any(vapply(split(fold.id, group), function(z) length(unique(z)) != 1L, logical(1))))
        stop("fold.id must keep each group in one fold, with at least two folds.")
    pred <- baseline <- rep(NA_real_, length(y)); tuning <- list()
    for (fold in sort(unique(fold.id))) {
        train <- fold.id != fold; test <- !train
        inner <- .field.folds(group[train], inner.folds, seed + 1L)
        tuned <- .field.tune(operator, vertex[train], y[train], group[train], lambdas, inner)
        fit <- fit.masked.graph.field(operator, vertex[train], y[train],
            tuned$lambda, .field.weights(group[train]))
        pred[test] <- fit$fitted[vertex[test]]
        baseline[test] <- sum(.field.weights(group[train]) * y[train])
        tuning[[as.character(fold)]] <- tuned
        if (!is.null(progress)) progress(list(fold=fold, lambda=tuned$lambda,
            boundary=tuned$boundary, heldout.coverage=mean(is.finite(pred[test]))))
    }
    final.folds <- .field.folds(group, inner.folds, seed + 2L)
    final.tuning <- .field.tune(operator, vertex, y, group, lambdas, final.folds)
    ok <- is.finite(pred)
    mse <- .field.score(y, pred, group)
    baseline.mse <- if (any(ok)) .field.score(y[ok], baseline[ok], group[ok]) else NA_real_
    list(predicted = pred, baseline = baseline, observed = y, vertex = vertex,
         group = group, fold.id = fold.id, mse = mse, rmse = sqrt(mse),
         baseline.mse = baseline.mse,
         relative.mse = if (is.finite(baseline.mse) && baseline.mse > 0) mse / baseline.mse else NA_real_,
         coverage = mean(ok), participant.coverage = mean(tapply(ok, group, all)),
         tuning = tuning, final.tuning = final.tuning,
         final = fit.masked.graph.field(operator, vertex, y, final.tuning$lambda,
                                        .field.weights(group)))
}

#' Compare Masked Graph Fields with Independent Graph Workers
#'
#' Runs one candidate graph per process, with dynamically scheduled jobs and no
#' nested fold parallelism. On Unix, forked workers inherit the loaded package
#' sources. Run from Rscript, not a GUI or multithreaded parent process.
#' @param operators Named list of metric graph low-pass operators.
#' @param studies Named list of argument lists containing vertex, y and group.
#' @param workers Number of graph workers; choose using available memory and CPU.
#'   Defaults to serial execution. Windows currently requires one worker.
#' @param ... Additional arguments passed to cv.masked.graph.field.
#' @return Named graph results, each containing named study results or an explicit
#'   error record. Retains per-observation predictions for common-support scoring.
#' @usage compare.masked.graph.fields(operators, studies, workers = 1L, ...)
#' @examples
#' op <- metric.graph.lowpass.operator(list(2L, c(1L, 3L), 2L),
#'                                    list(1, c(1, 1), 1))
#' studies <- list(example = list(vertex = rep(c(1, 3), 6), y = 1:12,
#'                                group = letters[1:12]))
#' comparison <- compare.masked.graph.fields(list(path = op), studies,
#'     workers = 1, lambdas = c(0.01, 1), outer.folds = 3, inner.folds = 2)
#' comparison$path$example$coverage
#' @export compare.masked.graph.fields
compare.masked.graph.fields <- function(operators, studies, workers = 1L, ...) {
    if (!length(operators) || is.null(names(operators)) || anyDuplicated(names(operators)) ||
        !length(studies) || is.null(names(studies)) || anyDuplicated(names(studies)))
        stop("Supply uniquely named operators and studies.")
    if (length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != as.integer(workers))
        stop("workers must be a positive integer.")
    workers <- min(workers, length(operators))
    if (.Platform$OS.type == "windows" && workers > 1L)
        stop("Parallel graph workers currently require Unix; use workers = 1 on Windows.")
    orders <- vapply(operators, function(op) nrow(.field.laplacian(op)), integer(1))
    if (length(unique(orders)) != 1L) stop("All graphs must have the same vertex count and ordering.")
    args <- list(...)
    job <- function(j) {
        lapply(studies, function(study) tryCatch(
            do.call(cv.masked.graph.field, c(list(operator = operators[[j]]), study, args)),
            error = function(e) list(error = conditionMessage(e))))
    }
    out <- if (workers == 1L) lapply(seq_along(operators), job) else
        parallel::mclapply(seq_along(operators), job, mc.cores = workers,
                          mc.preschedule = FALSE, mc.set.seed = FALSE)
    names(out) <- names(operators)
    out
}


#' Evaluate Nested Selection Among Candidate Graphs
#'
#' Chooses the graph using inner-fold loss only, then assembles its outer-fold
#' predictions. Only graphs with full inner-fold prediction coverage compete;
#' unsupported outer predictions remain NA. Candidate vertex identities and
#' ordering must agree. Final graph selection uses full-data inner CV.
#' @param comparison Output from compare.masked.graph.fields.
#' @param study Name of the study used to choose graphs (for example HMP).
#' @return Selected graph per outer fold, outer predictions and coverage,
#'   participant-balanced MSE, final graph name, and its final fitted field.
#' @usage select.masked.graph.field(comparison, study)
#' @examples
#' op <- metric.graph.lowpass.operator(list(2L, c(1L, 3L), 2L),
#'                                    list(1, c(1, 1), 1))
#' studies <- list(example = list(vertex = rep(c(1, 3), 6), y = 1:12,
#'                                group = letters[1:12]))
#' comparison <- compare.masked.graph.fields(list(path = op), studies,
#'     lambdas = c(0.01, 1), outer.folds = 3, inner.folds = 2)
#' selected <- select.masked.graph.field(comparison, "example")
#' selected$final.graph
#' selected$coverage
#' @export select.masked.graph.field
select.masked.graph.field <- function(comparison, study) {
    fits <- lapply(comparison, function(x) x[[study]])
    valid <- vapply(fits, function(x) !is.null(x$predicted) && is.null(x$error), logical(1))
    fits <- fits[valid]
    if (!length(fits)) stop("No successful graph fits for this study.")
    reference <- fits[[1L]]
    for (fit in fits) for (key in c("observed", "group", "vertex", "fold.id"))
        if (!identical(fit[[key]], reference[[key]]))
            stop("Graph fits must have identical observations, order and folds.")
    score <- function(tune) {
        allowed <- is.finite(tune$scores) & tune$coverage == 1
        if (any(allowed)) min(tune$scores[allowed]) else Inf
    }
    pred <- rep(NA_real_, length(reference$observed)); selected <- list()
    for (fold in sort(unique(reference$fold.id))) {
        scores <- vapply(fits, function(x) score(x$tuning[[as.character(fold)]]), numeric(1))
        if (!any(is.finite(scores))) stop("No graph has full inner-fold coverage; inspect support explicitly.")
        winner <- which.min(scores)
        selected[[as.character(fold)]] <- names(fits)[winner]
        test <- reference$fold.id == fold
        pred[test] <- fits[[winner]]$predicted[test]
    }
    scores <- vapply(fits, function(x) score(x$final.tuning), numeric(1))
    if (!any(is.finite(scores))) stop("No final candidate has full inner-fold coverage.")
    winner <- which.min(scores)
    list(selected.by.fold = selected, predicted = pred, coverage = mean(is.finite(pred)),
         mse = .field.score(reference$observed, pred, reference$group),
         final.graph = names(fits)[winner], final = fits[[winner]]$final)
}

#' Summarize Measurement Support for a Graph Field
#'
#' Reports graph-hop distance to the nearest observed vertex and the number of
#' distinct measured participants at each vertex and its immediate neighbors.
#' These are descriptive support measures, not uncertainty intervals or a
#' guarantee of transportability. Vertex order follows the supplied operator.
#' @param operator An unnormalized metric.graph.lowpass.operator.
#' @param vertex Integer vertex indices of observed responses.
#' @param group Nonmissing participant identifiers for the observations.
#' @return A data frame with nearest.observation.hops (Inf if unreachable),
#'   local.participants and component.supported for each vertex.
#' @examples
#' op <- metric.graph.lowpass.operator(list(2L, c(1L, 3L), 2L),
#'                                    list(1, c(1, 1), 1))
#' support.masked.graph.field(op, vertex = c(1, 3), group = c("a", "b"))
#' @export
support.masked.graph.field <- function(operator, vertex, group) {
    n <- nrow(.field.laplacian(operator))
    .field.observations(vertex, rep(0, length(vertex)), n)
    if (length(group) != length(vertex) || anyNA(group) || any(!nzchar(as.character(group))))
        stop("Supply one nonmissing group identifier per observation.")
    adj <- operator$graph$adj.list
    distance <- rep(Inf, n); sources <- unique(vertex)
    queue <- integer(n); queue[seq_along(sources)] <- sources
    head <- 1L; tail <- length(sources); distance[sources] <- 0
    while (head <= tail) {
        v <- queue[head]; head <- head + 1L
        new <- unique(adj[[v]][is.infinite(distance[adj[[v]]])])
        if (length(new)) {
            distance[new] <- distance[v] + 1
            queue[tail + seq_along(new)] <- new; tail <- tail + length(new)
        }
    }
    people <- split(as.character(group), factor(vertex, levels=seq_len(n)))
    counts <- vapply(seq_len(n), function(v)
        length(unique(unlist(people[unique(c(v, adj[[v]]))], use.names=FALSE))), integer(1))
    data.frame(nearest.observation.hops=distance, local.participants=counts,
               component.supported=is.finite(distance))
}
