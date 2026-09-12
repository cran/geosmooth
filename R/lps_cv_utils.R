#' Build a Cluster-Respecting Fold Assignment
#'
#' Assigns whole clusters to cross-validation folds so that no cluster is
#' split across folds: every observation of a cluster receives the same fold
#' id, and therefore no training set ever shares a cluster with the fold held
#' out against it.  Folds are built by a deterministic size-balanced greedy
#' rule: clusters are taken largest first (ties by first appearance in
#' \code{cluster.id}) and each is placed into the currently smallest fold
#' (ties by lowest fold index).  With \code{v} equal to the number of
#' clusters this reduces to leave-cluster-out.
#'
#' The assignment is fully deterministic when \code{shuffle.seed} is
#' \code{NULL}.  Supplying \code{shuffle.seed} applies one seeded permutation
#' to the cluster order before the greedy pass (the seed is consumed
#' immediately via \code{set.seed}), giving a reproducible randomized
#' variant; callers following the LPS evidence conventions should record the
#' seed they pass.
#'
#' @param cluster.id Vector of cluster labels (factor, character, or
#'   integer-like), one per observation, no missing values.
#' @param v Number of folds; an integer between 2 and the number of distinct
#'   clusters.
#' @param shuffle.seed Optional integer seed for a randomized cluster order;
#'   \code{NULL} (default) keeps the deterministic order.
#' @return An integer fold-id vector of length \code{length(cluster.id)} with
#'   values in \code{1:v}; every fold is nonempty and every cluster maps to
#'   exactly one fold.
#' @seealso [lps.nested.cv()] for nested cross-validation that can consume
#'   grouped folds at both the outer and inner level.
#' @examples
#' cluster <- rep(letters[1:6], each = 2)
#' lps.grouped.foldid(cluster, v = 3L)
#' @export
lps.grouped.foldid <- function(cluster.id, v = 5L, shuffle.seed = NULL) {
    if (!length(cluster.id) || anyNA(cluster.id)) {
        stop("'cluster.id' must be a non-empty vector without missing values.",
             call. = FALSE)
    }
    cluster.chr <- as.character(cluster.id)
    cluster.names <- unique(cluster.chr)
    n.clusters <- length(cluster.names)
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) ||
        v != round(v)) {
        stop("'v' must be a single whole number.", call. = FALSE)
    }
    v <- as.integer(v)
    if (v < 2L || v > n.clusters) {
        stop("'v' must be an integer between 2 and the number of distinct ",
             "clusters (", n.clusters, ").", call. = FALSE)
    }
    sizes <- as.integer(table(factor(cluster.chr, levels = cluster.names)))
    order.idx <- seq_len(n.clusters)
    if (!is.null(shuffle.seed)) {
        set.seed(shuffle.seed)
        order.idx <- sample.int(n.clusters)
    }
    # Stable sort: size descending, ties by (possibly shuffled) order.
    take <- order.idx[order(-sizes[order.idx])]
    fold.of.cluster <- integer(n.clusters)
    fold.load <- integer(v)
    for (cl in take) {
        target <- which.min(fold.load)
        fold.of.cluster[[cl]] <- target
        fold.load[[target]] <- fold.load[[target]] + sizes[[cl]]
    }
    foldid <- fold.of.cluster[match(cluster.chr, cluster.names)]
    as.integer(foldid)
}

#' Nested Cross-Validation for LPS with Explicit, Recorded Folds
#'
#' Runs outer-fold nested cross-validation around [fit.lps()]: for each outer
#' fold, candidate selection is one ordinary \code{fit.lps} call on the
#' inner-training rows only, with an explicit inner fold id and
#' \code{X.eval} set to the held-out outer-test rows, so the held-out fold
#' never participates in inner selection.  The pooled outer-test error is the
#' nested generalization estimate.  The same call also computes the
#' \emph{selected-min} arm — an ordinary \code{fit.lps} on all rows using the
#' \emph{same} \code{outer.foldid} — so optimism comparisons are paired on
#' one fold assignment by construction.
#'
#' The existing \code{fit.lps} behavior is untouched: this utility consumes
#' the public API with explicit \code{foldid} only.  It currently supports
#' \code{outcome.family = "gaussian"} (the default) and scores with RMSE.
#'
#' Inner folds are built per outer fold and fully recorded in the return
#' value: \code{"round.robin"} assigns \code{rep_len(1:inner.folds, n)} over
#' the inner-training rows in row order (deterministic; if
#' \code{inner.shuffle.seed} is supplied, fold \code{k} — by position — uses
#' \code{set.seed(inner.shuffle.seed + k)} and permutes the balanced
#' assignment, giving a reproducible randomized variant);
#' \code{"grouped"} builds cluster-respecting inner folds with
#' [lps.grouped.foldid()] on the inner-training rows of \code{cluster.id}
#' (passing the same per-fold offset seed when \code{inner.shuffle.seed} is
#' supplied).
#'
#' @param X Numeric matrix of observations (rows) used for training.
#' @param y Numeric response vector with \code{length(y) == nrow(X)}.
#' @param outer.foldid Positive integer vector of length \code{nrow(X)}
#'   assigning rows to outer folds (at least two distinct folds).
#' @param fit.args Named list of additional arguments passed to every
#'   \code{fit.lps} call (grids, backend, design settings, ...).  Must not
#'   contain \code{X}, \code{y}, \code{foldid}, or \code{X.eval}: the fold
#'   plumbing is owned by this function so both arms see the same folds.
#' @param inner.folds Number of inner selection folds (integer >= 2).
#' @param cluster.id Optional cluster labels (length \code{nrow(X)}),
#'   required for \code{inner.foldid.method = "grouped"}.
#' @param inner.foldid.method Inner fold construction: \code{"round.robin"}
#'   (default) or \code{"grouped"} (cluster-respecting).
#' @param inner.shuffle.seed Optional integer; per-fold offset seed for the
#'   randomized variants described above.  Recorded in the return value.
#' @return A list of class \code{"lps_nested_cv"}:
#'   \describe{
#'     \item{\code{nested.rmse}}{pooled RMSE of the outer-test predictions
#'       over all rows with finite predictions.}
#'     \item{\code{n.missing.predictions}}{count of non-finite outer-test
#'       predictions (e.g. \code{unstable.action = "na"} fallbacks).}
#'     \item{\code{folds}}{one row per outer fold: fold label, test size,
#'       selected \code{support.size} / \code{degree} / \code{kernel} /
#'       \code{bandwidth.multiplier}, the inner selected-min CV score, and
#'       the fold's outer-test RMSE.}
#'     \item{\code{predictions}}{length-\code{nrow(X)} vector of outer-test
#'       predictions (each row predicted by the fold that held it out).}
#'     \item{\code{selected.min}}{the paired selected-min arm on the same
#'       \code{outer.foldid}: \code{selected} (the \code{fit.lps} selected
#'       row), \code{cv.score} (its observed CV RMSE), \code{foldid}, and
#'       \code{fit} (the full-data \code{fit.lps} object, usable for
#'       deployment on an external test set).}
#'     \item{\code{outer.foldid}, \code{test.index}, \code{train.index},
#'       \code{inner.foldid}, \code{inner.foldid.used},
#'       \code{inner.cv.table}}{complete fold/index telemetry: the realized
#'       inner fold ids as constructed and as recorded by each inner
#'       \code{fit.lps} (\code{$foldid}), per-fold index sets, and each
#'       fold's full inner CV candidate table.}
#'     \item{\code{outer.cluster.whole}}{\code{TRUE}/\code{FALSE} whether
#'       every cluster lies wholly inside one outer fold (\code{NA} when
#'       \code{cluster.id} is absent).}
#'     \item{\code{inner.folds}, \code{inner.foldid.method},
#'       \code{inner.shuffle.seed}, \code{fit.args}, \code{call}}{the
#'       recorded configuration.}
#'   }
#' @seealso [lps.grouped.foldid()] for the grouped fold constructor.
#' @examples
#' X <- matrix(seq(0, 1, length.out = 18), ncol = 1)
#' y <- sin(2 * pi * X[, 1])
#' nested <- lps.nested.cv(
#'   X, y, outer.foldid = rep(1:3, length.out = 18), inner.folds = 2L,
#'   fit.args = list(support.grid = 6L, degree.grid = 1L,
#'                   kernel.grid = "gaussian", backend = "R")
#' )
#' nested$rmse
#' @export
lps.nested.cv <- function(X, y, outer.foldid,
                          fit.args = list(),
                          inner.folds = 5L,
                          cluster.id = NULL,
                          inner.foldid.method = c("round.robin", "grouped"),
                          inner.shuffle.seed = NULL) {
    X <- as.matrix(X)
    y <- as.numeric(y)
    n <- nrow(X)
    if (length(y) != n) {
        stop("'y' must have length nrow(X).", call. = FALSE)
    }
    if (!is.numeric(outer.foldid) || length(outer.foldid) != n ||
        anyNA(outer.foldid) || any(outer.foldid != as.integer(outer.foldid)) ||
        any(outer.foldid < 1L)) {
        stop("'outer.foldid' must be a positive integer vector of length ",
             "nrow(X).", call. = FALSE)
    }
    outer.foldid <- as.integer(outer.foldid)
    fold.labels <- sort(unique(outer.foldid))
    if (length(fold.labels) < 2L) {
        stop("'outer.foldid' must define at least two outer folds.",
             call. = FALSE)
    }
    if (!is.list(fit.args)) {
        stop("'fit.args' must be a named list.", call. = FALSE)
    }
    if (length(fit.args) &&
        (is.null(names(fit.args)) || any(!nzchar(names(fit.args))))) {
        stop("'fit.args' must be a named list.", call. = FALSE)
    }
    reserved <- intersect(names(fit.args), c("X", "y", "foldid", "X.eval"))
    if (length(reserved)) {
        stop("'fit.args' must not contain ",
             paste(sQuote(reserved), collapse = ", "),
             ": fold and data plumbing is owned by lps.nested.cv so both ",
             "arms of a comparison see the same folds.", call. = FALSE)
    }
    outcome.family <- fit.args[["outcome.family"]]
    if (!is.null(outcome.family) && !identical(outcome.family, "gaussian")) {
        stop("lps.nested.cv currently supports outcome.family = 'gaussian' ",
             "only.", call. = FALSE)
    }
    if (length(inner.folds) != 1L || !is.numeric(inner.folds) ||
        !is.finite(inner.folds) || inner.folds != round(inner.folds)) {
        stop("'inner.folds' must be a single whole number.", call. = FALSE)
    }
    inner.folds <- as.integer(inner.folds)
    if (inner.folds < 2L) {
        stop("'inner.folds' must be an integer >= 2.", call. = FALSE)
    }
    inner.foldid.method <- match.arg(inner.foldid.method)
    if (identical(inner.foldid.method, "grouped") && is.null(cluster.id)) {
        stop("inner.foldid.method = 'grouped' requires 'cluster.id'.",
             call. = FALSE)
    }
    if (!is.null(cluster.id)) {
        if (length(cluster.id) != n || anyNA(cluster.id)) {
            stop("'cluster.id' must have length nrow(X) without missing ",
                 "values.", call. = FALSE)
        }
        cluster.id <- as.character(cluster.id)
    }

    rmse <- function(pred, obs) {
        ok <- is.finite(pred) & is.finite(obs)
        if (!any(ok)) return(NA_real_)
        sqrt(mean((pred[ok] - obs[ok])^2))
    }
    selected.scalar <- function(selected, name, default) {
        value <- selected[[name]]
        if (is.null(value) || !length(value)) default else value[[1L]]
    }

    n.folds <- length(fold.labels)
    predictions <- rep(NA_real_, n)
    test.index <- vector("list", n.folds)
    train.index <- vector("list", n.folds)
    inner.foldid <- vector("list", n.folds)
    inner.foldid.used <- vector("list", n.folds)
    inner.cv.table <- vector("list", n.folds)
    fold.rows <- vector("list", n.folds)
    names(test.index) <- names(train.index) <- names(inner.foldid) <-
        names(inner.foldid.used) <- names(inner.cv.table) <-
        as.character(fold.labels)

    for (pos in seq_len(n.folds)) {
        label <- fold.labels[[pos]]
        test.idx <- which(outer.foldid == label)
        train.idx <- which(outer.foldid != label)
        n.train <- length(train.idx)
        if (n.train < inner.folds) {
            stop("Outer fold ", label, " leaves ", n.train,
                 " training rows, fewer than inner.folds = ", inner.folds,
                 ".", call. = FALSE)
        }
        if (identical(inner.foldid.method, "round.robin")) {
            assignment <- rep_len(seq_len(inner.folds), n.train)
            if (!is.null(inner.shuffle.seed)) {
                set.seed(inner.shuffle.seed + pos)
                assignment <- sample(assignment)
            }
            fold.inner <- as.integer(assignment)
        } else {
            seed.k <- if (is.null(inner.shuffle.seed)) {
                NULL
            } else {
                inner.shuffle.seed + pos
            }
            fold.inner <- lps.grouped.foldid(
                cluster.id = cluster.id[train.idx],
                v = inner.folds,
                shuffle.seed = seed.k
            )
        }
        fit.k <- do.call(fit.lps, c(
            list(
                X = X[train.idx, , drop = FALSE],
                y = y[train.idx],
                foldid = fold.inner,
                X.eval = X[test.idx, , drop = FALSE]
            ),
            fit.args
        ))
        predictions[test.idx] <- fit.k$fitted.values
        test.index[[pos]] <- test.idx
        train.index[[pos]] <- train.idx
        inner.foldid[[pos]] <- fold.inner
        inner.foldid.used[[pos]] <- fit.k$foldid
        inner.cv.table[[pos]] <- fit.k$cv.table
        fold.rows[[pos]] <- data.frame(
            fold = label,
            n.test = length(test.idx),
            selected.support.size = as.integer(
                selected.scalar(fit.k$selected, "support.size", NA_integer_)
            ),
            selected.degree = as.integer(
                selected.scalar(fit.k$selected, "degree", NA_integer_)
            ),
            selected.kernel = as.character(
                selected.scalar(fit.k$selected, "kernel", NA_character_)
            ),
            selected.bandwidth.multiplier = as.numeric(
                selected.scalar(fit.k$selected, "bandwidth.multiplier", 1)
            ),
            inner.cv.rmse = as.numeric(
                selected.scalar(fit.k$selected, "cv.rmse.observed", NA_real_)
            ),
            test.rmse = rmse(predictions[test.idx], y[test.idx]),
            stringsAsFactors = FALSE
        )
    }

    fit.full <- do.call(fit.lps, c(
        list(X = X, y = y, foldid = outer.foldid),
        fit.args
    ))

    outer.cluster.whole <- if (is.null(cluster.id)) {
        NA
    } else {
        all(vapply(
            split(outer.foldid, cluster.id),
            function(f) length(unique(f)) == 1L,
            logical(1L)
        ))
    }

    out <- list(
        nested.rmse = rmse(predictions, y),
        n.missing.predictions = sum(!is.finite(predictions)),
        folds = do.call(rbind, fold.rows),
        predictions = predictions,
        selected.min = list(
            selected = fit.full$selected,
            cv.score = selected.scalar(fit.full$selected,
                                       "cv.rmse.observed", NA_real_),
            foldid = fit.full$foldid,
            fit = fit.full
        ),
        outer.foldid = outer.foldid,
        test.index = test.index,
        train.index = train.index,
        inner.foldid = inner.foldid,
        inner.foldid.used = inner.foldid.used,
        inner.cv.table = inner.cv.table,
        outer.cluster.whole = outer.cluster.whole,
        inner.folds = inner.folds,
        inner.foldid.method = inner.foldid.method,
        inner.shuffle.seed = inner.shuffle.seed,
        fit.args = fit.args,
        call = match.call()
    )
    class(out) <- c("lps_nested_cv", "list")
    out
}
