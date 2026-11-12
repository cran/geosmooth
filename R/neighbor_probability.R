#' Nearest-Person and Nearest-Case Binary Probability Estimators
#'
#' Estimates the probability of a binary outcome at a point from the outcomes
#' of its nearest neighbors, with one of three estimators. The ordinary
#' nearest-neighbor estimator uses the \code{k} nearest people. The two
#' nearest-case estimators instead extend the neighborhood until it contains
#' \code{k} cases. All three are computed by the same single-threaded C++
#' code, which handles one query point at a time and never forms the matrix
#' of all pairwise distances.
#'
#' @param X Numeric matrix with one training record per row. Supply either
#'   \code{X} or \code{distances}.
#' @param y Binary responses, with \code{NA} for unanswered outcomes.
#' @param method \code{"knn"}, \code{"case_raw"}, or \code{"case_corrected"}.
#' @param k Fixed positive integer neighborhood size. Counts people for
#'   \code{"knn"} and cases for the other methods. Corrected sizes must be at
#'   least two. Omit for cross-validation using \code{foldid}.
#' @param k.grid Optional candidate sizes for cross-validation. If omitted,
#'   use every size up to the minimum feasible training count across folds.
#' @param distances Optional square, nonnegative training dissimilarity matrix
#'   in record order. Positive infinity denotes an unreachable record.
#'   No graph or powered path distances are constructed internally.
#' @param metric Coordinate distance: \code{"euclidean"}, \code{"chebyshev"}, or
#'   \code{"linf.simplex"}. The last computes exact intrinsic Euclidean-length
#'   paths on the nonnegative max-normalized surface, using dgraphs' C++ solver.
#'   Paths can cross intermediate faces; this is not Chebyshev distance.
#' @param normalization Row normalization: \code{"none"}, \code{"sum"}, or
#'   \code{"max"}. The latter two require nonnegative rows with positive mass.
#'   Supply already-computed distances with \code{normalization = "none"}.
#'   Normalization changes coordinates, not the distance definition: Euclidean
#'   distances on max-normalized profiles are ambient chords, not intrinsic
#'   paths on the cube surface. Use \code{metric = "linf.simplex"} for those
#'   paths, with \code{normalization = "max"}. Alternatively, \code{"none"}
#'   requires nonnegative rows whose maxima equal one within
#'   \code{sqrt(.Machine$double.eps)}; admitted rows are divided by their maxima.
#'   Sum normalization is rejected for intrinsic coordinate input.
#' @param search Neighbor search: \code{"auto"} uses ANN for Euclidean
#'   coordinates and brute force otherwise; \code{"ann"} explicitly requests
#'   the ANN kd-tree; \code{"brute"} is the exhaustive reference backend.
#'   ANN requires ambient Euclidean coordinates, including normalized
#'   coordinates. Intrinsic surface search always exhausts eligible candidates
#'   using one query's workspace, with no all-pairs distance matrix.
#' @param eps Nonnegative ANN approximation tolerance, default zero (exact).
#'   Positive values require ANN. The distance approximation guarantee does
#'   not bound probability error or case stopping ranks after exclusions.
#' @param id Unique nonmissing record IDs. Defaults to row numbers as strings.
#'   Supply persistent IDs to obtain ordering invariant to row permutations.
#' @param group Nonmissing participant IDs, defaulting to \code{id}.
#' @param foldid Optional fold vector or matrix (one repeat per column).
#'   Every participant must stay in one fold per repeat. Missing outcomes are
#'   excluded from scoring, but their groups still obey the split.
#' @param selection Minimum Brier error, largest size within one standard error
#'   (\code{"one.se"}), or smallest such size (\code{"locality.one.se"}).
#'
#' @details
#' Rank answering training records by distance, then by the SHA-256 hash of
#' \code{paste0("nearest-cases-v1|", id)}. Hash collisions use ID order. Keep only
#' the closest answering record per participant, before inspecting its outcome.
#' Exclude the query's entire participant when its group is supplied.
#'
#' If \eqn{T_k} is the rank including the kth case, the estimates are
#' \eqn{K^{-1}\sum_{j=1}^K y_{(j)}}, \eqn{k/T_k}, and
#' \eqn{(k-1)/(T_k-1)}. No clipping, pseudocount or automatic baseline fallback
#' is applied. Insufficient reachable people or cases return \code{NA}.
#' The number of cases reachable can depend on the query when a participant's
#' visits have discordant outcomes.
#'
#' ANN expands candidate neighborhoods until the requested people or cases
#' are available. Exact search completes the boundary before applying the ID
#' tie rule. Approximate search sorts retrieved candidates by exact distance;
#' for each retrieved participant it also checks all answering visits, keeping
#' the closest. It can miss other participants and tied records. Approximate
#' boundary diagnostics describe retrieved participants, not a complete shell.
#' Available people/cases diagnostics still describe the full training set.
#' Intrinsic surface search prepares sparse training profiles once per
#' prediction call or CV split and streams a complete eligible distance row
#' for each query. It uses dgraphs' shared intrinsic solver directly, without
#' constructing a sample graph, limiting paths to endpoint faces, or truncating
#' candidate neighbors. All three estimators share the same participant and
#' tie rules. Per-pair worst-case work is quadratic in the number of features
#' present in either profile. Distance workspace is linear in training records.
#'
#' Trees are reused across queries and sizes, and across participants within
#' each CV split. They are rebuilt for each prediction call and are not stored
#' in the serializable fit. CV retains one tree or sparse intrinsic profile
#' cache per split and repeat; a new prediction call rebuilds its cache.
#' Native caches are released on return or error from prediction or tuning.
#'
#' Tuning averages losses within participant, then across repeats, then across
#' participants. Its standard error is the standard deviation of those
#' participant losses divided by the square root of the number of people.
#' Candidates with any unsupported validation prediction are unselectable.
#' Numerical ties within \code{1e-14} favor larger sizes for minimum-error
#' selection. The tuning score is not an independent performance estimate;
#' use external participant-held outer folds for assessment.
#'
#' This implementation matches the sPTB study's prediction arithmetic and
#' minimum-error tuning when there is one record per participant, the same
#' distances, IDs and folds. The repeated-record policy is an explicit extension.
#' At a tied distance the prefix can split a tied shell, as in that study.
#'
#' @return A \code{neighbor.probability.fit} with selected \code{k},
#'   \code{method}, \code{backend = "C++"}, candidate scores, folds, and the
#'   \code{metric}, \code{distance.backend}, and training data needed for
#'   prediction. No training predictions are computed
#'   during fitting; \code{predict(fit)} computes leave-participant-out values.
#' @examples
#' X <- matrix(0:11, ncol = 1)
#' y <- rep(c(0, 1, 1), 4)
#' fit <- fit.neighbor.probability(X, y, method = "case_corrected", k = 2)
#' predict(fit, newdata = matrix(c(1.5, 5.5), ncol = 1))
#' tuned <- fit.neighbor.probability(X, y, method = "knn",
#'                                   foldid = rep(1:3, 4), k.grid = 1:4)
#' tuned$cv.table
#' # Intrinsic surface neighbors can differ from ambient Euclidean neighbors.
#' profiles <- rbind(c(0, 1, 0, 0), c(1, 1, 1, 1))
#' surface <- fit.neighbor.probability(profiles, c(1, 0), k = 1,
#'     metric = "linf.simplex", normalization = "max")
#' predict(surface, rbind(c(1, 0, 0, 0)))
#' @export
fit.neighbor.probability <- function(
    X = NULL, y, method = c("knn", "case_raw", "case_corrected"),
    k = NULL, k.grid = NULL, distances = NULL,
    metric = c("euclidean", "chebyshev", "linf.simplex"),
    normalization = c("none", "sum", "max"),
    id = NULL, group = NULL, foldid = NULL,
    selection = c("min", "one.se", "locality.one.se"),
    search = c("auto", "ann", "brute"), eps = 0) {

    method <- match.arg(method)
    metric <- match.arg(metric)
    normalization <- match.arg(normalization)
    selection <- match.arg(selection)
    search <- match.arg(search)
    if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps < 0 ||
        !is.finite((1 + eps)^2))
        stop("'eps' must be a finite nonnegative ANN tolerance.")
    if (search == "auto")
        search <- if (!is.null(X) && metric == "euclidean") "ann" else "brute"
    if (search == "ann" && (is.null(X) || metric != "euclidean"))
        stop("ANN search requires Euclidean coordinates.")
    if (search != "ann" && eps > 0)
        stop("Positive 'eps' requires ANN search.")
    if (!(is.numeric(y) || is.logical(y)) || !is.null(dim(y)) || !length(y) ||
        any(is.nan(y)) || any(!is.na(y) & (!is.finite(y) | !y %in% c(0, 1))))
        stop("'y' must be a binary vector with NA for missing responses.")
    y <- as.numeric(y)
    n <- length(y)
    if (all(is.na(y))) stop("At least one observed response is required.")
    if (is.null(X) == is.null(distances))
        stop("Supply exactly one of 'X' and 'distances'.")
    if (!is.null(X)) {
        X <- .neighbor.coordinates(X, normalization, metric)
        if (nrow(X) != n) stop("nrow(X) must equal length(y).")
    } else {
        if (normalization != "none")
            stop("Precomputed distances require normalization = 'none'.")
        distances <- .neighbor.distances(distances, n)
        if (nrow(distances) != n) stop("Training distances must be square.")
    }
    if (is.null(id)) id <- as.character(seq_len(n))
    id <- .neighbor.ids(id, n, "id")
    if (anyDuplicated(id)) stop("'id' must contain unique record IDs.")
    if (is.null(group)) group <- id
    group <- .neighbor.ids(group, n, "group")
    hashes <- vapply(id, function(s) digest::digest(
        paste0("nearest-cases-v1|", s), algo = "sha256", serialize = FALSE), "")
    order.id <- order(hashes, id, method = "radix")
    priority <- integer(n)
    priority[order.id] <- seq_len(n) - 1L
    object <- list(X = X, distances = distances, y = y, id = id, group = group,
                   group.code = match(group, unique(group)),
                   priority = priority, method = method, metric = metric,
                   normalization = normalization, backend = "C++",
                   search = search, eps = eps,
                   distance.backend = if (is.null(X)) "precomputed" else
                       if (metric == "linf.simplex") "dgraphs" else "geosmooth",
                   tie.rule = "nearest-cases-v1| + SHA-256; ID collision tie",
                   selection = selection, k = NULL, cv.table = NULL,
                   foldid = NULL, call = match.call())
    class(object) <- "neighbor.probability.fit"
    if (!is.null(k)) {
        k <- .neighbor.sizes(k, method)
        if (length(k) != 1L) stop("Fixed 'k' must be scalar; use predict(k = ...) for a path.")
        if (!is.null(k.grid) || !is.null(foldid))
            stop("Choose fixed 'k' or cross-validation, not both.")
        object$k <- k
        return(object)
    }
    if (is.null(foldid)) stop("Supply fixed 'k' or participant-safe 'foldid'.")
    folds <- .neighbor.folds(foldid, group)
    object$foldid <- folds
    observed <- which(!is.na(y))
    splits <- lapply(seq_len(ncol(folds)), function(r) {
        out <- lapply(unique(folds[, r]), function(f) {
            observed[folds[observed, r] != f]
        })
        names(out) <- unique(folds[, r])
        out
    })
    if (is.null(k.grid)) {
        cap <- min(vapply(unlist(splits, recursive = FALSE), function(t) {
            if (method == "knn") length(unique(group[t]))
            else length(unique(group[t[y[t] == 1]]))
        }, 0L))
        lo <- if (method == "case_corrected") 2L else 1L
        if (cap < lo) stop("Insufficient training people/cases in at least one fold.")
        k.grid <- seq.int(lo, cap)
    }
    sizes <- .neighbor.sizes(k.grid, method)
    indices <- list()
    on.exit({
        for (repeat.indices in indices) for (index in repeat.indices)
            .neighbor.release(index, object$metric)
    }, add = TRUE)
    for (r in seq_along(splits)) {
        indices[[r]] <- stats::setNames(vector("list", length(splits[[r]])), names(splits[[r]]))
        for (f in names(splits[[r]]))
            indices[[r]][f] <- list(.neighbor.index(object, splits[[r]][[f]]))
    }
    total <- total2 <- numeric(length(sizes))
    valid <- rep(TRUE, length(sizes))
    people <- unique(group[observed])
    # Stream participant loss vectors; do not store an n-by-candidate CV matrix.
    for (h in people) {
        q <- observed[group[observed] == h]
        loss <- numeric(length(sizes))
        for (r in seq_len(ncol(folds))) {
            train <- splits[[r]][[folds[q[1L], r]]]
            pred <- .neighbor.at.rows(object, q, train, sizes, FALSE,
                indices[[r]][[folds[q[1L], r]]])$probability
            loss <- loss + colMeans((pred - y[q])^2) / ncol(folds)
        }
        valid <- valid & is.finite(loss)
        loss[!is.finite(loss)] <- 0
        total <- total + loss
        total2 <- total2 + loss^2
    }
    np <- length(people)
    score <- total / np
    se <- if (np > 1L) sqrt(pmax(0, (total2 - total^2 / np) / (np - 1L)) / np) else
        rep(Inf, length(sizes))
    score[!valid] <- Inf
    se[!valid] <- Inf
    if (!any(is.finite(score))) stop("No candidate has complete finite validation predictions.")
    best <- max(which(score <= min(score) + 1e-14))
    eligible <- if (selection == "min") best else
        which(is.finite(score) & score <= score[best] + se[best] + 1e-14)
    chosen <- if (selection == "locality.one.se") min(eligible) else max(eligible)
    object$k <- sizes[chosen]
    object$cv.table <- data.frame(k = sizes, brier = score, se = se,
                                  valid = valid, selected = seq_along(sizes) == chosen)
    object$cv.people <- np
    object
}

#' Predict Nearest-Person Binary Probabilities
#'
#' @param object A fitted \code{\link{fit.neighbor.probability}} object.
#' @param newdata Query coordinate matrix. With no query data, predict at the
#'   stored training rows while excluding each query's whole participant.
#' @param distances Query-by-training distance matrix for a precomputed fit.
#'   Columns must retain the original training record order.
#' @param query.group Optional participant IDs for query rows. Matching training
#'   participants are excluded. Omit only for genuinely new participants.
#' @param k One or more requested sizes; defaults to the fitted size.
#' @param type Return a probability matrix or a list including matrices of
#'   neighbor counts, case counts, radii and boundary-tie diagnostics.
#' @param ... Reserved; additional arguments are rejected.
#' @return A query-by-size matrix, or a diagnostic list with that matrix as
#'   \code{probability}. Unsupported estimates remain \code{NA}.
#'   \code{boundary.ties} counts distances within \code{1e-12} of the cutoff;
#'   ordering itself uses exact distance ties. Diagnoses count distinct people.
#'   With positive \code{eps}, boundary counts cover retrieved participants
#'   only (\code{boundary.complete = FALSE}). \code{search.candidates} gives
#'   the final number of candidate records, including visit refinement.
#' @examples
#' X <- matrix(0:5, ncol = 1)
#' f <- fit.neighbor.probability(X, c(1, 0, 1, 0, 0, 1), k = 2)
#' predict(f) # leave-participant-out, not resubstitution
#' predict(f, matrix(2.5, ncol = 1), type = "diagnostics")
#' @export
predict.neighbor.probability.fit <- function(
    object, newdata = NULL, distances = NULL, query.group = NULL,
    k = object$k, type = c("response", "diagnostics"), ...) {
    if (length(list(...))) stop("Unused prediction arguments.")
    type <- match.arg(type)
    sizes <- .neighbor.sizes(k, object$method)
    precomputed <- is.null(object$X)
    training.query <- is.null(newdata) && is.null(distances)
    if (training.query) {
        query <- object$X
        distance <- object$distances
        nq <- length(object$y)
        if (is.null(query.group)) query.group <- object$group
    } else if (precomputed) {
        if (!is.null(newdata) || is.null(distances))
            stop("A precomputed fit requires query-by-training 'distances'.")
        distance <- .neighbor.distances(distances, length(object$y))
        nq <- nrow(distance)
        query <- NULL
    } else {
        if (!is.null(distances) || is.null(newdata))
            stop("A coordinate fit requires 'newdata'.")
        query <- .neighbor.coordinates(newdata, object$normalization, object$metric)
        if (ncol(query) != ncol(object$X)) stop("Query feature count differs from training.")
        nq <- nrow(query)
        distance <- NULL
    }
    group.code <- integer(nq)
    if (!is.null(query.group)) {
        query.group <- .neighbor.ids(query.group, nq, "query.group")
        group.code <- match(query.group, unique(object$group), nomatch = 0L)
    }
    result <- .neighbor.call(object, query, distance, group.code,
                             which(!is.na(object$y)), sizes, type == "diagnostics")
    colnames(result$probability) <- as.character(sizes)
    if (training.query) rownames(result$probability) <- object$id
    if (type == "response") return(result$probability)
    result$status <- ifelse(is.finite(result$probability), "ok",
                            if (object$method == "knn") "insufficient_people" else "insufficient_cases")
    result$k <- sizes
    result$method <- object$method
    result$metric <- object$metric
    result$distance.backend <- object$distance.backend
    result$search <- if (is.null(object$search)) "brute" else object$search
    result$eps <- if (is.null(object$eps)) 0 else object$eps
    result$boundary.complete <- result$eps == 0
    result
}

.neighbor.at.rows <- function(object, rows, train, sizes, diagnostics, index = NULL) {
    .neighbor.call(object,
                   if (!is.null(object$X)) object$X[rows, , drop = FALSE] else NULL,
                   if (!is.null(object$distances)) object$distances[rows, , drop = FALSE] else NULL,
                   object$group.code[rows], train, sizes, diagnostics, index)
}

.neighbor.index <- function(object, train) {
    if (is.null(object$X) || !length(train)) return(NULL)
    if (identical(object$metric, "linf.simplex"))
        return(rcpp_neighbor_intrinsic_index(object$X, as.integer(train - 1L)))
    if (!identical(object$search, "ann")) return(NULL)
    rcpp_neighbor_ann_index(object$X, as.integer(train - 1L), object$priority,
                            as.integer(object$group.code), object$y)
}

.neighbor.release <- function(index, metric) {
    if (!is.null(index)) rcpp_neighbor_release_index(index, identical(metric, "linf.simplex"))
    invisible(NULL)
}

.neighbor.call <- function(object, query, distance, query.group, train, sizes,
                           diagnostics, index = NULL) {
    if (is.null(index)) {
        index <- .neighbor.index(object, train)
        on.exit(.neighbor.release(index, object$metric), add = TRUE)
    }
    empty <- matrix(numeric(), 0L, 0L)
    rcpp_neighbor_probability(
        if (is.null(object$X)) empty else object$X,
        if (is.null(query)) empty else query,
        if (is.null(distance)) empty else distance,
        object$y, as.integer(train - 1L), as.integer(object$group.code),
        as.integer(query.group), as.integer(object$priority), sizes,
        match(object$method, c("knn", "case_raw", "case_corrected")) - 1L,
        match(object$metric, c("euclidean", "chebyshev", "linf.simplex")) - 1L,
        is.null(object$X), diagnostics, index,
        if (is.null(object$eps)) 0 else object$eps)
}

.neighbor.ids <- function(x, n, name) {
    if (!is.atomic(x) || !is.null(dim(x)) || length(x) != n || anyNA(x))
        stop("'", name, "' must contain one nonmissing ID per row.")
    x <- enc2utf8(as.character(x))
    if (any(!nzchar(x))) stop("'", name, "' cannot contain empty IDs.")
    x
}

.neighbor.coordinates <- function(X, normalization, metric = "euclidean") {
    if (!is.matrix(X) || !is.numeric(X) || ncol(X) < 1L || nrow(X) < 1L ||
        any(!is.finite(X))) stop("Coordinates must be a finite numeric matrix.")
    storage.mode(X) <- "double"
    if (metric == "linf.simplex") {
        if (normalization == "sum")
            stop("Intrinsic surface distances require normalization = 'max' or already max-normalized coordinates.")
        if (any(X < 0)) stop("Intrinsic surface coordinates must be nonnegative.")
        scale <- apply(X, 1L, max)
        if (any(scale <= 0)) stop("Intrinsic surface rows must have positive maxima.")
        if (normalization == "none" && any(abs(scale - 1) > sqrt(.Machine$double.eps)))
            stop("Intrinsic surface coordinates must be max-normalized; use normalization = 'max'.")
        return(X / scale)
    }
    if (normalization != "none") {
        if (any(X < 0)) stop("Row normalization requires nonnegative coordinates.")
        scale <- if (normalization == "sum") rowSums(X) else apply(X, 1L, max)
        if (any(!is.finite(scale) | scale <= 0))
            stop("Normalized rows must have finite positive mass.")
        X <- X / scale
    }
    X
}

.neighbor.distances <- function(D, n) {
    if (!is.matrix(D) || !is.numeric(D) || ncol(D) != n || nrow(D) < 1L ||
        anyNA(D) || any(D < 0))
        stop("Distances must be a nonnegative numeric query-by-training matrix; Inf is allowed.")
    storage.mode(D) <- "double"
    D
}

.neighbor.sizes <- function(k, method) {
    lo <- if (method == "case_corrected") 2L else 1L
    if (!is.numeric(k) || !length(k) || any(!is.finite(k)) ||
        any(k != floor(k) | k < lo | k > .Machine$integer.max))
        stop("'k' must contain integers >= ", lo, ".")
    sort(unique(as.integer(k)))
}

.neighbor.folds <- function(folds, group) {
    if (!is.atomic(folds) || anyNA(folds)) stop("Fold IDs must be nonmissing.")
    if (is.null(dim(folds))) folds <- matrix(folds, ncol = 1L)
    if (!is.matrix(folds) || nrow(folds) != length(group) || ncol(folds) < 1L)
        stop("'foldid' must have one row per training record.")
    folds <- matrix(as.character(folds), nrow(folds), ncol(folds))
    if (any(!nzchar(folds))) stop("Fold IDs cannot be empty.")
    for (r in seq_len(ncol(folds))) {
        if (length(unique(folds[,r])) < 2L) stop("Each repeat requires at least two folds.")
        if (any(vapply(split(folds[,r], group), function(v) length(unique(v)) != 1L, TRUE)))
            stop("A participant cannot cross folds.")
    }
    folds
}
