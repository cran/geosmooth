#' Fit Prediction-Synchronized Local Polynomial Smoothing
#'
#' Experimental R reference implementation of prediction-synchronized local
#' polynomial smoothing (PS-LPS).  For fixed local polynomial parameters, PS-LPS
#' fits one local polynomial chart per anchor and adds a quadratic penalty that
#' synchronizes chart predictions on overlap points.
#'
#' @param X Numeric covariate matrix.
#' @param y Numeric response vector.
#' @param foldid Optional integer cross-validation fold assignment.
#' @param support.size Optional single neighborhood size for the fixed local
#'   setup path. Use \code{support.grid} for CV selection over neighborhood
#'   sizes.
#' @param degree Single local polynomial degree for the fixed local setup path.
#'   Use \code{degree.grid} for CV selection over degrees.
#' @param kernel Single kernel name for the fixed local setup path. Use
#'   \code{kernel.grid} for CV selection over kernels.
#' @param chart.dim Chart dimension for the local PCA charts. A scalar fixes
#'   one dimension for all anchors; an integer vector supplies one local chart
#'   dimension per anchor. The special value \code{"auto"} estimates one global
#'   chart dimension from observed \code{X} only. The experimental special
#'   value \code{"local.auto"} estimates one local chart dimension per anchor,
#'   using the same local-PCA dimension rule as \code{\link{fit.lps}}.
#' @param support.grid Candidate neighborhood sizes. When supplied, PS-LPS
#'   selects over support size, degree, kernel, and synchronization strength by
#'   materialized-fold CV. If both \code{support.size} and \code{support.grid}
#'   are absent, defaults to \code{c(10L, 15L, 20L)}, matching
#'   \code{\link{fit.lps}}.
#' @param degree.grid Candidate local polynomial degrees for grid selection.
#'   Defaults to the scalar \code{degree}.
#' @param kernel.grid Candidate kernels for grid selection. Defaults to the
#'   scalar \code{kernel}.
#' @param auto.chart.support.metric Support system used when
#'   \code{chart.dim = "auto"} or \code{"local.auto"}. Because PS-LPS uses
#'   coordinate supports, \code{"operator"} is equivalent to
#'   \code{"coordinates"}.
#' @param auto.chart.selection.metric Which auto chart-dimension diagnostic to
#'   select when both diagnostics are requested.
#' @param chart.dim.grid Optional candidate chart dimensions for experimental
#'   local-candidate selection. Numeric grids can be evaluated with
#'   \code{selection.strategy = "sparse_kd"}.
#' @param selection.strategy Candidate-selection strategy. \code{"grid"}
#'   evaluates the requested candidate grid. \code{"sparse_kd"} evaluates a
#'   sparse support-size by chart-dimension skeleton when
#'   \code{chart.dim.grid} is supplied. \code{"plateau_kd"} uses the same
#'   geometry-only support-size and chart-dimension plateau rule as
#'   \code{\link{fit.lps}}.
#' @param chart.dim.max Optional explicit maximum chart dimension for the
#'   sparse coupled candidate family.
#' @param design.margin Nonnegative integer feasibility margin used to screen
#'   local polynomial design size before candidate evaluation.
#' @param lambda.sync.grid Candidate synchronization strengths.
#' @param lambda.sync.search Lambda-search policy.  \code{"grid"} evaluates the
#'   supplied grid exactly.  \code{"guarded"} uses an experimental guarded
#'   coarse-to-refine search with boundary expansion.
#' @param lambda.sync.selection Lambda-selection mode. \code{"cv"} performs the
#'   usual internal fold-weighted lambda selection. \code{"fixed"} treats the
#'   single value in \code{lambda.sync.grid} as preselected and skips internal
#'   lambda CV; this is intended for outer selection loops that already score
#'   scalar candidates.
#' @param local.candidate.search Local-candidate search policy when
#'   \code{support.grid}, \code{degree.grid}, or \code{kernel.grid} contains
#'   more than one candidate. \code{"screened"} is the routine default: it first
#'   ranks candidates by ordinary LPS materialized-fold CV, then runs PS-LPS only
#'   on a screened subset plus guard candidates. \code{"full"} evaluates PS-LPS
#'   lambda search for every local candidate and is the exact audit/reference
#'   path. \code{"subgrid"} skips the ordinary-LPS screening pass and evaluates
#'   only a deterministic support/kernel guard subgrid; this is intended for
#'   high-dimensional preflight runs where the screening pass is itself too
#'   expensive.
#' @param local.candidate.search.control Optional list controlling screened
#'   local-candidate search. Supported fields are \code{top.n} (default
#'   \code{8}), \code{max.candidates} (default \code{12}),
#'   \code{neighbor.radius} (default \code{1}), and
#'   \code{guard.support.quantiles} (default \code{c(0, 0.5, 1)}).
#' @param lambda.sync.search.control Optional list controlling guarded search.
#'   Supported fields are \code{coarse.size} (default \code{5}),
#'   \code{refine.radius} (default \code{2}), \code{rel.tol} (default
#'   \code{0.002}), \code{boundary.guard.rel.tol} (default \code{0.01}),
#'   \code{boundary.expand} (default \code{TRUE}), \code{boundary.factor}
#'   (default \code{3}), \code{max.boundary.expansions} (default \code{2}),
#'   and \code{max.candidates} (default \code{25}).  Boundary expansion may
#'   evaluate positive candidates outside the supplied \code{lambda.sync.grid}.
#'   \code{max.candidates} is a global cap on distinct evaluated candidates; a
#'   very small cap can prevent local refinement or boundary expansion.
#' @param lambda.diagnostics Lambda-diagnostic evaluation policy.
#'   \code{"all"} computes synchronization-energy and local-GCV diagnostics for
#'   every evaluated lambda and is the audit/reference mode. \code{"selected"}
#'   computes those diagnostics only for the selected lambda, which is faster
#'   for routine selection.
#' @param lambda.ridge Nonnegative scale-relative ridge used in the chart
#'   coefficient solve. Use \code{0} for the unregularized least-squares model.
#' @param design.basis Local polynomial design backend. See
#'   \code{\link{fit.lps}}. In PS-LPS, \code{"weighted.qr.drop"} drops
#'   numerically dependent columns separately in each anchor chart before the
#'   synchronized system is assembled, and
#'   \code{"orthogonal.polynomial.drop"} builds each synchronized chart in a
#'   weighted-orthogonal polynomial basis.
#' @param design.drop.tol Relative QR tolerance used by
#'   \code{design.basis = "weighted.qr.drop"} or
#'   \code{design.basis = "orthogonal.polynomial.drop"}.
#' @param ridge.multiplier.grid Optional nonnegative ridge multipliers for
#'   adaptive scale-relative ridge selection. When supplied, the solver uses the
#'   smallest multiplier whose penalized normal equations pass
#'   \code{ridge.condition.max}. If \code{NULL}, \code{lambda.ridge} is used as
#'   the single multiplier for backward compatibility.
#' @param ridge.condition.max Maximum allowed condition number for adaptive
#'   ridge selection. Use \code{Inf} to disable the condition-number guard.
#' @param sync.neighbor.size Number of nearby anchor pairs considered for
#'   synchronization from each anchor support.
#' @param overlap.weight Overlap weighting rule.
#' @param chart.activation Optional sparse-response chart activation rule.
#'   \code{"none"} preserves ordinary PS-LPS. \code{"subject.od"} is intended
#'   for subject-occupation density workflows and deactivates local charts with
#'   no subject mass, too few positive subject-visited support points, or only
#'   fringe occupation.
#' @param chart.activation.response Optional nonnegative vector used only for
#'   \code{chart.activation = "subject.od"} to decide whether each chart is
#'   active. When omitted, \code{y} is used.
#' @param chart.activation.control List controlling sparse chart activation;
#'   see \code{\link{fit.lps}}.
#' @param ps.lps.geometry.cache Optional precomputed geometry cache used by
#'   occupation-density cross-validation. This is an implementation detail;
#'   ordinary callers should leave it as \code{NULL}.
#' @param ps.lps.local.pca.supports Optional precomputed local-PCA supports used
#'   when activation-specific frames must be rebuilt. This is an implementation
#'   detail; ordinary callers should leave it as \code{NULL}.
#' @param cv.folds Number of folds when \code{foldid} is absent.
#' @param cv.seed Fold seed when \code{foldid} is absent.
#' @return A list with fitted values, selected lambda, CV table, diagnostics,
#'   and fitted chart coefficients.
#' @examples
#' X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
#' fit.ps.lps(
#'   X, sin(2 * pi * X[, 1]), support.size = 7L, degree = 1L,
#'   chart.dim = 1L, lambda.sync.grid = 0,
#'   lambda.sync.selection = "fixed", cv.folds = 2L
#' )
#' @export
fit.ps.lps <- function(
    X, y, foldid = NULL,
    support.size = NULL,
    degree = 2L,
    kernel = "gaussian",
    chart.dim = NULL,
    support.grid = NULL,
    degree.grid = NULL,
    kernel.grid = NULL,
    auto.chart.support.metric = c("coordinates", "operator", "both"),
    auto.chart.selection.metric = c("coordinates", "operator"),
    chart.dim.grid = NULL,
    selection.strategy = c("grid", "sparse_kd", "plateau_kd"),
    chart.dim.max = NULL,
    design.margin = 2L,
    lambda.sync.grid = c(0, 1e-3, 1e-2, 1e-1, 1, 10),
    lambda.sync.search = c("grid", "guarded"),
    lambda.sync.selection = c("cv", "fixed"),
    local.candidate.search = c("screened", "full", "subgrid"),
    local.candidate.search.control = list(),
    lambda.sync.search.control = list(),
    lambda.diagnostics = c("all", "selected"),
    lambda.ridge = 1e-8,
    design.basis = c("monomial", "weighted.qr", "weighted.qr.drop",
                     "orthogonal.polynomial.drop"),
    design.drop.tol = sqrt(.Machine$double.eps),
    ridge.multiplier.grid = NULL,
    ridge.condition.max = Inf,
    sync.neighbor.size = NULL,
    overlap.weight = c("normalized.product", "product"),
    chart.activation = c("none", "subject.od"),
    chart.activation.response = NULL,
    chart.activation.control = list(),
    ps.lps.geometry.cache = NULL,
    ps.lps.local.pca.supports = NULL,
    cv.folds = 5L,
    cv.seed = 1L) {

    X <- as.matrix(X)
    y <- as.numeric(y)
    if (!is.numeric(X) || !length(X) || any(!is.finite(X))) {
        stop("'X' must be a finite numeric matrix.", call. = FALSE)
    }
    if (length(y) != nrow(X) || any(!is.finite(y))) {
        stop("'y' must be finite and have length nrow(X).", call. = FALSE)
    }
    auto.chart.support.metric <- match.arg(auto.chart.support.metric)
    auto.chart.selection.metric <- match.arg(auto.chart.selection.metric)
    selection.strategy <- .coupled.kd.selection.strategy(selection.strategy)
    lambda.sync.search <- match.arg(lambda.sync.search)
    lambda.sync.selection <- match.arg(lambda.sync.selection)
    local.candidate.search <- match.arg(local.candidate.search)
    if (!is.null(chart.dim.grid) &&
        identical(local.candidate.search, "screened")) {
        local.candidate.search <- "full"
    }
    lambda.diagnostics <- match.arg(lambda.diagnostics)
    overlap.weight <- match.arg(overlap.weight)
    chart.activation <- match.arg(chart.activation)
    chart.activation.info <- .klp.prepare.chart.activation(
        chart.activation = chart.activation,
        chart.activation.response = chart.activation.response,
        fallback.response = y,
        n = nrow(X),
        control = chart.activation.control
    )
    design.basis <- match.arg(design.basis)
    design.drop.tol <- .klp.validate.nonnegative.scalar(
        design.drop.tol,
        "design.drop.tol"
    )
    ridge.condition.max <- .klp.validate.positive.scalar(
        ridge.condition.max,
        "ridge.condition.max",
        allow.infinite = TRUE
    )
    lambda.sync.grid <- sort(unique(as.numeric(lambda.sync.grid)))
    if (!length(lambda.sync.grid) || any(!is.finite(lambda.sync.grid)) ||
        any(lambda.sync.grid < 0)) {
        stop("'lambda.sync.grid' must contain finite nonnegative values.",
             call. = FALSE)
    }
    if (identical(lambda.sync.selection, "fixed") &&
        length(lambda.sync.grid) != 1L) {
        stop("'lambda.sync.selection = \"fixed\"' requires exactly one ",
             "lambda in 'lambda.sync.grid'.", call. = FALSE)
    }
    lambda.ridge <- as.numeric(lambda.ridge[[1L]])
    if (!is.finite(lambda.ridge) || lambda.ridge < 0) {
        stop("'lambda.ridge' must be a finite nonnegative scalar.",
             call. = FALSE)
    }
    if (is.null(ridge.multiplier.grid)) {
        ridge.multiplier.grid <- lambda.ridge
    }
    ridge.multiplier.grid <- .klp.clean.ridge.multiplier.grid(
        ridge.multiplier.grid
    )
    design.margin <- as.integer(design.margin)
    if (length(design.margin) != 1L || !is.finite(design.margin) ||
        design.margin < 0L) {
        stop("'design.margin' must be a nonnegative integer scalar.",
             call. = FALSE)
    }
    if (is.null(chart.dim) && is.null(chart.dim.grid)) {
        stop("Supply 'chart.dim' or 'chart.dim.grid' for fit.ps.lps().",
             call. = FALSE)
    }
    foldid <- .klp.prepare.foldid(nrow(X), foldid, cv.folds, cv.seed)
    local.spec <- .ps.lps.resolve.local.candidate.spec(
        X = X,
        support.size = support.size,
        support.grid = support.grid,
        degree = degree,
        degree.grid = degree.grid,
        kernel = kernel,
        kernel.grid = kernel.grid,
        chart.dim = chart.dim,
        chart.dim.grid = chart.dim.grid,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
        selection.strategy = selection.strategy,
        chart.dim.max = chart.dim.max,
        design.margin = design.margin,
        n = nrow(X)
    )
    local.grid <- local.spec$local.grid
    if (nrow(local.grid) > 1L || !is.null(chart.dim.grid)) {
        return(.ps.lps.fit.local.grid(
            X = X,
            y = y,
            foldid = foldid,
            local.grid = local.grid,
            chart.dim = chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric,
            lambda.sync.grid = lambda.sync.grid,
            lambda.sync.search = lambda.sync.search,
            lambda.sync.selection = lambda.sync.selection,
            local.candidate.search = local.candidate.search,
            local.candidate.search.control = local.candidate.search.control,
            lambda.sync.search.control = lambda.sync.search.control,
            lambda.diagnostics = lambda.diagnostics,
            lambda.ridge = lambda.ridge,
            design.basis = design.basis,
            design.drop.tol = design.drop.tol,
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max,
            sync.neighbor.size = sync.neighbor.size,
            overlap.weight = overlap.weight,
            cv.folds = cv.folds,
            cv.seed = cv.seed,
            chart.activation = chart.activation,
            chart.activation.response = chart.activation.response,
            chart.activation.control = chart.activation.control,
            selection.strategy = selection.strategy,
            coupled.kd.selection = local.spec$telemetry,
            coupled.kd.candidate.plan = local.spec$coupled.plan
        ))
    }
    support.size <- local.grid$support.size[[1L]]
    degree <- local.grid$degree[[1L]]
    kernel <- local.grid$kernel[[1L]]
    sync.neighbor.size <- .ps.lps.resolve.sync.neighbor.size(
        sync.neighbor.size,
        support.size
    )
    timing.start <- proc.time()
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])
    if (is.null(ps.lps.geometry.cache)) {
        chart.dim.info <- .ps.lps.resolve.chart.dim(
            X = X,
            support.size = support.size,
            degree = degree,
            chart.dim = chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
    } else {
        chart.dim.info <- ps.lps.geometry.cache$chart.dim.info
    }
    chart.dim.by.anchor <- chart.dim.info$chart.dim.by.anchor
    if (is.null(ps.lps.geometry.cache)) {
        t.frames <- proc.time()
        frames <- .ps.lps.prepare.frames(
            X = X,
            y = y,
            support.size = support.size,
            degree = degree,
            kernel = kernel,
            chart.dim.by.anchor = chart.dim.by.anchor,
            design.basis = design.basis,
            design.drop.tol = design.drop.tol,
            local.pca.supports = ps.lps.local.pca.supports,
            chart.activation.info = chart.activation.info
        )
        phase.frames.sec <- elapsed(t.frames)
        t.sync.rows <- proc.time()
        sync.rows <- .ps.lps.prepare.sync.rows(
            frames = frames,
            sync.neighbor.size = sync.neighbor.size,
            overlap.weight = overlap.weight
        )
        phase.sync.rows.sec <- elapsed(t.sync.rows)
        geometry.cache.used <- FALSE
    } else {
        frames <- ps.lps.geometry.cache$frames
        sync.rows <- ps.lps.geometry.cache$sync.rows
        phase.frames.sec <- 0
        phase.sync.rows.sec <- 0
        geometry.cache.used <- TRUE
    }
    folds <- sort(unique(foldid))
    has.positive.sync <- any(lambda.sync.grid > 0)
    needs.fold.component.cache <-
        has.positive.sync && identical(lambda.sync.selection, "cv")
    fold.component.caches <- vector("list", length(folds))
    names(fold.component.caches) <- as.character(folds)
    full.component.cache <- NULL
    phase.system.cache.sec <- 0
    phase.fold.component.cache.sec <- 0
    phase.full.component.cache.sec <- 0
    if (has.positive.sync) {
        t.system <- proc.time()
        system.cache <- .ps.lps.prepare.system.cache(frames, sync.rows)
        phase.system.cache.sec <- elapsed(t.system)
        if (needs.fold.component.cache) {
            t.fold.cache <- proc.time()
            for (fold in folds) {
                response.weights <- as.numeric(foldid != fold)
                fold.component.caches[[as.character(fold)]] <-
                    .ps.lps.prepare.component.cache(
                        cache = system.cache,
                        y = y,
                        response.weights = response.weights
                    )
            }
            phase.fold.component.cache.sec <- elapsed(t.fold.cache)
        }
        t.full.cache <- proc.time()
        full.component.cache <- .ps.lps.prepare.component.cache(
            cache = system.cache,
            y = y,
            response.weights = rep(1, length(y))
        )
        phase.full.component.cache.sec <- elapsed(t.full.cache)
    }

    evaluate.lambda <- function(lambda) {
        t.eval <- proc.time()
        pred <- rep(NA_real_, length(y))
        fold.solve.elapsed <- 0
        for (fold in folds) {
            response.weights <- as.numeric(foldid != fold)
            t.fold.solve <- proc.time()
            fit.fold <- if (lambda > 0) {
                .ps.lps.solve.component.cached(
                    component.cache = fold.component.caches[[as.character(fold)]],
                    lambda.sync = lambda,
                    lambda.ridge = lambda.ridge,
                    ridge.multiplier.grid = ridge.multiplier.grid,
                    ridge.condition.max = ridge.condition.max
                )
            } else {
                .ps.lps.solve(
                    frames = frames,
                    y = y,
                    response.weights = response.weights,
                    lambda.sync = lambda,
                    lambda.ridge = lambda.ridge,
                    ridge.multiplier.grid = ridge.multiplier.grid,
                    ridge.condition.max = ridge.condition.max,
                    sync.rows = sync.rows
                )
            }
            fold.solve.elapsed <- fold.solve.elapsed + elapsed(t.fold.solve)
            pred[foldid == fold] <- fit.fold$fitted.values[foldid == fold]
        }
        cv.rmse <- .klp.rmse(pred, y)
        fit.diag <- NULL
        diag.elapsed <- 0
        if (identical(lambda.diagnostics, "all")) {
            t.diag <- proc.time()
            fit.diag <- if (lambda > 0) {
                .ps.lps.solve.component.cached(
                    component.cache = full.component.cache,
                    lambda.sync = lambda,
                    lambda.ridge = lambda.ridge,
                    ridge.multiplier.grid = ridge.multiplier.grid,
                    ridge.condition.max = ridge.condition.max,
                    coefficients.only = TRUE
                )
            } else {
                .ps.lps.solve(
                    frames = frames,
                    y = y,
                    response.weights = rep(1, length(y)),
                    lambda.sync = lambda,
                    lambda.ridge = lambda.ridge,
                    ridge.multiplier.grid = ridge.multiplier.grid,
                    ridge.condition.max = ridge.condition.max,
                    sync.rows = sync.rows,
                    coefficients.only = TRUE
                )
            }
            diag.elapsed <- elapsed(t.diag)
        }
        data.frame(
            lambda.sync = lambda,
            lambda.ridge = lambda.ridge,
            cv.rmse.observed = cv.rmse,
            total.local.gcv.ps = fit.diag$total.local.gcv.ps %||% NA_real_,
            sync.energy = fit.diag$sync.energy %||% NA_real_,
            mean.sync.squared.disagreement =
                fit.diag$mean.sync.squared.disagreement %||% NA_real_,
            ridge.median = fit.diag$ridge.median %||% NA_real_,
            ridge.max = fit.diag$ridge.max %||% NA_real_,
            ridge.multiplier.selected =
                fit.diag$ridge.multiplier.selected %||% NA_real_,
            ridge.condition = fit.diag$ridge.condition %||% NA_real_,
            ridge.status = fit.diag$ridge.status %||% NA_character_,
            evaluation.elapsed.sec = elapsed(t.eval),
            fold.solve.elapsed.sec = fold.solve.elapsed,
            diagnostic.elapsed.sec = diag.elapsed,
            stringsAsFactors = FALSE
        )
    }

    t.lambda.search <- proc.time()
    if (identical(lambda.sync.selection, "fixed")) {
        selected.lambda <- lambda.sync.grid[[1L]]
        cv.table <- data.frame(
            lambda.sync = selected.lambda,
            lambda.ridge = lambda.ridge,
            cv.rmse.observed = NA_real_,
            total.local.gcv.ps = NA_real_,
            sync.energy = NA_real_,
            mean.sync.squared.disagreement = NA_real_,
            ridge.median = NA_real_,
            ridge.max = NA_real_,
            ridge.multiplier.selected = NA_real_,
            ridge.condition = NA_real_,
            ridge.status = NA_character_,
            evaluation.elapsed.sec = 0,
            fold.solve.elapsed.sec = 0,
            diagnostic.elapsed.sec = 0,
            stringsAsFactors = FALSE
        )
        selected <- cv.table
        search.telemetry <- data.frame(
            stage = "fixed",
            lambda.sync = selected.lambda,
            boundary = "none",
            expansion = 0L,
            selected.after.stage = TRUE,
            stringsAsFactors = FALSE
        )
    } else if (identical(lambda.sync.search, "guarded")) {
        search.out <- .ps.lps.search.lambda.sync(
            evaluate = evaluate.lambda,
            lambda.grid = lambda.sync.grid,
            control = lambda.sync.search.control
        )
        cv.table <- search.out$evaluated
        selected <- search.out$selected
        search.telemetry <- search.out$telemetry
    } else {
        cv.table <- do.call(rbind, lapply(lambda.sync.grid, evaluate.lambda))
        cv.table <- cv.table[order(cv.table$lambda.sync), , drop = FALSE]
        selected <- .ps.lps.select.lambda.table(cv.table, rel.tol = 0)
        search.telemetry <- .ps.lps.grid.search.telemetry(lambda.sync.grid)
    }
    phase.lambda.search.sec <- elapsed(t.lambda.search)

    t.final <- proc.time()
    final <- if (selected$lambda.sync[[1L]] > 0) {
        .ps.lps.solve.component.cached(
            component.cache = full.component.cache,
            lambda.sync = selected$lambda.sync[[1L]],
            lambda.ridge = lambda.ridge,
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max
        )
    } else {
        .ps.lps.solve(
            frames = frames,
            y = y,
            response.weights = rep(1, length(y)),
            lambda.sync = selected$lambda.sync[[1L]],
            lambda.ridge = lambda.ridge,
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max,
            sync.rows = sync.rows
        )
    }
    phase.final.solve.sec <- elapsed(t.final)
    selected.lambda <- selected$lambda.sync[[1L]]
    cv.table <- .ps.lps.fill.selected.lambda.diagnostics(
        cv.table = cv.table,
        selected.lambda = selected.lambda,
        final = final
    )
    selected <- cv.table[
        which(cv.table$lambda.sync == selected.lambda)[[1L]],
        ,
        drop = FALSE
    ]
    timing <- list(
        total.elapsed.sec = elapsed(timing.start),
        phase_frames_sec = phase.frames.sec,
        phase_sync_rows_sec = phase.sync.rows.sec,
        phase_system_cache_sec = phase.system.cache.sec,
        phase_fold_component_cache_sec = phase.fold.component.cache.sec,
        phase_full_component_cache_sec = phase.full.component.cache.sec,
        phase_lambda_search_sec = phase.lambda.search.sec,
        phase_final_solve_sec = phase.final.solve.sec,
        evaluated_lambda_count = nrow(cv.table),
        unique_lambda_count = length(unique(cv.table$lambda.sync)),
        boundary_expansion_count = sum(grepl("^boundary_expand_",
                                             search.telemetry$stage))
    )
    out <- c(
        list(
            method.id = "ps_lps",
            method.label = "PS-LPS",
            X = X,
            y = y,
            support.size = support.size,
            support.grid = support.size,
            degree = degree,
            degree.grid = degree,
            kernel = kernel,
            kernel.grid = kernel,
            requested.chart.dim = chart.dim,
            chart.dim = if (length(unique(chart.dim.by.anchor)) == 1L) {
                unique(chart.dim.by.anchor)
            } else {
                stats::median(chart.dim.by.anchor)
            },
            chart.dim.mode = chart.dim.info$chart.dim.mode,
            auto.chart.dim = chart.dim.info$auto.chart.dim,
            diagnostics = list(
                chart.dim = .local.chart.dimension.telemetry(
                    chart.dim.info = chart.dim.info,
                    chart.dim.by.anchor = chart.dim.by.anchor,
                    source.path = "fit.ps.lps.frame.chart_dim_resolution"
                )
            ),
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric,
            lambda.sync.grid = lambda.sync.grid,
            lambda.sync.search = lambda.sync.search,
            lambda.sync.selection = lambda.sync.selection,
            selection.strategy = selection.strategy,
            lambda.diagnostics = lambda.diagnostics,
            local.candidate.search = local.candidate.search,
            local.candidate.search.control = local.candidate.search.control,
            lambda.sync.search.telemetry = search.telemetry,
            lambda.ridge = lambda.ridge,
            design.basis = design.basis,
            design.drop.tol = design.drop.tol,
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max,
            selected = selected,
            cv.table = cv.table,
            foldid = foldid,
            sync.neighbor.size = sync.neighbor.size,
            overlap.weight = overlap.weight,
            cache.backend = if (has.positive.sync) {
                "component"
            } else {
                "independent"
            },
            geometry.cache.used = geometry.cache.used,
            chart.activation = chart.activation,
            chart.activation.control = chart.activation.info$control,
            chart.activation.diagnostics =
                attr(frames, "chart.activation.diagnostics"),
            ps.lps.timing = timing,
            frame.design.summary = .ps.lps.frame.design.summary(frames)
        ),
        final
    )
    out$diagnostics$chart.activation <- .klp.chart.activation.summary(
        attr(frames, "chart.activation.diagnostics"),
        chart.activation.info
    )
    out$diagnostics$coupled.kd.selection <- local.spec$telemetry
    out$coupled.kd.candidate.plan <- local.spec$coupled.plan
    out$selection.strategy <- selection.strategy
    class(out) <- c("ps_lps", "list")
    out
}

.ps.lps.resolve.local.grid <- function(support.size = NULL,
                                       support.grid = NULL,
                                       degree = 2L,
                                       degree.grid = NULL,
                                       kernel = "gaussian",
                                       kernel.grid = NULL,
                                       n) {
    if (!is.null(support.size) && !is.null(support.grid)) {
        stop("Use either 'support.size' or 'support.grid', not both.",
             call. = FALSE)
    }
    if (is.null(support.grid)) {
        support.grid <- if (is.null(support.size)) {
            c(10L, 15L, 20L)
        } else {
            support.size
        }
    }
    if (is.null(degree.grid)) degree.grid <- degree
    if (is.null(kernel.grid)) kernel.grid <- kernel
    support.grid <- .klp.clean.support.grid(support.grid, n)
    degree.grid <- .klp.clean.degree.grid(degree.grid)
    kernel.grid <- .klp.clean.kernel.grid(kernel.grid)
    cand <- expand.grid(
        support.size = support.grid,
        degree = degree.grid,
        kernel = kernel.grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    cand[order(cand$support.size, cand$degree, cand$kernel), , drop = FALSE]
}

.ps.lps.resolve.local.candidate.spec <- function(X,
                                                 support.size = NULL,
                                                 support.grid = NULL,
                                                 degree = 2L,
                                                 degree.grid = NULL,
                                                 kernel = "gaussian",
                                                 kernel.grid = NULL,
                                                 chart.dim = NULL,
                                                 chart.dim.grid = NULL,
                                                 auto.chart.support.metric,
                                                 auto.chart.selection.metric,
                                                 selection.strategy = "grid",
                                                 chart.dim.max = NULL,
                                                 design.margin = 2L,
                                                 n) {
    if (is.null(chart.dim.grid)) {
        return(list(
            local.grid = .ps.lps.resolve.local.grid(
                support.size = support.size,
                support.grid = support.grid,
                degree = degree,
                degree.grid = degree.grid,
                kernel = kernel,
                kernel.grid = kernel.grid,
                n = n
            ),
            telemetry = NULL,
            coupled.plan = NULL
        ))
    }
    support.grid <- if (is.null(support.grid)) {
        if (is.null(support.size)) c(10L, 15L, 20L) else support.size
    } else {
        if (!is.null(support.size)) {
            stop("Use either 'support.size' or 'support.grid', not both.",
                 call. = FALSE)
        }
        support.grid
    }
    if (is.null(degree.grid)) degree.grid <- degree
    if (is.null(kernel.grid)) kernel.grid <- kernel
    candidate.spec <- .coupled.kd.lps.candidate.spec(
        X = X,
        support.grid = support.grid,
        degree.grid = degree.grid,
        kernel.grid = kernel.grid,
        bandwidth.multiplier.grid = 1,
        chart.dim = chart.dim,
        chart.dim.grid = chart.dim.grid,
        coordinate.method = "local.pca",
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
        selection.strategy = selection.strategy,
        chart.dim.max = chart.dim.max,
        design.margin = design.margin,
        reuse.type = "weighted"
    )
    local.grid <- candidate.spec$candidates
    keep <- intersect(
        c("candidate.id", "support.size", "degree", "kernel", "chart.dim",
          "chart.dim.rank", "chart.dim.source", "chart.dim.raw",
          "chart.dim.clipped", "chart.dim.seed.clipped", "chart.dim.max",
          "design.ncol", "design.margin", "reuse.key",
          "reuse.chart.dim.max"),
        names(local.grid)
    )
    local.grid <- local.grid[, keep, drop = FALSE]
    local.grid$local.candidate.id <- seq_len(nrow(local.grid))
    list(
        local.grid = local.grid,
        telemetry = candidate.spec$telemetry,
        coupled.plan = candidate.spec$coupled.plan
    )
}

.ps.lps.resolve.sync.neighbor.size <- function(sync.neighbor.size,
                                               support.size) {
    if (is.null(sync.neighbor.size)) {
        sync.neighbor.size <- min(8L, as.integer(support.size) - 1L)
    }
    sync.neighbor.size <- as.integer(sync.neighbor.size[[1L]])
    if (!is.finite(sync.neighbor.size) || sync.neighbor.size < 1L ||
        sync.neighbor.size >= support.size) {
        stop("'sync.neighbor.size' must be an integer in [1, support.size - 1].",
             call. = FALSE)
    }
    sync.neighbor.size
}

.ps.lps.local.search.control <- function(control) {
    defaults <- list(
        top.n = 8L,
        max.candidates = 12L,
        neighbor.radius = 1L,
        guard.support.quantiles = c(0, 0.5, 1)
    )
    if (is.null(control)) control <- list()
    if (!is.list(control)) {
        stop("'local.candidate.search.control' must be a list.",
             call. = FALSE)
    }
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) {
        stop("Unknown local candidate search control field(s): ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- utils::modifyList(defaults, control)
    out$top.n <- as.integer(out$top.n[[1L]])
    out$max.candidates <- as.integer(out$max.candidates[[1L]])
    out$neighbor.radius <- as.integer(out$neighbor.radius[[1L]])
    out$guard.support.quantiles <-
        as.numeric(out$guard.support.quantiles)
    if (!is.finite(out$top.n) || out$top.n < 1L) {
        stop("'top.n' must be a positive integer.", call. = FALSE)
    }
    if (!is.finite(out$max.candidates) || out$max.candidates < 1L) {
        stop("'max.candidates' must be a positive integer.", call. = FALSE)
    }
    if (!is.finite(out$neighbor.radius) || out$neighbor.radius < 0L) {
        stop("'neighbor.radius' must be a nonnegative integer.",
             call. = FALSE)
    }
    if (!length(out$guard.support.quantiles) ||
        any(!is.finite(out$guard.support.quantiles)) ||
        any(out$guard.support.quantiles < 0 | out$guard.support.quantiles > 1)) {
        stop("'guard.support.quantiles' must contain values in [0, 1].",
             call. = FALSE)
    }
    out
}

.ps.lps.local.grid.key <- function(support.size, degree, kernel) {
    paste(as.integer(support.size), as.integer(degree),
          as.character(kernel), sep = "\r")
}

.ps.lps.local.grid.chart.dim <- function(row, chart.dim) {
    if ("chart.dim" %in% names(row) && !is.na(row$chart.dim[[1L]])) {
        return(.local.chart.decode.chart.dim(row$chart.dim[[1L]]))
    }
    chart.dim
}

.ps.lps.screen.local.grid <- function(
    X, y, foldid, local.grid, chart.dim, auto.chart.support.metric,
    auto.chart.selection.metric, local.candidate.search,
    local.candidate.search.control, cv.folds, cv.seed,
    design.basis = "monomial",
    design.drop.tol = sqrt(.Machine$double.eps),
    ridge.multiplier.grid = 0,
    ridge.condition.max = Inf) {

    out <- local.grid
    out$local.candidate.id <- seq_len(nrow(out))
    out$screening.cv.rmse.observed <- NA_real_
    out$screening.rank <- NA_integer_
    out$screening.reason <- "full_search"
    if (identical(local.candidate.search, "full") || nrow(local.grid) <= 1L) {
        return(list(
            active.ids = seq_len(nrow(local.grid)),
            table = out,
            control = list(),
            lps.screen = NULL
        ))
    }

    has.candidate.chart.dim <- "chart.dim" %in% names(local.grid)
    ctl <- .ps.lps.local.search.control(local.candidate.search.control)
    if (identical(local.candidate.search, "subgrid")) {
        guard.supports <- sort(unique(local.grid$support.size))
        guard.supports <- as.integer(stats::quantile(
            guard.supports,
            probs = unique(ctl$guard.support.quantiles),
            type = 1,
            names = FALSE
        ))
        guard <- out[out$support.size %in% guard.supports, , drop = FALSE]
        guard$kernel.order <- match(guard$kernel, unique(local.grid$kernel))
        guard <- guard[order(guard$support.size, guard$degree,
                             guard$kernel.order), , drop = FALSE]
        active.ids <- head(guard$local.candidate.id,
                           min(ctl$max.candidates, nrow(guard)))
        out$screening.reason <- "screened_out"
        out$screening.reason[out$local.candidate.id %in% active.ids] <-
            "subgrid_guard"
        out$screening.rank[out$local.candidate.id %in% active.ids] <-
            seq_along(active.ids)
        return(list(
            active.ids = as.integer(active.ids),
            table = out,
            control = ctl,
            lps.screen = NULL
        ))
    }

    if (isTRUE(has.candidate.chart.dim)) {
        stop("'local.candidate.search = \"screened\"' is not yet supported ",
             "with candidate-specific 'chart.dim'. Use ",
             "'local.candidate.search = \"full\"' or \"subgrid\".",
             call. = FALSE)
    }
    if (!(identical(chart.dim, "auto") ||
          .klp.is.local.auto.chart.dim(chart.dim) ||
          (length(chart.dim) == 1L && is.numeric(chart.dim)))) {
        stop("'local.candidate.search = \"screened\"' supports scalar, ",
             "'auto', or 'local.auto' chart dimensions.", call. = FALSE)
    }
    lps.screen <- fit.lps(
        X = X,
        y = y,
        foldid = foldid,
        support.grid = sort(unique(local.grid$support.size)),
        degree.grid = sort(unique(local.grid$degree)),
        kernel.grid = unique(local.grid$kernel),
        cv.folds = cv.folds,
        cv.seed = cv.seed,
        coordinate.method = "local.pca",
        chart.dim = chart.dim,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
        design.basis = design.basis,
        design.drop.tol = design.drop.tol,
        ridge.multiplier.grid = ridge.multiplier.grid,
        ridge.condition.max = ridge.condition.max,
        backend = "auto"
    )
    screen.tab <- lps.screen$cv.table
    screen.tab$screening.key <- .ps.lps.local.grid.key(
        screen.tab$support.size,
        screen.tab$degree,
        screen.tab$kernel
    )
    out$screening.key <- .ps.lps.local.grid.key(out$support.size, out$degree,
                                                out$kernel)
    m <- match(out$screening.key, screen.tab$screening.key)
    out$screening.cv.rmse.observed <- screen.tab$cv.rmse.observed[m]
    finite <- is.finite(out$screening.cv.rmse.observed)
    ranked <- order(out$screening.cv.rmse.observed[finite],
                    out$support.size[finite],
                    out$degree[finite],
                    out$kernel[finite])
    finite.ids <- which(finite)
    ranked.ids <- finite.ids[ranked]
    out$screening.rank[ranked.ids] <- seq_along(ranked.ids)

    selected.key <- .ps.lps.local.grid.key(
        lps.screen$selected$support.size[[1L]],
        lps.screen$selected$degree[[1L]],
        lps.screen$selected$kernel[[1L]]
    )
    selected.id <- out$local.candidate.id[out$screening.key == selected.key]
    top.ids <- head(ranked.ids, min(ctl$top.n, length(ranked.ids)))
    seed.ids <- unique(c(selected.id, top.ids))

    neighbor.ids <- integer(0)
    if (ctl$neighbor.radius > 0L && length(seed.ids)) {
        for (id in seed.ids) {
            row <- out[id, , drop = FALSE]
            neighbor.ids <- c(neighbor.ids, out$local.candidate.id[
                out$degree == row$degree[[1L]] &
                    out$kernel == row$kernel[[1L]] &
                    abs(out$support.size - row$support.size[[1L]]) <=
                    ctl$neighbor.radius
            ])
        }
    }

    guard.supports <- sort(unique(local.grid$support.size))
    guard.supports <- as.integer(stats::quantile(
        guard.supports,
        probs = unique(ctl$guard.support.quantiles),
        type = 1,
        names = FALSE
    ))
    guard.ids <- out$local.candidate.id[
        out$support.size %in% guard.supports
    ]

    reason <- rep("not_selected", nrow(out))
    reason[selected.id] <- "lps_selected"
    reason[top.ids] <- ifelse(reason[top.ids] == "lps_selected",
                              "lps_selected_top", "top_lps_cv")
    reason[neighbor.ids] <- ifelse(
        reason[neighbor.ids] %in% c("lps_selected", "lps_selected_top",
                                    "top_lps_cv"),
        paste(reason[neighbor.ids], "neighbor", sep = "+"),
        "neighbor"
    )
    reason[guard.ids] <- ifelse(
        reason[guard.ids] != "not_selected",
        paste(reason[guard.ids], "guard", sep = "+"),
        "guard"
    )

    active.ids <- unique(c(selected.id, top.ids, neighbor.ids, guard.ids))
    priority <- data.frame(
        id = active.ids,
        selected = ifelse(active.ids %in% selected.id, 0L, 1L),
        top = ifelse(active.ids %in% top.ids, 0L, 1L),
        neighbor = ifelse(active.ids %in% neighbor.ids, 0L, 1L),
        guard = ifelse(active.ids %in% guard.ids, 0L, 1L),
        rank = out$screening.rank[active.ids],
        stringsAsFactors = FALSE
    )
    priority$rank[!is.finite(priority$rank)] <- Inf
    priority <- priority[order(priority$selected, priority$top,
                               priority$neighbor, priority$guard,
                               priority$rank, priority$id), , drop = FALSE]
    active.ids <- head(priority$id, min(ctl$max.candidates,
                                       nrow(priority)))
    out$screening.reason <- reason
    out$screening.reason[!out$local.candidate.id %in% active.ids] <-
        "screened_out"
    out$screening.key <- NULL
    list(
        active.ids = as.integer(active.ids),
        table = out,
        control = ctl,
        lps.screen = lps.screen
    )
}

.ps.lps.fit.local.grid <- function(
    X, y, foldid, local.grid, chart.dim, auto.chart.support.metric,
    auto.chart.selection.metric, lambda.sync.grid, lambda.sync.search,
    lambda.sync.selection,
    local.candidate.search, local.candidate.search.control,
    lambda.sync.search.control, lambda.diagnostics, lambda.ridge, design.basis,
    design.drop.tol, ridge.multiplier.grid, ridge.condition.max,
    sync.neighbor.size, overlap.weight, cv.folds, cv.seed,
    chart.activation = "none",
    chart.activation.response = NULL,
    chart.activation.control = list(),
    selection.strategy = "grid",
    coupled.kd.selection = NULL,
    coupled.kd.candidate.plan = NULL) {

    local.grid.timing.start <- proc.time()
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])
    t.screen <- proc.time()
    screen <- .ps.lps.screen.local.grid(
        X = X,
        y = y,
        foldid = foldid,
        local.grid = local.grid,
        chart.dim = chart.dim,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric,
        local.candidate.search = local.candidate.search,
        local.candidate.search.control = local.candidate.search.control,
        cv.folds = cv.folds,
        cv.seed = cv.seed,
        design.basis = design.basis,
        design.drop.tol = design.drop.tol,
        ridge.multiplier.grid = ridge.multiplier.grid,
        ridge.condition.max = ridge.condition.max
    )
    phase.screening.sec <- elapsed(t.screen)
    active.ids <- screen$active.ids
    fits <- vector("list", nrow(local.grid))
    candidate.rows <- vector("list", nrow(local.grid))
    lambda.rows <- vector("list", nrow(local.grid))
    weighted.reuse.plan <- if ("chart.dim" %in% names(local.grid)) {
        .coupled.kd.reuse.plan(local.grid, reuse.type = "weighted")
    } else {
        data.frame()
    }
    weighted.support.cache <- new.env(parent = emptyenv())
    for (ii in seq_len(nrow(local.grid))) {
        row <- local.grid[ii, , drop = FALSE]
        row.chart.dim <- .ps.lps.local.grid.chart.dim(row, chart.dim)
        row.chart.dim.label <- if ("chart.dim" %in% names(row)) {
            as.character(row$chart.dim[[1L]])
        } else {
            NA_character_
        }
        screen.row <- screen$table[screen$table$local.candidate.id == ii, ,
                                   drop = FALSE]
        if (!ii %in% active.ids) {
            candidate.rows[[ii]] <- data.frame(
                local.candidate.id = ii,
                support.size = row$support.size[[1L]],
                degree = row$degree[[1L]],
                kernel = row$kernel[[1L]],
                chart.dim = row.chart.dim.label,
                chart.dim.mode = NA_character_,
                selected.lambda.sync = NA_real_,
                selected.lambda.ridge = NA_real_,
                selected.cv.rmse.observed = NA_real_,
                selected.total.local.gcv.ps = NA_real_,
                selected.sync.energy = NA_real_,
                selected.mean.sync.squared.disagreement = NA_real_,
                evaluated.lambda.count = 0L,
                unique.lambda.count = 0L,
                boundary.expansion.count = 0L,
                local.candidate.elapsed.sec = NA_real_,
                lambda.search.elapsed.sec = NA_real_,
                frame.prep.elapsed.sec = NA_real_,
                system.cache.elapsed.sec = NA_real_,
                fold.component.cache.elapsed.sec = NA_real_,
                final.solve.elapsed.sec = NA_real_,
                local.candidate.status = "screened_out",
                screening.cv.rmse.observed =
                    screen.row$screening.cv.rmse.observed[[1L]],
                screening.rank = screen.row$screening.rank[[1L]],
                screening.reason = screen.row$screening.reason[[1L]],
                stringsAsFactors = FALSE
            )
            next
        }
        t.candidate <- proc.time()
        geometry.cache <- NULL
        if ("chart.dim" %in% names(row)) {
            candidate.sync.size <- .ps.lps.resolve.sync.neighbor.size(
                sync.neighbor.size,
                row$support.size[[1L]]
            )
            local.pca.supports <- .coupled.kd.shared.local.pca.supports(
                X = X,
                support.size = row$support.size[[1L]],
                chart.dim = row.chart.dim,
                kernel = row$kernel[[1L]],
                reuse.plan = weighted.reuse.plan,
                cache.env = weighted.support.cache
            )
            geometry.cache <- .ps.lps.prepare.geometry.cache(
                X = X,
                support.size = row$support.size[[1L]],
                degree = row$degree[[1L]],
                kernel = row$kernel[[1L]],
                chart.dim = row.chart.dim,
                auto.chart.support.metric = auto.chart.support.metric,
                auto.chart.selection.metric = auto.chart.selection.metric,
                sync.neighbor.size = candidate.sync.size,
                overlap.weight = overlap.weight,
                design.basis = design.basis,
                design.drop.tol = design.drop.tol,
                local.pca.supports = local.pca.supports,
                chart.activation = chart.activation,
                chart.activation.response = chart.activation.response,
                chart.activation.control = chart.activation.control
            )
        }
        fit <- fit.ps.lps(
            X = X,
            y = y,
            foldid = foldid,
            support.size = row$support.size[[1L]],
            degree = row$degree[[1L]],
            kernel = row$kernel[[1L]],
            chart.dim = row.chart.dim,
            support.grid = NULL,
            degree.grid = NULL,
            kernel.grid = NULL,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric,
            lambda.sync.grid = lambda.sync.grid,
            lambda.sync.search = lambda.sync.search,
            lambda.sync.selection = lambda.sync.selection,
            local.candidate.search = "full",
            local.candidate.search.control = list(),
            lambda.sync.search.control = lambda.sync.search.control,
            lambda.diagnostics = lambda.diagnostics,
            lambda.ridge = lambda.ridge,
            design.basis = design.basis,
            design.drop.tol = design.drop.tol,
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max,
            sync.neighbor.size = sync.neighbor.size,
            overlap.weight = overlap.weight,
            chart.activation = chart.activation,
            chart.activation.response = chart.activation.response,
            chart.activation.control = chart.activation.control,
            ps.lps.geometry.cache = geometry.cache,
            cv.folds = cv.folds,
            cv.seed = cv.seed
        )
        candidate.elapsed <- elapsed(t.candidate)
        fits[[ii]] <- fit
        selected <- fit$selected[1L, , drop = FALSE]
        timing <- fit$ps.lps.timing %||% list()
        candidate.rows[[ii]] <- data.frame(
            local.candidate.id = ii,
            support.size = fit$support.size,
            degree = fit$degree,
            kernel = fit$kernel,
            chart.dim = fit$chart.dim,
            chart.dim.mode = fit$chart.dim.mode,
            selected.lambda.sync = selected$lambda.sync[[1L]],
            selected.lambda.ridge = selected$lambda.ridge[[1L]],
            selected.cv.rmse.observed = selected$cv.rmse.observed[[1L]],
            selected.total.local.gcv.ps = selected$total.local.gcv.ps[[1L]],
            selected.sync.energy = selected$sync.energy[[1L]],
            selected.mean.sync.squared.disagreement =
                selected$mean.sync.squared.disagreement[[1L]],
            evaluated.lambda.count = nrow(fit$cv.table),
            unique.lambda.count =
                length(unique(fit$cv.table$lambda.sync)),
            boundary.expansion.count =
                sum(grepl("^boundary_expand_",
                          fit$lambda.sync.search.telemetry$stage)),
            local.candidate.elapsed.sec = candidate.elapsed,
            lambda.search.elapsed.sec =
                as.numeric(timing$phase_lambda_search_sec %||% NA_real_),
            frame.prep.elapsed.sec =
                as.numeric(timing$phase_frames_sec %||% NA_real_),
            system.cache.elapsed.sec =
                as.numeric(timing$phase_system_cache_sec %||% NA_real_),
            fold.component.cache.elapsed.sec =
                as.numeric(timing$phase_fold_component_cache_sec %||% NA_real_),
            final.solve.elapsed.sec =
                as.numeric(timing$phase_final_solve_sec %||% NA_real_),
            local.candidate.status = "evaluated",
            screening.cv.rmse.observed =
                screen.row$screening.cv.rmse.observed[[1L]],
            screening.rank = screen.row$screening.rank[[1L]],
            screening.reason = screen.row$screening.reason[[1L]],
            stringsAsFactors = FALSE
        )
        tab <- fit$cv.table
        tab$local.candidate.id <- ii
        tab$support.size <- fit$support.size
        tab$degree <- fit$degree
        tab$kernel <- fit$kernel
        tab$chart.dim <- fit$chart.dim
        tab$chart.dim.mode <- fit$chart.dim.mode
        lambda.rows[[ii]] <- tab[, c(
            "local.candidate.id", "support.size", "degree", "kernel",
            "chart.dim", "chart.dim.mode",
            setdiff(names(tab), c("local.candidate.id", "support.size",
                                  "degree", "kernel", "chart.dim",
                                  "chart.dim.mode"))
        ), drop = FALSE]
    }
    candidate.table <- do.call(rbind, candidate.rows)
    ok <- is.finite(candidate.table$selected.cv.rmse.observed)
    if (!any(ok)) {
        stop("No PS-LPS local candidate produced a finite CV score.",
             call. = FALSE)
    }
    best.order <- order(
        candidate.table$selected.cv.rmse.observed[ok],
        candidate.table$selected.lambda.sync[ok],
        candidate.table$support.size[ok],
        candidate.table$degree[ok],
        candidate.table$kernel[ok]
    )
    ok.idx <- which(ok)
    best.idx <- ok.idx[best.order[[1L]]]
    best.fit <- fits[[best.idx]]
    best.fit$support.grid <- sort(unique(local.grid$support.size))
    best.fit$degree.grid <- sort(unique(local.grid$degree))
    best.fit$kernel.grid <- sort(unique(local.grid$kernel))
    best.fit$local.candidate.table <- candidate.table
    best.fit$lambda.cv.table <- do.call(rbind, lambda.rows)
    best.fit$ps.lps.local.grid.timing <- list(
        total.elapsed.sec = elapsed(local.grid.timing.start),
        phase_screening_sec = phase.screening.sec,
        phase_candidate_loop_sec = sum(
            candidate.table$local.candidate.elapsed.sec[
                is.finite(candidate.table$local.candidate.elapsed.sec)
            ],
            na.rm = TRUE
        ),
        planned_local_candidate_count = nrow(local.grid),
        evaluated_local_candidate_count = sum(
            candidate.table$local.candidate.status == "evaluated",
            na.rm = TRUE
        ),
        total_lambda_rows = nrow(best.fit$lambda.cv.table),
        total_unique_lambda_rows = length(unique(
            best.fit$lambda.cv.table$lambda.sync
        ))
    )
    best.fit$local.candidate.search <- local.candidate.search
    best.fit$selection.strategy <- selection.strategy
    best.fit$local.candidate.search.control <- screen$control
    best.fit$local.candidate.screen.table <- screen$table
    best.fit$local.candidate.screen.lps.selected <-
        if (is.null(screen$lps.screen)) NULL else screen$lps.screen$selected
    best.fit$selected.local.candidate <-
        candidate.table[best.idx, , drop = FALSE]
    best.fit$selection.contract <-
        if (!is.null(coupled.kd.selection) &&
            isTRUE(coupled.kd.selection$coupled.chart.dim.search)) {
            "sparse_kd_coupled_support_chart_dim_with_lambda_sync"
        } else if (identical(local.candidate.search, "screened")) {
            "screened_lps_cv_then_materialized_fold_cv_with_lambda_sync"
        } else if (identical(local.candidate.search, "subgrid")) {
            "subgrid_then_materialized_fold_cv_with_lambda_sync"
        } else {
            "materialized_fold_cv_over_support_kernel_degree_with_lambda_sync"
        }
    best.fit$selected <- cbind(
        candidate.table[best.idx, c("local.candidate.id", "support.size",
                                    "degree", "kernel", "chart.dim",
                                    "chart.dim.mode"), drop = FALSE],
        best.fit$selected
    )
    best.fit$diagnostics$coupled.kd.selection <- coupled.kd.selection
    best.fit$coupled.kd.candidate.plan <- coupled.kd.candidate.plan
    best.fit
}

.ps.lps.select.lambda.table <- function(x, rel.tol = 0) {
    if (!nrow(x)) stop("'x' must contain at least one lambda row.", call. = FALSE)
    rel.tol <- as.numeric(rel.tol[[1L]])
    if (!is.finite(rel.tol) || rel.tol < 0) {
        stop("'rel.tol' must be a finite nonnegative scalar.", call. = FALSE)
    }
    best <- min(x$cv.rmse.observed, na.rm = TRUE)
    if (!is.finite(best)) {
        stop("No finite CV RMSE values are available.", call. = FALSE)
    }
    cutoff <- best * (1 + rel.tol)
    tied <- x[is.finite(x$cv.rmse.observed) &
                  x$cv.rmse.observed <= cutoff, , drop = FALSE]
    tied[order(tied$lambda.sync), , drop = FALSE][1L, , drop = FALSE]
}

.ps.lps.fill.selected.lambda.diagnostics <- function(cv.table,
                                                     selected.lambda,
                                                     final) {
    if (!nrow(cv.table)) {
        return(cv.table)
    }
    selected.lambda <- as.numeric(selected.lambda[[1L]])
    idx <- which(cv.table$lambda.sync == selected.lambda)
    if (!length(idx)) {
        return(cv.table)
    }
    idx <- idx[[1L]]
    replacements <- list(
        total.local.gcv.ps = final$total.local.gcv.ps %||% NA_real_,
        sync.energy = final$sync.energy %||% NA_real_,
        mean.sync.squared.disagreement =
            final$mean.sync.squared.disagreement %||% NA_real_,
        ridge.median = final$ridge.median %||% NA_real_,
        ridge.max = final$ridge.max %||% NA_real_,
        ridge.multiplier.selected =
            final$ridge.multiplier.selected %||% NA_real_,
        ridge.condition = final$ridge.condition %||% NA_real_,
        ridge.status = final$ridge.status %||% NA_character_
    )
    for (nm in names(replacements)) {
        if (!nm %in% names(cv.table)) {
            cv.table[[nm]] <- if (is.character(replacements[[nm]])) {
                NA_character_
            } else {
                NA_real_
            }
        }
        cv.table[[nm]][idx] <- replacements[[nm]]
    }
    cv.table
}

.ps.lps.raw.best.lambda.table <- function(x) {
    if (!nrow(x)) stop("'x' must contain at least one lambda row.", call. = FALSE)
    x[order(x$cv.rmse.observed, x$lambda.sync), , drop = FALSE][1L, ,
                                                                drop = FALSE]
}

.ps.lps.search.control <- function(control) {
    defaults <- list(
        coarse.size = 5L,
        refine.radius = 2L,
        rel.tol = 0.002,
        boundary.guard.rel.tol = 0.01,
        boundary.expand = TRUE,
        boundary.factor = 3,
        max.boundary.expansions = 2L,
        max.candidates = 25L
    )
    if (is.null(control)) control <- list()
    if (!is.list(control)) {
        stop("'lambda.sync.search.control' must be a list.", call. = FALSE)
    }
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown)) {
        stop("Unknown lambda search control field(s): ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- utils::modifyList(defaults, control)
    out$coarse.size <- as.integer(out$coarse.size[[1L]])
    out$refine.radius <- as.integer(out$refine.radius[[1L]])
    out$max.boundary.expansions <-
        as.integer(out$max.boundary.expansions[[1L]])
    out$max.candidates <- as.integer(out$max.candidates[[1L]])
    out$rel.tol <- as.numeric(out$rel.tol[[1L]])
    out$boundary.guard.rel.tol <-
        as.numeric(out$boundary.guard.rel.tol[[1L]])
    out$boundary.factor <- as.numeric(out$boundary.factor[[1L]])
    out$boundary.expand <- isTRUE(out$boundary.expand)
    if (!is.finite(out$coarse.size) || out$coarse.size < 2L) {
        stop("'coarse.size' must be an integer at least 2.", call. = FALSE)
    }
    if (!is.finite(out$refine.radius) || out$refine.radius < 1L) {
        stop("'refine.radius' must be a positive integer.", call. = FALSE)
    }
    if (!is.finite(out$max.boundary.expansions) ||
        out$max.boundary.expansions < 0L) {
        stop("'max.boundary.expansions' must be a nonnegative integer.",
             call. = FALSE)
    }
    if (!is.finite(out$max.candidates) || out$max.candidates < 1L) {
        stop("'max.candidates' must be a positive integer.", call. = FALSE)
    }
    if (!is.finite(out$rel.tol) || out$rel.tol < 0) {
        stop("'rel.tol' must be a finite nonnegative scalar.", call. = FALSE)
    }
    if (!is.finite(out$boundary.guard.rel.tol) ||
        out$boundary.guard.rel.tol < 0) {
        stop("'boundary.guard.rel.tol' must be a finite nonnegative scalar.",
             call. = FALSE)
    }
    if (!is.finite(out$boundary.factor) || out$boundary.factor <= 1) {
        stop("'boundary.factor' must be a finite scalar larger than 1.",
             call. = FALSE)
    }
    out
}

.ps.lps.grid.search.telemetry <- function(lambda.grid) {
    lambda.grid <- sort(unique(as.numeric(lambda.grid)))
    data.frame(
        stage = "full_grid",
        lambda.sync = lambda.grid,
        boundary = "none",
        expansion = 0L,
        selected.after.stage = NA_real_,
        stringsAsFactors = FALSE
    )
}

.ps.lps.search.lambda.sync <- function(evaluate, lambda.grid, control = list()) {
    if (!is.function(evaluate)) {
        stop("'evaluate' must be a function.", call. = FALSE)
    }
    ctl <- .ps.lps.search.control(control)
    lambda.grid <- sort(unique(as.numeric(lambda.grid)))
    if (!length(lambda.grid) || any(!is.finite(lambda.grid)) ||
        any(lambda.grid < 0)) {
        stop("'lambda.grid' must contain finite nonnegative values.",
             call. = FALSE)
    }
    positive <- lambda.grid[lambda.grid > 0]
    if (!length(positive)) {
        row <- evaluate(0)
        return(list(
            evaluated = row,
            selected = row,
            telemetry = data.frame(
                stage = "zero_only",
                lambda.sync = 0,
                boundary = "zero",
                expansion = 0L,
                selected.after.stage = 0,
                stringsAsFactors = FALSE
            )
        ))
    }

    eval.env <- new.env(parent = emptyenv())
    evaluated.count <- function() length(ls(eval.env, all.names = TRUE))
    eval.one <- function(lambda, stage, boundary = "none", expansion = 0L) {
        lambda <- as.numeric(lambda[[1L]])
        key <- sprintf("%.17g", lambda)
        if (!exists(key, envir = eval.env, inherits = FALSE)) {
            assign(key, evaluate(lambda), envir = eval.env)
        }
        row <- get(key, envir = eval.env, inherits = FALSE)
        row$stage <- stage
        row$boundary <- boundary
        row$expansion <- as.integer(expansion)
        row
    }
    collect.evaluated <- function() {
        keys <- ls(eval.env, all.names = TRUE)
        if (!length(keys)) {
            return(data.frame())
        }
        out <- do.call(rbind, lapply(keys, function(key) {
            get(key, envir = eval.env, inherits = FALSE)
        }))
        out[order(out$lambda.sync), , drop = FALSE]
    }
    append.telemetry <- function(rows, selected, telemetry) {
        if (!nrow(rows)) return(telemetry)
        add <- data.frame(
            stage = rows$stage,
            lambda.sync = rows$lambda.sync,
            boundary = rows$boundary,
            expansion = rows$expansion,
            selected.after.stage = selected$lambda.sync[[1L]],
            stringsAsFactors = FALSE
        )
        rbind(telemetry, add)
    }

    if (0 %in% lambda.grid) {
        eval.set <- 0
    } else {
        eval.set <- numeric(0)
    }
    remaining <- max(0L, ctl$max.candidates - length(eval.set))
    coarse.n <- min(length(positive), ctl$coarse.size, remaining)
    if (coarse.n > 0L) {
        coarse.idx <- unique(round(seq(1, length(positive),
                                       length.out = coarse.n)))
        coarse <- positive[coarse.idx]
        eval.set <- sort(unique(c(eval.set, coarse)))
    }
    telemetry <- data.frame(
        stage = character(),
        lambda.sync = numeric(),
        boundary = character(),
        expansion = integer(),
        selected.after.stage = numeric(),
        stringsAsFactors = FALSE
    )
    stage.rows <- do.call(rbind, lapply(eval.set, eval.one,
                                        stage = "coarse"))
    selected <- .ps.lps.select.lambda.table(
        collect.evaluated(),
        rel.tol = ctl$rel.tol
    )
    telemetry <- append.telemetry(stage.rows, selected, telemetry)

    coarse.raw.best <- .ps.lps.raw.best.lambda.table(collect.evaluated())
    refine.anchor <- if (coarse.raw.best$lambda.sync[[1L]] > 0) {
        coarse.raw.best$lambda.sync[[1L]]
    } else {
        selected$lambda.sync[[1L]]
    }
    if (refine.anchor > 0) {
        best.idx <- which(positive == refine.anchor)
        if (length(best.idx)) {
            lo <- max(1L, best.idx - ctl$refine.radius)
            hi <- min(length(positive), best.idx + ctl$refine.radius)
            refine <- positive[lo:hi]
        } else {
            log.grid <- log10(positive)
            nearest <- which.min(abs(log.grid - log10(refine.anchor)))
            lo <- max(1L, nearest - ctl$refine.radius)
            hi <- min(length(positive), nearest + ctl$refine.radius)
            refine <- positive[lo:hi]
        }
        already <- collect.evaluated()$lambda.sync
        refine <- refine[!refine %in% already]
        remaining <- max(0L, ctl$max.candidates - evaluated.count())
        refine <- head(refine, remaining)
        if (length(refine)) {
            stage.rows <- do.call(rbind, lapply(refine, eval.one,
                                                stage = "refine"))
            selected <- .ps.lps.select.lambda.table(
                collect.evaluated(),
                rel.tol = ctl$rel.tol
            )
            telemetry <- append.telemetry(stage.rows, selected, telemetry)
        }
    }

    expansion <- 0L
    while (ctl$boundary.expand && expansion < ctl$max.boundary.expansions &&
           evaluated.count() < ctl$max.candidates &&
           any(collect.evaluated()$lambda.sync > 0)) {
        evaluated <- collect.evaluated()
        positive.evaluated <- evaluated$lambda.sync[evaluated$lambda.sync > 0]
        if (!length(positive.evaluated)) break
        raw.best <- .ps.lps.raw.best.lambda.table(evaluated)
        best.cv <- raw.best$cv.rmse.observed[[1L]]
        left.lambda <- min(positive.evaluated)
        right.lambda <- max(positive.evaluated)
        left.row <- evaluated[evaluated$lambda.sync == left.lambda, ,
                              drop = FALSE][1L, , drop = FALSE]
        right.row <- evaluated[evaluated$lambda.sync == right.lambda, ,
                               drop = FALSE][1L, , drop = FALSE]
        near.left <- is.finite(left.row$cv.rmse.observed[[1L]]) &&
            left.row$cv.rmse.observed[[1L]] <=
            best.cv * (1 + ctl$boundary.guard.rel.tol)
        near.right <- is.finite(right.row$cv.rmse.observed[[1L]]) &&
            right.row$cv.rmse.observed[[1L]] <=
            best.cv * (1 + ctl$boundary.guard.rel.tol)
        boundary <- "none"
        candidate <- NA_real_
        if (raw.best$lambda.sync[[1L]] > 0 &&
            raw.best$lambda.sync[[1L]] == left.lambda) {
            boundary <- "left"
            candidate <- raw.best$lambda.sync[[1L]] / ctl$boundary.factor
        } else if (raw.best$lambda.sync[[1L]] > 0 &&
                   raw.best$lambda.sync[[1L]] == right.lambda) {
            boundary <- "right"
            candidate <- raw.best$lambda.sync[[1L]] * ctl$boundary.factor
        } else if (selected$lambda.sync[[1L]] > 0 &&
                   selected$lambda.sync[[1L]] == left.lambda) {
            boundary <- "left"
            candidate <- selected$lambda.sync[[1L]] / ctl$boundary.factor
        } else if (selected$lambda.sync[[1L]] > 0 &&
                   selected$lambda.sync[[1L]] == right.lambda) {
            boundary <- "right"
            candidate <- selected$lambda.sync[[1L]] * ctl$boundary.factor
        } else if (near.right) {
            boundary <- "right"
            candidate <- right.lambda * ctl$boundary.factor
        } else if (near.left) {
            boundary <- "left"
            candidate <- left.lambda / ctl$boundary.factor
        }
        if (!is.finite(candidate) || boundary == "none" || candidate <= 0 ||
            candidate %in% evaluated$lambda.sync) {
            break
        }
        expansion <- expansion + 1L
        stage.rows <- eval.one(
            candidate,
            stage = paste0("boundary_expand_", boundary),
            boundary = boundary,
            expansion = expansion
        )
        selected <- .ps.lps.select.lambda.table(
            collect.evaluated(),
            rel.tol = ctl$rel.tol
        )
        telemetry <- append.telemetry(stage.rows, selected, telemetry)
    }

    evaluated <- collect.evaluated()
    list(
        evaluated = evaluated,
        selected = selected,
        telemetry = telemetry
    )
}

.ps.lps.prepare.chart.dim <- function(chart.dim, n, p) {
    if (length(chart.dim) == 1L) {
        out <- rep(as.integer(chart.dim), n)
    } else {
        out <- as.integer(chart.dim)
    }
    if (length(out) != n || any(!is.finite(out)) || any(out < 1L) ||
        any(out > p)) {
        stop("'chart.dim' must be a scalar or length-n integer vector in [1, ncol(X)].",
             call. = FALSE)
    }
    out
}

.ps.lps.resolve.chart.dim <- function(
    X, support.size, degree, chart.dim, auto.chart.support.metric,
    auto.chart.selection.metric) {

    if (identical(chart.dim, "auto")) {
        info <- .klp.resolve.chart.dim(
            X = X,
            support.size = support.size,
            degree = degree,
            coordinate.method = "local.pca",
            chart.dim = chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        dim <- as.integer(info$chart.dim)
        return(list(
            chart.dim.by.anchor = rep(dim, nrow(X)),
            chart.dim.mode = "global.auto",
            auto.chart.dim = dim,
            auto.chart.dim.diagnostics = info$diagnostics
        ))
    }
    if (.klp.is.local.auto.chart.dim(chart.dim)) {
        info <- .klp.resolve.chart.dim(
            X = X,
            support.size = support.size,
            degree = degree,
            coordinate.method = "local.pca",
            chart.dim = chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric
        )
        pred.info <- .klp.resolve.prediction.chart.dim(
            X.train = X,
            X.eval = X,
            support.size = support.size,
            degree = degree,
            coordinate.method = "local.pca",
            chart.dim = chart.dim,
            auto.chart.support.metric = auto.chart.support.metric,
            auto.chart.selection.metric = auto.chart.selection.metric,
            summary.dim = info$chart.dim
        )
        dims <- .ps.lps.prepare.chart.dim(
            chart.dim = pred.info$chart.dim.by.eval,
            n = nrow(X),
            p = ncol(X)
        )
        return(list(
            chart.dim.by.anchor = dims,
            chart.dim.mode = "local.auto",
            auto.chart.dim = as.integer(info$chart.dim),
            auto.chart.dim.diagnostics = info$diagnostics
        ))
    }
    list(
        chart.dim.by.anchor = .ps.lps.prepare.chart.dim(
            chart.dim,
            nrow(X),
            ncol(X)
        ),
        chart.dim.mode = "fixed",
        auto.chart.dim = NA_integer_,
        auto.chart.dim.diagnostics = NULL
    )
}

.ps.lps.prepare.geometry.cache <- function(X, support.size, degree, kernel,
                                           chart.dim,
                                           auto.chart.support.metric,
                                           auto.chart.selection.metric,
                                           sync.neighbor.size,
                                           overlap.weight,
                                           design.basis,
                                           design.drop.tol,
                                           local.pca.supports = NULL,
                                           chart.activation = "none",
                                           chart.activation.response = NULL,
                                           chart.activation.control = list()) {
    chart.activation.info <- .klp.prepare.chart.activation(
        chart.activation = chart.activation,
        chart.activation.response = chart.activation.response,
        fallback.response = rep(0, nrow(X)),
        n = nrow(X),
        control = chart.activation.control
    )
    chart.dim.info <- .ps.lps.resolve.chart.dim(
        X = X,
        support.size = support.size,
        degree = degree,
        chart.dim = chart.dim,
        auto.chart.support.metric = auto.chart.support.metric,
        auto.chart.selection.metric = auto.chart.selection.metric
    )
    frames <- .ps.lps.prepare.frames(
        X = X,
        y = rep(0, nrow(X)),
        support.size = support.size,
        degree = degree,
        kernel = kernel,
        chart.dim.by.anchor = chart.dim.info$chart.dim.by.anchor,
        design.basis = design.basis,
        design.drop.tol = design.drop.tol,
        local.pca.supports = local.pca.supports,
        chart.activation.info = chart.activation.info
    )
    sync.rows <- .ps.lps.prepare.sync.rows(
        frames = frames,
        sync.neighbor.size = sync.neighbor.size,
        overlap.weight = overlap.weight
    )
    list(
        chart.dim.info = chart.dim.info,
        local.pca.supports = local.pca.supports,
        chart.activation.info = chart.activation.info,
        chart.activation.diagnostics =
            attr(frames, "chart.activation.diagnostics"),
        frames = frames,
        sync.rows = sync.rows
    )
}

.ps.lps.prepare.frames <- function(X, y, support.size, degree, kernel,
                                   chart.dim.by.anchor,
                                   design.basis = "monomial",
                                   design.drop.tol =
                                       sqrt(.Machine$double.eps),
                                   local.pca.supports = NULL,
                                   chart.activation.info = NULL) {
    n <- nrow(X)
    frames <- vector("list", n)
    activation.rows <- vector("list", n)
    native.local <- local.pca.supports
    if (!is.null(native.local) && length(native.local) != n) {
        stop("'local.pca.supports' must contain one entry per anchor.",
             call. = FALSE)
    }
    if (is.null(native.local)) {
        native.local <- tryCatch(
            rcpp_ps_lps_local_pca_supports(
                X = as.matrix(X),
                support_size = as.integer(support.size),
                chart_dim_by_anchor = as.integer(chart.dim.by.anchor),
                kernel = as.character(kernel[[1L]])
            ),
            error = function(e) NULL
        )
    }
    for (ii in seq_len(n)) {
        d <- chart.dim.by.anchor[[ii]]
        if (!is.null(native.local)) {
            local <- native.local[[ii]]
            idx <- local$index
            dist <- local$distances
            weights <- local$weights
            coords <- as.matrix(local$coordinates)
            if (ncol(coords) < d) {
                stop("Precomputed local PCA coordinates have fewer columns ",
                     "than the requested chart dimension.", call. = FALSE)
            }
            coords <- coords[, seq_len(d), drop = FALSE]
        } else {
            ordered <- .klp.local.order(
                X.train = X,
                center = X[ii, , drop = TRUE],
                support.size = support.size
            )
            idx <- ordered$index
            dist <- ordered$distances
            weights <- .klp.kernel.weights(dist, kernel)
            if (!any(weights > 0)) weights[] <- 1
            coords <- .klp.local.coordinates(
                X.support = X[idx, , drop = FALSE],
                center = X[ii, , drop = TRUE],
                coordinate.method = "local.pca",
                chart.dim = d,
                local.chart.method = "pca",
                weights = weights,
                return.chart = FALSE
            )
        }
        activation <- .klp.chart.activation.row(
            info = chart.activation.info,
            anchor = ii,
            support.index = idx,
            weights = weights
        )
        if (!isTRUE(activation$active)) {
            row.by.point <- seq_along(idx)
            names(row.by.point) <- as.character(idx)
            frames[[ii]] <- list(
                anchor = ii,
                index = idx,
                distances = dist,
                weights = rep(0, length(idx)),
                design = matrix(0, nrow = length(idx), ncol = 1L),
                anchor.design = matrix(0, nrow = 1L, ncol = 1L),
                chart.dim = d,
                design.basis = design.basis,
                solver.design.basis = "inactive_zero",
                design.drop.tol = design.drop.tol,
                design.columns.original = 1L,
                design.columns.kept = integer(0),
                q = 1L,
                rank = 0L,
                row.by.point = row.by.point,
                active = FALSE,
                inactive.reason = activation$reason
            )
            activation.rows[[ii]] <- activation$row
            next
        }
        raw.design <- .local.polynomial.design.matrix(coords, degree)
        raw.anchor.design <- .local.polynomial.design.matrix(
            matrix(0, nrow = 1L, ncol = d),
            degree
        )
        keep <- seq_len(ncol(raw.design))
        solver.design.basis <- design.basis
        if (identical(design.basis, "orthogonal.polynomial.drop")) {
            transformed <- .klp.orthogonal.polynomial.transform(
                design = raw.design,
                weights = weights,
                prediction.rows = raw.anchor.design,
                design.drop.tol = design.drop.tol
            )
            if (isTRUE(transformed$ok)) {
                design <- transformed$design
                anchor.design <- transformed$prediction.rows
                keep <- seq_len(ncol(design))
                solver.design.basis <- "orthogonal.polynomial.transformed"
            } else {
                design <- raw.design[, 1L, drop = FALSE]
                anchor.design <- raw.anchor.design[, 1L, drop = FALSE]
                keep <- 1L
                solver.design.basis <- "monomial"
            }
        } else if (identical(design.basis, "weighted.qr.drop")) {
            keep <- .klp.weighted.qr.keep.columns(
                raw.design,
                weights,
                design.drop.tol
            )
            if (!length(keep)) keep <- 1L
            design <- raw.design[, keep, drop = FALSE]
            anchor.design <- raw.anchor.design[, keep, drop = FALSE]
        } else {
            design <- raw.design[, keep, drop = FALSE]
            anchor.design <- raw.anchor.design[, keep, drop = FALSE]
        }
        row.names <- as.character(idx)
        row.by.point <- seq_along(idx)
        names(row.by.point) <- row.names
        xw <- design * sqrt(weights)
        rank <- qr(xw, LAPACK = TRUE)$rank
        frames[[ii]] <- list(
            anchor = ii,
            index = idx,
            distances = dist,
            weights = weights,
            design = design,
            anchor.design = anchor.design,
            chart.dim = d,
            design.basis = design.basis,
            solver.design.basis = solver.design.basis,
            design.drop.tol = design.drop.tol,
            design.columns.original = ncol(raw.design),
            design.columns.kept = keep,
            q = ncol(design),
            rank = as.integer(rank),
            row.by.point = row.by.point,
            active = TRUE,
            inactive.reason = "active"
        )
        activation.rows[[ii]] <- activation$row
    }
    q <- vapply(frames, `[[`, integer(1L), "q")
    offsets <- cumsum(c(0L, head(q, -1L)))
    for (ii in seq_along(frames)) frames[[ii]]$offset <- offsets[[ii]]
    attr(frames, "ncoef") <- sum(q)
    attr(frames, "chart.activation.diagnostics") <- do.call(
        rbind,
        activation.rows
    )
    frames
}

.ps.lps.frame.design.summary <- function(frames) {
    data.frame(
        anchor = vapply(frames, `[[`, integer(1L), "anchor"),
        chart.dim = vapply(frames, `[[`, integer(1L), "chart.dim"),
        design.columns.original =
            vapply(frames, `[[`, integer(1L), "design.columns.original"),
        design.columns.kept = vapply(frames, `[[`, integer(1L), "q"),
        design.rank = vapply(frames, `[[`, integer(1L), "rank"),
        active = vapply(frames, function(fr) isTRUE(fr$active %||% TRUE),
                        logical(1L)),
        inactive.reason = vapply(frames, function(fr) {
            as.character(fr$inactive.reason %||% "active")
        }, character(1L)),
        stringsAsFactors = FALSE
    )
}

.ps.lps.prepare.sync.rows <- function(frames, sync.neighbor.size,
                                      overlap.weight) {
    any.inactive <- any(vapply(
        frames,
        function(fr) isFALSE(fr$active %||% TRUE),
        logical(1L)
    ))
    native.rows <- if (isTRUE(any.inactive)) {
        NULL
    } else {
        tryCatch(
            rcpp_ps_lps_prepare_sync_rows(
                frames = frames,
                sync_neighbor_size = as.integer(sync.neighbor.size),
                overlap_weight = as.character(overlap.weight[[1L]])
            ),
            error = function(e) NULL
        )
    }
    if (!is.null(native.rows)) {
        return(native.rows)
    }
    n <- length(frames)
    pair.keys <- character(0)
    for (ii in seq_len(n)) {
        nbr <- frames[[ii]]$index[frames[[ii]]$index != ii]
        nbr <- head(nbr, max(0L, as.integer(sync.neighbor.size)))
        if (!length(nbr)) next
        lo <- pmin(ii, nbr)
        hi <- pmax(ii, nbr)
        pair.keys <- c(pair.keys, paste(lo, hi, sep = "_"))
    }
    pair.keys <- unique(pair.keys)
    out <- vector("list", length(pair.keys))
    keep <- logical(length(pair.keys))
    for (pp in seq_along(pair.keys)) {
        pair <- as.integer(strsplit(pair.keys[[pp]], "_", fixed = TRUE)[[1L]])
        ii <- pair[[1L]]
        jj <- pair[[2L]]
        if (isFALSE(frames[[ii]]$active %||% TRUE) ||
            isFALSE(frames[[jj]]$active %||% TRUE)) {
            next
        }
        overlap <- intersect(frames[[ii]]$index, frames[[jj]]$index)
        if (!length(overlap)) next
        row.i <- frames[[ii]]$row.by.point[as.character(overlap)]
        row.j <- frames[[jj]]$row.by.point[as.character(overlap)]
        prod <- frames[[ii]]$weights[row.i] * frames[[jj]]$weights[row.j]
        if (identical(overlap.weight, "normalized.product")) {
            omega <- length(overlap) * prod / (sum(prod) + sqrt(.Machine$double.eps))
        } else {
            omega <- prod
        }
        ok <- is.finite(omega) & omega > 0
        if (!any(ok)) next
        out[[pp]] <- list(
            i = ii,
            j = jj,
            point = overlap[ok],
            row.i = as.integer(row.i[ok]),
            row.j = as.integer(row.j[ok]),
            omega = as.numeric(omega[ok])
        )
        keep[[pp]] <- TRUE
    }
    out[keep]
}

.ps.lps.prepare.system.cache <- function(frames, sync.rows) {
    data.blocks <- lapply(frames, function(fr) {
        qseq <- seq_len(fr$q)
        list(
            index = fr$index,
            weights = fr$weights,
            design = fr$design,
            q = fr$q,
            cols = fr$offset + qseq
        )
    })
    sync.blocks <- vector("list", length(sync.rows))
    if (length(sync.rows)) {
        for (ss in seq_along(sync.rows)) {
            sr <- sync.rows[[ss]]
            fi <- frames[[sr$i]]
            fj <- frames[[sr$j]]
            qi <- seq_len(fi$q)
            qj <- seq_len(fj$q)
            sync.blocks[[ss]] <- list(
                i = sr$i,
                j = sr$j,
                point = sr$point,
                omega.sqrt = sqrt(sr$omega),
                cols.i = fi$offset + qi,
                cols.j = fj$offset + qj,
                values.i = fi$design[sr$row.i, qi, drop = FALSE],
                values.j = fj$design[sr$row.j, qj, drop = FALSE],
                q.i = fi$q,
                q.j = fj$q
            )
        }
    }
    structure(
        list(
            frames = frames,
            sync.rows = sync.rows,
            data.blocks = data.blocks,
            sync.blocks = sync.blocks,
            ncoef = attr(frames, "ncoef")
        ),
        class = "ps_lps_system_cache"
    )
}

.ps.lps.prepare.normal.cache <- function(cache, y, response.weights,
                                         lambda.sync) {
    if (!inherits(cache, "ps_lps_system_cache")) {
        stop("'cache' must be a PS-LPS system cache.", call. = FALSE)
    }
    lambda.sync <- as.numeric(lambda.sync[[1L]])
    if (!is.finite(lambda.sync) || lambda.sync <= 0) {
        stop("'lambda.sync' must be positive for a PS-LPS normal cache.",
             call. = FALSE)
    }
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])
    t.native <- proc.time()
    assembled <- rcpp_ps_lps_assemble_cached_system(
        cache = cache,
        y = y,
        response_weights = response.weights,
        lambda_sync = lambda.sync
    )
    phase.native.sec <- elapsed(t.native)
    t.sparse <- proc.time()
    A <- Matrix::sparseMatrix(
        i = assembled$rows,
        j = assembled$cols,
        x = assembled$vals,
        dims = c(assembled$nrow, assembled$ncol)
    )
    phase.sparse.sec <- elapsed(t.sparse)
    t.cross <- proc.time()
    cross <- Matrix::crossprod(A)
    phase.cross.sec <- elapsed(t.cross)
    t.rhs <- proc.time()
    rhs.cross <- Matrix::crossprod(A, assembled$rhs)
    phase.rhs.sec <- elapsed(t.rhs)
    scale <- max(Matrix::diag(cross), na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    timings <- list(
        phase_count_sec = NA_real_,
        phase_fill_triplets_sec = phase.native.sec,
        phase_sparse_matrix_sec = phase.sparse.sec,
        phase_assembly_sec = phase.native.sec + phase.sparse.sec,
        phase_crossprod_sec = phase.cross.sec,
        phase_rhs_crossprod_sec = phase.rhs.sec,
        phase_ridge_normal_sec = NA_real_,
        phase_solve_sec = NA_real_,
        phase_fallback_solve_sec = 0,
        phase_diagnostics_sec = NA_real_,
        phase_fitted_sec = NA_real_,
        phase_normal_cache_sec = phase.native.sec + phase.sparse.sec +
            phase.cross.sec + phase.rhs.sec
    )
    structure(
        list(
            frames = cache$frames,
            sync.rows = cache$sync.rows,
            y = y,
            response.weights = response.weights,
            lambda.sync = lambda.sync,
            cross = cross,
            rhs.cross = rhs.cross,
            scale = scale,
            n.system.rows = assembled$nrow,
            n.system.cols = assembled$ncol,
            n.system.nnz = assembled$nnz,
            solve.phase.timings = timings
        ),
        class = "ps_lps_normal_cache"
    )
}

.ps.lps.prepare.component.cache <- function(cache, y, response.weights) {
    if (!inherits(cache, "ps_lps_system_cache")) {
        stop("'cache' must be a PS-LPS system cache.", call. = FALSE)
    }
    y <- as.numeric(y)
    response.weights <- as.numeric(response.weights)
    if (length(y) != length(response.weights) ||
        any(!is.finite(y)) || any(!is.finite(response.weights))) {
        stop("'y' and 'response.weights' must be finite vectors of equal length.",
             call. = FALSE)
    }
    ncoef <- cache$ncoef
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])

    t.data.native <- proc.time()
    data.assembled <- rcpp_ps_lps_assemble_cached_system(
        cache = cache,
        y = y,
        response_weights = response.weights,
        lambda_sync = 0
    )
    phase.data.native.sec <- elapsed(t.data.native)
    t.data.sparse <- proc.time()
    A.data <- Matrix::sparseMatrix(
        i = data.assembled$rows,
        j = data.assembled$cols,
        x = data.assembled$vals,
        dims = c(data.assembled$nrow, data.assembled$ncol)
    )
    phase.data.sparse.sec <- elapsed(t.data.sparse)
    t.data.cross <- proc.time()
    cross.data <- Matrix::crossprod(A.data)
    phase.data.cross.sec <- elapsed(t.data.cross)
    t.data.rhs <- proc.time()
    rhs.data <- Matrix::crossprod(A.data, data.assembled$rhs)
    phase.data.rhs.sec <- elapsed(t.data.rhs)

    sync.nrow <- 0L
    sync.nnz <- 0L
    cross.sync <- Matrix::sparseMatrix(
        i = integer(),
        j = integer(),
        x = numeric(),
        dims = c(ncoef, ncoef)
    )
    phase.sync.native.sec <- 0
    phase.sync.sparse.sec <- 0
    phase.sync.cross.sec <- 0
    if (length(cache$sync.blocks)) {
        t.sync.native <- proc.time()
        sync.assembled <- rcpp_ps_lps_assemble_cached_system(
            cache = cache,
            y = y,
            response_weights = rep(0, length(y)),
            lambda_sync = 1
        )
        phase.sync.native.sec <- elapsed(t.sync.native)
        sync.nrow <- sync.assembled$sync_nrow
        sync.nnz <- sync.assembled$nnz
        t.sync.sparse <- proc.time()
        A.sync <- Matrix::sparseMatrix(
            i = sync.assembled$rows,
            j = sync.assembled$cols,
            x = sync.assembled$vals,
            dims = c(sync.assembled$nrow, sync.assembled$ncol)
        )
        phase.sync.sparse.sec <- elapsed(t.sync.sparse)
        t.sync.cross <- proc.time()
        cross.sync <- Matrix::crossprod(A.sync)
        phase.sync.cross.sec <- elapsed(t.sync.cross)
    }

    timings <- list(
        phase_data_native_sec = phase.data.native.sec,
        phase_data_sparse_matrix_sec = phase.data.sparse.sec,
        phase_data_crossprod_sec = phase.data.cross.sec,
        phase_data_rhs_crossprod_sec = phase.data.rhs.sec,
        phase_sync_native_sec = phase.sync.native.sec,
        phase_sync_sparse_matrix_sec = phase.sync.sparse.sec,
        phase_sync_crossprod_sec = phase.sync.cross.sec,
        phase_component_cache_sec = phase.data.native.sec +
            phase.data.sparse.sec + phase.data.cross.sec +
            phase.data.rhs.sec + phase.sync.native.sec +
            phase.sync.sparse.sec + phase.sync.cross.sec
    )

    structure(
        list(
            frames = cache$frames,
            sync.rows = cache$sync.rows,
            y = y,
            response.weights = response.weights,
            cross.data = cross.data,
            cross.sync = cross.sync,
            rhs.data = rhs.data,
            n.system.cols = ncoef,
            data.nrow = data.assembled$data_nrow,
            data.nnz = data.assembled$nnz,
            sync.nrow = sync.nrow,
            sync.nnz = sync.nnz,
            solve.phase.timings = timings
        ),
        class = "ps_lps_component_cache"
    )
}

.ps.lps.component.normal.cache <- function(component.cache, lambda.sync) {
    if (!inherits(component.cache, "ps_lps_component_cache")) {
        stop("'component.cache' must be a PS-LPS component cache.",
             call. = FALSE)
    }
    lambda.sync <- as.numeric(lambda.sync[[1L]])
    if (!is.finite(lambda.sync) || lambda.sync <= 0) {
        stop("'lambda.sync' must be positive for a PS-LPS component normal cache.",
             call. = FALSE)
    }
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])
    t.combine <- proc.time()
    cross <- component.cache$cross.data +
        lambda.sync * component.cache$cross.sync
    scale <- max(Matrix::diag(cross), na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    phase.combine.sec <- elapsed(t.combine)
    timings <- c(
        component.cache$solve.phase.timings,
        list(
            phase_component_combine_sec = phase.combine.sec,
            phase_count_sec = NA_real_,
            phase_fill_triplets_sec = NA_real_,
            phase_sparse_matrix_sec = NA_real_,
            phase_assembly_sec = NA_real_,
            phase_crossprod_sec = NA_real_,
            phase_rhs_crossprod_sec = NA_real_,
            phase_ridge_normal_sec = NA_real_,
            phase_solve_sec = NA_real_,
            phase_fallback_solve_sec = 0,
            phase_diagnostics_sec = NA_real_,
            phase_fitted_sec = NA_real_,
            phase_normal_cache_sec = phase.combine.sec
        )
    )
    structure(
        list(
            frames = component.cache$frames,
            sync.rows = component.cache$sync.rows,
            y = component.cache$y,
            response.weights = component.cache$response.weights,
            lambda.sync = lambda.sync,
            cross = cross,
            rhs.cross = component.cache$rhs.data,
            scale = scale,
            n.system.rows = component.cache$data.nrow + component.cache$sync.nrow,
            n.system.cols = component.cache$n.system.cols,
            n.system.nnz = component.cache$data.nnz + component.cache$sync.nnz,
            solve.phase.timings = timings
        ),
        class = "ps_lps_normal_cache"
    )
}

.ps.lps.matrix.condition <- function(M) {
    rc <- tryCatch(
        rcond(as.matrix(M)),
        error = function(e) NA_real_
    )
    rc <- as.numeric(rc[[1L]])
    if (!is.finite(rc) || rc <= 0) return(Inf)
    1 / rc
}

.ps.lps.choose.global.ridge <- function(cross, scale,
                                        ridge.multiplier.grid,
                                        ridge.condition.max,
                                        sparse = TRUE) {
    ridge.multiplier.grid <- .klp.clean.ridge.multiplier.grid(
        ridge.multiplier.grid
    )
    scale <- as.numeric(scale[[1L]])
    if (!is.finite(scale) || scale <= 0) scale <- 1
    ncoef <- ncol(cross)
    if (!is.finite(ridge.condition.max)) {
        rho <- ridge.multiplier.grid[[1L]]
        ridge <- rho * scale
        normal <- if (isTRUE(sparse)) {
            cross + Matrix::Diagonal(ncoef, x = ridge)
        } else {
            cross + diag(ridge, ncoef)
        }
        return(list(
            ridge.multiplier = rho,
            ridge = ridge,
            normal = normal,
            condition = NA_real_,
            status = "unguarded"
        ))
    }
    for (rho in ridge.multiplier.grid) {
        ridge <- rho * scale
        normal <- if (isTRUE(sparse)) {
            cross + Matrix::Diagonal(ncoef, x = ridge)
        } else {
            cross + diag(ridge, ncoef)
        }
        cond <- .ps.lps.matrix.condition(normal)
        if (!is.finite(ridge.condition.max) ||
            (is.finite(cond) && cond <= ridge.condition.max)) {
            return(list(
                ridge.multiplier = rho,
                ridge = ridge,
                normal = normal,
                condition = cond,
                status = "ok"
            ))
        }
    }
    rho <- utils::tail(ridge.multiplier.grid, 1L)
    ridge <- rho * scale
    normal <- if (isTRUE(sparse)) {
        cross + Matrix::Diagonal(ncoef, x = ridge)
    } else {
        cross + diag(ridge, ncoef)
    }
    list(
        ridge.multiplier = rho,
        ridge = ridge,
        normal = normal,
        condition = .ps.lps.matrix.condition(normal),
        status = "condition_not_met"
    )
}

.ps.lps.solve.component.cached <- function(component.cache, lambda.sync,
                                          lambda.ridge = 1e-8,
                                          ridge.multiplier.grid = lambda.ridge,
                                          ridge.condition.max = Inf,
                                          coefficients.only = FALSE) {
    normal.cache <- .ps.lps.component.normal.cache(
        component.cache = component.cache,
        lambda.sync = lambda.sync
    )
    .ps.lps.solve.normal.cached(
        normal.cache = normal.cache,
        lambda.ridge = lambda.ridge,
        ridge.multiplier.grid = ridge.multiplier.grid,
        ridge.condition.max = ridge.condition.max,
        coefficients.only = coefficients.only
    )
}

.ps.lps.solve.cached <- function(cache, y, response.weights, lambda.sync,
                                 lambda.ridge = 1e-8,
                                 ridge.multiplier.grid = lambda.ridge,
                                 ridge.condition.max = Inf,
                                 coefficients.only = FALSE) {
    if (!inherits(cache, "ps_lps_system_cache")) {
        stop("'cache' must be a PS-LPS system cache.", call. = FALSE)
    }
    frames <- cache$frames
    sync.rows <- cache$sync.rows
    n <- length(frames)
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])
    if (lambda.sync == 0) {
        solve.start <- proc.time()
        independent <- .ps.lps.solve.independent(
            frames = frames,
            y = y,
            response.weights = response.weights,
            lambda.ridge = lambda.ridge,
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max
        )
        diagnostics <- .ps.lps.diagnostics(
            frames = frames,
            y = y,
            beta = independent$coefficients,
            sync.rows = sync.rows,
            ridge.values = independent$ridge.values
        )
        diagnostics$ridge.multiplier.selected <-
            independent$ridge.multiplier.selected
        diagnostics$ridge.condition <- independent$ridge.condition
        diagnostics$ridge.status <- independent$ridge.status
        timings <- list(
            phase_count_sec = NA_real_,
            phase_fill_triplets_sec = NA_real_,
            phase_sparse_matrix_sec = NA_real_,
            phase_assembly_sec = NA_real_,
            phase_crossprod_sec = NA_real_,
            phase_rhs_crossprod_sec = NA_real_,
            phase_ridge_normal_sec = NA_real_,
            phase_solve_sec = elapsed(solve.start),
            phase_fallback_solve_sec = 0,
            phase_diagnostics_sec = NA_real_,
            phase_fitted_sec = NA_real_
        )
        diagnostics <- c(diagnostics, list(solve.phase.timings = timings))
        if (coefficients.only) return(diagnostics)
        return(c(independent, diagnostics))
    }
    normal.cache <- .ps.lps.prepare.normal.cache(
        cache = cache,
        y = y,
        response.weights = response.weights,
        lambda.sync = lambda.sync
    )
    .ps.lps.solve.normal.cached(
        normal.cache = normal.cache,
        lambda.ridge = lambda.ridge,
        ridge.multiplier.grid = ridge.multiplier.grid,
        ridge.condition.max = ridge.condition.max,
        coefficients.only = coefficients.only
    )
}

.ps.lps.solve.normal.cached <- function(normal.cache, lambda.ridge = 1e-8,
                                        ridge.multiplier.grid = lambda.ridge,
                                        ridge.condition.max = Inf,
                                        coefficients.only = FALSE) {
    if (!inherits(normal.cache, "ps_lps_normal_cache")) {
        stop("'normal.cache' must be a PS-LPS normal cache.", call. = FALSE)
    }
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])
    frames <- normal.cache$frames
    sync.rows <- normal.cache$sync.rows
    y <- normal.cache$y
    cross <- normal.cache$cross
    rhs.cross <- normal.cache$rhs.cross
    ncoef <- normal.cache$n.system.cols
    n <- length(frames)
    t.ridge <- proc.time()
    scale <- normal.cache$scale
    ridge.choice <- .ps.lps.choose.global.ridge(
        cross = cross,
        scale = scale,
        ridge.multiplier.grid = ridge.multiplier.grid,
        ridge.condition.max = ridge.condition.max
    )
    ridge <- ridge.choice$ridge
    normal <- ridge.choice$normal
    phase.ridge.sec <- elapsed(t.ridge)
    t.solve <- proc.time()
    fallback.sec <- 0
    beta <- tryCatch(
        as.numeric(Matrix::solve(normal, rhs.cross)),
        error = function(e) rep(NA_real_, ncoef)
    )
    phase.solve.sec <- elapsed(t.solve)
    if (any(!is.finite(beta))) {
        t.fallback <- proc.time()
        ridge <- sqrt(.Machine$double.eps) * scale
        normal <- cross + Matrix::Diagonal(ncoef, x = ridge)
        beta <- as.numeric(Matrix::solve(normal, rhs.cross))
        fallback.sec <- elapsed(t.fallback)
        ridge.choice$status <- "fallback_ridge"
        ridge.choice$ridge.multiplier <- sqrt(.Machine$double.eps)
        ridge.choice$condition <- .ps.lps.matrix.condition(normal)
    }
    t.diagnostics <- proc.time()
    diagnostics <- .ps.lps.diagnostics(
        frames = frames,
        y = y,
        beta = beta,
        sync.rows = sync.rows,
        ridge.values = ridge
    )
    phase.diagnostics.sec <- elapsed(t.diagnostics)
    timings <- normal.cache$solve.phase.timings
    timings$phase_ridge_normal_sec <- phase.ridge.sec
    timings$phase_solve_sec <- phase.solve.sec
    timings$phase_fallback_solve_sec <- fallback.sec
    timings$phase_diagnostics_sec <- phase.diagnostics.sec
    timings$phase_fitted_sec <- NA_real_
    diagnostics <- c(diagnostics, list(solve.phase.timings = timings))
    diagnostics$ridge.multiplier.selected <- ridge.choice$ridge.multiplier
    diagnostics$ridge.condition <- ridge.choice$condition
    diagnostics$ridge.status <- ridge.choice$status
    if (coefficients.only) return(diagnostics)
    t.fitted <- proc.time()
    fitted <- vapply(seq_len(n), function(ii) {
        fr <- frames[[ii]]
        idx <- fr$offset + seq_len(fr$q)
        sum(fr$anchor.design[1L, ] * beta[idx])
    }, numeric(1L))
    diagnostics$solve.phase.timings$phase_fitted_sec <- elapsed(t.fitted)
    out <- c(
        list(
            coefficients = beta,
            fitted.values = fitted,
            n.system.rows = normal.cache$n.system.rows,
            n.system.cols = ncoef,
            lambda.ridge = lambda.ridge
        ),
        diagnostics
    )
    out$ridge.multiplier.selected <- ridge.choice$ridge.multiplier
    out$ridge.condition <- ridge.choice$condition
    out$ridge.status <- ridge.choice$status
    out
}

.ps.lps.solve <- function(frames, y, response.weights, lambda.sync, sync.rows,
                          lambda.ridge = 1e-8,
                          ridge.multiplier.grid = lambda.ridge,
                          ridge.condition.max = Inf,
                          coefficients.only = FALSE) {
    elapsed <- function(start) unname((proc.time() - start)[["elapsed"]])
    if (lambda.sync == 0) {
        solve.start <- proc.time()
        independent <- .ps.lps.solve.independent(
            frames = frames,
            y = y,
            response.weights = response.weights,
            lambda.ridge = lambda.ridge,
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max
        )
        diagnostics <- .ps.lps.diagnostics(
            frames = frames,
            y = y,
            beta = independent$coefficients,
            sync.rows = sync.rows,
            ridge.values = independent$ridge.values
        )
        diagnostics$ridge.multiplier.selected <-
            independent$ridge.multiplier.selected
        diagnostics$ridge.condition <- independent$ridge.condition
        diagnostics$ridge.status <- independent$ridge.status
        timings <- list(
            phase_count_sec = NA_real_,
            phase_fill_triplets_sec = NA_real_,
            phase_sparse_matrix_sec = NA_real_,
            phase_assembly_sec = NA_real_,
            phase_crossprod_sec = NA_real_,
            phase_rhs_crossprod_sec = NA_real_,
            phase_ridge_normal_sec = NA_real_,
            phase_solve_sec = elapsed(solve.start),
            phase_fallback_solve_sec = 0,
            phase_diagnostics_sec = NA_real_,
            phase_fitted_sec = NA_real_
        )
        diagnostics <- c(diagnostics, list(solve.phase.timings = timings))
        if (coefficients.only) return(diagnostics)
        return(c(independent, diagnostics))
    }
    n <- length(frames)
    ncoef <- attr(frames, "ncoef")
    t.count <- proc.time()
    data.nrows <- 0L
    data.nnz <- 0L
    for (ii in seq_len(n)) {
        fr <- frames[[ii]]
        rw <- response.weights[fr$index]
        ww <- rw * fr$weights
        ok <- is.finite(rw) & is.finite(ww) & ww > 0
        if (any(ok)) {
            n.ok <- sum(ok)
            data.nrows <- data.nrows + n.ok
            data.nnz <- data.nnz + n.ok * fr$q
        }
    }
    sync.nrows <- 0L
    sync.nnz <- 0L
    if (lambda.sync > 0 && length(sync.rows)) {
        for (sr in sync.rows) {
            n.sr <- length(sr$point)
            if (!n.sr) next
            fi <- frames[[sr$i]]
            fj <- frames[[sr$j]]
            sync.nrows <- sync.nrows + n.sr
            sync.nnz <- sync.nnz + n.sr * (fi$q + fj$q)
        }
    }
    phase.count.sec <- elapsed(t.count)
    rr <- data.nrows + sync.nrows
    nnz <- data.nnz + sync.nnz
    if (!rr || !nnz) stop("PS-LPS system has no rows.", call. = FALSE)
    t.fill <- proc.time()
    rows <- integer(nnz)
    cols <- integer(nnz)
    vals <- numeric(nnz)
    rhs <- numeric(rr)
    row.pos <- 0L
    nz.pos <- 0L
    for (ii in seq_len(n)) {
        fr <- frames[[ii]]
        qseq <- seq_len(fr$q)
        rw <- response.weights[fr$index]
        ww <- rw * fr$weights
        ok <- is.finite(rw) & is.finite(ww) & ww > 0
        if (!any(ok)) next
        for (aa in which(ok)) {
            point <- fr$index[[aa]]
            sw <- sqrt(ww[[aa]])
            row.pos <- row.pos + 1L
            next.nz <- nz.pos + fr$q
            idx <- seq.int(nz.pos + 1L, next.nz)
            rows[idx] <- row.pos
            cols[idx] <- fr$offset + qseq
            vals[idx] <- sw * fr$design[aa, qseq]
            rhs[[row.pos]] <- sw * y[[point]]
            nz.pos <- next.nz
        }
    }
    if (lambda.sync > 0 && length(sync.rows)) {
        for (sr in sync.rows) {
            n.sr <- length(sr$point)
            if (!n.sr) next
            fi <- frames[[sr$i]]
            fj <- frames[[sr$j]]
            qi <- seq_len(fi$q)
            qj <- seq_len(fj$q)
            scale <- sqrt(lambda.sync * sr$omega)
            for (aa in seq_len(n.sr)) {
                row.pos <- row.pos + 1L
                len.i <- fi$q
                len.j <- fj$q
                next.nz.i <- nz.pos + len.i
                idx.i <- seq.int(nz.pos + 1L, next.nz.i)
                rows[idx.i] <- row.pos
                cols[idx.i] <- fi$offset + qi
                vals[idx.i] <- scale[[aa]] * fi$design[sr$row.i[[aa]], qi]
                next.nz.j <- next.nz.i + len.j
                idx.j <- seq.int(next.nz.i + 1L, next.nz.j)
                rows[idx.j] <- row.pos
                cols[idx.j] <- fj$offset + qj
                vals[idx.j] <- -scale[[aa]] * fj$design[sr$row.j[[aa]], qj]
                nz.pos <- next.nz.j
            }
        }
    }
    if (row.pos != rr || nz.pos != nnz) {
        stop("Internal PS-LPS sparse assembly count mismatch.", call. = FALSE)
    }
    phase.fill.sec <- elapsed(t.fill)
    t.sparse <- proc.time()
    A <- Matrix::sparseMatrix(
        i = rows,
        j = cols,
        x = vals,
        dims = c(rr, ncoef)
    )
    phase.sparse.sec <- elapsed(t.sparse)
    t.cross <- proc.time()
    cross <- Matrix::crossprod(A)
    phase.cross.sec <- elapsed(t.cross)
    t.ridge <- proc.time()
    scale <- max(Matrix::diag(cross), na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    ridge.choice <- .ps.lps.choose.global.ridge(
        cross = cross,
        scale = scale,
        ridge.multiplier.grid = ridge.multiplier.grid,
        ridge.condition.max = ridge.condition.max
    )
    ridge <- ridge.choice$ridge
    normal <- ridge.choice$normal
    phase.ridge.sec <- elapsed(t.ridge)
    t.rhs <- proc.time()
    rhs.cross <- Matrix::crossprod(A, rhs)
    phase.rhs.sec <- elapsed(t.rhs)
    t.solve <- proc.time()
    fallback.sec <- 0
    beta <- tryCatch(
        as.numeric(Matrix::solve(normal, rhs.cross)),
        error = function(e) rep(NA_real_, ncoef)
    )
    phase.solve.sec <- elapsed(t.solve)
    if (any(!is.finite(beta))) {
        t.fallback <- proc.time()
        ridge <- sqrt(.Machine$double.eps) * scale
        normal <- cross + Matrix::Diagonal(ncoef, x = ridge)
        beta <- as.numeric(Matrix::solve(normal, rhs.cross))
        fallback.sec <- elapsed(t.fallback)
        ridge.choice$status <- "fallback_ridge"
        ridge.choice$ridge.multiplier <- sqrt(.Machine$double.eps)
        ridge.choice$condition <- .ps.lps.matrix.condition(normal)
    }
    t.diagnostics <- proc.time()
    diagnostics <- .ps.lps.diagnostics(
        frames = frames,
        y = y,
        beta = beta,
        sync.rows = sync.rows,
        ridge.values = ridge
    )
    phase.diagnostics.sec <- elapsed(t.diagnostics)
    timings <- list(
        phase_count_sec = phase.count.sec,
        phase_fill_triplets_sec = phase.fill.sec,
        phase_sparse_matrix_sec = phase.sparse.sec,
        phase_assembly_sec = phase.count.sec + phase.fill.sec + phase.sparse.sec,
        phase_crossprod_sec = phase.cross.sec,
        phase_rhs_crossprod_sec = phase.rhs.sec,
        phase_ridge_normal_sec = phase.ridge.sec,
        phase_solve_sec = phase.solve.sec,
        phase_fallback_solve_sec = fallback.sec,
        phase_diagnostics_sec = phase.diagnostics.sec,
        phase_fitted_sec = NA_real_
    )
    diagnostics <- c(diagnostics, list(solve.phase.timings = timings))
    diagnostics$ridge.multiplier.selected <- ridge.choice$ridge.multiplier
    diagnostics$ridge.condition <- ridge.choice$condition
    diagnostics$ridge.status <- ridge.choice$status
    if (coefficients.only) return(diagnostics)
    t.fitted <- proc.time()
    fitted <- vapply(seq_len(n), function(ii) {
        fr <- frames[[ii]]
        idx <- fr$offset + seq_len(fr$q)
        sum(fr$anchor.design[1L, ] * beta[idx])
    }, numeric(1L))
    diagnostics$solve.phase.timings$phase_fitted_sec <- elapsed(t.fitted)
    out <- c(
        list(
            coefficients = beta,
            fitted.values = fitted,
            n.system.rows = rr,
            n.system.cols = ncoef,
            lambda.ridge = lambda.ridge
        ),
        diagnostics
    )
    out$ridge.multiplier.selected <- ridge.choice$ridge.multiplier
    out$ridge.condition <- ridge.choice$condition
    out$ridge.status <- ridge.choice$status
    out
}

.ps.lps.solve.independent <- function(frames, y, response.weights,
                                      lambda.ridge = 0,
                                      ridge.multiplier.grid = lambda.ridge,
                                      ridge.condition.max = Inf) {
    n <- length(frames)
    ncoef <- attr(frames, "ncoef")
    beta <- rep(0, ncoef)
    ridge.values <- numeric(n)
    ridge.multiplier.values <- rep(NA_real_, n)
    ridge.condition.values <- rep(NA_real_, n)
    ridge.status.values <- rep(NA_character_, n)
    for (ii in seq_len(n)) {
        fr <- frames[[ii]]
        qseq <- seq_len(fr$q)
        w <- fr$weights * response.weights[fr$index]
        ok <- is.finite(y[fr$index]) & is.finite(w) & w > 0
        if (!any(ok)) {
            beta[fr$offset + qseq] <- 0
            ridge.values[[ii]] <- 0
            ridge.multiplier.values[[ii]] <- 0
            ridge.status.values[[ii]] <- "empty"
            next
        }
        design <- fr$design[ok, , drop = FALSE]
        yy <- y[fr$index][ok]
        ww <- w[ok]
        solved <- .klp.solve.local.wls(
            design = design,
            y = yy,
            weights = ww,
            design.basis = fr$solver.design.basis %||%
                fr$design.basis %||% "monomial",
            design.drop.tol = fr$design.drop.tol %||%
                sqrt(.Machine$double.eps),
            ridge.multiplier.grid = ridge.multiplier.grid,
            ridge.condition.max = ridge.condition.max
        )
        if (!is.null(solved) && isTRUE(solved$ok) &&
            length(solved$coefficients) && all(is.finite(solved$coefficients))) {
            coef <- rep(0, fr$q)
            keep <- solved$kept.columns %||% seq_along(solved$coefficients)
            keep <- keep[keep >= 1L & keep <= fr$q]
            coef[keep] <- solved$coefficients[seq_along(keep)]
            beta[fr$offset + qseq] <- coef
            ridge.values[[ii]] <- solved$ridge.lambda %||% 0
            ridge.multiplier.values[[ii]] <- solved$ridge.multiplier %||% 0
            ridge.condition.values[[ii]] <- solved$condition %||% NA_real_
            ridge.status.values[[ii]] <- solved$status %||% "ok"
        } else {
            legacy.unguarded <- identical(fr$design.basis %||% "monomial",
                                           "monomial") &&
                length(ridge.multiplier.grid) == 1L &&
                ridge.multiplier.grid[[1L]] == 0 &&
                !is.finite(ridge.condition.max)
            if (isTRUE(legacy.unguarded)) {
                coef <- c(stats::weighted.mean(yy, ww),
                          rep(0, fr$q - 1L))
                beta[fr$offset + qseq] <- coef
                ridge.values[[ii]] <- 0
                ridge.multiplier.values[[ii]] <- 0
                ridge.condition.values[[ii]] <- solved$condition %||% NA_real_
                ridge.status.values[[ii]] <- paste0(
                    "legacy_mean_fallback_",
                    solved$status %||% "local_solve_failed"
                )
            } else {
                beta[fr$offset + qseq] <- NA_real_
                ridge.values[[ii]] <- NA_real_
                ridge.multiplier.values[[ii]] <- NA_real_
                ridge.condition.values[[ii]] <- solved$condition %||% NA_real_
                ridge.status.values[[ii]] <- paste0(
                    "unstable_",
                    solved$status %||% "local_solve_failed"
                )
            }
        }
    }
    fitted <- vapply(seq_len(n), function(ii) {
        fr <- frames[[ii]]
        idx <- fr$offset + seq_len(fr$q)
        sum(fr$anchor.design[1L, ] * beta[idx])
    }, numeric(1L))
    list(
        coefficients = beta,
        fitted.values = fitted,
        n.system.rows = sum(vapply(frames, function(fr) length(fr$index),
                                   integer(1L))),
        n.system.cols = ncoef,
        lambda.ridge = lambda.ridge,
        ridge.values = ridge.values,
        ridge.multiplier.selected = {
            mult <- ridge.multiplier.values[is.finite(ridge.multiplier.values)]
            if (length(mult)) max(mult) else NA_real_
        },
        ridge.condition = {
            cond <- ridge.condition.values[is.finite(ridge.condition.values)]
            if (length(cond)) max(cond) else NA_real_
        },
        ridge.status = paste(unique(stats::na.omit(ridge.status.values)),
                             collapse = ";")
    )
}

.ps.lps.diagnostics <- function(frames, y, beta, sync.rows,
                                ridge.values = numeric(0)) {
    n <- length(frames)
    local.gcv <- numeric(n)
    local.rss <- numeric(n)
    df.ratio <- numeric(n)
    for (ii in seq_len(n)) {
        fr <- frames[[ii]]
        idx <- fr$offset + seq_len(fr$q)
        pred <- as.numeric(fr$design %*% beta[idx])
        resid <- y[fr$index] - pred
        rss <- sum(fr$weights * resid^2)
        denom <- 1 - fr$rank / length(fr$index)
        local.rss[[ii]] <- rss
        df.ratio[[ii]] <- fr$rank / length(fr$index)
        local.gcv[[ii]] <- if (denom > 0) {
            (rss / length(fr$index)) / denom^2
        } else {
            Inf
        }
    }
    sync.energy <- 0
    sync.weight <- 0
    if (length(sync.rows)) {
        for (sr in sync.rows) {
            fi <- frames[[sr$i]]
            fj <- frames[[sr$j]]
            bi <- beta[fi$offset + seq_len(fi$q)]
            bj <- beta[fj$offset + seq_len(fj$q)]
            pred.i <- as.numeric(fi$design[sr$row.i, , drop = FALSE] %*% bi)
            pred.j <- as.numeric(fj$design[sr$row.j, , drop = FALSE] %*% bj)
            diff <- pred.i - pred.j
            sync.energy <- sync.energy + 0.5 * sum(sr$omega * diff^2)
            sync.weight <- sync.weight + sum(sr$omega)
        }
    }
    list(
        total.local.gcv.ps = sum(local.gcv[is.finite(local.gcv)]),
        mean.local.gcv.ps = mean(local.gcv[is.finite(local.gcv)]),
        local.gcv.ps = local.gcv,
        local.rss.ps = local.rss,
        local.df.ratio = df.ratio,
        sync.energy = sync.energy,
        mean.sync.squared.disagreement = if (sync.weight > 0) {
            2 * sync.energy / sync.weight
        } else {
            0
        },
        mean.sync.disagreement = if (sync.weight > 0) {
            2 * sync.energy / sync.weight
        } else {
            0
        },
        ridge.min = if (length(ridge.values)) min(ridge.values) else 0,
        ridge.median = if (length(ridge.values)) stats::median(ridge.values) else 0,
        ridge.max = if (length(ridge.values)) max(ridge.values) else 0
    )
}
