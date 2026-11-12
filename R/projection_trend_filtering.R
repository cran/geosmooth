#' Projection Pursuit with Order-2 Trend Filtering
#'
#' Fit a sum of piecewise-quadratic curves on fixed or learned linear
#' projections. Gaussian models use squared loss; binomial models add the
#' curves on the logit scale.
#'
#' @param X Finite numeric predictor matrix or data frame.
#' @param y Finite numeric response. For `family="binomial"`, values must lie
#'   in `[0,1]`; proportions can be accompanied by trial weights.
#' @param lambda Nonnegative scalar or `Inf`, on the mean-loss scale.
#' @param family `"gaussian"` or `"binomial"`.
#' @param directions Matrix with one direction per row and one column per
#'   original predictor. Directions refer to the preprocessed coordinates.
#'   Entries for dropped constant columns must be zero. In learned mode these
#'   are additional candidates; in fixed mode these are the fitted directions.
#' @param direction.search `"fixed"` fits supplied directions (coordinate
#'   directions when omitted); `"learned"` selects and refines directions
#'   with pattern search. `"gn"` uses multistart damped Gauss--Newton
#'   proposals; `"pls"` uses one covariance start per component with the
#'   same Gauss--Newton refinement. The latter two modes construct their
#'   own starts and do not accept supplied directions.
#' @param nterms Number of supplied directions in fixed mode; maximum number
#'   of components in learned mode. Defaults to the supplied row count, all
#'   nonconstant coordinates in fixed mode, or two in learned mode. Zero fits
#'   an intercept only.
#' @param weights Positive finite observation weights, default one.
#' @param offset Known additive predictor, default zero.
#' @param curve.ridge Nonnegative component ridge coefficient on the mean-loss
#'   scale; default `0.001`. Positive ridge stabilizes component decomposition.
#' @param preprocessing `"scale"` centers and scales each nonconstant column;
#'   `"isotropic"` centers and uses a common scale, preserving rotations and
#'   distance ratios; `"none"` leaves retained columns unchanged. All
#'   transformations use training weights only. The default is `"isotropic"`
#'   for `"gn"` and `"pls"`, and `"scale"` otherwise.
#' @param seed Integer random seed. Local randomness leaves the caller's RNG
#'   state unchanged.
#' @param control Named list. `max.sweeps=200`, `tolerance=1e-7`,
#'   `gap.tolerance=1e-8`, and `objective.tolerance=1e-10` govern joint fitting.
#'   `inner=list()` passes controls to the univariate solver.
#'   `joint.solver="backfit"` retains alternating component updates;
#'   `joint.solver="dense"` attempts a simultaneous dense solve when there
#'   are exactly two components, with guarded fallback to backfitting.
#'   `joint.max.workspace.bytes=1024^3` bounds the estimated dense workspace.
#'   The joint certificate and objective safeguard apply to both solvers.
#'   Search controls are `starts=2`, `sparse.directions=8`,
#'   `dense.directions=8`, `refinement.passes=2`, `refinement.proposals=16`,
#'   `steps=c(0.25,0.1,0.04,0.01)`, and `improvement.tolerance=1e-9`.
#'   Zero refinement passes select from a fixed candidate bank.
#'   The `"gn"` and `"pls"` modes instead use `gn.max.steps=8` accepted
#'   proposals and `gn.max.fits=40` native fits per start, with
#'   `gn.shared.fits=400` shared fits for previously selected components.
#'   `stream.key=character()` can identify a dataset or representation;
#'   local start streams also depend on the seed and component number.
#' @param backend Univariate solver: `"auto"`, `"dense"`, or `"banded"`.
#' @param warn Warn when joint fitting or a terminal search failure makes the
#'   returned fit ineligible. Failed candidates are retained in `search.history`.
#'
#' @details For training weight mass W, the objective is
#'   \deqn{W^{-1}\sum_i w_i\ell(y_i,\eta_i) +
#'   \lambda\sum_j\|D_j^{(3)}h_j\|_1 +
#'   \frac{\rho}{2W}\sum_{ij}w_i h_j(u_j^T z_i)^2.}
#'   Each component has training weighted mean zero; directions have unit
#'   norm and a canonical sign. Projection ranges and the falling-factorial
#'   basis follow \code{\link{fit.gaussian.trend.filtering}} and
#'   \code{\link{fit.logistic.trend.filtering}}. Native lambda and curve ridge
#'   are each multiplied by W exactly once. Gaussian loss includes a factor
#'   of one half. `lambda=Inf` restricts components to quadratics; it does not
#'   select a constant model.
#'
#'   Fixed-direction fitting is convex. Backfitting updates a curve with the
#'   other components supplied as an offset and checks joint stationarity,
#'   dual feasibility, complementarity and the original objective after each
#'   sweep. The `"learned"` mode uses greedy additions and derivative-free
#'   tangent proposals, followed by joint refitting. The search is nonconvex;
#'   its recorded stopping status is not a global optimality certificate.
#'   Gauss--Newton modes use one-sided derivatives of the falling-factorial
#'   curve, solve in the tangent space of the unit direction, and accept only
#'   certified joint-objective improvements. Binomial curvature is floored
#'   at `1e-6` without changing the likelihood score. All curve fits use the
#'   original responses and offsets, including for binomial outcomes.
#'   Use `preprocessing="isotropic"` to preserve the geometry of supplied
#'   score coordinates. PCA and variance-threshold tuning are not performed
#'   automatically by this function.
#'   In `"gn"` mode the ordered starts are coordinate axes, weighted ridge,
#'   weighted covariance, sliced inverse regression, and twice as many random
#'   directions as retained coordinates. Sign-equivalent duplicates retain
#'   their first origin. SIR uses two outcome classes for the first Bernoulli
#'   component with unit trial weights; otherwise it slices the centered
#'   working response by weight mass, keeping exact response ties together.
#'   Both supervised seeds and fitted curves eliminate a weighted intercept.
#'   A start must obtain a completed joint certificate before it can be
#'   selected. A per-start limit can retain that certificate; exhaustion of
#'   the candidate or shared allowance makes the candidate ineligible.
#'   With q retained coordinates and requested K components, the default
#'   native-fit caps are `40*K*(3*q+3)+400` for `"gn"` and `40*K+400`
#'   for `"pls"`. An intercept-only fit makes no native curve-fitting calls.
#'   Directions with fewer than four distinct projections are ineligible.
#'   Exact ties are retained. Distinct projected locations that collapse
#'   during floating-point normalization produce an explicit failed fit;
#'   locations are never rounded, merged or jittered to repair a candidate.
#'   One-class binomial responses have no finite
#'   unpenalized intercept and produce an error.
#'
#' @return A `projection.trend.filtering.fit` containing centered `terms`,
#'   `intercept`, canonical `directions` in retained predictor coordinates,
#'   component fits, saved `preprocessing`, `checks`, `converged`, `status`,
#'   `objective.history`, and a candidate-level `search.history`.
#'   `checks$converged` reports joint numerical optimality; `converged` also
#'   excludes terminal search failures. Search termination is recorded
#'   separately in `search.status`. Ordinary predictions are NA for failed fits.
#'   Gauss--Newton modes additionally return `numerical.certified`,
#'   `cv.eligible`, `counters`, `curve.diagnostics`, and `search.banks`.
#'   An interrupted search may preserve a numerical certificate while being
#'   ineligible for CV. Requested and effective component counts can differ.
#'   In `counters`, `candidate.fits` counts returned attempts, while
#'   `completed.candidate.fits` counts eligible, certified results.
#' @examples
#' x <- seq(-1, 1, length.out = 24)
#' X <- cbind(x, cos(3*x))
#' fit <- fit.projection.trend.filtering(X, sin(2*x), lambda = 0.002,
#'     directions = matrix(c(1, 0), 1), direction.search = "fixed")
#' predict(fit, X[1:3, , drop = FALSE])
#' @export
fit.projection.trend.filtering <- function(
    X, y, lambda, family = c("gaussian", "binomial"), directions = NULL,
    direction.search = c("learned", "fixed", "gn", "pls"), nterms = NULL,
    weights = NULL, offset = NULL, curve.ridge = 0.001,
    preprocessing = c("scale", "isotropic", "none"), seed = 1L,
    control = list(), backend = c("auto", "dense", "banded"), warn = TRUE) {
    call <- match.call()
    family <- match.arg(family); direction.search <- match.arg(direction.search)
    if (missing(preprocessing) && direction.search %in% c("gn", "pls")) preprocessing <- "isotropic"
    preprocessing <- match.arg(preprocessing); backend <- match.arg(backend)
    control <- .ptf.control(control)
    if (direction.search %in% c("gn", "pls") && !is.null(directions))
        stop("gn and pls construct their own starts; omit directions.", call. = FALSE)
    .ptf.scalar(lambda, "lambda", infinite = TRUE)
    .ptf.scalar(curve.ridge, "curve.ridge")
    .ptf.scalar(seed, "seed", integer = TRUE)
    .ptf.logical(warn, "warn")
    problem <- .ptf.problem(X, y, weights, offset, family, preprocessing,
                            lambda, curve.ridge, control, backend)
    p <- ncol(problem$Z)
    if (!is.null(directions)) {
        directions <- .ptf.directions(directions, problem$preprocessing)
    }
    if (is.null(nterms)) nterms <- if (!is.null(directions)) nrow(directions) else
        if (direction.search == "fixed") p else 2L
    .ptf.scalar(nterms, "nterms", integer = TRUE)
    if (nterms > 0 && !p) stop("No nonconstant predictors.", call. = FALSE)
    if (direction.search == "fixed") {
        if (is.null(directions)) directions <- if (nterms == 0) matrix(0, 0, p) else diag(p)
        if (nrow(directions) != nterms)
            stop("nterms must equal the number of fixed directions.", call. = FALSE)
        for (j in seq_len(nterms)) .ptf.projection(problem, directions[j, ])
        state <- .ptf.backfit(problem, .ptf.initial(problem, directions))
        search <- list(history = data.frame(), status = "fixed_directions",
                       starts = data.frame(), solver.calls = state$solver.calls)
    } else {
        searched <- if (direction.search %in% c("gn", "pls"))
            .ptf.gn.search(problem, nterms, seed, direction.search) else
            .ptf.with.seed(seed, .ptf.search(problem, directions, nterms))
        state <- searched$state; search <- searched$search
    }
    out <- c(list(call = call, family = family, order = 2L, X = problem$X,
        y = y, weights = problem$weights, offset = problem$offset,
        offset.used = !is.null(offset), lambda = lambda, curve.ridge = curve.ridge,
        weight.mass = sum(problem$weights), preprocessing = problem$preprocessing,
        direction.search = direction.search, requested.nterms = nterms, seed = seed,
        control = control, backend = backend), state)
    out$nterms <- ncol(state$terms)
    out$search.history <- search$history; out$search.status <- search$status
    out$search.starts <- search$starts; out$solver.calls <- search$solver.calls
    if (!is.null(search$banks)) out$search.banks <- search$banks
    class(out) <- "projection.trend.filtering.fit"
    if (warn && !out$converged)
        warning("Projection fit failed convergence checks: ", out$status,
                ". Inspect checks and search.history.", call. = FALSE)
    out
}

.ptf.scalar <- function(x, name, integer = FALSE, positive = FALSE, infinite = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        (!infinite && !is.finite(x)) || x < 0 || (positive && x == 0) ||
        (integer && (!is.finite(x) || x != floor(x) || x > .Machine$integer.max)))
        stop(name, " must be a ", if (positive) "positive" else "nonnegative",
             if (integer) " integer." else if (infinite) " scalar or Inf." else " finite scalar.",
             call. = FALSE)
    invisible(x)
}

.ptf.logical <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
        stop(name, " must be TRUE or FALSE.", call. = FALSE)
}

.ptf.control <- function(control) {
    defaults <- list(max.sweeps = 200L, tolerance = 1e-7, gap.tolerance = 1e-8,
        objective.tolerance = 1e-10, inner = list(), starts = 2L,
        sparse.directions = 8L, dense.directions = 8L,
        refinement.passes = 2L, refinement.proposals = 16L,
        steps = c(.25, .1, .04, .01), improvement.tolerance = 1e-9,
        gn.max.steps = 8L, gn.max.fits = 40L, gn.shared.fits = 400L,
        stream.key = character(), joint.solver = "backfit",
        joint.max.workspace.bytes = 1024^3)
    if (!is.list(control) || (length(control) && (is.null(names(control)) ||
        any(!names(control) %in% names(defaults)) || anyDuplicated(names(control)))))
        stop("Unknown or unnamed control fields.", call. = FALSE)
    out <- utils::modifyList(defaults, control)
    for (nm in c("max.sweeps", "starts", "sparse.directions", "dense.directions",
                 "refinement.passes", "refinement.proposals",
                 "gn.max.steps", "gn.max.fits", "gn.shared.fits"))
        .ptf.scalar(out[[nm]], nm, integer = TRUE, positive = nm %in% c("max.sweeps", "starts"))
    for (nm in c("tolerance", "gap.tolerance", "objective.tolerance", "improvement.tolerance"))
        .ptf.scalar(out[[nm]], nm, positive = TRUE)
    if (!is.numeric(out$steps) || !length(out$steps) || any(!is.finite(out$steps)) || any(out$steps <= 0))
        stop("steps must be a positive finite numeric vector.", call. = FALSE)
    if (!is.list(out$inner)) stop("inner must be a list.", call. = FALSE)
    if (!is.character(out$stream.key) || anyNA(out$stream.key))
        stop("stream.key must be a nonmissing character vector.", call. = FALSE)
    if (!is.character(out$joint.solver) || length(out$joint.solver) != 1L ||
        is.na(out$joint.solver) || !out$joint.solver %in% c("backfit", "dense"))
        stop("joint.solver must be backfit or dense.", call. = FALSE)
    .ptf.scalar(out$joint.max.workspace.bytes, "joint.max.workspace.bytes", positive = TRUE)
    # Validate even when an intercept-only fit never calls the native solver.
    allowed <- c("tolerance", "gap.tolerance", "max.iterations", "max.workspace.bytes")
    if (length(out$inner) && (is.null(names(out$inner)) ||
        any(!names(out$inner) %in% allowed) || anyDuplicated(names(out$inner))))
        stop("Unknown or unnamed inner control fields.", call. = FALSE)
    for (nm in names(out$inner)) .ptf.scalar(out$inner[[nm]], nm,
        integer = nm == "max.iterations", positive = TRUE)
    out
}

.ptf.matrix <- function(X, name = "X") {
    if (is.data.frame(X) && all(vapply(X, is.numeric, logical(1)))) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)))
        stop(name, " must be a finite numeric matrix or numeric data frame.", call. = FALSE)
    if (!is.null(colnames(X)) && (anyNA(colnames(X)) || anyDuplicated(colnames(X))))
        stop(name, " must have unique, nonmissing column names.", call. = FALSE)
    X
}

.ptf.problem <- function(X, y, weights, offset, family, preprocessing, lambda, rho, control, backend) {
    X <- .ptf.matrix(X); n <- nrow(X)
    if (n < 4L) stop("X needs at least four rows.", call. = FALSE)
    if (!is.numeric(y) || !is.null(dim(y)) || length(y) != n || any(!is.finite(y)) ||
        (family == "binomial" && any(y < 0 | y > 1)))
        stop("y must be finite, match X rows, and lie in [0,1] for binomial fits.", call. = FALSE)
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(offset)) offset <- numeric(n)
    for (nm in c("weights", "offset")) {
        v <- get(nm)
        if (!is.numeric(v) || !is.null(dim(v)) || length(v) != n || any(!is.finite(v)))
            stop(nm, " must be finite with one value per row.", call. = FALSE)
    }
    if (any(weights <= 0) || !is.finite(sum(weights)))
        stop("weights must be positive with finite sum.", call. = FALSE)
    if (family == "binomial" && (all(y == 0) || all(y == 1)))
        stop("No finite unpenalized intercept: y contains only one outcome.", call. = FALSE)
    w <- weights/sum(weights)
    retained <- which(vapply(seq_len(ncol(X)), function(j) length(unique(X[, j])) > 1L, logical(1)))
    Z <- X[, retained, drop = FALSE]
    center <- if (preprocessing == "none") numeric(ncol(Z)) else colSums(Z*w)
    Z <- sweep(Z, 2, center)
    scale <- if (preprocessing == "none") rep(1, ncol(Z)) else sqrt(colSums(Z^2*w))
    if (preprocessing == "isotropic" && length(scale)) scale[] <- sqrt(mean(scale^2))
    if (any(!is.finite(center)) || any(!is.finite(scale)) || any(scale <= 0))
        stop("Nonfinite or degenerate preprocessing; rescale X.", call. = FALSE)
    Z <- sweep(Z, 2, scale, "/")
    if (any(!is.finite(Z))) stop("Nonfinite transformed predictors.", call. = FALSE)
    list(X = X, Z = Z, y = y, weights = weights, w = w, offset = offset,
         family = family, lambda = lambda, rho = rho, control = control, backend = backend,
         preprocessing = list(method = preprocessing, columns = colnames(X),
             ncol = ncol(X), retained = retained, center = center, scale = scale))
}

.ptf.canonical <- function(u) {
    if (!length(u)) stop("Direction must be finite and nonzero.", call. = FALSE)
    size <- max(abs(u))
    if (!is.finite(size) || size == 0)
        stop("Direction must be finite and nonzero.", call. = FALSE)
    u <- u/size; u <- u/sqrt(sum(u^2))
    if (u[which.max(abs(u))] < 0) u <- -u
    u
}

.ptf.directions <- function(U, prep) {
    U <- .ptf.matrix(U, "directions")
    if (ncol(U) != prep$ncol) stop("directions must match original predictor columns.", call. = FALSE)
    if (any(U[, setdiff(seq_len(prep$ncol), prep$retained), drop = FALSE] != 0))
        stop("Directions must be zero on constant predictor columns.", call. = FALSE)
    U <- U[, prep$retained, drop = FALSE]
    for (j in seq_len(nrow(U))) U[j, ] <- .ptf.canonical(U[j, ])
    if (nrow(.ptf.unique.directions(U)) != nrow(U))
        stop("Duplicate directions are not allowed.", call. = FALSE)
    U
}

.ptf.unique.directions <- function(U) {
    keep <- integer()
    for (j in seq_len(nrow(U))) if (!length(keep) ||
        all(vapply(keep, function(k) max(abs(U[k, ]-U[j, ])) > 1e-12, logical(1))))
        keep <- c(keep, j)
    U[keep, , drop = FALSE]
}

.ptf.projection <- function(problem, u) {
    t <- .ptf.project(problem$X, u, problem$preprocessing)
    if (any(!is.finite(t)) || length(unique(t)) < 4L || !is.finite(diff(range(t))))
        stop("A direction needs at least four distinct finite projected values.", call. = FALSE)
    t
}

.ptf.project <- function(X, u, prep) {
    # Explicit products avoid fused multiply-add changing exact cancellations
    # in symmetric designs. Center after summing, so centering separate columns
    # cannot split exact projection ties. Do not round, merge or jitter values.
    X <- X[, prep$retained, drop = FALSE]
    if (length(prep$scale) && all(prep$scale == prep$scale[1L]))
        return((rowSums(sweep(X, 2L, u, "*"))-sum(prep$center*u))/prep$scale[1L])
    coefficient <- u/prep$scale
    rowSums(sweep(X, 2L, coefficient, "*"))-sum(prep$center*coefficient)
}

.ptf.intercept <- function(problem, offset) {
    if (problem$family == "gaussian") return(sum(problem$w*(problem$y-offset)))
    target <- sum(problem$w*problem$y)
    b <- stats::qlogis(target)
    score <- function(a) sum(problem$w*(stats::plogis(offset+a)-problem$y))
    stats::uniroot(score, c(b-max(offset)-1, b-min(offset)+1), tol = 1e-12)$root
}

.ptf.initial <- function(problem, U) {
    list(intercept = .ptf.intercept(problem, problem$offset), directions = U,
         terms = matrix(0, nrow(problem$Z), nrow(U)), components = vector("list", nrow(U)),
         objective.history = numeric(), solver.calls = 0L)
}

.ptf.loss <- function(y, eta, family) {
    if (family == "gaussian") (y-eta)^2/2 else
        (1-y)*pmax(eta, 0) + y*pmax(-eta, 0) + log1p(exp(-abs(eta)))
}

.ptf.diagnose <- function(problem, state) {
    .ptf.observed(problem$observer, "certificate", list(),
        .ptf.diagnose.core(problem, state), identity)
}

.ptf.diagnose.core <- function(problem, state) {
    eta <- problem$offset + state$intercept + rowSums(state$terms)
    residual <- problem$w * (if (problem$family == "gaussian") eta-problem$y else
        stats::plogis(eta)-problem$y)
    stationarity <- abs(sum(residual)); roughness <- complementarity <- dual <- 0
    ready <- TRUE
    for (j in seq_along(state$components)) {
        component <- state$components[[j]]
        if (is.null(component)) { ready <- FALSE; next }
        fit <- component$fit; d <- fit$design; f <- fit$fits[[1L]]
        beta <- f$coefficients[-1L]; v <- f$dual[-1L]
        penalty <- c(0, 0, 1/d$scales[-c(1L, 2L)])
        roughness <- roughness + sum(penalty*abs(beta))
        dual <- max(dual, abs(v)-1)
        if (is.finite(problem$lambda)) complementarity <- complementarity +
            problem$lambda*sum(penalty*(abs(beta)-v*beta))
        u <- (fit$x-d$lo)/d$span
        grad <- residual + problem$rho*problem$w*state$terms[, j]
        # Scan columns without materializing an n-by-knots matrix. This also
        # checks stale component duals against the current joint residual.
        cols <- if (is.infinite(problem$lambda)) 1:2 else seq_along(beta)
        for (k in cols) {
            h <- if (k == 1L) u else if (k == 2L) u^2 else
                .5*(u-d$knots[k-1L])*(u-d$knots[k])*(u > d$knots[k])
            g <- sum((h-d$center[k])*grad)/d$scales[k]
            if (is.finite(problem$lambda)) g <- g + problem$lambda*penalty[k]*v[k]
            stationarity <- max(stationarity, abs(g))
        }
        ready <- ready && isTRUE(f$converged)
    }
    loss <- sum(problem$w*.ptf.loss(problem$y, eta, problem$family))
    ridge <- problem$rho/2*sum(state$terms^2*problem$w)
    l1 <- if (is.infinite(problem$lambda)) {
        if (roughness == 0) 0 else Inf
    } else problem$lambda*roughness
    objective <- loss + ridge + l1
    finite <- all(is.finite(c(eta, objective, stationarity, complementarity, dual)))
    list(objective = objective, loss = loss, l1 = l1, ridge = ridge,
         stationarity = if (ready) stationarity else Inf,
         complementarity = complementarity, dual.feasibility = max(0, dual),
         converged = ready && finite && stationarity <= problem$control$tolerance &&
             abs(complementarity) <= problem$control$gap.tolerance && dual <= 1e-12)
}

.ptf.update <- function(problem, state, j, u = state$directions[j, ]) {
    .ptf.observed(problem$observer, "component", list(term = j, direction = u),
        .ptf.update.core(problem, state, j, u),
        function(value) list(status = "complete", diagnostics = value$components[[j]]$fit$diagnostics))
}

.ptf.update.core <- function(problem, state, j, u = state$directions[j, ]) {
    x <- .ptf.projection(problem, u)
    offset <- problem$offset + state$intercept + rowSums(state$terms)-state$terms[, j]
    fit <- .fit.univariate.trend.filtering(x, problem$y,
        problem$lambda*sum(problem$weights), weights = problem$weights, offset = offset,
        curve.ridge = problem$rho*sum(problem$weights), family = problem$family,
        control = problem$control$inner, backend = problem$backend, warn = FALSE,
        before.fit = problem$before.fit, observer = problem$observer)
    if (!isTRUE(fit$diagnostics$converged))
        stop("component_solver_", fit$diagnostics$status, call. = FALSE)
    curve <- as.numeric(stats::predict(fit, x, offset = numeric(length(x)), type = "link"))
    center <- sum(problem$w*curve)
    state$intercept <- state$intercept + center
    state$terms[, j] <- curve-center
    state$directions[j, ] <- u
    state$components[[j]] <- list(fit = fit, center = center)
    state
}

.ptf.backfit <- function(problem, state) {
    .ptf.observed(problem$observer, "joint", list(checks = state$checks),
        if (identical(problem$control$joint.solver, "dense") && length(state$components) == 2L)
            .ptf.backfit.dense(problem, state) else .ptf.backfit.core(problem, state), function(value)
            list(status = value$status, certified = value$converged,
                 checks = value$checks, completed.sweeps = value$sweeps,
                 attempted.updates = value$solver.calls))
}

.ptf.backfit.core <- function(problem, state) {
    check <- .ptf.diagnose(problem, state)
    history <- check$objective; calls <- 0L; status <- "iteration_limit"; sweeps <- 0L
    if (check$converged) status <- "converged"
    else for (s in seq_len(problem$control$max.sweeps)) {
        if (!is.null(problem$counter)) problem$counter$sweeps <- problem$counter$sweeps+1
        sweep.id <- .ptf.observe.begin(problem$observer, "sweep", list(index = s))
        previous <- state
        updated <- 0L
        failure <- NULL
        for (j in seq_along(state$components)) {
            calls <- calls + 1L
            candidate <- tryCatch(.ptf.update(problem, state, j), error = .ptf.numerical.error)
            if (inherits(candidate, "error")) { failure <- conditionMessage(candidate); break }
            state <- candidate
            updated <- updated+1L
        }
        if (!is.null(failure)) {
            state <- previous; status <- failure
            .ptf.observe.end(problem$observer, sweep.id,
                list(status = "failed", reason = failure, completed.updates = updated,
                     rollback = TRUE))
            break
        }
        check.new <- .ptf.diagnose(problem, state)
        if (!is.finite(check.new$objective) || check.new$objective > utils::tail(history, 1L) +
            problem$control$objective.tolerance*max(1, abs(utils::tail(history, 1L)))) {
            state <- previous; status <- "objective_increase"
            .ptf.observe.end(problem$observer, sweep.id,
                list(status = status, completed.updates = updated, checks = check.new, rollback = TRUE))
            break
        }
        sweeps <- s; check <- check.new; history <- c(history, check$objective)
        .ptf.observe.end(problem$observer, sweep.id,
            list(status = "complete", completed.updates = updated, checks = check, rollback = FALSE))
        if (check$converged) { status <- "converged"; break }
    }
    state$checks <- .ptf.diagnose(problem, state)
    state$converged <- identical(status, "converged") && state$checks$converged
    state$status <- status; state$sweeps <- sweeps
    state$objective.history <- history; state$solver.calls <- calls
    state
}

#' Predict from Projection Trend Filtering
#' @param object A `projection.trend.filtering.fit` or a completed
#'   `projection.trend.filtering.cv`.
#' @param newdata Predictor matrix or numeric data frame; defaults to training
#'   rows. Named columns are aligned to training names when both are supplied.
#' @param offset Known offset, required for new rows if used in fitting.
#' @param type `"response"`, `"link"`, or `"terms"`. Terms exclude the common
#'   intercept and offset, which are returned as attributes.
#' @param allow.unconverged Explicitly allow predictions from a failed fit.
#' @param ... Reserved for generic compatibility.
#' @return A vector, or a matrix for `type="terms"`. The `outside.range`
#'   attribute is a row-by-component matrix identifying constant-boundary
#'   extension of each univariate curve.
#' @examples
#' X <- matrix(seq(0, 1, length.out = 12), ncol = 1)
#' y <- sin(3*X[,1])
#' fit <- fit.projection.trend.filtering(X, y, lambda=0.001,
#'     direction.search="fixed")
#' predict(fit, matrix(c(0.2, 0.8), ncol=1))
#' head(fitted(fit))
#' head(residuals(fit))
#' print(fit)
#' summary(fit)
#' cv <- cv.projection.trend.filtering(X, y, foldid=rep(1:3, 4),
#'     lambda.grid=Inf, nterms.grid=1, direction.search="fixed")
#' predict(cv, X[1:2, , drop=FALSE])
#' @export
predict.projection.trend.filtering.fit <- function(object, newdata = NULL,
    offset = NULL, type = c("response", "link", "terms"), allow.unconverged = FALSE, ...) {
    type <- match.arg(type); .ptf.logical(allow.unconverged, "allow.unconverged")
    training <- is.null(newdata)
    if (training) { newdata <- object$X; if (is.null(offset)) offset <- object$offset }
    X <- .ptf.matrix(newdata, "newdata"); d <- object$preprocessing
    if (ncol(X) != d$ncol) stop("newdata must match training columns.", call. = FALSE)
    if (!is.null(d$columns) && !is.null(colnames(X))) {
        if (!setequal(colnames(X), d$columns)) stop("newdata names must match training columns.", call. = FALSE)
        X <- X[, match(d$columns, colnames(X)), drop = FALSE]
    }
    n <- nrow(X)
    if (is.null(offset)) {
        if (object$offset.used && !training) stop("Supply offset for newdata.", call. = FALSE)
        offset <- numeric(n)
    }
    if (!is.numeric(offset) || !is.null(dim(offset)) || length(offset) != n || any(!is.finite(offset)))
        stop("offset must be finite with one value per newdata row.", call. = FALSE)
    Z <- sweep(sweep(X[, d$retained, drop = FALSE], 2, d$center), 2, d$scale, "/")
    if (any(!is.finite(Z))) stop("Nonfinite transformed newdata.", call. = FALSE)
    terms <- matrix(0, n, object$nterms); outside <- matrix(FALSE, n, object$nterms)
    for (j in seq_len(object$nterms)) {
        comp <- object$components[[j]]
        if (is.null(comp)) { terms[, j] <- NA_real_; outside[, j] <- NA; next }
        t <- .ptf.project(X, object$directions[j, ], d)
        v <- stats::predict(comp$fit, t, offset = numeric(n), type = "link",
                            allow.unconverged = allow.unconverged)
        terms[, j] <- as.numeric(v)-comp$center
        outside[, j] <- attr(v, "outside.range")
    }
    colnames(terms) <- colnames(outside) <- if (object$nterms) paste0("term", seq_len(object$nterms)) else character()
    eta <- offset + object$intercept + rowSums(terms)
    ans <- switch(type, terms = terms, link = eta,
        response = if (object$family == "binomial") stats::plogis(eta) else eta)
    if (!.ptf.eligible(object) && !allow.unconverged) ans[] <- NA_real_
    attr(ans, "outside.range") <- outside
    if (type == "terms") { attr(ans, "intercept") <- object$intercept; attr(ans, "offset") <- offset }
    ans
}

#' @rdname predict.projection.trend.filtering.fit
#' @export
fitted.projection.trend.filtering.fit <- function(object, ...) stats::predict(object, ...)
#' @rdname predict.projection.trend.filtering.fit
#' @export
residuals.projection.trend.filtering.fit <- function(object, ...) object$y-stats::fitted(object, ...)
#' @rdname predict.projection.trend.filtering.fit
#' @param x A fitted projection trend-filtering model.
#' @method print projection.trend.filtering.fit
#' @export
print.projection.trend.filtering.fit <- function(x, ...) {
    cat("Order-2 projection trend filtering (", x$family, ")\n", sep = "")
    cat("Components:", x$nterms, " Fit:", x$status, " Search:", x$search.status, "\n")
    cat("Mean-loss objective:", format(x$checks$objective), "\n")
    invisible(x)
}
#' @rdname predict.projection.trend.filtering.fit
#' @method summary projection.trend.filtering.fit
#' @export
summary.projection.trend.filtering.fit <- function(object, ...) {
    structure(list(family = object$family, nterms = object$nterms,
        lambda = object$lambda, curve.ridge = object$curve.ridge, checks = object$checks,
        status = object$status, search.status = object$search.status,
        solver.calls = object$solver.calls, directions = object$directions), class = "list")
}
