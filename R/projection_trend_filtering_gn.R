.ptf.native.budget <- function(K, q, method, per.start = 40, shared = 400, total = NULL) {
    e <- new.env(parent = emptyenv())
    e$per.start <- per.start; e$shared.limit <- shared
    starts <- if (method == "pls") 1 else 3*q+3
    e$total.limit <- if (is.null(total)) if (K == 0) 0 else per.start*K*starts+shared else total
    e$native <- e$start.used <- e$shared.used <- e$sweeps <- e$evaluations <- 0
    e$mode <- "start"; e$blocked <- NULL; e$global.exhausted <- FALSE
    e
}

.ptf.budget.enter <- function(budget, mode = "start") {
    budget$mode <- mode; budget$blocked <- NULL
    if (mode == "start") budget$start.used <- 0
    invisible(budget)
}

.ptf.reserve.native <- function(budget) {
    scope <- if (budget$native >= budget$total.limit) "candidate" else
        if (budget$mode == "shared" && budget$shared.used >= budget$shared.limit) "shared" else
        if (budget$mode == "start" && budget$start.used >= budget$per.start) "start" else NULL
    if (!is.null(scope)) {
        budget$blocked <- scope
        if (scope != "start") budget$global.exhausted <- TRUE
        stop(structure(list(message = "budget_exhausted", call = NULL, scope = scope),
                       class = c("projection_fit_budget", "error", "condition")))
    }
    budget$native <- budget$native+1
    if (budget$mode == "shared") budget$shared.used <- budget$shared.used+1 else
        budget$start.used <- budget$start.used+1
    invisible(NULL)
}

.ptf.addition.outcome <- function(completed = 0L, admissible = 0L, failures = 0L,
                                   interrupted = 0L, global.exhausted = FALSE,
                                   selected.failed = FALSE, improving = FALSE) {
    status <- if (global.exhausted) "budget_exhausted" else
        if (selected.failed) "addition_refit_failed" else
        if (completed > 0) if (improving) "accepted" else "no_improving_addition" else
        if (!admissible) "no_admissible_direction" else
        if (interrupted > 0) "budget_exhausted" else "all_candidates_failed"
    list(status = status, eligible = status %in%
        c("accepted", "no_improving_addition", "no_admissible_direction"),
        failures = failures, interrupted = interrupted)
}

.ptf.eligible <- function(fit) {
    isTRUE(fit$converged) && (is.null(fit$cv.eligible) || isTRUE(fit$cv.eligible))
}

.ptf.gn.trial <- function(problem, state, term, u) {
    problem$counter$evaluations <- problem$counter$evaluations+1
    trial <- tryCatch(.ptf.update(problem, state, term, u), error = .ptf.numerical.error)
    if (inherits(trial, "error")) return(list(state = NULL, status = conditionMessage(trial)))
    joint <- .ptf.backfit(problem, trial)
    list(state = if (joint$converged) joint else NULL, status = joint$status)
}

.ptf.gn.refine <- function(problem, state, term, budget, record, start) {
    accepted <- 0L; final.mu <- NA_real_; stopping <- "step_limit"
    for (iteration in seq_len(problem$control$gn.max.steps)) {
        fit <- state$components[[term]]$fit; u <- state$directions[term, ]
        g <- tryCatch(.ptf.curve.slopes(fit), error = .ptf.numerical.error)
        if (inherits(g, "error") || any(!is.finite(g))) {
            stopping <- "nonfinite_slope"
            record(start, "refinement", term, stopping, state, u = u)
            break
        }
        flat.threshold <- 1e-10*max(1, abs(state$terms[, term]))
        flat.norm <- max(abs(g))*fit$design$span
        if (all(g == 0) || flat.norm <= flat.threshold || length(u) == 1L) {
            stopping <- if (all(g == 0)) "flat_curve" else if (flat.norm <= flat.threshold)
                "near_flat" else "rank_one"
            record(start, "refinement", term, stopping, state, u = u,
                   slope.norm = flat.norm, slope.threshold = flat.threshold)
            break
        }
        eta <- problem$offset+state$intercept+rowSums(state$terms)
        proposal <- tryCatch({
            eq <- .ptf.gn.equations(problem$Z, problem$y, eta, g, problem$w, problem$family)
            list(eq = eq, inc = .ptf.gn.increment(eq$A, eq$b, u))
        }, error = .ptf.numerical.error)
        if (inherits(proposal, "error")) {
            stopping <- "proposal_failed"
            record(start, "refinement", term, stopping, state, u = u)
            break
        }
        eq <- proposal$eq; inc <- proposal$inc
        if (inc$status != "proposal") {
            stopping <- inc$status
            record(start, "refinement", term, stopping, state, u = u,
                   score.norm = inc$tangent.score, score.threshold = inc$score.threshold)
            break
        }
        improved <- FALSE
        for (escalation in 0:2) {
            final.mu <- inc$mu*10^escalation
            delta <- .ptf.gn.increment(eq$A, eq$b, u, final.mu)$delta
            for (fraction in c(1, .5, .25, .125)) {
                proposal <- .ptf.retract(u, delta, fraction)
                other <- state$directions[-term, , drop = FALSE]
                geometry <- tryCatch({
                    if (.ptf.duplicate(proposal, other)) stop("duplicate_direction")
                    t <- .ptf.projection(problem, proposal)
                    if (any(diff((sort(unique(t))-min(t))/diff(range(t))) <= 0))
                        stop("collapsed_projection")
                    TRUE
                }, error = .ptf.numerical.error)
                if (inherits(geometry, "error")) {
                    record(start, "refinement", term, conditionMessage(geometry), u = proposal,
                           mu = final.mu, fraction = fraction)
                    next
                }
                trial <- .ptf.gn.trial(problem, state, term, proposal)
                take <- !is.null(trial$state) &&
                    .ptf.improves(problem, trial$state$checks$objective, state$checks$objective)
                record(start, "refinement", term,
                    if (take) "accepted" else if (is.null(trial$state)) trial$status else "rejected",
                    trial$state, accepted = take, u = proposal, mu = final.mu,
                    fraction = fraction, floored = eq$floored,
                    slope.norm = flat.norm, slope.threshold = flat.threshold,
                    score.norm = inc$tangent.score, score.threshold = inc$score.threshold)
                if (take) { state <- trial$state; improved <- TRUE; accepted <- accepted+1L; break }
                if (!is.null(budget$blocked)) break
            }
            if (improved || !is.null(budget$blocked)) break
        }
        if (!is.null(budget$blocked)) { stopping <- "budget_exhausted"; break }
        if (!improved) { stopping <- "rejected"; break }
    }
    list(state = state, status = stopping, accepted = accepted, mu = final.mu)
}

.ptf.improves <- function(problem, value, previous) {
    is.finite(value) && value < previous-problem$control$improvement.tolerance*max(1, abs(previous))
}

.ptf.gn.search <- function(problem, nterms, seed, method = "gn", budget = NULL) {
    p <- ncol(problem$Z)
    if (is.null(budget)) budget <- .ptf.native.budget(nterms, p, method,
        problem$control$gn.max.fits, problem$control$gn.shared.fits)
    problem$counter <- budget
    problem$before.fit <- function() .ptf.reserve.native(budget)
    state <- .ptf.backfit(problem, .ptf.initial(problem, matrix(0, 0, p)))
    path <- state$checks$objective; rows <- starts <- banks <- list()
    status <- if (nterms) "component_limit" else "intercept_only"; eligible <- TRUE
    record <- function(start, stage, term, status, fit = NULL, accepted = FALSE,
                       u = NULL, mu = NA_real_, fraction = NA_real_, floored = NA_integer_,
                       slope.norm = NA_real_, slope.threshold = NA_real_,
                       score.norm = NA_real_, score.threshold = NA_real_) {
        rows[[length(rows)+1L]] <<- data.frame(start = start, stage = stage, term = term,
            status = status, objective = if (is.null(fit)) NA_real_ else fit$checks$objective,
            accepted = accepted, native.fits = budget$native, mu = mu, fraction = fraction,
            curvature.floored = floored, slope.norm = slope.norm, slope.threshold = slope.threshold,
            score.norm = score.norm, score.threshold = score.threshold, direction = I(list(u)))
    }
    for (term in seq_len(nterms)) {
        bank <- .ptf.gn.starts(problem, state, term, seed, method, problem$control$stream.key)
        banks[[term]] <- bank
        completed <- admissible <- failures <- interrupted <- 0L
        best <- NULL
        for (s in seq_len(nrow(bank$directions))) {
            .ptf.budget.enter(budget)
            u <- bank$directions[s, ]; before <- budget$native
            geometry <- tryCatch({
                if (.ptf.duplicate(u, state$directions)) stop("duplicate_direction")
                t <- .ptf.projection(problem, u)
                if (any(diff((sort(unique(t))-min(t))/diff(range(t))) <= 0))
                    stop("collapsed_projection")
                TRUE
            }, error = .ptf.numerical.error)
            candidate <- NULL; proposals <- 0L; mu <- NA_real_
            if (inherits(geometry, "error")) result <- conditionMessage(geometry) else {
                admissible <- admissible+1L
                expanded <- state
                expanded$directions <- rbind(expanded$directions, u)
                expanded$terms <- cbind(expanded$terms, 0)
                expanded$components <- c(expanded$components, list(NULL))
                budget$evaluations <- budget$evaluations+1
                initial <- tryCatch(.ptf.update(problem, expanded, term, u), error = .ptf.numerical.error)
                if (inherits(initial, "error")) result <- conditionMessage(initial) else
                    if (all(initial$terms[, term] == 0)) {
                        completed <- completed+1L; result <- "flat_start"
                    } else {
                        joint <- .ptf.backfit(problem, initial)
                        result <- joint$status
                        if (joint$converged) {
                            completed <- completed+1L
                            refined <- .ptf.gn.refine(problem, joint, term, budget, record, s)
                            candidate <- refined$state; result <- refined$status
                            proposals <- refined$accepted; mu <- refined$mu
                        }
                    }
                if (!is.null(budget$blocked)) interrupted <- interrupted+1L else
                    if (is.null(candidate) && result != "flat_start") failures <- failures+1L
            }
            if (!is.null(candidate) && (is.null(best) ||
                .ptf.improves(problem, candidate$checks$objective, best$checks$objective))) best <- candidate
            starts[[length(starts)+1L]] <- data.frame(term = term, start = s,
                origin = bank$origins[[s]][1L], status = result, native.fits = budget$native-before,
                certified = !is.null(candidate), accepted.proposals = proposals, mu = mu,
                angle = if (is.null(candidate)) NA_real_ else
                    acos(min(1, abs(sum(u*candidate$directions[term, ]))))*180/pi,
                objective = if (is.null(candidate)) NA_real_ else candidate$checks$objective)
            record(s, "start", term, result, candidate, u = u)
            if (budget$global.exhausted) break
        }
        improving <- !is.null(best) && .ptf.improves(problem, best$checks$objective, state$checks$objective)
        decision <- .ptf.addition.outcome(completed, admissible, failures, interrupted,
                                         budget$global.exhausted, improving = improving)
        if (decision$status != "accepted") {
            status <- decision$status; eligible <- decision$eligible; break
        }
        .ptf.budget.enter(budget, "shared")
        check <- .ptf.diagnose(problem, best)
        if (!check$converged) best <- .ptf.backfit(problem, best) else best$checks <- check
        decision <- .ptf.addition.outcome(global.exhausted = budget$global.exhausted,
            selected.failed = !best$converged ||
                !.ptf.improves(problem, best$checks$objective, state$checks$objective),
            completed = 1L, improving = TRUE)
        if (!decision$eligible) { status <- decision$status; eligible <- FALSE; break }
        state <- best; path <- c(path, state$checks$objective)
        record(NA_integer_, "addition", term, "accepted", state, TRUE, state$directions[term, ])
        if (term > 1L) for (j in seq_len(term-1L)) {
            refined <- .ptf.gn.refine(problem, state, j, budget, record, NA_integer_)
            state <- refined$state; path <- c(path, state$checks$objective)
            if (budget$global.exhausted) break
        }
        if (budget$global.exhausted) { status <- "budget_exhausted"; eligible <- FALSE; break }
    }
    state$checks <- .ptf.diagnose(problem, state)
    state$directions <- unname(state$directions)
    state$numerical.certified <- state$checks$converged
    state$cv.eligible <- eligible && state$checks$converged
    state$converged <- state$cv.eligible
    state$status <- if (state$converged) "converged" else status
    state$objective.history <- path
    state$curve.diagnostics <- lapply(state$components, .ptf.curve.diagnostics)
    state$counters <- list(native.fits = budget$native, joint.sweeps = budget$sweeps,
        trial.evaluations = budget$evaluations, candidate.fits = 1L,
        completed.candidate.fits = as.integer(state$cv.eligible),
        candidate.cap = budget$total.limit, shared.used = budget$shared.used,
        shared.cap = budget$shared.limit)
    list(state = state, search = list(history = if (length(rows)) do.call(rbind, rows) else data.frame(),
        starts = if (length(starts)) do.call(rbind, starts) else data.frame(), status = status,
        solver.calls = budget$native, banks = banks))
}

.ptf.weighted.anchor <- function(problem, bank, before.fit = NULL) {
    values <- rep(NA_real_, nrow(bank)); rows <- vector("list", nrow(bank)); calls <- 0L
    mass <- sum(problem$weights)
    enter <- function() { if (!is.null(before.fit)) before.fit(); calls <<- calls+1L }
    for (j in seq_len(nrow(bank))) {
        result <- tryCatch({
            x <- .ptf.projection(problem, bank[j, ])
            fit <- function(lambda) .fit.univariate.trend.filtering(x, problem$y,
                lambda*mass, weights = problem$weights, offset = problem$offset,
                curve.ridge = problem$rho*mass, family = problem$family,
                backend = problem$backend, control = problem$control$inner,
                warn = FALSE, before.fit = enter)
            null <- fit(Inf)
            if (!null$diagnostics$converged) stop("quadratic_anchor_failed")
            eta <- as.numeric(predict(null, type = "link")); h <- eta-problem$offset
            h <- h-sum(problem$w*h)
            grad <- problem$w*((if (problem$family == "gaussian") eta else stats::plogis(eta))-
                problem$y+problem$rho*h)
            d <- null$design; H <- .logistic.tf.basis((x-d$lo)/d$span, d$knots)
            g <- as.vector(crossprod(sweep(H, 2, d$center), grad))
            anchor <- max(1e-10, 1.01*max(abs(g[-c(1L, 2L)])))
            ok <- FALSE
            for (extension in 0:2) {
                finite <- fit(anchor)
                ok <- finite$diagnostics$converged &&
                    max(abs(as.numeric(predict(finite, type = "link"))-eta)) < 1e-6
                if (ok) break
                anchor <- anchor*10
            }
            if (!ok) stop("finite_anchor_failed")
            list(value = anchor, status = "certified")
        }, error = function(e) {
            if (inherits(e, "projection_fit_budget")) stop(e)
            .ptf.numerical.error(e)
            list(value = NA_real_, status = conditionMessage(e))
        })
        values[j] <- result$value
        rows[[j]] <- data.frame(direction = j, value = result$value, status = result$status)
    }
    if (!any(is.finite(values))) stop("No certified projection penalty anchor.")
    list(value = max(values, na.rm = TRUE), records = do.call(rbind, rows), native.fits = calls)
}
