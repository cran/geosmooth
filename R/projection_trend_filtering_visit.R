# Evaluated-trial Gauss--Newton refinement with optional practical stopping.
.ptf.gn.visit <- function(problem, state, term, limit, start, pass, addition,
                           evaluate, record, improves, accepted.limit = Inf) {
    .ptf.scalar(accepted.limit, "accepted.limit", positive = TRUE, infinite = TRUE)
    if (is.finite(accepted.limit) && accepted.limit != floor(accepted.limit))
        stop("accepted.limit must be an integer or Inf.", call. = FALSE)
    o <- problem$observer
    visit.id <- .ptf.observe.begin(o, "visit", list(start = start, pass = pass,
        addition = addition, term = term, limit = limit))
    rule <- if (is.null(problem$practical)) .ptf.practical(0) else problem$practical
    rule$streak <- 0L; rule$stopped <- FALSE
    used <- generated <- accepted <- calls <- 0L; accepted.rows <- integer()
    stopping <- "proposal_limit"; events <- list(); objectives <- numeric()
    while (used < limit) {
        fit <- state$components[[term]]$fit; u <- state$directions[term, ]
        g <- tryCatch(.ptf.observed(o, "slopes", list(),
            .ptf.curve.slopes(fit)), error = .ptf.numerical.error)
        if (inherits(g, "error") || any(!is.finite(g))) { stopping <- "nonfinite_slope"; break }
        flat.norm <- max(abs(g))*fit$design$span
        flat.threshold <- 1e-10*max(1, abs(state$terms[, term]))
        if (all(g == 0) || flat.norm <= flat.threshold || length(u) == 1L) {
            stopping <- if (all(g == 0)) "flat_curve" else if (flat.norm <= flat.threshold)
                "near_flat" else "rank_one"
            break
        }
        eta <- problem$offset+state$intercept+rowSums(state$terms)
        calc <- tryCatch(.ptf.observed(o, "proposal_calculation", list(), {
            eq <- .ptf.gn.equations(problem$Z, problem$y, eta, g, problem$w, problem$family)
            list(eq = eq, inc = .ptf.gn.increment(eq$A, eq$b, u))
        }), error = .ptf.numerical.error)
        if (inherits(calc, "error")) { stopping <- "proposal_failed"; break }
        eq <- calc$eq; inc <- calc$inc
        if (inc$status != "proposal") { stopping <- inc$status; break }
        improved <- FALSE
        for (escalation in 0:2) {
            mu <- inc$mu*10^escalation
            delta <- .ptf.gn.increment(eq$A, eq$b, u, mu)$delta
            for (fraction in c(1, .5, .25, .125)) {
                if (used >= limit) break
                proposal <- .ptf.retract(u, delta, fraction)
                generated <- generated+1L
                # Duplicate proposals consume no evaluated-trial slot.
                if (.ptf.duplicate(proposal, state$directions[-term, , drop = FALSE])) next
                used <- used+1L
                before <- state$checks$objective; streak.before <- rule$streak
                proposal.id <- .ptf.observe.begin(o, "proposal", list(slot = used,
                    before.direction = u, direction = proposal,
                    angle.degrees = acos(min(1, abs(sum(u*proposal))))*180/pi,
                    mu = mu, fraction = fraction, escalation = escalation,
                    curvature.floored = eq$floored, slope.norm = flat.norm,
                    slope.threshold = flat.threshold, tangent.score = inc$tangent.score))
                joint.called <- FALSE
                trial <- evaluate(state, proposal, term, start, "refinement")
                status <- "component_rejected"; take <- FALSE
                if (is.null(trial)) status <- "component_failed" else
                if (improves(trial$state$checks$objective, state$checks$objective)) {
                    joint.called <- TRUE
                    joint <- .ptf.backfit(problem, trial$state)
                    calls <- calls+joint$solver.calls
                    if (!joint$converged) {
                        record(start, "joint_refinement", term, joint$status, joint$checks$objective)
                        status <- "joint_failed"
                    } else if (improves(joint$checks$objective, state$checks$objective)) {
                        state <- joint; accepted.rows <- c(accepted.rows, trial$row)
                        objectives <- c(objectives, state$checks$objective)
                        accepted <- accepted+1L; take <- improved <- TRUE; status <- "accepted"
                    } else status <- "joint_rejected"
                }
                rule <- .ptf.practical.update(rule, take, before, state$checks$objective)
                if (take) .ptf.observe.snapshot(problem, state, "accepted_refinement")
                .ptf.observe.end(o, proposal.id, list(status = status, accepted = take,
                    joint.called = joint.called,
                    component.objective = if (is.null(trial)) NA_real_ else trial$state$checks$objective,
                    before.objective = before, retained.objective = state$checks$objective,
                    nominal.tolerance = problem$control$improvement.tolerance*max(1, abs(before)),
                    gain = rule$gain, null.loss = rule$scale, epsilon = rule$epsilon,
                    relative.gain = if (rule$scale > 0) rule$gain/rule$scale else NA_real_,
                    practical.threshold = rule$epsilon*rule$scale, small = rule$small,
                    streak.before = streak.before, streak.after = rule$streak,
                    practical.stop = rule$stopped,
                    retained.snapshot = if (is.null(o)) NULL else o$last.snapshot))
                events[[length(events)+1L]] <- data.frame(slot = used, step = accepted+!take,
                    escalation = escalation, mu = mu, fraction = fraction, status = status,
                    accepted = take, curvature.floored = eq$floored, slope.norm = flat.norm,
                    tangent.score = inc$tangent.score)
                if (take) break
            }
            if (improved || used >= limit) break
        }
        if (rule$stopped) { stopping <- "practical_gain_stop"; break }
        if (accepted >= accepted.limit && used < limit) {
            stopping <- "accepted_depth_stop"; break
        }
        if (!improved) { stopping <- if (used >= limit) "proposal_limit" else "rejected"; break }
    }
    id <- length(problem$counter$visits)+1L
    problem$counter$visits[[id]] <- data.frame(start = start, addition = addition,
        pass = pass, term = term, generated = generated, evaluated = used,
        accepted = accepted, stopping = stopping,
        proposal.limit.reached = used >= limit, unused.slots = limit-used,
        null.loss = rule$scale, epsilon = rule$epsilon, streak = rule$streak)
    problem$counter$events[[id]] <- if (length(events)) do.call(rbind, events) else data.frame()
    if (is.finite(accepted.limit)) {
        problem$counter$visits[[id]]$accepted.limit <- accepted.limit
        problem$counter$visits[[id]]$depth.reached <- accepted >= accepted.limit
    }
    .ptf.observe.end(o, visit.id, list(status = stopping, evaluated = used,
        generated = generated, accepted = accepted, proposal.limit.reached = used >= limit,
        unused.slots = limit-used, streak = rule$streak))
    list(state = state, solver.calls = calls, accepted.rows = accepted.rows, objectives = objectives)
}
