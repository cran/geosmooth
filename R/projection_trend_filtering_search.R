.ptf.with.seed <- function(seed, code) {
    present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (present) previous <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(if (present) assign(".Random.seed", previous, envir = .GlobalEnv) else
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
    force(code)
}

.ptf.bank <- function(problem, supplied) {
    p <- ncol(problem$Z); control <- problem$control
    U <- diag(p)
    if (!is.null(supplied)) U <- rbind(U, supplied)
    for (k in seq_len(control$sparse.directions)) {
        u <- numeric(p); size <- min(p, 2L+(k %% 2L))
        ix <- sample.int(p, size); u[ix] <- stats::rnorm(size)
        U <- rbind(U, .ptf.canonical(u))
    }
    for (k in seq_len(control$dense.directions))
        U <- rbind(U, .ptf.canonical(stats::rnorm(p)))
    .ptf.unique.directions(U)
}

.ptf.duplicate <- function(u, U) {
    nrow(U) > 0L && any(vapply(seq_len(nrow(U)), function(j)
        max(abs(u-U[j, ])) <= 1e-12, logical(1)))
}

.ptf.search <- function(problem, supplied, nterms) {
    control <- problem$control; p <- ncol(problem$Z)
    null <- .ptf.backfit(problem, .ptf.initial(problem, matrix(0, 0L, p)))
    if (nterms == 0L) return(list(state = null, search = list(history = data.frame(),
        status = "intercept_only", starts = data.frame(), solver.calls = 0L)))
    rows <- list(); starts <- list(); calls <- 0L; best <- NULL
    record <- function(start, stage, term, status, objective = NA_real_, accepted = FALSE,
                       direction = NULL) {
        rows[[length(rows)+1L]] <<- data.frame(start = start, stage = stage, term = term,
            status = status, objective = objective, accepted = accepted,
            direction = I(list(direction)))
        length(rows)
    }
    evaluate <- function(state, u, j, start, stage) {
        calls <<- calls+1L
        candidate <- tryCatch(.ptf.update(problem, state, j, u), error = .ptf.numerical.error)
        if (inherits(candidate, "error")) {
            record(start, stage, j, conditionMessage(candidate), direction = u)
            return(NULL)
        }
        candidate$checks <- .ptf.diagnose(problem, candidate)
        id <- record(start, stage, j, "evaluated", candidate$checks$objective, direction = u)
        if (!is.finite(candidate$checks$objective)) return(NULL)
        list(state = candidate, row = id)
    }
    improves <- function(a, b) is.finite(a) && a < b-control$improvement.tolerance*max(1,abs(b))
    for (start in seq_len(control$starts)) {
        bank <- .ptf.bank(problem, supplied); state <- null
        path <- state$checks$objective; search.status <- "component_limit"
        for (j in seq_len(nterms)) {
            expanded <- state
            expanded$directions <- rbind(expanded$directions, numeric(p))
            expanded$terms <- cbind(expanded$terms, 0)
            expanded$components <- c(expanded$components, list(NULL))
            candidate.best <- NULL; objective <- state$checks$objective
            evaluated <- 0L; attempted <- 0L
            for (k in seq_len(nrow(bank))) {
                u <- bank[k, ]
                if (.ptf.duplicate(u, state$directions)) next
                attempted <- attempted+1L
                trial <- evaluate(expanded, u, j, start, "addition")
                if (is.null(trial)) next
                evaluated <- evaluated+1L
                if (improves(trial$state$checks$objective, objective)) {
                    candidate.best <- trial; objective <- trial$state$checks$objective
                }
            }
            if (is.null(candidate.best)) {
                search.status <- if (attempted && !evaluated) "all_candidates_failed" else "no_improving_addition"
                if (search.status == "all_candidates_failed") {
                    state$converged <- FALSE; state$status <- search.status
                }
                break
            }
            joint <- .ptf.backfit(problem, candidate.best$state)
            calls <- calls+joint$solver.calls
            if (!joint$converged || !improves(joint$checks$objective, state$checks$objective)) {
                record(start,"joint_addition",j,joint$status,joint$checks$objective)
                search.status <- "addition_refit_failed"
                state$converged <- FALSE; state$status <- search.status
                break
            }
            rows[[candidate.best$row]]$accepted <- TRUE
            state <- joint; path <- c(path, state$checks$objective)
            for (pass in seq_len(control$refinement.passes)) {
                accepted.pass <- FALSE
                for (term in seq_len(nrow(state$directions))) {
                    proposals <- 0L
                    for (step in control$steps) {
                        if (proposals >= control$refinement.proposals) break
                        quota <- min(control$refinement.proposals-proposals,
                                     ceiling(control$refinement.proposals/length(control$steps)))
                        step.proposals <- 0L
                        # Coordinate tangent polls and random tangents also allow
                        # a sparse starting support to gain or lose coordinates.
                        size <- max(1L, ceiling(quota/4L))
                        tangents <- rbind(diag(p)[sample.int(p,min(p,size)),,drop=FALSE],
                                         matrix(stats::rnorm(p*size),size,p))
                        for (k in seq_len(nrow(tangents))) {
                            if (step.proposals >= quota) break
                            for (sign in c(-1, 1)) {
                                if (step.proposals >= quota) break
                                u <- state$directions[term, ]
                                tangent <- tangents[k, ]-sum(tangents[k, ]*u)*u
                                length.tangent <- sqrt(sum(tangent^2))
                                if (length.tangent < 1e-12) next
                                proposal <- .ptf.canonical(u+sign*step*tangent/length.tangent)
                                if (.ptf.duplicate(proposal,state$directions[-term,,drop=FALSE])) next
                                proposals <- proposals+1L
                                step.proposals <- step.proposals+1L
                                trial <- evaluate(state,proposal,term,start,"refinement")
                                if (is.null(trial) || !improves(trial$state$checks$objective,state$checks$objective)) next
                                joint <- .ptf.backfit(problem,trial$state)
                                calls <- calls+joint$solver.calls
                                if (!joint$converged) {
                                    record(start,"joint_refinement",term,joint$status,joint$checks$objective)
                                    next
                                }
                                if (improves(joint$checks$objective,state$checks$objective)) {
                                    rows[[trial$row]]$accepted <- TRUE
                                    state <- joint; path <- c(path,state$checks$objective)
                                    accepted.pass <- TRUE
                                }
                            }
                        }
                    }
                }
                if (!accepted.pass) break
            }
        }
        state$objective.history <- path
        starts[[start]] <- data.frame(start = start, objective = state$checks$objective,
            nterms = ncol(state$terms), converged = state$converged,
            fit.status = state$status, search.status = search.status)
        if (is.null(best) || (state$converged && !best$converged) ||
            (identical(state$converged,best$converged) && state$checks$objective < best$checks$objective)) {
            best <- state; best.status <- search.status
        }
    }
    list(state = best, search = list(history = if (length(rows)) do.call(rbind,rows) else data.frame(),
         status = best.status, starts = do.call(rbind,starts), solver.calls = calls))
}
