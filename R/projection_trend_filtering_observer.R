# Optional internal observation. Numerical code never uses observer output.
.ptf.infrastructure.error <- function(message, resource = FALSE) {
    structure(list(message = message, call = NULL),
        class = c(if (resource) "projection_resource_error" else "projection_observer_error",
                  "projection_infrastructure_error", "error", "condition"))
}

.ptf.numerical.error <- function(e) {
    if (inherits(e, "projection_infrastructure_error")) stop(e)
    e
}

.ptf.observer.guard <- function(code) {
    tryCatch(force(code), error = function(e) {
        if (inherits(e, "projection_infrastructure_error")) stop(e)
        stop(.ptf.infrastructure.error(conditionMessage(e)))
    })
}

.ptf.observer <- function(sink, metadata = list()) {
    stopifnot(is.function(sink), is.list(metadata))
    o <- new.env(parent = emptyenv())
    o$sink <- sink; o$metadata <- metadata; o$sequence <- 0L
    o$stack <- list(); o$native <- 0L; o$diagnostic.seconds <- 0
    o$serialization.seconds <- 0
    o$last.snapshot <- NULL
    o
}

.ptf.observer.clock <- function() .ptf.observer.guard(unname(proc.time()[3L]))

.ptf.observe.diagnostic <- function(o, code) {
    if (is.null(o)) return(force(code))
    started <- .ptf.observer.clock(); previous <- o$diagnostic.seconds
    value <- .ptf.observer.guard(force(code))
    o$diagnostic.seconds <- o$diagnostic.seconds + max(0,
        .ptf.observer.clock()-started-(o$diagnostic.seconds-previous))
    value
}

.ptf.observe.emit <- function(o, record) {
    if (is.null(o)) return(invisible(NULL))
    .ptf.observer.guard({
        started <- .ptf.observer.clock()
        previous <- o$serialization.seconds
        o$sink(c(list(schema = 1L, metadata = o$metadata), record))
        elapsed <- .ptf.observer.clock()-started
        o$diagnostic.seconds <- o$diagnostic.seconds + elapsed
        # Count the union of sink spans, including any nested serialization once.
        o$serialization.seconds <- o$serialization.seconds + max(0,
            elapsed-(o$serialization.seconds-previous))
    })
    invisible(NULL)
}

.ptf.observe.begin <- function(o, kind, data = list()) {
    if (is.null(o)) return(NULL)
    .ptf.observer.guard({
        o$sequence <- o$sequence+1L
        id <- o$sequence
        parent <- if (length(o$stack)) tail(o$stack, 1L)[[1L]]$id else NA_integer_
        frame <- list(id = id, parent = parent, kind = kind,
            start = .ptf.observer.clock(), diagnostic.start = o$diagnostic.seconds,
            serialization.start = o$serialization.seconds,
            child.elapsed = 0, child.diagnostic = 0, native.start = o$native)
        .ptf.observe.emit(o, list(id = id, parent = parent, kind = kind,
                                event = "begin", data = data))
        o$stack[[length(o$stack)+1L]] <- frame
        id
    })
}

.ptf.observe.end <- function(o, id, data = list()) {
    if (is.null(o)) return(invisible(NULL))
    .ptf.observer.guard({
        n <- length(o$stack)
        if (!n || !identical(o$stack[[n]]$id, id)) stop("Observer span nesting mismatch.")
        f <- o$stack[[n]]
        elapsed <- .ptf.observer.clock()-f$start
        diagnostic <- o$diagnostic.seconds-f$diagnostic.start
        timing <- list(inclusive.elapsed = elapsed, diagnostic.elapsed = diagnostic,
            serialization.elapsed = o$serialization.seconds-f$serialization.start,
            exclusive.fitting.elapsed = max(0, elapsed-f$child.elapsed-
                (diagnostic-f$child.diagnostic)))
        .ptf.observe.emit(o, list(id = id, parent = f$parent, kind = f$kind,
            event = "end", native.entries = o$native-f$native.start,
            timing = timing, data = data))
        # Charge closing serialization to the parent's diagnostic time.
        total <- .ptf.observer.clock()-f$start
        total.diagnostic <- o$diagnostic.seconds-f$diagnostic.start
        o$stack[[n]] <- NULL
        if (n > 1L) {
            o$stack[[n-1L]]$child.elapsed <- o$stack[[n-1L]]$child.elapsed+total
            o$stack[[n-1L]]$child.diagnostic <- o$stack[[n-1L]]$child.diagnostic+total.diagnostic
        }
    })
    invisible(NULL)
}

.ptf.observed <- function(o, kind, data, code, summarize = function(value) list(status = "complete")) {
    if (is.null(o)) return(force(code))
    id <- .ptf.observe.begin(o, kind, data)
    value <- tryCatch(force(code), error = function(e) {
        .ptf.numerical.error(e)
        .ptf.observe.end(o, id, list(status = "error", message = conditionMessage(e)))
        stop(e)
    })
    detail <- .ptf.observe.diagnostic(o, summarize(value))
    .ptf.observe.end(o, id, detail)
    value
}

.ptf.observe.snapshot <- function(problem, state, reason) {
    o <- problem$observer
    if (is.null(o)) return(invisible(NULL))
    .ptf.observer.guard({
        id <- .ptf.observe.begin(o, "snapshot", list(reason = reason))
        snapshot <- .ptf.observe.diagnostic(o, list(reason = reason,
            state = state[intersect(names(state), c("intercept", "directions", "terms",
                "components", "checks", "converged", "status"))],
            problem = problem[c("X", "family", "preprocessing", "weights", "w", "offset",
                "y", "lambda", "rho", "backend")],
            link = problem$offset+state$intercept+rowSums(state$terms)))
        .ptf.observe.end(o, id, snapshot)
        o$last.snapshot <- id
    })
    invisible(NULL)
}

.ptf.native.call <- function(o, before.fit, metadata, code) {
    if (is.null(o)) {
        if (!is.null(before.fit)) before.fit()
        return(force(code))
    }
    id <- .ptf.observe.begin(o, "native", metadata)
    started <- NULL
    value <- tryCatch({
        .ptf.observed(o, "native_reservation", list(),
            if (!is.null(before.fit)) before.fit())
        o$native <- o$native+1L
        started <- .ptf.observer.clock()
        result <- force(code)
        elapsed <- .ptf.observer.clock()-started
        result
    }, error = function(e) {
        .ptf.numerical.error(e)
        elapsed <- if (is.null(started)) NA_real_ else .ptf.observer.clock()-started
        .ptf.observe.end(o, id, list(status = "error", message = conditionMessage(e),
            solve.elapsed = elapsed, entered = !is.null(started)))
        stop(e)
    })
    detail <- .ptf.observe.diagnostic(o, lapply(value, function(f)
        f[intersect(names(f), c("status", "converged", "iterations", "barrier.stages",
            "backtracks", "max.damping", "initialization.iterations", "active.knots"))]))
    .ptf.observe.end(o, id, list(status = "returned", diagnostics = detail,
                               solve.elapsed = elapsed, entered = TRUE))
    value
}

# Practical stopping does not decide numerical acceptance.
.ptf.practical <- function(null.loss, epsilon = 0, patience = 2L) {
    .ptf.scalar(null.loss, "null.loss")
    .ptf.scalar(epsilon, "epsilon")
    .ptf.scalar(patience, "patience", integer = TRUE, positive = TRUE)
    list(scale = null.loss, epsilon = epsilon, patience = patience, streak = 0L,
         stopped = FALSE, gain = NA_real_, small = FALSE)
}

.ptf.practical.update <- function(rule, accepted, previous, current) {
    rule$gain <- NA_real_; rule$small <- FALSE
    if (!isTRUE(accepted)) return(rule)
    rule$gain <- previous-current
    rule$small <- rule$scale > 0 && rule$epsilon > 0 &&
        is.finite(rule$gain) && rule$gain > 0 && rule$gain <= rule$epsilon*rule$scale
    rule$streak <- if (rule$small) rule$streak+1L else 0L
    rule$stopped <- rule$streak >= rule$patience
    rule
}
