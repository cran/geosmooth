test_that("accepted-depth caps count adopted joint updates and preserve slot precedence", {
    skip <- 0L
    testthat::local_mocked_bindings(
        .ptf.curve.slopes = function(...) rep(1, 6),
        .ptf.gn.equations = function(...) list(A = diag(2), b = c(0, 1), floored = 0),
        .ptf.gn.increment = function(...) list(status = "proposal", mu = .01,
                                              delta = c(0, .1), tangent.score = 1),
        .ptf.backfit = function(problem, state) { state$solver.calls <- 0L; state },
        .ptf.duplicate = function(...) {
            if (skip > 0L) { skip <<- skip-1L; TRUE } else FALSE
        },
        .package = "geosmooth")
    initial <- list(directions = matrix(c(1, 0), 1), terms = matrix(1, 6, 1),
        components = list(list(fit = list(design = list(span = 1)))),
        checks = list(objective = 1), converged = TRUE, status = "converged")
    run <- function(depth, reject = 0L, fail.joint = FALSE) {
        counter <- new.env(); counter$visits <- counter$events <- list()
        p <- list(control = list(improvement.tolerance = 1e-9), counter = counter,
            offset = rep(0, 6), Z = matrix(0, 6, 2), y = rep(0, 6), w = rep(1/6, 6),
            family = "gaussian")
        i <- 0L
        evaluate <- function(state, u, ...) {
            i <<- i+1L
            if (i <= reject) return(NULL)
            state$directions[1, ] <- u
            state$checks$objective <- state$checks$objective-.001
            if (fail.joint) state$converged <- FALSE
            list(state = state, row = i)
        }
        result <- geosmooth:::.ptf.gn.visit(p, initial, 1, 8, 1, 1, 1,
            evaluate, function(...) NULL, function(a, b) a < b-1e-9, accepted.limit = depth)
        list(result = result, visits = counter$visits, events = counter$events)
    }
    a <- run(Inf); b <- run(8)
    b$visits[[1]]$accepted.limit <- b$visits[[1]]$depth.reached <- NULL
    expect_identical(a, b)
    expect_identical(a$visits[[1]]$stopping, "proposal_limit")
    for (reject in c(0L, 2L, 7L)) {
        a <- run(1, reject)
        expect_equal(a$visits[[1]]$accepted, 1)
        expect_equal(a$visits[[1]]$evaluated, reject+1)
        expect_identical(a$visits[[1]]$stopping,
            if (reject == 7L) "proposal_limit" else "accepted_depth_stop")
        expect_equal(a$result$state$checks$objective, .999)
        expect_true(a$result$state$converged)
    }
    expect_equal(run(1, 8)$visits[[1]]$accepted, 0)
    failed <- run(1, fail.joint = TRUE)
    expect_equal(failed$visits[[1]]$accepted, 0)
    expect_identical(failed$result$state, initial)
    expect_equal(failed$visits[[1]]$evaluated, 8)
    # Each invocation starts a fresh visit rather than sharing an accepted counter.
    expect_identical(run(1), run(1))
    skip <- 2L
    duplicates <- run(1)
    expect_equal(duplicates$visits[[1]]$generated, 3)
    expect_equal(duplicates$visits[[1]]$evaluated, 1)
    expect_equal(duplicates$visits[[1]]$accepted, 1)
    expect_error(run(0), "accepted.limit")
    expect_error(run(1.5), "integer")
})
