joint_fixture_problem <- function(family = "gaussian", lambda = Inf, X = NULL,
                                  y = NULL, control = list()) {
    if (is.null(X)) X <- cbind(seq(-1, 1, length.out = 24), cos(seq_len(24)))
    n <- nrow(X)
    eta <- .25+.6*X[, 1]+.25*X[, 2]^2+.2*pmax(X[, 1]-.1, 0)^2
    if (is.null(y)) y <- if (family == "gaussian") eta else plogis(eta)
    problem <- geosmooth:::.ptf.problem(X, y, rep(c(1, 2, 4), length.out = n),
        .15*sin(seq_len(n)), family, "isotropic", lambda, .001,
        geosmooth:::.ptf.control(utils::modifyList(list(joint.solver = "dense"), control)), "banded")
    sink <- getOption("geosmooth.fixture.observer")
    if (is.function(sink)) problem$observer <- geosmooth:::.ptf.observer(sink)
    problem
}

joint_fixture_record <- function(label, data = list()) {
    callback <- getOption("geosmooth.fixture.record")
    if (is.function(callback)) callback(label, data)
    invisible(NULL)
}

joint_fixture_close <- function(actual, reference, tolerance) {
    expect_lte(max(abs(actual-reference)), tolerance*max(1, abs(reference)))
}

test_that("dense joint Gaussian and binomial fits reconstruct and certify the original objective", {
    for (family in c("gaussian", "binomial")) for (lambda in c(0, .001, .01, Inf)) {
        label <- paste(family, "joint", lambda)
        joint_fixture_record(paste("START", label))
        p <- joint_fixture_problem(family, lambda)
        initial <- geosmooth:::.ptf.initial(p, diag(2))
        g <- geosmooth:::.ptf.joint.geometry(p, initial)
        fit <- geosmooth:::.ptf.backfit(p, initial)
        expect_true(fit$converged, info = label)
        expect_false(fit$joint.work$fallback, info = label)
        expect_equal(fit$joint.work$direct.entries, 1)
        expect_true(geosmooth:::.ptf.diagnose(p, fit)$converged)
        if (isTRUE(fit$joint.work$fallback) || !isTRUE(fit$converged)) next
        native <- fit$joint.diagnostics$native
        eta <- p$offset+fit$intercept+rowSums(fit$terms)
        joint_fixture_close(eta, native$linear.predictors, 1e-10)
        joint_fixture_close(fit$checks$objective, native$objective/sum(p$weights), 1e-10)
        joint_fixture_close(colSums(fit$terms*p$w), c(0, 0), 1e-10)
        for (j in 1:2) {
            component <- fit$components[[j]]; f <- component$fit; d <- f$design
            expect_identical(f$backend, "dense")
            beta <- f$fits[[1]]$coefficients[-1]/d$scales
            t <- (f$x-d$lo)/d$span
            # Explicit polynomial evaluation, independent of the matrix constructor.
            raw <- beta[1]*(t-d$center[1])+beta[2]*(t^2-d$center[2])
            for (k in 3:length(beta)) raw <- raw+beta[k]*(
                .5*(t-d$knots[k-1])*(t-d$knots[k])*(t > d$knots[k])-d$center[k])
            joint_fixture_close(raw-component$center, fit$terms[, j], 1e-10)
            joint_fixture_close(as.numeric(predict(f, f$x, offset = rep(0, length(t)),
                type = "link"))-component$center, fit$terms[, j], 1e-10)
            mid <- d$lo+d$span*(head(d$knots, -1)+tail(d$knots, -1))/2
            step <- 1e-6*d$span
            slope <- (as.numeric(predict(f, mid+step, offset = rep(0, length(mid)), type = "link"))-
                as.numeric(predict(f, mid-step, offset = rep(0, length(mid)), type = "link")))/(2*step)
            joint_fixture_close(geosmooth:::.ptf.curve.slopes(f, mid), slope, 2e-3)
            outside <- c(d$lo-d$span, d$lo, d$lo+d$span, d$lo+2*d$span)
            pred <- as.numeric(predict(f, outside, offset = rep(0, 4), type = "link"))
            expect_identical(pred[1], pred[2]); expect_identical(pred[3], pred[4])
            expect_equal(f$fits[[1]]$dual[-1], native$dual[g$indices[[j]]], tolerance = 0)
        }
        b <- native$coefficients; mass <- sum(p$weights)
        smooth <- function(b) {
            z <- as.vector(g$A %*% b+p$offset)
            loss <- if (family == "gaussian") (z-p$y)^2/2 else
                pmax(z, 0)-p$y*z+log1p(exp(-abs(z)))
            sum(p$w*loss)+sum((g$R %*% b)^2)/(2*mass)
        }
        residual <- if (family == "gaussian") eta-p$y else plogis(eta)-p$y
        gradient <- as.vector(crossprod(g$A, p$w*residual)+crossprod(g$R, g$R %*% b)/mass)
        selected <- if (is.infinite(lambda)) which(g$penalty == 0) else seq_along(b)
        finite.difference <- vapply(selected, function(k) {
            step <- 1e-6; plus <- minus <- b; plus[k] <- plus[k]+step; minus[k] <- minus[k]-step
            (smooth(plus)-smooth(minus))/(2*step)
        }, numeric(1))
        joint_fixture_close(gradient[selected], finite.difference, 2e-6)
        kkt <- gradient[selected]+if (is.finite(lambda))
            lambda*g$penalty[selected]*native$dual[selected] else 0
        expect_lte(max(abs(kkt)), 1e-7)
        p$control$joint.solver <- "backfit"
        reference <- geosmooth:::.ptf.backfit(p, initial)
        differences <- list()
        if (reference$converged) {
            reference.eta <- p$offset+reference$intercept+rowSums(reference$terms)
            differences <- list(link = max(abs(eta-reference.eta)),
                terms = max(abs(fit$terms-reference$terms)),
                objective = abs(fit$checks$objective-reference$checks$objective))
            joint_fixture_close(eta, reference.eta, 2e-4)
            joint_fixture_close(fit$terms, reference$terms, 5e-4)
            joint_fixture_close(fit$checks$objective, reference$checks$objective, 1e-8)
        }
        joint_fixture_record(paste("END", label), list(checks = fit$checks,
            direct = fit$joint.work, native.iterations = native$iterations,
            reference.certified = reference$converged, reference.status = reference$status,
            differences = differences))
    }
})

test_that("Gaussian joint limits agree with independent quadratic solutions", {
    for (lambda in c(0, Inf)) {
        p <- joint_fixture_problem(lambda = lambda)
        initial <- geosmooth:::.ptf.initial(p, diag(2))
        fit <- geosmooth:::.ptf.backfit(p, initial)
        expect_true(fit$converged); expect_false(fit$joint.work$fallback)
        if (lambda == 0) {
            b <- sum(p$w*(p$y-p$offset))
            h <- (p$y-p$offset-b)/(2+p$rho)
            joint_fixture_close(fit$intercept, b, 1e-6)
            joint_fixture_close(fit$intercept+rowSums(fit$terms), b+2*h, 1e-6)
            joint_fixture_close(fit$terms, cbind(h, h), 5e-4)
            joint_fixture_record("END gaussian zero-penalty analytic solution", list(
                link.error = max(abs(fit$intercept+rowSums(fit$terms)-b-2*h)),
                component.error = max(abs(fit$terms-cbind(h,h)))))
        } else {
            z <- lapply(1:2, function(j) {
                t <- geosmooth:::.ptf.projection(p, diag(2)[j, ])
                H <- cbind(t, t^2); sweep(H, 2, colSums(H*p$w))
            })
            A <- cbind(1, z[[1]], z[[2]])
            ridge <- matrix(0, 5, 5)
            for (j in 1:2) ridge[2*j+(0:1), 2*j+(0:1)] <-
                p$rho*crossprod(z[[j]], z[[j]]*p$w)
            b <- solve(crossprod(A, A*p$w)+ridge, crossprod(A, (p$y-p$offset)*p$w))
            joint_fixture_close(fit$intercept+rowSums(fit$terms), as.vector(A %*% b), 1e-6)
        }
    }
})

test_that("screened structured components convert to consistent dense geometry", {
    for (family in c("gaussian", "binomial")) {
        p <- joint_fixture_problem(family, .001)
        state <- geosmooth:::.ptf.initial(p, diag(2))
        state <- geosmooth:::.ptf.update(p, state, 1)
        state <- geosmooth:::.ptf.update(p, state, 2)
        g <- geosmooth:::.ptf.joint.geometry(p, state)
        for (block in g$blocks) {
            expect_false(is.null(block$geometry.difference))
            expect_lte(block$geometry.difference$center, 1e-10)
            expect_lte(block$geometry.difference$scales, 1e-10)
        }
        fit <- geosmooth:::.ptf.backfit(p, state)
        expect_true(fit$converged); expect_false(fit$joint.work$fallback)
        joint_fixture_record(paste("END", family, "screened geometry"), list(
            geometry.difference = lapply(g$blocks, `[[`, "geometry.difference"), checks = fit$checks))
    }
})

test_that("failed native and fallback entries are counted without certifying a truncated fit", {
    for (family in c("gaussian", "binomial")) {
        p <- joint_fixture_problem(family, Inf, control = list(inner = list(max.iterations = 1L)))
        entries <- 0L; p$before.fit <- function() entries <<- entries+1L
        events <- list()
        p$observer <- geosmooth:::.ptf.observer(function(e) events[[length(events)+1L]] <<- e)
        state <- geosmooth:::.ptf.initial(p, diag(2))
        fit <- geosmooth:::.ptf.backfit(p, state)
        expect_false(fit$converged); expect_true(fit$joint.work$fallback)
        expect_equal(fit$joint.work$direct.entries, 1)
        expect_equal(fit$solver.calls, 2); expect_equal(entries, 2)
        expect_equal(ncol(fit$terms), 2)
        expect_identical(fit$terms, state$terms)
        expect_length(p$observer$stack, 0)
        calls <- Filter(function(e) e$kind == "native" && e$event == "begin", events)
        expect_identical(vapply(calls, function(e) e$data$call.type, ""), c("joint", "univariate"))
        expect_equal(p$observer$native, entries)
    }
})

test_that("joint geometry handles ties and close directions without changing the loss", {
    tied <- as.matrix(expand.grid(seq(-1, 1, length.out = 6), seq(-1, 1, length.out = 4)))
    for (family in c("gaussian", "binomial")) for (kind in c("ties", "close")) {
        p <- joint_fixture_problem(family, .001, if (kind == "ties") tied else NULL)
        U <- if (kind == "ties") diag(2) else rbind(c(1, 0), geosmooth:::.ptf.canonical(c(1, .02)))
        fit <- geosmooth:::.ptf.backfit(p, geosmooth:::.ptf.initial(p, U))
        expect_true(fit$converged); expect_false(fit$joint.work$fallback)
        expect_lte(fit$checks$stationarity, 1e-7)
        expect_lte(abs(fit$checks$complementarity), 1e-8)
        p2 <- p; p2$weights <- p$weights*7
        scaled <- geosmooth:::.ptf.backfit(p2, geosmooth:::.ptf.initial(p2, U))
        joint_fixture_close(scaled$terms, fit$terms, 5e-4)
        joint_fixture_close(scaled$checks$objective, fit$checks$objective, 1e-8)
        joint_fixture_record(paste("END", family, kind), list(checks = fit$checks, direct = fit$joint.work))
    }
    p <- joint_fixture_problem()
    initial <- geosmooth:::.ptf.initial(p, diag(2))
    p$X[, 1] <- rep(1:3, 8)
    expect_error(geosmooth:::.ptf.joint.geometry(p, initial), "four distinct")
    p$X[, 1] <- c(-1e308, -1e-300, 0, rep(1, 21))
    p$preprocessing <- list(retained = 1:2, center = c(0, 0), scale = c(1, 1))
    expect_error(geosmooth:::.ptf.joint.geometry(p, initial), "collapsed")
})

test_that("binomial single-class errors differ from ridge-regularized separation", {
    X <- cbind(seq(-1, 1, length.out = 24), cos(seq_len(24)))
    for (y in list(rep(0, 24), rep(1, 24)))
        expect_error(joint_fixture_problem("binomial", Inf, X, y), "finite unpenalized intercept")
    for (y in list(as.numeric(X[, 1] > 0), rep(c(0, 1, 1, 0), 6))) {
        p <- joint_fixture_problem("binomial", Inf, X, y)
        fit <- geosmooth:::.ptf.backfit(p, geosmooth:::.ptf.initial(p, diag(2)))
        expect_true(fit$converged); expect_false(fit$joint.work$fallback)
        expect_true(all(is.finite(fit$terms)))
        joint_fixture_record("END binomial mixed outcomes", list(checks = fit$checks))
    }
})

test_that("joint policy leaves zero and one components unchanged and reuses certified states", {
    for (family in c("gaussian", "binomial")) {
        p <- joint_fixture_problem(family)
        for (K in 0:1) {
            U <- diag(2)[seq_len(K), , drop = FALSE]
            initial <- geosmooth:::.ptf.initial(p, U)
            a <- geosmooth:::.ptf.backfit(p, initial)
            q <- p; q$control$joint.solver <- "backfit"
            expect_identical(a, geosmooth:::.ptf.backfit(q, initial))
        }
        fit <- geosmooth:::.ptf.backfit(p, geosmooth:::.ptf.initial(p, diag(2)))
        entries <- 0L; p$before.fit <- function() entries <<- entries+1L
        reused <- geosmooth:::.ptf.backfit(p, fit)
        expect_equal(entries, 0); expect_equal(reused$solver.calls, 0)
        expect_equal(reused$joint.work$direct.entries, 0)
        expect_identical(reused$joint.work$direct.status, "incoming_certified")
        expect_identical(reused$terms, fit$terms)
    }
    expect_error(geosmooth:::.ptf.control(list(joint.solver = "other")), "joint.solver")
    expect_error(geosmooth:::.ptf.control(list(joint.max.workspace.bytes = 0)), "workspace")
})

test_that("direct failure and objective guards fall back from untouched input with full accounting", {
    p <- joint_fixture_problem()
    initial <- geosmooth:::.ptf.initial(p, diag(2))
    original <- geosmooth:::.ptf.backfit.core
    incoming <- serialize(initial, NULL)
    fallback.calls <- 0L; entries <- 0L
    p$before.fit <- function() entries <<- entries+1L
    reason <- "certificate"
    testthat::local_mocked_bindings(
        .ptf.joint.reconstruct = function(problem, state, geometry, native) {
            if (reason == "certificate") stop(geosmooth:::.ptf.joint.failure("injected_certificate"))
            if (reason == "unexpected") stop("unexpected programming error")
            state$checks <- list(objective = geosmooth:::.ptf.diagnose(problem, state)$objective+1)
            state
        },
        .ptf.backfit.core = function(problem, state) {
            fallback.calls <<- fallback.calls+1L
            expect_identical(serialize(state, NULL), incoming)
            original(problem, state)
        }, .package = "geosmooth")
    for (why in c("certificate", "objective")) {
        reason <- why; entries <- 0L
        fit <- geosmooth:::.ptf.backfit(p, initial)
        expect_true(fit$converged); expect_true(fit$joint.work$fallback)
        expect_equal(fit$joint.work$direct.entries, 1)
        expect_equal(fit$solver.calls, entries)
        if (why == "objective") expect_identical(fit$joint.work$direct.status, "objective_increase")
    }
    expect_equal(fallback.calls, 2)
    reason <- "unexpected"
    expect_error(geosmooth:::.ptf.backfit(p, initial), "unexpected programming")
    expect_equal(fallback.calls, 2)
    p$control$joint.max.workspace.bytes <- 1
    entries <- 0L
    fit <- geosmooth:::.ptf.backfit(p, initial)
    expect_true(fit$joint.work$fallback); expect_equal(fit$joint.work$direct.entries, 0)
    expect_identical(fit$joint.work$failure.kind, "workspace")
    expect_equal(fit$solver.calls, entries)
})

test_that("joint infrastructure failures propagate and entered work stays observable", {
    for (where in c("native_begin", "native_end", "joint_reconstruction_begin")) {
        p <- joint_fixture_problem()
        events <- list()
        p$observer <- geosmooth:::.ptf.observer(function(e) {
            events[[length(events)+1L]] <<- e
            if (paste(e$kind, e$event, sep = "_") == where)
                stop(geosmooth:::.ptf.infrastructure.error("injected resource failure", resource = TRUE))
        })
        expect_error(geosmooth:::.ptf.backfit(p, geosmooth:::.ptf.initial(p, diag(2))),
                     class = "projection_resource_error")
        expect_equal(p$observer$native, if (where == "native_begin") 0 else 1)
        expect_false(any(vapply(events, function(e) e$kind == "joint_fallback", logical(1))))
    }
})
