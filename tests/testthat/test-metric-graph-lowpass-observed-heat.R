heat_test_graph <- function(n = 5L) {
    adj <- lapply(seq_len(n), function(i) as.integer(c(i - 1L, i + 1L)[c(i > 1L, i < n)]))
    list(adj = adj, w = lapply(adj, function(x) rep(1, length(x))))
}
heat_test_basis <- function(n = 5L) {
    g <- heat_test_graph(n)
    metric.graph.lowpass.basis(g$adj, g$w, weight.type = "conductance", spectrum = "full")
}
heat_test_fit <- function(y, times = c(0.01, 0.1, 1), ...) {
    g <- heat_test_graph(length(y))
    fit.metric.graph.lowpass(g$adj, g$w, y, weight.type = "conductance",
                            spectrum = "full", eta.grid = times, ...)
}

# An independent definition: exponentiate +tL, add the observation mask,
# and solve the penalized least-squares normal equations.
heat_penalty_reference <- function(L, y, t) {
    ev <- eigen(L, symmetric = TRUE)
    Q <- ev$vectors %*% diag(expm1(t * pmax(ev$values, 0))) %*% t(ev$vectors)
    O <- which(!is.na(y))
    R <- diag(length(y))[O, , drop = FALSE]
    B <- solve(crossprod(R) + Q, t(R))
    S <- R %*% B
    f <- as.vector(B %*% y[O])
    q <- length(O)
    list(f = f, df = sum(diag(S)),
         gcv = q * sum((y[O] - f[O])^2) / (q - sum(diag(S)))^2)
}

test_that("direct conductances retain magnitudes, zeros and metadata", {
    adj <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
    w <- list(c(2, 0), c(2, 3), c(0, 3))
    op <- metric.graph.lowpass.operator(adj, w, weight.type = "conductance")
    expect_equal(as.matrix(op$laplacian$matrix),
                 matrix(c(2, -2, 0, -2, 5, -3, 0, -3, 3), 3))
    expect_equal(op$weight.type, "conductance")
    expect_equal(op$conductance$rule, "identity")
    expect_true(all(is.na(op$edge.table$length)))
    expect_equal(op$graph$adj.list, list(2L, c(1L, 3L), 2L))
    tiny <- metric.graph.lowpass.operator(list(2L, 1L), list(1e-310, 1e-310),
                                          weight.type = "conductance")
    expect_identical(tiny$conductance$values, 1e-310)
    expect_equal(tiny$degree, rep(1e-310, 2))
    w[[1]][1] <- 1e-20
    expect_error(metric.graph.lowpass.operator(adj, w, weight.type = "conductance"), "mismatch|agree")
    expect_error(metric.graph.lowpass.operator(list(2L, 1L), list(0, 1e-20),
                                               weight.type = "conductance"), "agree")
    expect_error(metric.graph.lowpass.operator(list(1L, integer()), list(1, numeric()),
                                               weight.type = "conductance"), "self-loops")
    expect_error(metric.graph.lowpass.operator(list(2L, 1L), list(-1, -1),
                                               weight.type = "conductance"), "nonnegative")
})

test_that("full spectrum overrides auto threshold and verifies all eigenpairs", {
    g <- heat_test_graph(205L)
    b <- metric.graph.lowpass.basis(g$adj, g$w, spectrum = "full",
                                    dense.eigen.threshold = 10L, n.eigenpairs = 2L)
    expect_equal(dim(b$spectral$eigenvectors), c(205L, 205L))
    expect_equal(b$spectral$backend, "dense")
    expect_true(b$spectral$is.complete)
    expect_lt(b$spectral$verification$orthogonality.error, 1e-8)
    expect_lt(b$spectral$verification$relative.residual, 1e-8)
    expect_error(metric.graph.lowpass.basis(g$adj, g$w, spectrum = "full",
                                             eigen.solver = "sparse"), "requires auto or dense")
    expect_error(metric.graph.lowpass.basis(g$adj, g$w, n.eigenpairs = 205L,
                                             eigen.solver = "sparse"), "n-2")
    automatic <- metric.graph.lowpass.basis(g$adj, g$w, n.eigenpairs = 205L,
                                           dense.eigen.threshold = 10L)
    expect_true(automatic$spectral$is.complete)
})

test_that("fully observed heat GCV matches conventional dense calculation", {
    y <- c(1, 0, 4, -2, 3)
    fit <- heat_test_fit(y, c(0, 0.04, 0.2, 1))
    L <- as.matrix(fit$laplacian$matrix)
    for (j in 2:4) {
        ref <- heat_penalty_reference(L, y, fit$gcv$eta.grid[j])
        expect_equal(fit$gcv$gcv.scores[j], ref$gcv, tolerance = 1e-10)
        expect_equal(fit$gcv$effective.df.grid[j], ref$df, tolerance = 1e-10)
    }
    expect_true(is.na(fit$gcv$gcv.scores[1]))
    expect_equal(fit$gcv$status[1], "identity_undefined_gcv")
    expect_identical(fit$gcv$identity$fitted.values, as.double(y))
    expect_gt(fit$gcv$eta.optimal, 0)
    expect_error(heat_test_fit(y, 0), "No valid GCV")
    tiny <- heat_test_fit(y, c(1e-200, 1e-100, 1e-20))
    limit <- length(y) * sum((L %*% y)^2) / sum(diag(L))^2
    expect_equal(tiny$gcv$gcv.scores, rep(limit, 3), tolerance = 1e-10)
    expect_equal(heat_test_fit(rep(7, 5))$fitted.values, rep(7, 5), tolerance = 1e-12)
})

test_that("missing heat fit matches penalized normal equations and observed GCV", {
    y <- c(1, NA, 4, NA, 2)
    fit <- heat_test_fit(y, c(0, 0.04, 0.2, 1))
    L <- as.matrix(fit$laplacian$matrix)
    b <- heat_test_basis()
    path <- apply.metric.graph.lowpass.path(b, y, c(0, 0.04, 0.2, 1))
    for (j in 2:4) {
        ref <- heat_penalty_reference(L, y, fit$gcv$eta.grid[j])
        expect_equal(path$fitted.values[, j], ref$f, tolerance = 1e-10)
        expect_equal(fit$gcv$gcv.scores[j], ref$gcv, tolerance = 1e-10)
        expect_equal(fit$gcv$effective.df.grid[j], ref$df, tolerance = 1e-10)
    }
    expect_equal(fit$gcv$n.observed, 3L)
    expect_true(all(is.na(fit$residuals[c(2, 4)])))
    expect_true(all(is.finite(fit$fitted.values)))
    expect_equal(path$fitted.values[, 1], y)
    expect_equal(heat_test_fit(c(7, NA, 7, NA, 7))$fitted.values, rep(7, 5), tolerance = 1e-12)
    tiny <- heat_test_fit(y, c(1e-200, 1e-100, 1e-20))
    O <- which(!is.na(y)); U <- which(is.na(y))
    K <- L[O, O] - L[O, U] %*% solve(L[U, U], L[U, O])
    limit <- length(O) * sum((K %*% y[O])^2) / sum(diag(K))^2
    expect_equal(tiny$gcv$gcv.scores, rep(limit, 3), tolerance = 1e-10)
    expect_equal(tiny$fitted.values[U], as.vector(-solve(L[U, U], L[U, O]) %*% y[O]), tolerance = 1e-10)
    huge <- heat_test_fit(y, 1e6)
    expect_equal(huge$fitted.values, rep(mean(y[O]), 5), tolerance = 1e-10)
})

test_that("each feature mask and refit uses its own observed fitting operator", {
    Y <- cbind(first = c(1, NA, 4, NA, 2), second = c(NA, 3, 2, 1, NA), full = 1:5)
    fit <- heat_test_fit(Y[, 1])
    refitted <- refit(fit, Y, per.column.gcv = TRUE, eta.grid = c(.01, .1, 1))
    for (j in 1:3) {
        separate <- heat_test_fit(Y[, j])
        expect_equal(refitted$fitted.values[, j], separate$fitted.values, ignore_attr = TRUE)
        expect_equal(refitted$gcv.scores[j], separate$gcv$gcv.optimal, ignore_attr = TRUE)
    }
    expect_equal(unname(refitted$n.observed), c(3, 3, 5))
    fixed <- refit(fit, Y)
    for (j in 1:3) {
        ref <- heat_penalty_reference(as.matrix(fit$laplacian$matrix), Y[, j], fit$spectral$eta.optimal)
        expect_equal(fixed$fitted.values[, j], ref$f, ignore_attr = TRUE, tolerance = 1e-10)
    }
    path <- apply.metric.graph.lowpass.path(heat_test_basis(), Y, c(.01, .1, 1))
    expect_equal(dim(path$fitted.values), c(5, 3, 3))
    expect_equal(path$fitted.values[, , 3],
                 apply.metric.graph.lowpass.path(heat_test_basis(), Y[, 3], c(.01, .1, 1))$fitted.values,
                 ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("unidentified components and unsupported missing models fail explicitly", {
    adj <- list(2L, 1L, 4L, 3L, integer())
    weights <- list(1, 1, 1, 1, numeric())
    b <- metric.graph.lowpass.basis(adj, weights, weight.type = "conductance", spectrum = "full")
    expect_equal(b$spectral$nullity, 3L)
    expect_error(apply.metric.graph.lowpass.path(b, c(1, NA, NA, NA, 4), .1), "component needs")
    p <- apply.metric.graph.lowpass.path(b, c(1, NA, 3, NA, 4), .1)
    expect_equal(as.vector(p$fitted.values), c(1, 1, 3, 3, 4), tolerance = 1e-12)
    expect_true(is.na(p$gcv.scores))
    expect_error(heat_test_fit(c(1, NA, NA, NA, NA)), "No valid GCV")
    g <- heat_test_graph()
    partial <- metric.graph.lowpass.basis(g$adj, g$w, n.eigenpairs = 2L)
    expect_error(apply.metric.graph.lowpass.path(partial, c(1, NA, 2, 3, 4), .1), "complete spectrum")
    expect_error(heat_test_fit(c(1, NA, 2, 3, 4), filter.type = "tikhonov"), "heat_kernel")
    expect_error(heat_test_fit(c(1, NaN, 2, 3, 4)), "NaN/Inf")
    expect_error(heat_test_fit(rep(NA_real_, 5)), "observed value")
})


test_that("normalized direct conductances avoid degree-product underflow and overflow", {
    for (w in c(1e-310, 1e300)) {
        op <- metric.graph.lowpass.operator(list(2L, 1L), list(w, w),
                                            weight.type = "conductance",
                                            laplacian.type = "symmetric.normalized")
        expect_equal(as.matrix(op$laplacian$matrix), matrix(c(1, -1, -1, 1), 2),
                     tolerance = 1e-12)
    }
})
