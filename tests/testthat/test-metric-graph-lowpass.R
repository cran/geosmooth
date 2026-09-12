make_path_graph_lengths <- function(lengths) {
  n <- length(lengths) + 1L
  adj <- vector("list", n)
  w <- vector("list", n)
  add_edge <- function(i, j, ell) {
    adj[[i]] <<- c(adj[[i]], j)
    w[[i]] <<- c(w[[i]], ell)
    adj[[j]] <<- c(adj[[j]], i)
    w[[j]] <<- c(w[[j]], ell)
  }
  for (i in seq_along(lengths)) add_edge(i, i + 1L, lengths[i])
  list(adj.list = lapply(adj, as.integer), weight.list = lapply(w, as.double))
}

laplacian_dense <- function(op) {
  as.matrix(op$laplacian$matrix)
}

expect_metric_fit_finite <- function(fit) {
  expect_true(all(is.finite(fit$fitted.values)))
  expect_true(all(is.finite(fit$residuals)))
  expect_true(all(is.finite(fit$gcv$eta.optimal)))
  expect_true(all(is.finite(fit$gcv$gcv.scores)))
  expect_true(all(is.finite(fit$gcv$gcv.optimal)))
  expect_true(all(is.finite(fit$gcv$effective.df)))
  expect_true(all(is.finite(fit$spectral$eigenvalues)))
  expect_true(all(is.finite(fit$spectral$eigenvectors)))
  expect_true(all(is.finite(fit$spectral$filtered.eigenvalues)))
}

expect_metric_refit_finite <- function(refit) {
  expect_true(all(is.finite(refit$fitted.values)))
  expect_true(all(is.finite(refit$residuals)))
  if (isTRUE(refit$per.column.gcv)) {
    expect_true(all(is.finite(refit$eta.optimal)))
    expect_true(all(is.finite(refit$gcv.scores)))
    expect_true(all(is.finite(refit$effective.df)))
  } else {
    expect_true(all(is.finite(refit$eta.used)))
  }
}

test_that("metric operator validates undirected positive length graphs", {
  graph <- make_path_graph_lengths(c(1, 2))
  graph.bad <- graph
  graph.bad$adj.list[[2]] <- as.integer(1L)
  graph.bad$weight.list[[2]] <- 1

  expect_error(
    metric.graph.lowpass.operator(graph.bad$adj.list, graph.bad$weight.list),
    "no reciprocal entry"
  )

  graph.bad <- graph
  graph.bad$weight.list[[1]][1] <- -1
  expect_error(
    metric.graph.lowpass.operator(graph.bad$adj.list, graph.bad$weight.list),
    "non-positive"
  )
})

test_that("metric conductance transforms match hand calculations", {
  graph <- make_path_graph_lengths(c(2, 4))
  eps <- 1e-8

  op.inv <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "inverse.length.power",
    conductance.epsilon = eps,
    conductance.alpha = 2
  )
  expect_equal(op.inv$edge.table$conductance, (c(2, 4) + eps)^(-2), tolerance = 1e-12)

  op.exp <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "exp.length",
    conductance.sigma = 2
  )
  expect_equal(op.exp$edge.table$conductance, exp(-c(2, 4) / 2), tolerance = 1e-12)

  op.exp2 <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "exp.length.squared",
    conductance.sigma = 2
  )
  expect_equal(op.exp2$edge.table$conductance, exp(-(c(2, 4)^2) / 4), tolerance = 1e-12)
})

test_that("self-tuned Gaussian uses local incident edge scales", {
  graph <- make_path_graph_lengths(c(2, 4))
  op <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "self.tuned.gaussian",
    conductance.local.k = 1L
  )
  expect_equal(op$conductance$local.scales, c(2, 2, 4), tolerance = 1e-12)
  expect_equal(
    op$edge.table$conductance,
    c(exp(-(2^2) / (2 * 2)), exp(-(4^2) / (2 * 4))),
    tolerance = 1e-12
  )
})

test_that("unnormalized weighted Laplacian assembly matches R reference", {
  graph <- make_path_graph_lengths(c(2, 4))
  op <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "inverse.length.power",
    conductance.epsilon = 1e-8
  )
  c12 <- 1 / (2 + 1e-8)
  c23 <- 1 / (4 + 1e-8)
  L.ref <- matrix(c(
     c12, -c12,     0,
    -c12, c12+c23, -c23,
        0, -c23,    c23
  ), nrow = 3, byrow = TRUE)
  L <- laplacian_dense(op)
  expect_equal(L, L.ref, tolerance = 1e-12)
  expect_equal(rowSums(L), rep(0, 3), tolerance = 1e-12)
  expect_true(isTRUE(all.equal(L, t(L), tolerance = 1e-12)))
  expect_true(min(eigen(L, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
})

test_that("symmetric normalized Laplacian assembly matches R reference", {
  graph <- make_path_graph_lengths(c(2, 4))
  op <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    conductance.rule = "inverse.length.power",
    conductance.epsilon = 1e-8,
    laplacian.type = "symmetric.normalized"
  )
  c12 <- 1 / (2 + 1e-8)
  c23 <- 1 / (4 + 1e-8)
  degree <- c(c12, c12 + c23, c23)
  L.ref <- diag(1, 3)
  L.ref[1, 2] <- L.ref[2, 1] <- -c12 / sqrt(degree[1] * degree[2])
  L.ref[2, 3] <- L.ref[3, 2] <- -c23 / sqrt(degree[2] * degree[3])

  L <- laplacian_dense(op)
  expect_equal(L, L.ref, tolerance = 1e-12)
  expect_equal(as.vector(L %*% sqrt(degree)), rep(0, 3), tolerance = 1e-10)
  ev <- eigen(L, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(ev) > -1e-10)
  expect_true(max(ev) < 2 + 1e-10)
})

test_that("C++ dense spectral path matches R dense eigen reference", {
  graph <- make_path_graph_lengths(c(1, 2, 3, 4))
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list,
    y = sin(seq_len(5)),
    n.eigenpairs = 5L,
    eigen.solver = "dense",
    conductance.rule = "inverse.length.power",
    conductance.epsilon = 1e-8
  )
  L <- laplacian_dense(fit$operator)
  ev <- eigen(L, symmetric = TRUE)
  expect_equal(fit$spectral$eigenvalues, sort(ev$values), tolerance = 1e-9)

  U.ref <- ev$vectors[, order(ev$values), drop = FALSE]
  cors <- abs(diag(crossprod(fit$spectral$eigenvectors, U.ref)))
  expect_true(all(cors > 1 - 1e-7))
})

test_that("metric low-pass returns finite fit diagnostics across conductance rules", {
  graph <- make_path_graph_lengths(c(0.2, 0.5, 0.9, 1.7, 2.8, 4.6))
  y <- sin(seq_len(7) / 2) + cos(seq_len(7) / 5)
  eta.grid <- c(0.001, 0.01, 0.1, 1, 10)
  specs <- list(
    list(conductance.rule = "inverse.length.power", conductance.alpha = 1.5),
    list(conductance.rule = "exp.length", conductance.sigma.rule = "median"),
    list(conductance.rule = "exp.length.squared", conductance.sigma.rule = "edge.quantile"),
    list(conductance.rule = "self.tuned.gaussian", conductance.local.k = 2L)
  )

  for (spec in specs) {
    args <- c(
      list(
        adj.list = graph$adj.list,
        weight.list = graph$weight.list,
        y = y,
        n.eigenpairs = 7L,
        eigen.solver = "dense",
        eta.grid = eta.grid
      ),
      spec
    )
    fit <- do.call(fit.metric.graph.lowpass, args)
    expect_metric_fit_finite(fit)
  }
})

test_that("sparse spectral path handles inverse-square metric conductances", {
  lengths <- exp(seq(log(1e-4), log(2e-2), length.out = 79))
  graph <- make_path_graph_lengths(lengths)
  y <- sin(seq_len(80) / 7)
  eta.grid <- c(1e-8, 1e-6, 1e-4)

  fit.sparse <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    conductance.rule = "inverse.length.power",
    conductance.alpha = 2,
    n.eigenpairs = 20L,
    eigen.solver = "sparse",
    dense.fallback = "never",
    eta.grid = eta.grid
  )
  fit.dense <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    conductance.rule = "inverse.length.power",
    conductance.alpha = 2,
    n.eigenpairs = 20L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )

  expect_equal(fit.sparse$spectral$backend, "sparse.shift")
  expect_equal(fit.sparse$spectral$eigenvalues, fit.dense$spectral$eigenvalues, tolerance = 1e-6)
  expect_equal(fit.sparse$fitted.values, fit.dense$fitted.values, tolerance = 1e-8)
  expect_metric_fit_finite(fit.sparse)
  expect_metric_fit_finite(fit.dense)
})

test_that("metric low-pass finite guards hold under aggressive conductance scaling", {
  lengths <- exp(seq(log(1e-7), log(1e1), length.out = 99))
  graph <- make_path_graph_lengths(lengths)
  y <- sin(seq_len(100) / 4) + 0.25 * cos(seq_len(100) / 9)
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    conductance.rule = "inverse.length.power",
    conductance.alpha = 3,
    conductance.epsilon = 1e-12,
    n.eigenpairs = 25L,
    eigen.solver = "sparse",
    dense.fallback = "never",
    eta.grid = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  )
  expect_metric_fit_finite(fit)

  fixed.refit <- refit.metric.graph.lowpass(fit, rev(y))
  expect_metric_refit_finite(fixed.refit)

  Y <- cbind(reversed = rev(y), shifted = y + seq_along(y) / 100)
  gcv.refit <- refit.metric.graph.lowpass(
    fit, Y,
    per.column.gcv = TRUE,
    eta.grid = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  )
  expect_metric_refit_finite(gcv.refit)
})

test_that("symmetric normalized smoothing preserves sqrt-degree null vector", {
  graph <- make_path_graph_lengths(c(1, 2, 4))
  op <- metric.graph.lowpass.operator(
    graph$adj.list, graph$weight.list,
    laplacian.type = "symmetric.normalized"
  )
  y.null <- sqrt(op$degree)

  fit.null <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y.null,
    laplacian.type = "symmetric.normalized",
    n.eigenpairs = 4L,
    eigen.solver = "dense",
    eta.grid = 1e6
  )
  expect_equal(fit.null$fitted.values, y.null, tolerance = 1e-8)

  fit.constant <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, rep(1, 4),
    laplacian.type = "symmetric.normalized",
    n.eigenpairs = 4L,
    eigen.solver = "dense",
    eta.grid = 1e6
  )
  expect_gt(max(abs(fit.constant$fitted.values - 1)), 1e-3)
})

test_that("metric low-pass preserves constant responses and refit fixed eta", {
  graph <- make_path_graph_lengths(rep(1, 4))
  y <- rep(3, 5)
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    n.eigenpairs = 5L,
    eigen.solver = "dense",
    eta.grid = c(0.1, 1, 10)
  )
  expect_equal(fit$fitted.values, y, tolerance = 1e-10)

  y2 <- seq_len(5)
  refit <- refit.metric.graph.lowpass(fit, y2)
  V <- fit$spectral$eigenvectors
  f <- fit$spectral$filtered.eigenvalues
  y2.hat <- as.vector(V %*% (f * as.vector(crossprod(V, y2))))
  expect_equal(refit$fitted.values, y2.hat, tolerance = 1e-10)
  expect_metric_fit_finite(fit)
  expect_metric_refit_finite(refit)
})

test_that("effective degrees of freedom decrease with eta for heat kernel", {
  graph <- make_path_graph_lengths(rep(1, 5))
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list,
    y = sin(seq_len(6)),
    n.eigenpairs = 6L,
    eigen.solver = "dense",
    eta.grid = c(0.01, 0.1, 1, 10),
    filter.type = "heat_kernel"
  )
  weights <- compute.filter.weights.matrix(
    fit$spectral$eigenvalues,
    fit$gcv$eta.grid,
    fit$spectral$filter.type
  )
  expect_true(all(diff(colSums(weights)) <= 1e-12))
})

test_that("per-column GCV refit matches independent fits", {
  graph <- make_path_graph_lengths(rep(1, 7))
  Y <- cbind(
    sin(seq_len(8) / 2),
    cos(seq_len(8) / 3)
  )
  eta.grid <- c(0.01, 0.1, 1, 10)
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, Y[, 1],
    n.eigenpairs = 8L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )
  refit <- refit.metric.graph.lowpass(
    fit, Y,
    per.column.gcv = TRUE,
    eta.grid = eta.grid
  )
  fixed.refit <- refit.metric.graph.lowpass(fit, Y)

  fit1 <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, Y[, 1],
    n.eigenpairs = 8L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )
  fit2 <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, Y[, 2],
    n.eigenpairs = 8L,
    eigen.solver = "dense",
    eta.grid = eta.grid
  )

  expect_equal(refit$eta.optimal, c(fit1$gcv$eta.optimal, fit2$gcv$eta.optimal), tolerance = 1e-12)
  expect_equal(refit$fitted.values[, 1], fit1$fitted.values, tolerance = 1e-10)
  expect_equal(refit$fitted.values[, 2], fit2$fitted.values, tolerance = 1e-10)
  expect_metric_fit_finite(fit)
  expect_metric_refit_finite(refit)
  expect_metric_refit_finite(fixed.refit)
})

test_that("reusable basis preserves the established fit eigensystem", {
  graph <- make_path_graph_lengths(c(0.4, 0.7, 1.1, 1.6, 2.2))
  y <- sin(seq_len(6) / 3)
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    conductance.rule = "self.tuned.gaussian",
    conductance.local.k = 2L,
    n.eigenpairs = 6L,
    eigen.solver = "dense"
  )
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    conductance.rule = "self.tuned.gaussian",
    conductance.local.k = 2L,
    n.eigenpairs = 6L,
    eigen.solver = "dense",
    eta.grid = c(0.1, 1, 10)
  )

  expect_s3_class(basis, "metric.graph.lowpass.basis")
  expect_true(basis$spectral$is.complete)
  expect_equal(basis$spectral$eigenvalues, fit$spectral$eigenvalues,
               tolerance = 1e-12)
  expect_equal(
    abs(crossprod(basis$spectral$eigenvectors,
                  fit$spectral$eigenvectors)),
    diag(6),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(basis$laplacian$matrix),
    as.matrix(fit$laplacian$matrix),
    tolerance = 1e-12
  )
})

test_that("low-pass path reproduces the W1 spectral heat formula", {
  graph <- make_path_graph_lengths(c(0.3, 0.5, 0.8, 1.3, 2.1, 3.4))
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    conductance.rule = "self.tuned.gaussian",
    conductance.local.k = 2L,
    n.eigenpairs = 7L,
    eigen.solver = "dense"
  )
  y <- seq_len(7) / sum(seq_len(7))
  eta <- c(0, 0.05, 0.5, 5)
  path <- apply.metric.graph.lowpass.path(
    basis, y, eta,
    filter.type = "heat_kernel",
    unresolved.action = "error"
  )
  V <- basis$spectral$eigenvectors
  lambda <- basis$spectral$eigenvalues
  expected <- vapply(eta, function(value) {
    if (value == 0) return(y)
    as.numeric(V %*% (exp(-value * lambda) * crossprod(V, y)))
  }, numeric(length(y)))

  expect_s3_class(path, "metric.graph.lowpass.path")
  expect_equal(
    unname(path$fitted.values), unname(expected), tolerance = 1e-11
  )
  expect_true(all(path$resolution$spectrally.resolved))
  expect_true(path$resolution$exact.identity[[1L]])
  expect_equal(path$effective.df[[1L]], length(y))
})

test_that("low-pass path batches multiple responses without changing results", {
  graph <- make_path_graph_lengths(rep(1, 11))
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 12L,
    eigen.solver = "dense"
  )
  Y <- cbind(
    increasing = seq_len(12),
    oscillating = sin(seq_len(12)),
    decreasing = rev(seq_len(12))
  )
  eta <- c(0.01, 0.1, 1)
  batched <- apply.metric.graph.lowpass.path(
    basis, Y, eta, block.size = 2L
  )
  separate <- lapply(seq_len(ncol(Y)), function(j) {
    apply.metric.graph.lowpass.path(basis, Y[, j], eta)$fitted.values
  })

  expect_equal(dim(batched$fitted.values), c(12L, 3L, 3L))
  for (j in seq_along(separate)) {
    expect_equal(
      unname(batched$fitted.values[, , j]),
      unname(separate[[j]]),
      tolerance = 1e-11
    )
  }
})

test_that("heat-time grids distinguish W1 and guarded truncated rules", {
  graph <- make_path_graph_lengths(rep(1, 29))
  complete <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 30L,
    eigen.solver = "dense"
  )
  truncated <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 10L,
    eigen.solver = "dense"
  )
  w1 <- metric.graph.heat.eta.grid(
    complete, rule = "w1_inverse_spectrum", n.initial = 8L
  )
  guarded <- metric.graph.heat.eta.grid(
    truncated, rule = "spectral_guarded", n.initial = 8L,
    truncation.tol = 1e-6, equilibrium.tol = 1e-5
  )

  positive <- complete$spectral$eigenvalues[
    complete$spectral$eigenvalues > .Machine$double.eps
  ]
  expect_equal(unname(w1[[1L]]), 1 / max(positive), tolerance = 1e-12)
  expect_equal(unname(tail(w1, 1L)), 1 / min(positive), tolerance = 1e-10)
  expect_equal(attr(guarded, "rule"), "spectral_guarded")
  expect_false(attr(guarded, "basis.complete"))
  expect_error(
    metric.graph.heat.eta.grid(
      truncated, rule = "w1_inverse_spectrum", n.initial = 8L
    ),
    "requires a complete spectral basis"
  )
  expect_lte(
    exp(-guarded[[1L]] * attr(guarded, "lambda.cut")),
    1e-6 * (1 + 1e-10)
  )
})

test_that("guarded lower heat-time proposals are geometric and bounded", {
  graph <- make_path_graph_lengths(rep(1, 19))
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 20L,
    eigen.solver = "dense"
  )
  eta <- metric.graph.heat.eta.grid(
    basis, rule = "spectral_guarded", n.initial = 8L
  )
  extension <- metric.graph.heat.extend.lower(
    basis, eta,
    endpoint.status = "active",
    expansion.factor = 3,
    expansions.completed = 1L,
    identity.departure = 0.01
  )

  expect_s3_class(extension, "metric.graph.heat.lower.extension")
  expect_true(extension$admitted)
  expect_equal(extension$status, "admitted")
  expect_equal(extension$expansion.round, 2L)
  expect_equal(
    extension$proposed.eta,
    min(eta) / 3,
    tolerance = 1e-14
  )
  expect_equal(min(extension$eta.grid), extension$proposed.eta)
  expect_false(any(extension$eta.grid == 0))
  expect_true(extension$basis.complete)
})

test_that("guarded lower heat-time proposals respect stop states and identity floor", {
  graph <- make_path_graph_lengths(rep(1, 19))
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 20L,
    eigen.solver = "dense"
  )
  eta <- c(0.1, 0.3, 1)

  inactive <- metric.graph.heat.extend.lower(
    basis, eta, endpoint.status = "inactive"
  )
  capped <- metric.graph.heat.extend.lower(
    basis, eta, endpoint.status = "active",
    max.expansions = 2L, expansions.completed = 2L
  )
  floor.stop <- metric.graph.heat.extend.lower(
    basis, eta, endpoint.status = "active",
    identity.departure = 0.9
  )

  expect_false(inactive$admitted)
  expect_equal(inactive$status, "inactive_endpoint")
  expect_false(capped$admitted)
  expect_equal(capped$status, "expansion_cap_reached")
  expect_false(floor.stop$admitted)
  expect_equal(floor.stop$status, "identity_floor_reached")
  expect_equal(floor.stop$eta.grid, sort(eta))
})

test_that("guarded lower heat-time proposals reject unresolved truncation", {
  graph <- make_path_graph_lengths(rep(1, 39))
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 12L,
    eigen.solver = "dense"
  )
  eta <- c(1e-4, 1e-3, 1e-2)

  marked <- metric.graph.heat.extend.lower(
    basis, eta,
    endpoint.status = "active",
    identity.departure = 1e-6,
    truncation.tol = 1e-4,
    unresolved.action = "mark"
  )
  expect_false(marked$admitted)
  expect_equal(marked$status, "unresolved_truncated_basis")
  expect_false(marked$basis.complete)
  expect_gt(marked$omitted.attenuation.bound, marked$truncation.tol)
  expect_error(
    metric.graph.heat.extend.lower(
      basis, eta,
      endpoint.status = "active",
      identity.departure = 1e-6,
      truncation.tol = 1e-4,
      unresolved.action = "error"
    ),
    "not resolved by the truncated basis"
  )
})

test_that("main graph low-pass fit keeps fixed eta search as the default", {
  graph <- make_path_graph_lengths(rep(1, 19))
  y <- c(rep(0, 10), rep(1, 10))
  default <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    n.eigenpairs = 20L,
    eigen.solver = "dense",
    eta.grid = c(0.01, 1, 10)
  )
  explicit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    n.eigenpairs = 20L,
    eigen.solver = "dense",
    eta.grid = c(0.01, 1, 10),
    eta.search = "fixed"
  )

  expect_equal(default$fitted.values, explicit$fitted.values, tolerance = 0)
  expect_equal(default$gcv$eta.grid, c(0.01, 1, 10))
  expect_equal(default$gcv$search$mode, "fixed")
  expect_equal(default$gcv$search$status, "fixed_grid")
  expect_equal(default$gcv$search$expansions.completed, 0L)
})

test_that("main heat fit exposes guarded lower-time GCV search", {
  graph <- make_path_graph_lengths(rep(1, 19))
  y <- c(rep(0, 10), rep(1, 10))
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, y,
    n.eigenpairs = 20L,
    eigen.solver = "dense",
    eta.grid = c(0.01, 1, 10),
    eta.search = "guarded.gcv",
    eta.max.expansions = 3L
  )

  expect_equal(fit$gcv$search$mode, "guarded.gcv")
  expect_equal(fit$gcv$search$status, "interior_optimum")
  expect_equal(fit$gcv$search$expansions.completed, 1L)
  expect_equal(min(fit$gcv$eta.grid), 0.01 / 3, tolerance = 1e-14)
  expect_equal(nrow(fit$gcv$search$rounds), 2L)
  expect_true(fit$gcv$search$rounds$lower.endpoint.selected[[1L]])
  expect_false(fit$gcv$search$rounds$lower.endpoint.selected[[2L]])
  expect_equal(fit$gcv$search$rounds$extension.status[[1L]], "admitted")
  expect_equal(fit$parameters$eta.search, "guarded.gcv")
})

test_that("main guarded heat fit reports the identity floor and round cap", {
  graph <- make_path_graph_lengths(rep(1, 19))
  fit <- fit.metric.graph.lowpass(
    graph$adj.list, graph$weight.list, seq_len(20),
    n.eigenpairs = 20L,
    eigen.solver = "dense",
    eta.grid = c(0.01, 1, 10),
    eta.search = "guarded.gcv",
    eta.max.expansions = 2L,
    eta.identity.departure = 0.01
  )

  expect_equal(fit$gcv$search$status, "expansion_cap_reached")
  expect_equal(fit$gcv$search$expansions.completed, 2L)
  expect_equal(
    fit$gcv$search$rounds$extension.status,
    c("admitted", "admitted_at_identity_floor", "expansion_cap_reached")
  )
  expect_equal(nrow(fit$gcv$search$rounds), 3L)
})

test_that("main guarded eta search is heat-only and excludes identity", {
  graph <- make_path_graph_lengths(rep(1, 9))
  y <- seq_len(10)

  expect_error(
    fit.metric.graph.lowpass(
      graph$adj.list, graph$weight.list, y,
      n.eigenpairs = 10L,
      eigen.solver = "dense",
      filter.type = "tikhonov",
      eta.grid = c(0.1, 1),
      eta.search = "guarded.gcv"
    ),
    "requires filter.type = 'heat_kernel'"
  )
  expect_error(
    fit.metric.graph.lowpass(
      graph$adj.list, graph$weight.list, y,
      n.eigenpairs = 10L,
      eigen.solver = "dense",
      eta.grid = c(0, 0.1, 1),
      eta.search = "guarded.gcv"
    ),
    "positive competitive times"
  )
})

test_that("truncated paths guard weak smoothing and preserve exact eta zero", {
  graph <- make_path_graph_lengths(rep(1, 39))
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 12L,
    eigen.solver = "dense"
  )
  y <- sin(seq_len(40) / 3)

  expect_warning(
    path <- apply.metric.graph.lowpass.path(
      basis, y, c(0, 1e-8, 100),
      truncation.tol = 1e-4,
      unresolved.action = "warn"
    ),
    "exceed the truncated-basis attenuation tolerance"
  )
  expect_equal(path$fitted.values[, 1L], y, tolerance = 0)
  expect_true(path$resolution$exact.identity[[1L]])
  expect_false(path$resolution$spectrally.resolved[[2L]])
  expect_error(
    apply.metric.graph.lowpass.path(
      basis, y, 1e-8,
      truncation.tol = 1e-4,
      unresolved.action = "error"
    ),
    "increase n.eigenpairs"
  )
})

test_that("Butterworth unresolved advice does not reverse eta direction", {
  graph <- make_path_graph_lengths(rep(1, 39))
  basis <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 12L,
    eigen.solver = "dense"
  )
  y <- sin(seq_len(40) / 3)

  expect_warning(
    apply.metric.graph.lowpass.path(
      basis, y, eta.grid = 100,
      filter.type = "butterworth",
      truncation.tol = 1e-4,
      unresolved.action = "warn"
    ),
    "choose a filter parameter that more strongly attenuates omitted modes"
  )
})

test_that("guarded first-200-mode heat agrees with the complete heat path", {
  graph <- make_path_graph_lengths(rep(1, 239))
  full <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 240L,
    eigen.solver = "dense"
  )
  truncated <- metric.graph.lowpass.basis(
    graph$adj.list, graph$weight.list,
    n.eigenpairs = 200L,
    eigen.solver = "dense"
  )
  set.seed(20260724)
  y <- stats::rbinom(240L, 1L, 0.15)
  eta <- metric.graph.heat.eta.grid(
    truncated,
    rule = "spectral_guarded",
    n.initial = 6L,
    truncation.tol = 1e-6,
    equilibrium.tol = 1e-5
  )
  full.path <- apply.metric.graph.lowpass.path(full, y, eta)
  truncated.path <- apply.metric.graph.lowpass.path(
    truncated, y, eta,
    truncation.tol = 1e-6,
    unresolved.action = "error"
  )
  difference <- full.path$fitted.values - truncated.path$fitted.values

  expect_lt(sqrt(mean(difference^2)), 1e-7)
  expect_lt(max(abs(difference)), 2e-6)
})
