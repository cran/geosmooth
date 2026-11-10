test_that("GRIP and edge-KK adapters execute the requested current backends", {
  skip_if_not_installed("grip")
  subgraph <- list(adj.list = list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L),
                   weight.list = list(1, c(1, 2), c(2, 1), c(1, 3), 3))
  distances <- as.matrix(dist(c(0, 1, 3, 4, 7)))
  fit <- .transported.graph.hessian.grip.edge.kk.embedding(subgraph, distances, 2L)
  expect_identical(fit$backend.used, "grip+edge.kk")
  expect_true(is.na(fit$fallback.reason))
  expect_equal(dim(fit$coordinates), c(5L, 2L))
  expect_true(all(is.finite(fit$coordinates)))
  polished <- .transported.graph.hessian.mds.edge.kk.embedding(
    subgraph, distances, 2L, fallback.reason = NULL)
  expect_identical(polished$backend.used, "cmdscale+edge.kk")
  expect_true(is.na(polished$fallback.reason))
  expect_lt(polished$edge.stress, 1e-6)
})

test_that("an incompatible edge-KK signature cannot drop intended settings", {
  old <- list(name = "edge.kk", fun = function(coords, adj_list, weight_list) {
    stop("The incompatible function should not execute")
  })
  expect_error(.transported.graph.hessian.run.edge.kk(old, matrix(0, 2, 2),
    list(adj.list = list(2L, 1L), weight.list = list(1, 1)), 2L),
    "Incompatible grip edge.kk.*adj.list")
})
