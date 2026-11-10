test_that("occupation mixtures distinguish density from Bernoulli truth", {
  truth <- synthetic.truth.occupation.mixture(
    centers = rbind(c(-0.4, 0), c(0.5, 0.2)),
    covariances = list(diag(c(0.12, 0.2)), diag(c(0.18, 0.1))),
    weights = c(0.55, 0.45),
    gamma = 1.2,
    probability.maximum = 0.65)
  spec <- synthetic.spec(
    dgraphs::synthetic.quadform(2, 2, frame = "canonical"),
    dgraphs::synthetic.sampling.uniform.box(-1, 1),
    truth,
    synthetic.response.bernoulli(
      minimum.positive = 10, maximum.attempts = 5),
    recipe.id = "occupation-mixture-test")
  object <- materialize.synthetic(
    spec, n = 200, seed = 17, rng.policy = "named.stream.v1")
  expect_identical(object$truth.scope, "finite.design")
  expect_equal(max(object$truth), 0.65, tolerance = 1e-14)
  expect_true(all(object$truth >= 0 & object$truth <= 0.65))
  expect_true(all(object$response %in% c(0, 1)))
  expect_gte(sum(object$response), 10)
  density <- object$parameters$truth$analytic.density
  expect_length(density, object$n)
  expect_true(all(density > 0))
  expect_false(isTRUE(all.equal(density, object$truth)))
})

test_that("occupation-mixture component validates normalization inputs", {
  centers <- matrix(c(0, 0), nrow = 1)
  covariance <- list(diag(2))
  expect_error(
    synthetic.truth.occupation.mixture(
      centers, covariance, 1, probability.maximum = 1.1),
    "must not exceed")
  expect_error(
    synthetic.truth.occupation.mixture(
      centers, covariance, 1, normalization = "population"),
    "Only design.maximum")
})
