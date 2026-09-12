test_that("circle grid absorbs the historical noise-free convention", {
  n <- 9L
  radius <- 2.5
  angles <- seq(0, 2 * pi, length.out = n + 1L)[-1L]
  expected <- radius * cbind(cos(angles), sin(angles))

  spec <- synthetic.spec(
    geometry = synthetic.circle(radius = radius),
    sampling = synthetic.sampling.grid.interval(
      0, 2 * pi, endpoints = "exclude.lower"),
    truth = synthetic.truth.polynomial(c(b0 = 0)),
    response = synthetic.response.gaussian(0))
  observed <- materialize.synthetic(
    spec, n = n, seed = 17L, rng.policy = "legacy")

  expect_equal(observed$latent[, 1L], angles, tolerance = 0)
  expect_equal(observed$predictors, expected, tolerance = 0)
  expect_equal(observed$truth, rep(0, n), tolerance = 0)
})

test_that("random circle angles preserve the legacy base-R draw", {
  set.seed(73L)
  expected <- sort(stats::runif(12L, 0, 2 * pi))
  spec <- synthetic.spec(
    geometry = synthetic.circle(),
    sampling = synthetic.sampling.uniform.interval(
      0, 2 * pi, order = "ascending"),
    truth = synthetic.truth.polynomial(c(b0 = 0)),
    response = synthetic.response.gaussian(0))
  observed <- materialize.synthetic(
    spec, n = 12L, seed = 73L, rng.policy = "legacy")
  expect_equal(observed$latent[, 1L], expected, tolerance = 0)
})

test_that("trefoil grid absorbs the historical parameterization", {
  n <- 11L
  scale <- 1.3
  tt <- seq(0, 2 * pi, length.out = n)
  expected <- cbind(
    scale * sin(tt) + 2 * sin(2 * tt),
    scale * cos(tt) - 2 * cos(2 * tt),
    -scale * sin(3 * tt))
  spec <- synthetic.spec(
    geometry = synthetic.trefoil(scale = scale),
    sampling = synthetic.sampling.grid.interval(
      0, 2 * pi, endpoints = "include.both"),
    truth = synthetic.truth.polynomial(c(b0 = 0)),
    response = synthetic.response.gaussian(0))
  observed <- materialize.synthetic(
    spec, n = n, seed = 1L, rng.policy = "legacy")
  expect_equal(observed$predictors, expected, tolerance = 0)
})

test_that("curve domains and ambient frames are validated", {
  expect_error(
    embed.synthetic.geometry(synthetic.circle(), matrix(2 * pi + 1e-4)),
    "outside")
  geometry <- synthetic.circle(
    ambient.dim = 5L, frame = "random.orthonormal")
  spec <- synthetic.spec(
    geometry = geometry,
    sampling = synthetic.sampling.grid.interval(0, 2 * pi),
    truth = synthetic.truth.polynomial(c(b0 = 0)),
    response = synthetic.response.gaussian(0))
  data <- materialize.synthetic(spec, n = 8L, seed = 9L)
  expect_equal(dim(data$predictors), c(8L, 5L))
  expect_equal(
    sqrt(rowSums(
      sweep(data$predictors, 2L, geometry$parameters$offset)^2)),
    rep(1, 8L),
    tolerance = 1e-12)
})

test_that("synthetic dataset plot method is shared across families", {
  spec <- synthetic.spec(
    geometry = synthetic.circle(),
    sampling = synthetic.sampling.grid.interval(0, 2 * pi),
    truth = synthetic.truth.polynomial(c(b0 = 0)),
    response = synthetic.response.gaussian(0))
  data <- materialize.synthetic(spec, n = 8L, seed = 1L)
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  expect_invisible(plot(data))
})
