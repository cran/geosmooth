.expect.synthetic.rejected <- function(object, pattern) {
  expect_error(validate.synthetic.dataset(object), pattern)
  expect_error(synthetic.dataset.checksum(object), pattern)
}

test_that("ordinary identity and RNG corruptions are rejected", {
  object <- materialize.synthetic(
    synthetic.registry.spec("G2"),
    n = 40L, seed = 17L, rng.policy = "named.stream.v1")

  bad <- object
  bad$rng.streams$sampling[1L] <- bad$rng.streams$sampling[1L] + 1L
  .expect.synthetic.rejected(bad, "rng.streams")

  bad <- object
  bad$parameters$resolved.seeds$seed <- 18L
  .expect.synthetic.rejected(bad, "resolved.seeds")

  bad <- object
  bad$seed <- 18L
  .expect.synthetic.rejected(bad, "dataset.id")

  bad <- object
  bad$dataset.id <- sub("__seed17", "__seed18", bad$dataset.id, fixed = TRUE)
  .expect.synthetic.rejected(bad, "dataset.id")

  bad <- object
  bad$rng.policy <- "legacy"
  .expect.synthetic.rejected(bad, "dataset.id")

  bad <- object
  bad$dataset.id <- "forged-short-id"
  attr(bad, "frozen.instance") <- TRUE
  .expect.synthetic.rejected(bad, "dataset.id")
})

test_that("dimension, frame, predictor, and support corruptions are rejected", {
  object <- materialize.synthetic(
    synthetic.registry.spec("G2"),
    n = 40L, seed = 17L, rng.policy = "named.stream.v1")

  bad <- object
  bad$intrinsic.dim <- bad$intrinsic.dim + 1L
  .expect.synthetic.rejected(bad, "dimension and region")

  bad <- object
  bad$parameters$frame.matrix[] <- 0
  .expect.synthetic.rejected(bad, "not orthonormal")

  bad <- object
  bad$predictors[1L, 1L] <- bad$predictors[1L, 1L] + 1
  .expect.synthetic.rejected(bad, "do not reconstruct predictors")

  disk <- materialize.synthetic(
    synthetic.registry.spec("G3a"),
    n = 40L, seed = 23L, rng.policy = "legacy")
  disk$latent[1L, ] <- 100
  .expect.synthetic.rejected(disk, "sampling support")
})

test_that("stratified dimension names and regions are contractual", {
  object <- materialize.synthetic(
    synthetic.registry.spec("G4"),
    n = 40L, seed = 19L, rng.policy = "legacy")

  bad <- object
  names(bad$intrinsic.dim.by.region) <- rev(
    names(bad$intrinsic.dim.by.region))
  .expect.synthetic.rejected(bad, "dimension and region")

  bad <- object
  bad$observed.regions <- rev(bad$observed.regions)
  .expect.synthetic.rejected(bad, "Region declarations")
})

test_that("scientific comparison includes exact identity fields", {
  x <- materialize.synthetic(
    synthetic.registry.spec("G2"),
    n = 40L, seed = 17L, rng.policy = "named.stream.v1")
  y <- materialize.synthetic(
    synthetic.registry.spec("G2"),
    n = 40L, seed = 18L, rng.policy = "named.stream.v1")
  comparison <- compare.synthetic.dataset(x, y)
  expect_false(comparison$equal)
  expect_true(all(c(
    "dataset.id", "seed", "rng.streams") %in% comparison$mismatches))

  other.spec <- synthetic.registry.spec("G2", list(sigma = 0.2))
  z <- materialize.synthetic(
    other.spec, n = 40L, seed = 17L, rng.policy = "named.stream.v1")
  comparison <- compare.synthetic.dataset(x, z)
  expect_false(comparison$equal)
  expect_true(all(c(
    "dataset.id", "specification.sha256",
    "response.spec") %in% comparison$mismatches))
})

test_that("random-frame equivalence still requires valid frames", {
  x <- materialize.synthetic(
    synthetic.registry.spec("G2"),
    n = 40L, seed = 17L, rng.policy = "named.stream.v1")
  y <- x
  y$parameters$frame.matrix[, 1L] <-
    -y$parameters$frame.matrix[, 1L]
  y$predictors <- geosmooth:::.embed.synthetic.geometry(
    y$geometry.spec, y$latent,
    frame.matrix = y$parameters$frame.matrix)
  expect_identical(validate.synthetic.dataset(y), y)
  expect_true(compare.synthetic.dataset(x, y)$equal)

  y$parameters$frame.matrix[1L, 1L] <- 0
  expect_error(compare.synthetic.dataset(x, y), "not orthonormal")
})

test_that("named atomic vectors have a literal canonical checksum", {
  fixture <- utils::read.csv(
    test_path("..", "fixtures", "canonical_serializer_v1.csv"),
    stringsAsFactors = FALSE)
  coefficients <- c(b11 = 3, b0 = 1, b1 = 2)
  expect_identical(
    geosmooth:::.synthetic.sha256(coefficients),
    fixture$expected.sha256[
      fixture$case.id == "named.polynomial.coefficients.v1"])

  make.spec <- function(value) {
    synthetic.spec(
      dgraphs::synthetic.quadform(1L, 1L),
      dgraphs::synthetic.sampling.uniform.interval(0, 1),
      synthetic.truth.polynomial(value),
      synthetic.response.gaussian(0))
  }
  expect_identical(
    make.spec(c(b0 = 1, b1 = 2, b11 = 3))$specification.sha256,
    make.spec(c(b11 = 3, b0 = 1, b1 = 2))$specification.sha256)

  contractual.a <- list(
    intrinsic.dim.by.region = c(first = 1, second = 2))
  contractual.b <- list(
    intrinsic.dim.by.region = c(second = 2, first = 1))
  expect_false(identical(
    geosmooth:::.synthetic.sha256(contractual.a),
    geosmooth:::.synthetic.sha256(contractual.b)))
})

test_that("canonical checksums ignore the producing R version header", {
  payload <- list(
    intrinsic.dim = 1L,
    ambient.dim = 1L,
    forms = list(),
    offset = 0)
  raw <- serialize(
    geosmooth:::.synthetic.canonicalize(payload), NULL,
    ascii = FALSE, xdr = TRUE, version = 3)
  oldrel.raw <- release.raw <- raw
  oldrel.raw[7:10] <- as.raw(c(0, 4, 5, 2))
  release.raw[7:10] <- as.raw(c(0, 4, 6, 1))

  expect_identical(
    geosmooth:::.synthetic.normalize.xdr.v3.header(oldrel.raw),
    geosmooth:::.synthetic.normalize.xdr.v3.header(release.raw))
  expect_identical(
    geosmooth:::.synthetic.xdr.v3(payload)[7:10],
    as.raw(c(0, 4, 7, 0)))
})
