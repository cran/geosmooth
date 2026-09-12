test_that("synthetic specifications are draw-free and canonically hashable", {
  geometry <- synthetic.quadform(
    intrinsic.dim = 2, ambient.dim = 5,
    forms = list(diag(c(0.5, -0.25))),
    frame = "supplied",
    frame.matrix = rbind(diag(3), matrix(0, nrow = 2, ncol = 3)))
  spec <- synthetic.spec(
    geometry,
    synthetic.sampling.uniform.box(-1, 1),
    synthetic.truth.polynomial(c(b0 = 1, b1 = 2, b2 = -1)),
    synthetic.response.gaussian(0.1))
  expect_s3_class(spec, "synthetic_spec")
  expect_match(spec$specification.sha256, "^[0-9a-f]{64}$")
  expect_identical(
    spec$specification.sha256,
    synthetic.spec(
      geometry,
      synthetic.sampling.uniform.box(-1, 1),
      synthetic.truth.polynomial(c(b0 = 1, b1 = 2, b2 = -1)),
      synthetic.response.gaussian(0.1))$specification.sha256)
  expect_error(
    synthetic.spec(
      geometry,
      synthetic.sampling.uniform.box(-1, 1),
      structure(
        list(kind = "truth", family = "bad", version = 1L,
             parameters = list(fun = function(x) x)),
        class = c("synthetic_truth", "synthetic_component")),
      synthetic.response.gaussian(0.1)),
    "cannot contain functions")
})

test_that("materialized datasets enforce identity and dimension contracts", {
  spec <- synthetic.registry.spec("G3a", list(R = 2, truth = "linear"))
  ds <- materialize.synthetic(
    spec, n = 80, seed = 17, rng.policy = "named.stream.v1")
  expect_s3_class(ds, "synthetic_dataset")
  expect_identical(ds$intrinsic.dim, 2L)
  expect_identical(ds$ambient.dim, 3L)
  expect_identical(ds$codimension, 1L)
  expect_identical(nrow(ds$predictors), 80L)
  expect_match(ds$dataset.id, spec$specification.sha256, fixed = TRUE)
  expect_identical(validate.synthetic.dataset(ds), ds)
  expect_match(synthetic.dataset.checksum(ds), "^[0-9a-f]{64}$")

  altered <- ds
  altered$predictors[1, 1] <- altered$predictors[1, 1] + 1e-3
  expect_error(
    compare.synthetic.dataset(ds, altered),
    "do not reconstruct predictors")
  expect_error(
    synthetic.dataset.checksum(altered),
    "do not reconstruct predictors")
})

test_that("nullable stratified latent coordinates compare safely", {
  spec <- synthetic.registry.spec("G4")
  x <- materialize.synthetic(spec, n = 50, seed = 9, rng.policy = "legacy")
  y <- materialize.synthetic(spec, n = 50, seed = 9, rng.policy = "legacy")
  expect_true(anyNA(x$latent))
  expect_true(compare.synthetic.dataset(x, y)$equal)
})

test_that("named streams preserve caller RNG and isolate sampling from response", {
  set.seed(901)
  before <- .Random.seed
  base <- synthetic.registry.spec("G3a", list(sigma = 0))
  noisy <- synthetic.registry.spec("G3a", list(sigma = 0.7))
  x <- materialize.synthetic(
    base, n = 60, seed = 31, rng.policy = "named.stream.v1")
  expect_identical(.Random.seed, before)
  y <- materialize.synthetic(
    noisy, n = 60, seed = 31, rng.policy = "named.stream.v1")
  expect_identical(.Random.seed, before)
  expect_identical(x$latent, y$latent)
  expect_identical(x$predictors, y$predictors)
  expect_false(identical(x$response, y$response))
})

test_that("legacy derived seeds are preflighted before RNG state changes", {
  set.seed(902)
  before <- .Random.seed
  expect_error(
    materialize.synthetic(
      synthetic.registry.spec("G2"), n = 20,
      seed = .Machine$integer.max, rng.policy = "legacy"),
    "derived frame seed")
  expect_identical(.Random.seed, before)
  expect_error(
    materialize.synthetic(
      synthetic.registry.spec("G3a"), n = 20,
      seed = .Machine$integer.max, rng.policy = "legacy"),
    "derived response seed")
  expect_identical(.Random.seed, before)

  no.response.draw <- synthetic.registry.spec("G1", list(sigma = 0))
  expect_s3_class(
    materialize.synthetic(
      no.response.draw, n = 10, seed = .Machine$integer.max,
      rng.policy = "legacy"),
    "synthetic_dataset")
  expect_identical(.Random.seed, before)
})

test_that("G7 rejects an empty structural-zero subset explicitly", {
  condition <- tryCatch(
    synthetic.sampling.dirichlet.zeros(
      concentration = 1, zero.fraction = 0.5, zero.parts = integer()),
    error = identity)
  expect_s3_class(condition, "geosmooth_empty_zero_parts")
  expect_match(conditionMessage(condition), "at least one")
})

test_that("canonical frozen instances replay their committed checksums", {
  path <- geosmooth:::.synthetic.instance.path()
  expect_true(nzchar(path) && file.exists(path))
  registry <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_setequal(
    c("instance.id", "recipe.id", "n", "seed", "rng.policy",
      "checksum.id", "fixture.path", "status", "version",
      "content.sha256"),
    names(registry))
  expect_setequal(
    registry$recipe.id,
    synthetic.registry.ids()[grepl("^G", synthetic.registry.ids())])
  for (instance.id in registry$instance.id) {
    object <- materialize.synthetic.instance(instance.id)
    expected <- registry$content.sha256[
      registry$instance.id == instance.id]
    expect_identical(object$dataset.id, instance.id)
    verification.scope <- attr(
      object, "verification.scope", exact = TRUE)
    if (identical(verification.scope, "exact-environment")) {
      expect_identical(
        synthetic.dataset.checksum(object), expected,
        info = instance.id)
    } else {
      expect_identical(
        verification.scope, "scientific-parity", info = instance.id)
      fixture.path <- geosmooth:::.synthetic.registry.fixture.path(
        registry$fixture.path[registry$instance.id == instance.id])
      expect_true(
        compare.synthetic.dataset(object, readRDS(fixture.path))$equal,
        info = instance.id)
    }
  }
})
