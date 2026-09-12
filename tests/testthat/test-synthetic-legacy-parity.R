.legacy.parity.fixture.path <- function() {
  dev <- tryCatch(
    testthat::test_path("..", "fixtures", "g1_g7_legacy_parity_v1.rds"),
    error = function(e) "")
  if (nzchar(dev) && file.exists(dev)) return(dev)
  installed <- system.file(
    "testdata", "g1_g7_legacy_parity_v1.rds", package = "geosmooth")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  ""
}

.materialize.legacy.parity.case <- function(fixture) {
  recipe.id <- sub("^dgp\\.", "", fixture$fun)
  recipe.id <- switch(
    recipe.id,
    g1 = "G1", g2 = "G2", g3a = "G3a", g3b = "G3b",
    g3c = "G3c", g3d = "G3d", g4 = "G4", g5 = "G5",
    g6 = "G6", g7 = "G7")
  parameters <- fixture$args
  n <- parameters$n
  seed <- parameters$seed
  parameters$n <- NULL
  parameters$seed <- NULL
  if (identical(recipe.id, "G3c") && !is.null(parameters$c)) {
    parameters$pitch <- parameters$c
    parameters$c <- NULL
  }
  spec <- synthetic.registry.spec(recipe.id, parameters)
  materialize.synthetic(
    spec, n = n, seed = seed, rng.policy = "legacy")
}

.canonical.legacy.parity.field <- function(object, field) {
  switch(
    field,
    X = object$predictors,
    y = object$response,
    U = object$latent,
    d = object$intrinsic.dim,
    p = object$ambient.dim,
    sigma = switch(
      object$response.spec$family,
      gaussian = object$response.spec$parameters$parameters$sd,
      clustered.gaussian =
        object$response.spec$parameters$parameters$residual.sd,
      NA_real_),
    object[[field]])
}

test_that("G1--G7 registry recipes reproduce frozen scientific fields", {
  fixture.path <- .legacy.parity.fixture.path()
  skip_if(!nzchar(fixture.path), "legacy parity fixture is unavailable")
  fixtures <- readRDS(fixture.path)
  expect_identical(
    attr(fixtures, "source.revision"),
    "5ecf721eae2765d8cde01d7fb82b17ff8bde8599")

  for (case.id in names(fixtures)) {
    fixture <- fixtures[[case.id]]
    actual <- .materialize.legacy.parity.case(fixture)
    expected <- fixture$scientific
    for (field in names(expected)) {
      observed <- .canonical.legacy.parity.field(actual, field)
      if (is.numeric(expected[[field]])) {
        expect_equal(observed, expected[[field]], tolerance = 1e-14,
                     info = paste(case.id, field))
      } else {
        expect_identical(observed, expected[[field]],
                         info = paste(case.id, field))
      }
    }
  }
})
