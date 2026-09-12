.ssrhe.parity.fixture.path <- function() {
  dev <- tryCatch(
    testthat::test_path("..", "fixtures", "ssrhe_legacy_parity_v1.rds"),
    error = function(e) "")
  if (nzchar(dev) && file.exists(dev)) return(dev)
  ""
}

.expect.ssrhe.numeric.parity <- function(actual, expected, info) {
  expect_identical(dim(actual), dim(expected), info = info)
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)
  difference <- abs(actual - expected)
  bound <- 1e-14 + 1e-12 * abs(expected)
  expect_true(all(difference <= bound), info = info)
}

.ssrhe.canonical.field <- function(object, field) {
  switch(
    field,
    U = object$latent,
    X = object$predictors,
    y = object$response,
    object[[field]])
}

test_that("the SSRHE registry contains every S/V combination and surface lane", {
  ids <- synthetic.registry.ids()
  one.d <- ids[grepl("^S[0-9]{2}\\.V[1-3]$", ids)]
  flat <- ids[grepl("^ssrhe\\.flat\\.", ids)]
  quad <- ids[grepl("^ssrhe\\.quadform\\.", ids)]
  expect_length(one.d, 48L)
  expect_length(flat, 3L)
  expect_length(quad, 45L)
  expect_setequal(
    one.d,
    as.vector(outer(
      sprintf("S%02d", 1:16), paste0("V", 1:3), paste, sep = ".")))
})

test_that("legacy SSRHE seeds reproduce the documented formulas", {
  expect_identical(synthetic.registry.seed("S16.V3", replicate = 7),
                   289308L)
  expect_identical(
    synthetic.registry.seed("ssrhe.flat.d4.v1", n = 31, replicate = 2),
    287312L)
  expect_identical(
    synthetic.registry.seed(
      "ssrhe.quadform.d3.k1.c125.v1", n = 31, replicate = 2),
    294537L)
  expect_error(
    synthetic.registry.seed("ssrhe.flat.d2.v1"),
    "n must")
})

test_that("all 48 one-dimensional recipes preserve the frozen legacy grid", {
  fixture.path <- .ssrhe.parity.fixture.path()
  skip_if(!nzchar(fixture.path), "SSRHE parity fixture is unavailable")
  fixture <- readRDS(fixture.path)
  expect_identical(attr(fixture, "contract"), "ssrhe-legacy-parity-v1")
  expect_length(fixture$one.d, 720L)
  for (case.id in names(fixture$one.d)) {
    case <- fixture$one.d[[case.id]]
    actual <- materialize.synthetic(
      synthetic.registry.spec(case$recipe.id),
      n = case$n, seed = case$seed, rng.policy = "legacy")
    for (field in names(case$scientific)) {
      .expect.ssrhe.numeric.parity(
        .ssrhe.canonical.field(actual, field), case$scientific[[field]],
        paste(case.id, field))
    }
  }
})

test_that("flat and quadform recipes preserve dimensions 2--4", {
  fixture.path <- .ssrhe.parity.fixture.path()
  skip_if(!nzchar(fixture.path), "SSRHE parity fixture is unavailable")
  fixture <- readRDS(fixture.path)
  expect_length(fixture$surfaces, 126L)
  for (case.id in names(fixture$surfaces)) {
    case <- fixture$surfaces[[case.id]]
    actual <- materialize.synthetic(
      synthetic.registry.spec(case$recipe.id),
      n = case$n, seed = case$seed, rng.policy = "legacy")
    for (field in names(case$scientific)) {
      .expect.ssrhe.numeric.parity(
        .ssrhe.canonical.field(actual, field), case$scientific[[field]],
        paste(case.id, field))
    }
  }
})

test_that("one-dimensional variants retain their named mechanisms", {
  v1 <- synthetic.registry.spec("S01.V1")
  v2 <- synthetic.registry.spec("S01.V2")
  v3 <- synthetic.registry.spec("S01.V3")
  expect_identical(v1$sampling.spec$family, "uniform.interval")
  expect_identical(v2$sampling.spec$family, "truncated.normal")
  expect_identical(v3$sampling.spec$family, "gapped.uniform")
  expect_identical(v3$response.spec$family, "laplace.outlier")
  expect_identical(
    v3$response.spec$parameters$parameters$minimum.outliers, 1L)

  for (n in c(1L, 2L, 10L, 11L)) {
    object <- materialize.synthetic(
      v3, n = n,
      seed = synthetic.registry.seed("S01.V3", replicate = 1),
      rng.policy = "legacy")
    expect_identical(nrow(object$predictors), n)
    expect_true(all(diff(object$predictors[, 1]) >= 0))
  }
})
