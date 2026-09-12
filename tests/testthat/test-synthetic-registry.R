.read.synthetic.registry.table <- function(file) {
  path <- geosmooth:::.synthetic.registry.asset(file)
  expect_true(nzchar(path) && file.exists(path), info = file)
  utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE)
}

test_that("normalized component foreign keys resolve exactly once", {
  recipes <- .read.synthetic.registry.table("recipes.csv")
  table.by.field <- c(
    geometry.id = "geometries.csv",
    sampling.id = "samplings.csv",
    truth.id = "truths.csv",
    response.id = "responses.csv")
  kind.by.field <- c(
    geometry.id = "geometry", sampling.id = "sampling",
    truth.id = "truth", response.id = "response")

  expect_setequal(recipes$recipe.id, synthetic.registry.ids())
  expect_identical(anyDuplicated(recipes$recipe.id), 0L)
  expect_true(all(recipes$status == "active"))
  expect_true(all(recipes$version == 1L))
  expect_true(all(nzchar(recipes$compatibility.recipe)))
  expect_true(all(nzchar(recipes$metadata)))

  all.component.ids <- character()
  for (field in names(table.by.field)) {
    components <- .read.synthetic.registry.table(table.by.field[[field]])
    expect_identical(anyDuplicated(components$component.id), 0L, info = field)
    expect_true(all(components$kind == kind.by.field[[field]]), info = field)
    expect_true(all(components$version >= 1L), info = field)
    expect_true(all(nzchar(components$parameters)), info = field)
    for (id in recipes[[field]]) {
      expect_identical(sum(components$component.id == id), 1L, info = id)
      component <- geosmooth:::.synthetic.registry.component(
        id, kind.by.field[[field]])
      row <- components[components$component.id == id, , drop = FALSE]
      expect_identical(
        geosmooth:::.synthetic.sha256(component),
        row$component.sha256,
        info = id)
      expect_identical(
        geosmooth:::.synthetic.registry.json.text(component$parameters),
        row$parameters,
        info = id)
    }
    all.component.ids <- c(all.component.ids, components$component.id)
  }
  expect_identical(anyDuplicated(all.component.ids), 0L)
})

test_that("every registry recipe is reconstructed from normalized rows", {
  recipes <- .read.synthetic.registry.table("recipes.csv")
  for (i in seq_len(nrow(recipes))) {
    spec <- synthetic.registry.spec(recipes$recipe.id[i])
    expect_identical(
      spec$specification.sha256, recipes$specification.sha256[i],
      info = recipes$recipe.id[i])
    expect_identical(
      geosmooth:::.synthetic.registry.resolve(recipes$recipe.id[i]),
      spec,
      info = recipes$recipe.id[i])
  }
})

test_that("instance, checksum, and fixture foreign keys are complete", {
  recipes <- .read.synthetic.registry.table("recipes.csv")
  instances <- .read.synthetic.registry.table("instances.csv")
  checksums <- .read.synthetic.registry.table("checksums.csv")
  environments <- .read.synthetic.registry.table(
    "environment_fingerprints.csv")
  schema <- .read.synthetic.registry.table(
    "environment_fingerprint_schema.csv")
  required.scope <- c(
    "r.version", "platform", "architecture", "compiler",
    "os.name", "os.release", "os.version", "os.machine", "endianness",
    "rng.policy", "rng.kind", "normal.kind", "sample.kind", "rng.version",
    "blas.path", "blas.sha256", "lapack.path", "lapack.sha256",
    "lapack.version", "math.runtime", "geosmooth.version",
    "registry.schema.version", "evaluator.registry.version",
    "dependency.versions")

  expect_identical(anyDuplicated(instances$instance.id), 0L)
  expect_identical(anyDuplicated(instances$checksum.id), 0L)
  expect_identical(anyDuplicated(checksums$checksum.id), 0L)
  expect_identical(
    anyDuplicated(environments$environment.fingerprint.id), 0L)
  expect_true(all(instances$recipe.id %in% recipes$recipe.id))
  expect_true(all(instances$status == "active"))
  expect_true(all(instances$version == 1L))
  expect_identical(schema$field, required.scope)
  expect_true(all(schema$required))
  expect_true(all(required.scope %in% names(environments)))
  expect_false(any(required.scope %in% names(checksums)))
  expect_true(all(
    checksums$environment.fingerprint.id %in%
      environments$environment.fingerprint.id))
  expect_true(all(checksums$payload.contract == "synthetic-content-v1"))
  expect_true(all(checksums$serialization == "R-xdr-v3"))

  for (i in seq_len(nrow(instances))) {
    instance <- instances[i, , drop = FALSE]
    checksum <- checksums[
      checksums$checksum.id == instance$checksum.id, , drop = FALSE]
    expect_identical(nrow(checksum), 1L, info = instance$instance.id)
    expect_identical(
      checksum$instance.id, instance$instance.id,
      info = instance$instance.id)
    expect_identical(
      checksum$content.sha256, instance$content.sha256,
      info = instance$instance.id)
    expect_identical(
      sum(
        environments$environment.fingerprint.id ==
          checksum$environment.fingerprint.id),
      1L,
      info = instance$instance.id)
    fixture.path <- geosmooth:::.synthetic.registry.fixture.path(
      instance$fixture.path)
    expect_true(
      nzchar(fixture.path) && file.exists(fixture.path),
      info = instance$instance.id)
    fixture <- readRDS(fixture.path)
    expect_identical(
      validate.synthetic.dataset(fixture), fixture,
      info = instance$instance.id)
    object <- materialize.synthetic.instance(instance$instance.id)
    expect_identical(object$dataset.id, instance$instance.id)
    expect_identical(
      attr(object, "frozen.instance.id", exact = TRUE),
      instance$instance.id)
    if (geosmooth:::.synthetic.checksum.scope.matches(checksum)) {
      expect_identical(
        synthetic.dataset.checksum(object),
        instance$content.sha256,
        info = instance$instance.id)
      expect_identical(
        attr(object, "verification.scope", exact = TRUE),
        "exact-environment")
    } else {
      expect_true(
        compare.synthetic.dataset(object, fixture)$equal,
        info = instance$instance.id)
      expect_identical(
        attr(object, "verification.scope", exact = TRUE),
        "scientific-parity")
    }
  }

  for (field in required.scope) {
    value <- as.character(environments[[field]])
    expect_true(all(nzchar(value)), info = field)
    expect_false(any(grepl(
      "unknown|unavailable", value, ignore.case = TRUE)),
      info = field)
  }
  for (i in seq_len(nrow(environments))) {
    row <- environments[i, , drop = FALSE]
    fingerprint <-
      geosmooth:::.synthetic.environment.fingerprint.from.row(row)
    expect_identical(
      row$environment.fingerprint.id,
      geosmooth:::.synthetic.environment.fingerprint.id(fingerprint))
    expect_identical(
      row$fingerprint.sha256,
      geosmooth:::.synthetic.sha256(fingerprint))
  }
})

test_that("registry JSON is typed and matrix values are row-major", {
  value <- list(
    count = 2L,
    coefficients = c(b0 = 1, b1 = 2),
    frame = matrix(as.double(1:6), nrow = 2L, byrow = TRUE),
    missing = NA_integer_,
    absent = NULL)
  json <- geosmooth:::.synthetic.registry.json.text(value)
  expect_match(
    json,
    '"nrow":2,"ncol":3,"values":\\[\\[1,2,3\\],\\[4,5,6\\]\\]',
    perl = TRUE)
  expect_identical(
    geosmooth:::.synthetic.registry.json.value(json), value)
})

test_that("exact and scientific fallback verification are deterministic", {
  instances <- .read.synthetic.registry.table("instances.csv")
  checksums <- .read.synthetic.registry.table("checksums.csv")
  instance <- instances[1L, , drop = FALSE]
  checksum <- checksums[
    checksums$checksum.id == instance$checksum.id, , drop = FALSE]
  fixture.path <- geosmooth:::.synthetic.registry.fixture.path(
    instance$fixture.path)
  fixture <- readRDS(fixture.path)
  object <- materialize.synthetic.instance(instance$instance.id)
  current <- geosmooth:::.synthetic.environment.fingerprint(
    fixture$rng.policy)

  exact <- geosmooth:::.synthetic.verify.frozen.instance(
    object, fixture, instance, checksum, current.fingerprint = current)
  expected.scope <- if (geosmooth:::.synthetic.checksum.scope.matches(
      checksum, current = current)) {
    "exact-environment"
  } else {
    "scientific-parity"
  }
  expect_identical(
    attr(exact, "verification.scope", exact = TRUE),
    expected.scope)

  mismatched <- current
  mismatched$blas.sha256 <- paste0("mismatch-", mismatched$blas.sha256)
  expect_false(geosmooth:::.synthetic.checksum.scope.matches(
    checksum, current = mismatched))
  fallback <- geosmooth:::.synthetic.verify.frozen.instance(
    object, fixture, instance, checksum, current.fingerprint = mismatched)
  expect_identical(
    attr(fallback, "verification.scope", exact = TRUE),
    "scientific-parity")
})

test_that("ordinary API behavior has no process-wide verification bypass", {
  source.files <- Sys.glob(file.path(
    test_path("..", "..", "R"), "*.R"))
  source <- paste(
    unlist(lapply(source.files, readLines, warn = FALSE)),
    collapse = "\n")
  expect_false(grepl("geosmooth.registry.verify", source, fixed = TRUE))

  old <- options(geosmooth.registry.verify = FALSE)
  on.exit(options(old), add = TRUE)
  spec <- synthetic.registry.spec("G1")
  expect_identical(
    spec,
    geosmooth:::.synthetic.registry.resolve("G1"))
})

test_that("specification hashes do not depend on collation", {
  original <- Sys.getlocale("LC_COLLATE")
  on.exit(try(Sys.setlocale("LC_COLLATE", original), silent = TRUE), add = TRUE)
  expected <- synthetic.registry.spec("G2")$specification.sha256
  changed <- try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE)
  skip_if(inherits(changed, "try-error"), "C collation is unavailable")
  expect_identical(
    synthetic.registry.spec("G2")$specification.sha256, expected)
})

test_that("legacy parameter overrides remain content-bound but non-frozen", {
  registered <- synthetic.registry.spec("G3a")
  altered <- synthetic.registry.spec("G3a", list(R = 2))
  expect_false(identical(
    registered$specification.sha256, altered$specification.sha256))
  object <- materialize.synthetic(
    altered, n = 20, seed = 1, rng.policy = "legacy")
  expect_match(
    object$dataset.id, altered$specification.sha256, fixed = TRUE)
})
