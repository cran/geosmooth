.package.source.files <- function() {
  c(
    list.files("R", pattern = "\\.[Rr]$", full.names = TRUE),
    list.files("vignettes", pattern = "\\.[Rr](md)?$", full.names = TRUE))
}

test_that("installed synthetic code has no cross-repository source calls", {
  files <- .package.source.files()
  text <- unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl(
    "source\\s*\\(\\s*['\"](?:/Users/|\\.\\./)",
    text, perl = TRUE)))
})

test_that("legacy family disposition ledger is normalized and explicit", {
  path <- geosmooth:::.synthetic.registry.asset(
    "legacy_family_disposition.csv")
  expect_true(nzchar(path) && file.exists(path))
  ledger <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c(
    "family.id", "legacy.source", "legacy.symbols", "evidence", "decision",
    "canonical.owner", "implementation.status",
    "parity.or.retirement.record", "next.gate")
  expect_setequal(names(ledger), required)
  expect_identical(anyDuplicated(ledger$family.id), 0L)
  expect_true(all(nzchar(ledger$decision)))
  expect_true(all(nzchar(ledger$implementation.status)))
  expect_true(all(nzchar(ledger$next.gate)))
  expect_true(all(
    ledger$implementation.status %in% c(
      "implemented", "implemented-core", "implemented-noise-free",
      "implemented-g5-subset", "partial-private", "approved-pending",
      "pending-decision", "out-of-scope")))
})

test_that("retired compatibility APIs stay out of the package namespace", {
  retired <- c(
    "dgp.g1", "dgp.g2", "dgp.g3a", "dgp.g3b", "dgp.g3c", "dgp.g3d",
    "dgp.g4", "dgp.g5", "dgp.g6", "dgp.g7", "dgp.materialize",
    "dgp.content.sha256", "density.dependency.precheck",
    "fit.density.empirical", "fit.density.graph.random.walk",
    "graph.low.pass.filter", "print.dgp_dataset")
  namespace <- asNamespace("geosmooth")
  expect_false(any(retired %in% getNamespaceExports("geosmooth")))
  expect_false(any(vapply(
    retired, exists, logical(1), envir = namespace, inherits = FALSE)))

  expect_true(exists(
    ".fit.density.empirical", envir = namespace, inherits = FALSE))
  expect_true(exists(
    ".fit.density.graph.random.walk", envir = namespace, inherits = FALSE))
})

test_that("canonical APIs do not expose deprecated arguments or fields", {
  expect_false("sync.lambda.grid" %in% names(formals(fit.slpl.tf)))

  object <- materialize.synthetic(
    synthetic.registry.spec("G1"),
    n = 10L, seed = 1L, rng.policy = "named.stream.v1")
  aliases <- c("U", "Z", "X", "y", "d", "p", "sigma", "gtag", "params")
  expect_false(any(aliases %in% names(object)))
  expect_false(inherits(object, "dgp_dataset"))
})
