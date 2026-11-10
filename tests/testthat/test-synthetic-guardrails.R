.package.source.files <- function() {
  root <- geosmooth.test.source.root()
  if (is.null(root)) return(character())
  c(
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$",
               full.names = TRUE),
    list.files(file.path(root, "vignettes"), pattern = "\\.[Rr](md)?$",
               full.names = TRUE))
}

test_that("installed synthetic code has no cross-repository source calls", {
  files <- .package.source.files()
  skip_if(!length(files),
          "package source tree is unavailable in installed-package tests")
  expect_gt(length(files), 0L)
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

test_that("geometry re-exports and separate quadratic selection entry points are retired", {
  migrated <- c(
    "edge.lengths.synthetic.geometry", "embed.synthetic.geometry",
    "quadform.gradient", "quadform.metric", "synthetic.circle", "synthetic.helix",
    "synthetic.point.line.junction", "synthetic.quadform",
    "synthetic.sampling.clustered", "synthetic.sampling.dirichlet.zeros",
    "synthetic.sampling.gapped.uniform", "synthetic.sampling.grid.interval",
    "synthetic.sampling.truncated.normal", "synthetic.sampling.uniform.box",
    "synthetic.sampling.uniform.disk", "synthetic.sampling.uniform.interval",
    "synthetic.sampling.uniform.rectangle", "synthetic.simplex",
    "synthetic.sphere.cap", "synthetic.stratified", "synthetic.stratum.point",
    "synthetic.stratum.rectangle", "synthetic.stratum.segment",
    "synthetic.torus.patch", "synthetic.trefoil"
  )
  retired.fitters <- c("fit.ssrhe.hessian.regression.cv", "fit.ssrhe.hessian.regression.gcv")
  exports <- getNamespaceExports("geosmooth")
  expect_length(migrated, 25L)
  expect_false(any(c(migrated, retired.fitters) %in% exports))
  expect_true(all(migrated %in% getNamespaceExports("dgraphs")))
  expect_true("synthetic.sampling.stratified" %in% exports)
  expect_true("fit.ssrhe.hessian.regression" %in% exports)
  expect_false(any(vapply(retired.fitters, exists, logical(1),
                          envir = asNamespace("geosmooth"), inherits = FALSE)))
  # Exercise the user-visible dependency spelling, without attaching dgraphs.
  spec <- synthetic.spec(
    dgraphs::synthetic.circle(),
    dgraphs::synthetic.sampling.uniform.interval(0, 2 * pi),
    synthetic.truth.polynomial(c(b0 = 0)), synthetic.response.gaussian(sd = 0.1)
  )
  data <- materialize.synthetic(spec, n = 12L, seed = 31L)
  expect_s3_class(data, "synthetic_dataset")
  expect_silent(validate.synthetic.dataset(data))
})
