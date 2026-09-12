# Canonical materialized synthetic dataset contract.

.new.synthetic.dataset <- function(
    spec, n, seed, rng.policy, rng.streams, sample.out, frame.matrix,
    truth.out, response.out, parameters = list(), provenance = list()) {
  geometry <- spec$geometry.spec
  gp <- geometry$parameters
  if (is.null(sample.out$predictors)) {
    X <- .embed.synthetic.geometry(
      geometry, sample.out$latent, frame.matrix = frame.matrix)
  } else {
    X <- as.matrix(sample.out$predictors)
  }
  storage.mode(X) <- "double"
  latent <- sample.out$latent
  if (!is.null(latent)) {
    latent <- as.matrix(latent)
    storage.mode(latent) <- "double"
  }
  region <- sample.out$region
  if (!is.null(region)) region <- as.character(region)
  family <- geometry$family
  if (family == "g4.segment.rectangle") {
    intrinsic.dim <- NA_integer_
    intrinsic.dim.by.region <- gp$intrinsic.dim.by.region
    codimension <- NA_integer_
    codimension.by.region <- gp$codimension.by.region
    declared.regions <- gp$declared.regions
  } else if (family == "simplex") {
    intrinsic.dim <- NA_integer_
    intrinsic.dim.by.region <- c(
      interior = gp$parts - 1L,
      zero = gp$parts -
        length(spec$sampling.spec$parameters$zero.parts) - 1L)
    codimension <- NA_integer_
    codimension.by.region <- gp$parts - intrinsic.dim.by.region
    declared.regions <- c("interior", "zero")
  } else {
    intrinsic.dim <- gp$intrinsic.dim
    intrinsic.dim.by.region <- NULL
    codimension <- ncol(X) - intrinsic.dim
    codimension.by.region <- NULL
    declared.regions <- if (is.null(region)) NULL else unique(region)
  }
  observed.regions <- if (is.null(region)) NULL else {
    if (is.null(declared.regions)) unique(region)
    else declared.regions[declared.regions %in% unique(region)]
  }
  label <- .synthetic.dataset.label(spec$recipe.id)
  dataset.id <- paste0(
    label, "__spec", spec$specification.sha256,
    "__n", n, "__seed", seed, "__rng", rng.policy)
  truth.spec <- spec$truth.spec
  response.spec <- spec$response.spec
  truth.scope <- truth.spec$parameters$scope
  truth.estimand <- truth.spec$parameters$estimand
  object <- list(
    dataset.id = dataset.id,
    specification.sha256 = spec$specification.sha256,
    recipe.id = spec$recipe.id,
    registry.tag = spec$registry.tag,
    n = as.integer(n),
    intrinsic.dim = as.integer(intrinsic.dim),
    intrinsic.dim.by.region = intrinsic.dim.by.region,
    ambient.dim = as.integer(ncol(X)),
    codimension = as.integer(codimension),
    codimension.by.region = codimension.by.region,
    latent = latent,
    latent.mask = sample.out$latent.mask,
    predictors = X,
    truth = as.numeric(truth.out$value),
    truth.scope = truth.scope,
    truth.estimand = truth.estimand,
    response = as.numeric(response.out$value),
    region = region,
    declared.regions = declared.regions,
    observed.regions = observed.regions,
    geometry.spec = spec$geometry.spec,
    sampling.spec = spec$sampling.spec,
    truth.spec = truth.spec,
    response.spec = response.spec,
    seed = as.integer(seed),
    rng.policy = rng.policy,
    rng.streams = rng.streams,
    parameters = utils::modifyList(
      list(
        frame.matrix = frame.matrix,
        sampling = sample.out$parameters,
        truth = truth.out$parameters,
        response = response.out$parameters
      ),
      parameters),
    provenance = provenance
  )
  class(object) <- c("synthetic_dataset", "list")
  attr(object, "compatibility") <- spec$compatibility
  attr(object, "metadata") <- spec$metadata
  validate.synthetic.dataset(object)
  object
}

.synthetic.dataset.label <- function(recipe.id) {
  if (is.null(recipe.id)) "adhoc" else
    gsub("[^A-Za-z0-9_.-]", "-", recipe.id)
}

.synthetic.expected.dataset.id <- function(x) {
  paste0(
    .synthetic.dataset.label(x$recipe.id),
    "__spec", x$specification.sha256,
    "__n", x$n,
    "__seed", x$seed,
    "__rng", x$rng.policy)
}

.validate.synthetic.rng.metadata <- function(x) {
  random.frame <- identical(
    x$geometry.spec$parameters$frame, "random.orthonormal")
  response.random <- .synthetic.response.uses.rng(x$response.spec)
  compatibility <- attr(x, "compatibility", exact = TRUE)
  frame.seed <- compatibility$frame.seed %||% NULL
  seed.policy <- compatibility$seed.policy %||% "named.stream.v1"
  expected.seeds <- .validate.synthetic.seed.plan(
    x$seed,
    frame.seed = frame.seed,
    need.frame = x$rng.policy == "legacy" && random.frame,
    need.response = x$rng.policy == "legacy" && response.random &&
      !identical(seed.policy, "ssrhe.1d.v1"))
  if (!identical(x$parameters$resolved.seeds, expected.seeds)) {
    stop("parameters$resolved.seeds is inconsistent with the RNG policy.",
         call. = FALSE)
  }
  if (x$rng.policy == "named.stream.v1") {
    expected.names <- c(
      "sampling", "geometry.frame", "truth",
      "response", "response.auxiliary")
    valid.state <- function(state) {
      is.integer(state) && length(state) == 7L && !anyNA(state)
    }
    expected.streams <- .with.synthetic.rng.preserved(
      .synthetic.named.streams(x$seed))
    if (!is.list(x$rng.streams) ||
        !identical(names(x$rng.streams), expected.names) ||
        !all(vapply(x$rng.streams, valid.state, logical(1))) ||
        !identical(x$rng.streams, expected.streams)) {
      stop("rng.streams violates the named.stream.v1 state contract.",
           call. = FALSE)
    }
  } else {
    expected.names <- c("sampling", "geometry.frame", "response")
    scalar.seed.or.null <- function(value) {
      is.null(value) ||
        (is.integer(value) && length(value) == 1L && !is.na(value) &&
           value >= 0L)
    }
    if (!is.list(x$rng.streams) ||
        !identical(names(x$rng.streams), expected.names) ||
        !identical(x$rng.streams$sampling, x$seed) ||
        !all(vapply(x$rng.streams, scalar.seed.or.null, logical(1))) ||
        !identical(
          x$rng.streams$geometry.frame,
          expected.seeds$frame.seed %||% NULL) ||
        !identical(
          x$rng.streams$response,
          expected.seeds$response.seed %||% NULL)) {
      stop("rng.streams violates the legacy state contract.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.validate.synthetic.realized.support <- function(x) {
  geometry <- x$geometry.spec
  sampling <- x$sampling.spec
  gp <- geometry$parameters
  sp <- sampling$parameters
  family <- sampling$family
  .validate.synthetic.support(geometry, sampling)
  if (geometry$family == "g4.segment.rectangle") {
    A <- x$region == "A"
    B <- x$region == "B"
    if (!all(A | B) || !any(A) || !any(B) ||
        any(x$latent[A, 1L] < -1 | x$latent[A, 1L] > 0) ||
        any(!is.na(x$latent[A, 2L])) ||
        any(abs(x$predictors[A, 1L] - x$latent[A, 1L]) > 1e-12) ||
        any(x$latent[B, 1L] < 0 | x$latent[B, 1L] > 1) ||
        any(x$latent[B, 2L] < -0.5 | x$latent[B, 2L] > 0.5) ||
        any(abs(x$predictors[B, 1L] - x$latent[B, 1L]) > 1e-12) ||
        any(abs(x$predictors[B, 2L] - x$latent[B, 2L]) > 1e-12) ||
        any(abs(x$predictors[B, 3L]) > 1e-12)) {
      stop("Realized G4 rows violate their stratum geometry.", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  latent.source <- if (geometry$family == "simplex" &&
                         is.null(x$latent)) x$predictors else x$latent
  latent <- .validate.synthetic.latent(geometry, latent.source)
  if (family == "uniform.box") {
    lower <- rep(sp$lower, length.out = ncol(latent))
    upper <- rep(sp$upper, length.out = ncol(latent))
    valid <- rowSums(sweep(latent, 2L, lower, "<")) == 0L &
      rowSums(sweep(latent, 2L, upper, ">")) == 0L
  } else if (family == "uniform.disk") {
    valid <- rowSums(latent^2) <= sp$radius^2 * (1 + 1e-12)
  } else if (family %in% c(
      "uniform.interval", "grid.interval", "truncated.normal")) {
    valid <- latent[, 1L] >= sp$lower & latent[, 1L] <= sp$upper
  } else if (family == "uniform.rectangle") {
    valid <- rowSums(sweep(latent, 2L, sp$lower, "<")) == 0L &
      rowSums(sweep(latent, 2L, sp$upper, ">")) == 0L
  } else if (family == "gapped.uniform") {
    valid <- vapply(latent[, 1L], function(value) {
      any(value >= sp$intervals[, 1L] & value <= sp$intervals[, 2L])
    }, logical(1))
  } else {
    valid <- rep(TRUE, nrow(latent))
  }
  if (!all(valid)) {
    stop("Realized latent coordinates violate the sampling support.",
         call. = FALSE)
  }
  if (family == "dirichlet.zeros") {
    zero <- x$region == "zero"
    if (any(zero) &&
        any(abs(latent[zero, sp$zero.parts, drop = FALSE]) > 1e-12)) {
      stop("Structural-zero rows violate zero.parts.", call. = FALSE)
    }
  }
  frame <- x$parameters$frame.matrix
  if (is.null(gp$frame)) {
    if (!is.null(frame)) {
      stop("A frameless geometry cannot store frame.matrix.", call. = FALSE)
    }
    reconstructed <- .embed.synthetic.geometry(geometry, latent)
  } else {
    k <- gp$canonical.dim
    if (!is.matrix(frame) ||
        !identical(dim(frame), c(gp$ambient.dim, k)) ||
        any(!is.finite(frame)) ||
        !isTRUE(all.equal(
          crossprod(frame), diag(k), tolerance = 1e-10))) {
      stop("Stored frame.matrix has invalid dimensions or is not orthonormal.",
           call. = FALSE)
    }
    if (gp$frame == "canonical" &&
        !isTRUE(all.equal(frame, diag(k), tolerance = 0))) {
      stop("Stored canonical frame.matrix is inconsistent.", call. = FALSE)
    }
    if (gp$frame == "supplied" &&
        !isTRUE(all.equal(frame, gp$frame.matrix, tolerance = 0))) {
      stop("Stored supplied frame.matrix is inconsistent.", call. = FALSE)
    }
    reconstructed <- .embed.synthetic.geometry(
      geometry, latent, frame.matrix = frame)
  }
  if (!isTRUE(all.equal(
    unname(reconstructed), unname(x$predictors), tolerance = 1e-10))) {
    stop("Stored frame and latent coordinates do not reconstruct predictors.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a canonical synthetic dataset
#' @param x A `synthetic_dataset`.
#' @return `x`, invisibly.
#' @examples
#' x <- materialize.synthetic(synthetic.registry.spec("G1"), n = 20L, seed = 1L)
#' validate.synthetic.dataset(x)
#' @export
validate.synthetic.dataset <- function(x) {
  required <- c(
    "dataset.id", "specification.sha256", "recipe.id", "registry.tag", "n",
    "intrinsic.dim", "intrinsic.dim.by.region", "ambient.dim",
    "codimension", "codimension.by.region", "latent", "latent.mask",
    "predictors", "truth", "truth.scope", "truth.estimand", "response",
    "region", "declared.regions", "observed.regions", "geometry.spec",
    "sampling.spec", "truth.spec", "response.spec", "seed", "rng.policy",
    "rng.streams", "parameters", "provenance")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Missing synthetic dataset fields: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!inherits(x, "synthetic_dataset")) {
    stop("x must inherit from synthetic_dataset.", call. = FALSE)
  }
  if (!is.integer(x$n) || length(x$n) != 1L || is.na(x$n) || x$n < 1L ||
      !is.integer(x$ambient.dim) || length(x$ambient.dim) != 1L ||
      is.na(x$ambient.dim) || x$ambient.dim < 1L ||
      !is.integer(x$seed) || length(x$seed) != 1L ||
      is.na(x$seed) || x$seed < 0L) {
    stop("n, ambient.dim, and seed must be valid integer scalars.",
         call. = FALSE)
  }
  if (!is.matrix(x$predictors) || !identical(nrow(x$predictors), x$n) ||
      ncol(x$predictors) != x$ambient.dim ||
      any(!is.finite(x$predictors))) {
    stop("predictors violate the n by ambient.dim finite-matrix contract.",
         call. = FALSE)
  }
  if (length(x$truth) != x$n || length(x$response) != x$n ||
      any(!is.finite(x$truth)) || any(!is.finite(x$response))) {
    stop("truth and response must be finite length-n vectors.",
         call. = FALSE)
  }
  if (!x$truth.scope %in% c("population", "finite.design")) {
    stop("Unknown truth.scope.", call. = FALSE)
  }
  if (!is.character(x$rng.policy) || length(x$rng.policy) != 1L ||
      !x$rng.policy %in% c("named.stream.v1", "legacy")) {
    stop("rng.policy is missing or unknown.", call. = FALSE)
  }
  component.fields <- c(
    geometry.spec = "geometry", sampling.spec = "sampling",
    truth.spec = "truth", response.spec = "response")
  for (field in names(component.fields)) {
    component <- x[[field]]
    if (!.is.synthetic.component(component, component.fields[[field]]) ||
        !is.integer(component$version) || length(component$version) != 1L ||
        is.na(component$version) || component$version < 1L) {
      stop(field, " is not a valid versioned synthetic component.",
           call. = FALSE)
    }
  }
  expected.hash <- .synthetic.sha256(list(
    geometry.spec = x$geometry.spec,
    sampling.spec = x$sampling.spec,
    truth.spec = x$truth.spec,
    response.spec = x$response.spec,
    compatibility = attr(x, "compatibility", exact = TRUE),
    metadata = attr(x, "metadata", exact = TRUE)))
  if (!identical(expected.hash, x$specification.sha256) ||
      !grepl("^[0-9a-f]{64}$", x$specification.sha256)) {
    stop("dataset.id and specification.sha256 are inconsistent.",
         call. = FALSE)
  }
  frozen.id <- attr(x, "frozen.instance.id", exact = TRUE)
  if (is.null(frozen.id)) {
    if (!identical(x$dataset.id, .synthetic.expected.dataset.id(x))) {
      stop("Ordinary dataset.id does not exactly encode its canonical fields.",
           call. = FALSE)
    }
  } else {
    if (!is.character(frozen.id) || length(frozen.id) != 1L ||
        !identical(x$dataset.id, frozen.id)) {
      stop("Frozen instance identity metadata is inconsistent.",
           call. = FALSE)
    }
    row <- .synthetic.registry.row("instances.csv", "instance.id", frozen.id)
    spec <- synthetic.registry.spec(row$recipe.id)
    if (!identical(x$recipe.id, spec$recipe.id) ||
        !identical(x$specification.sha256, spec$specification.sha256) ||
        !identical(x$n, as.integer(row$n)) ||
        !identical(x$seed, as.integer(row$seed)) ||
        !identical(x$rng.policy, row$rng.policy)) {
      stop("Frozen dataset does not match its complete instance row.",
           call. = FALSE)
    }
  }
  geometry <- x$geometry.spec
  gp <- geometry$parameters
  if (geometry$family == "g4.segment.rectangle") {
    expected.intrinsic <- NA_integer_
    expected.intrinsic.by.region <- gp$intrinsic.dim.by.region
    expected.codimension <- NA_integer_
    expected.codimension.by.region <- gp$codimension.by.region
    expected.declared <- gp$declared.regions
  } else if (geometry$family == "simplex") {
    expected.intrinsic <- NA_integer_
    expected.intrinsic.by.region <- c(
      interior = gp$parts - 1L,
      zero = gp$parts -
        length(x$sampling.spec$parameters$zero.parts) - 1L)
    expected.codimension <- NA_integer_
    expected.codimension.by.region <-
      gp$parts - expected.intrinsic.by.region
    expected.declared <- c("interior", "zero")
  } else {
    expected.intrinsic <- as.integer(gp$intrinsic.dim)
    expected.intrinsic.by.region <- NULL
    expected.codimension <- as.integer(x$ambient.dim - expected.intrinsic)
    expected.codimension.by.region <- NULL
    expected.declared <- if (is.null(x$region)) NULL else unique(x$region)
  }
  if (!identical(x$ambient.dim, as.integer(gp$ambient.dim)) ||
      !identical(x$intrinsic.dim, expected.intrinsic) ||
      !identical(x$intrinsic.dim.by.region, expected.intrinsic.by.region) ||
      !identical(x$codimension, expected.codimension) ||
      !identical(x$codimension.by.region, expected.codimension.by.region) ||
      !identical(x$declared.regions, expected.declared)) {
    stop("Dataset dimension and region metadata are inconsistent.",
         call. = FALSE)
  }
  if (!identical(x$truth.scope, x$truth.spec$parameters$scope) ||
      !identical(x$truth.estimand, x$truth.spec$parameters$estimand)) {
    stop("Truth scope or estimand is inconsistent with truth.spec.",
         call. = FALSE)
  }
  if (!is.null(x$region)) {
    if (length(x$region) != x$n ||
        any(!x$region %in% x$declared.regions) ||
        !identical(x$observed.regions,
                   x$declared.regions[x$declared.regions %in%
                                        unique(x$region)])) {
      stop("Region declarations and observations are inconsistent.",
           call. = FALSE)
    }
  } else if (!is.null(x$declared.regions) ||
             !is.null(x$observed.regions)) {
    stop("Unstratified datasets cannot declare or observe regions.",
         call. = FALSE)
  }
  if (!is.null(x$latent)) {
    if (!is.matrix(x$latent) || nrow(x$latent) != x$n) {
      stop("latent must be NULL or an n-row matrix.", call. = FALSE)
    }
    if (is.null(x$latent.mask)) {
      if (any(!is.finite(x$latent))) {
        stop("Unmasked latent coordinates must be finite.", call. = FALSE)
      }
    } else {
      if (!identical(dim(x$latent.mask), dim(x$latent)) ||
          any(!is.finite(x$latent[x$latent.mask])) ||
          any(!is.na(x$latent[!x$latent.mask]))) {
        stop("latent and latent.mask are inconsistent.", call. = FALSE)
      }
    }
  } else if (!is.null(x$latent.mask)) {
    stop("latent.mask must be NULL when latent is NULL.", call. = FALSE)
  }
  .validate.synthetic.rng.metadata(x)
  .validate.synthetic.realized.support(x)
  invisible(x)
}

.synthetic.dataset.payload <- function(x) {
  list(
    contract.version = "synthetic-content-v1",
    dataset.id = x$dataset.id,
    specification.sha256 = x$specification.sha256,
    recipe.id = x$recipe.id,
    registry.tag = x$registry.tag,
    n = x$n,
    intrinsic.dim = x$intrinsic.dim,
    intrinsic.dim.by.region = x$intrinsic.dim.by.region,
    ambient.dim = x$ambient.dim,
    codimension = x$codimension,
    codimension.by.region = x$codimension.by.region,
    latent = x$latent,
    latent.mask = x$latent.mask,
    predictors = x$predictors,
    truth = x$truth,
    truth.scope = x$truth.scope,
    truth.estimand = x$truth.estimand,
    response = x$response,
    region = x$region,
    declared.regions = x$declared.regions,
    observed.regions = x$observed.regions,
    seed = x$seed,
    rng.policy = x$rng.policy,
    rng.streams = x$rng.streams,
    geometry.spec = x$geometry.spec,
    sampling.spec = x$sampling.spec,
    truth.spec = x$truth.spec,
    response.spec = x$response.spec,
    parameters = x$parameters
  )
}

#' Compute the canonical content checksum of a synthetic dataset
#' @param x A `synthetic_dataset`.
#' @return A lowercase SHA-256 string.
#' @examples
#' x <- materialize.synthetic(synthetic.registry.spec("G1"), n = 20L, seed = 1L)
#' synthetic.dataset.checksum(x)
#' @export
synthetic.dataset.checksum <- function(x) {
  validate.synthetic.dataset(x)
  .synthetic.sha256(.synthetic.dataset.payload(x))
}

#' Compare two synthetic datasets scientifically
#' @param x,y Synthetic datasets.
#' @param tolerance Absolute-plus-relative comparison tolerance.
#' @return A list with equality and mismatch details.
#' @examples
#' spec <- synthetic.registry.spec("G1")
#' x <- materialize.synthetic(spec, n = 20L, seed = 1L)
#' y <- materialize.synthetic(spec, n = 20L, seed = 1L)
#' compare.synthetic.dataset(x, y)
#' @export
compare.synthetic.dataset <- function(x, y, tolerance = c(1e-12, 1e-10)) {
  validate.synthetic.dataset(x)
  validate.synthetic.dataset(y)
  numeric.equal <- function(a, b) {
    if (!identical(dim(a), dim(b)) || length(a) != length(b) ||
        !identical(is.na(a), is.na(b))) {
      return(FALSE)
    }
    keep <- !is.na(a)
    all(abs(a[keep] - b[keep]) <=
          tolerance[1] + tolerance[2] * abs(b[keep]))
  }
  exact.fields <- c(
    "dataset.id", "specification.sha256", "recipe.id", "registry.tag",
    "n", "intrinsic.dim", "intrinsic.dim.by.region", "ambient.dim",
    "codimension", "codimension.by.region", "latent.mask",
    "truth.scope", "truth.estimand", "region", "declared.regions",
    "observed.regions", "seed", "rng.policy", "rng.streams",
    "geometry.spec", "sampling.spec", "truth.spec", "response.spec")
  mismatch <- exact.fields[
    !vapply(exact.fields, function(nm) identical(x[[nm]], y[[nm]]),
            logical(1))]
  numeric.fields <- c("latent", "truth", "response")
  mismatch <- c(
    mismatch,
    numeric.fields[
      !vapply(numeric.fields, function(nm) {
        if (is.null(x[[nm]]) || is.null(y[[nm]])) {
          return(identical(x[[nm]], y[[nm]]))
        }
        numeric.equal(x[[nm]], y[[nm]])
      }, logical(1))])
  random.frame <- identical(
    x$geometry.spec$parameters$frame, "random.orthonormal") &&
    identical(
      y$geometry.spec$parameters$frame, "random.orthonormal")
  orientation.invariant <- random.frame &&
    !identical(
      x$truth.spec$parameters$evaluation.coordinates, "predictors") &&
    !identical(
      y$truth.spec$parameters$evaluation.coordinates, "predictors")
  if (orientation.invariant) {
    dx <- as.matrix(stats::dist(x$predictors))^2
    dy <- as.matrix(stats::dist(y$predictors))^2
    if (!numeric.equal(dx, dy)) mismatch <- c(mismatch, "predictor.distances")
  } else if (!numeric.equal(x$predictors, y$predictors)) {
    mismatch <- c(mismatch, "predictors")
  }
  compare.parameters <- function(a, b, path = "parameters") {
    if (is.list(a) || is.list(b)) {
      if (!is.list(a) || !is.list(b) || !identical(names(a), names(b))) {
        return(path)
      }
      out <- character()
      for (name in names(a)) {
        if (orientation.invariant && path == "parameters" &&
            name == "frame.matrix") next
        out <- c(
          out,
          compare.parameters(
            a[[name]], b[[name]], paste0(path, "$", name)))
      }
      return(out)
    }
    if (is.numeric(a) || is.numeric(b)) {
      if (!is.numeric(a) || !is.numeric(b) || !numeric.equal(a, b)) {
        return(path)
      }
      return(character())
    }
    if (!identical(a, b)) path else character()
  }
  mismatch <- c(mismatch, compare.parameters(x$parameters, y$parameters))
  list(equal = !length(mismatch), mismatches = unique(mismatch))
}

#' Print a canonical synthetic dataset
#' @param x A `synthetic_dataset`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @method print synthetic_dataset
#' @export
print.synthetic_dataset <- function(x, ...) {
  cat("<synthetic_dataset> ", x$dataset.id, "\n", sep = "")
  cat("  n=", x$n, "  intrinsic=", if (is.na(x$intrinsic.dim)) "varies"
      else x$intrinsic.dim, "  ambient=", x$ambient.dim,
      "  truth=", x$truth.scope, "\n", sep = "")
  if (!is.null(x$observed.regions)) {
    cat("  observed regions: ", paste(x$observed.regions, collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a synthetic dataset as one metadata row
#' @param x A `synthetic_dataset`.
#' @param row.names Optional row names.
#' @param optional Ignored.
#' @param ... Ignored.
#' @return A one-row data frame.
#' @export
as.data.frame.synthetic_dataset <- function(
    x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    dataset.id = x$dataset.id,
    recipe.id = x$recipe.id %||% NA_character_,
    registry.tag = x$registry.tag %||% NA_character_,
    n = x$n,
    intrinsic.dim = x$intrinsic.dim,
    ambient.dim = x$ambient.dim,
    truth.scope = x$truth.scope,
    rng.policy = x$rng.policy,
    checksum = synthetic.dataset.checksum(x),
    stringsAsFactors = FALSE)
}

#' Plot a canonical synthetic dataset
#'
#' One-dimensional datasets show response and truth against the first latent
#' coordinate. Higher-dimensional datasets show the first two predictor
#' coordinates, colored by response or region.
#'
#' @param x A `synthetic_dataset`.
#' @param color Color mapping: response, truth, or region.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @method plot synthetic_dataset
#' @export
plot.synthetic_dataset <- function(
    x, color = c("response", "truth", "region"), ...) {
  validate.synthetic.dataset(x)
  color <- match.arg(color)
  if (!is.null(x$latent) && ncol(x$latent) == 1L) {
    graphics::plot(
      x$latent[, 1L], x$response,
      xlab = "latent coordinate", ylab = "response", ...)
    ord <- order(x$latent[, 1L], method = "radix")
    graphics::lines(
      x$latent[ord, 1L], x$truth[ord], col = "#D55E00", lwd = 2)
    return(invisible(x))
  }
  if (ncol(x$predictors) < 2L) {
    graphics::plot(
      seq_len(x$n), x$response,
      xlab = "observation", ylab = "response", ...)
    return(invisible(x))
  }
  group <- if (color == "region" && !is.null(x$region)) {
    as.integer(factor(x$region, levels = x$declared.regions))
  } else if (color == "truth") {
    x$truth
  } else {
    x$response
  }
  palette <- grDevices::hcl.colors(64L, "Viridis")
  if (color == "region" && !is.null(x$region)) {
    index <- pmax(1L, as.integer(group))
  } else {
    span <- range(group, finite = TRUE)
    index <- if (diff(span) == 0) rep.int(32L, length(group)) else
      1L + floor(63 * (group - span[1L]) / diff(span))
  }
  graphics::plot(
    x$predictors[, 1L], x$predictors[, 2L],
    col = palette[index], pch = 19,
    xlab = "predictor 1", ylab = "predictor 2", ...)
  invisible(x)
}
