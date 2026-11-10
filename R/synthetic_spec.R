# Synthetic data component and specification infrastructure.

.synthetic.scalar.character <- function(x, name, null.ok = FALSE) {
  if (null.ok && is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(name, " must be one nonempty character value.", call. = FALSE)
  }
  x
}

.synthetic.scalar.integer <- function(x, name, lower = NULL) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) ||
      x != trunc(x)) {
    stop(name, " must be one finite whole number.", call. = FALSE)
  }
  if (x < -.Machine$integer.max || x > .Machine$integer.max) {
    stop(name, " is outside R's representable integer seed range.",
         call. = FALSE)
  }
  value <- as.integer(x)
  if (!is.null(lower) && value < lower) {
    stop(name, " must be at least ", lower, ".", call. = FALSE)
  }
  value
}

.synthetic.scalar.double <- function(x, name, lower = NULL,
                                     strict = FALSE) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x)) {
    stop(name, " must be one finite numeric value.", call. = FALSE)
  }
  value <- as.double(x)
  if (!is.null(lower) &&
      if (strict) value <= lower else value < lower) {
    relation <- if (strict) "greater than" else "at least"
    stop(name, " must be ", relation, " ", lower, ".", call. = FALSE)
  }
  value
}

.synthetic.prohibited <- function(x) {
  if (is.function(x) || is.environment(x) || is.language(x) ||
      typeof(x) %in% c("externalptr", "weakref")) {
    return(TRUE)
  }
  if (is.list(x)) return(any(vapply(x, .synthetic.prohibited, logical(1))))
  FALSE
}

.synthetic.canonicalize <- function(x, field = NULL) {
  if (.synthetic.prohibited(x)) {
    stop("Synthetic specifications cannot contain functions, environments, ",
         "language objects, or external pointers.", call. = FALSE)
  }
  if (is.factor(x)) return(as.character(x))
  if (is.matrix(x) || is.array(x)) {
    dimnames(x) <- NULL
    return(x)
  }
  if (is.list(x)) {
    x <- unclass(x)
    if (!is.null(names(x))) {
      Encoding(names(x)) <- "UTF-8"
      x <- x[order(names(x), method = "radix")]
      return(stats::setNames(
        Map(
          function(value, name) .synthetic.canonicalize(value, name),
          x, names(x)),
        names(x)))
    }
    return(lapply(x, .synthetic.canonicalize, field = NULL))
  }
  if (is.character(x)) {
    Encoding(x) <- "UTF-8"
  }
  if (is.atomic(x) && !is.null(names(x))) {
    Encoding(names(x)) <- "UTF-8"
    contractual.fields <- c(
      "intrinsic.dim.by.region", "codimension.by.region")
    if (is.null(field) || !field %in% contractual.fields) {
      x <- x[order(names(x), method = "radix")]
    }
  }
  x
}

.synthetic.normalize.xdr.v3.header <- function(raw) {
  if (!is.raw(raw) || length(raw) < 10L ||
      !identical(raw[1:6], as.raw(c(0x58, 0x0a, 0, 0, 0, 3)))) {
    stop("Synthetic checksum payload is not XDR serialization version 3.",
         call. = FALSE)
  }
  # XDR serialization records the producing R version in bytes 7--10 even
  # when the serialized payload is otherwise byte-for-byte identical.  Fix
  # those non-payload bytes to the historical registry writer (R 4.7.0) so
  # existing identities remain valid and can be verified by older R releases.
  raw[7:10] <- as.raw(c(0, 4, 7, 0))
  raw
}

.synthetic.xdr.v3 <- function(x) {
  payload <- .synthetic.canonicalize(x)
  # A round trip removes incidental reference-sharing differences between
  # otherwise identical R objects before the contractual XDR serialization.
  normalized <- unserialize(
    serialize(payload, NULL, ascii = FALSE, xdr = TRUE, version = 3))
  .synthetic.normalize.xdr.v3.header(
    serialize(normalized, NULL, ascii = FALSE, xdr = TRUE, version = 3))
}

.synthetic.sha256 <- function(x) {
  digest::digest(
    .synthetic.xdr.v3(x), algo = "sha256", serialize = FALSE)
}

.new.synthetic.component <- function(kind, family, parameters,
                                     version = 1L, subclass) {
  kind <- .synthetic.scalar.character(kind, "kind")
  family <- .synthetic.scalar.character(family, "family")
  version <- .synthetic.scalar.integer(version, "version", 1L)
  if (.synthetic.prohibited(parameters)) {
    stop("Component parameters must be canonically serializable.",
         call. = FALSE)
  }
  structure(
    list(
      kind = kind,
      family = family,
      version = version,
      parameters = parameters
    ),
    class = c(subclass, paste0("synthetic_", kind), "synthetic_component")
  )
}

.is.synthetic.component <- function(x, kind) {
  inherits(x, paste0("synthetic_", kind)) &&
    identical(x$kind, kind)
}

#' Assemble a synthetic-data specification
#'
#' Combines draw-free geometry, sampling, truth, and response components into
#' one canonically hashable specification.
#'
#' @param geometry A synthetic geometry component.
#' @param sampling A synthetic sampling component.
#' @param truth A synthetic truth component.
#' @param response A synthetic response component.
#' @param recipe.id Optional readable recipe label.
#' @param registry.tag Optional stable experimental-plan tag.
#' @param compatibility Optional serializable legacy compatibility metadata.
#' @param metadata Optional serializable scientific metadata.
#' @return A `synthetic_spec`.
#' @examples
#' spec <- synthetic.spec(
#'   dgraphs::synthetic.circle(), dgraphs::synthetic.sampling.uniform.interval(0, 2 * pi),
#'   synthetic.truth.named("helix.sin.v1"), synthetic.response.gaussian(0.1)
#' )
#' spec$specification.sha256
#' @export
synthetic.spec <- function(geometry, sampling, truth, response,
                           recipe.id = NULL, registry.tag = NULL,
                           compatibility = NULL, metadata = list()) {
  if (!.is.synthetic.component(geometry, "geometry")) {
    stop("geometry must be a synthetic geometry component.", call. = FALSE)
  }
  if (!.is.synthetic.component(sampling, "sampling")) {
    stop("sampling must be a synthetic sampling component.", call. = FALSE)
  }
  if (!.is.synthetic.component(truth, "truth")) {
    stop("truth must be a synthetic truth component.", call. = FALSE)
  }
  if (!.is.synthetic.component(response, "response")) {
    stop("response must be a synthetic response component.", call. = FALSE)
  }
  recipe.id <- .synthetic.scalar.character(
    recipe.id, "recipe.id", null.ok = TRUE)
  registry.tag <- .synthetic.scalar.character(
    registry.tag, "registry.tag", null.ok = TRUE)
  if (.synthetic.prohibited(compatibility) ||
      .synthetic.prohibited(metadata)) {
    stop("compatibility and metadata must be canonically serializable.",
         call. = FALSE)
  }
  payload <- list(
    geometry.spec = geometry,
    sampling.spec = sampling,
    truth.spec = truth,
    response.spec = response,
    compatibility = compatibility,
    metadata = metadata
  )
  specification.sha256 <- .synthetic.sha256(payload)
  structure(
    c(payload, list(
      recipe.id = recipe.id,
      registry.tag = registry.tag,
      specification.sha256 = specification.sha256
    )),
    class = c("synthetic_spec", "list")
  )
}

.validate.synthetic.seed <- function(seed, name = "seed") {
  .synthetic.scalar.integer(seed, name, lower = 0L)
}

.derive.synthetic.seed <- function(seed, offset, name) {
  seed <- .validate.synthetic.seed(seed)
  offset <- .synthetic.scalar.integer(offset, paste0(name, " offset"))
  derived <- as.double(seed) + as.double(offset)
  .synthetic.scalar.integer(derived, name)
}

.validate.synthetic.seed.plan <- function(seed, frame.seed = NULL,
                                          need.frame = FALSE,
                                          need.response = FALSE) {
  seed <- .validate.synthetic.seed(seed)
  out <- list(seed = seed)
  if (need.frame) {
    out$frame.seed <- if (is.null(frame.seed)) {
      .derive.synthetic.seed(seed, 500000L, "derived frame seed")
    } else {
      .validate.synthetic.seed(frame.seed, "frame.seed")
    }
  }
  if (need.response) {
    out$response.seed <- .derive.synthetic.seed(
      seed, 900000L, "derived response seed")
  }
  out
}

.with.synthetic.rng.preserved <- function(code) {
  kinds <- RNGkind()
  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(kinds))
    if (had.seed) {
      assign(".Random.seed", old.seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(code)
}
