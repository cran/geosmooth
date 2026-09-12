# Canonical synthetic-data materialization and RNG policy implementation.

.synthetic.named.streams <- function(seed) {
  RNGversion(.synthetic.rng.version())
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  names <- c("sampling", "geometry.frame", "truth",
             "response", "response.auxiliary")
  states <- vector("list", length(names))
  states[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
  if (length(states) > 1L) {
    for (i in 2:length(states)) {
      states[[i]] <- parallel::nextRNGStream(states[[i - 1L]])
    }
  }
  stats::setNames(states, names)
}

.synthetic.state.for.attempt <- function(state, attempt) {
  if (attempt <= 1L) return(state)
  for (i in 2:attempt) state <- parallel::nextRNGSubStream(state)
  state
}

.with.synthetic.state <- function(state, code) {
  assign(".Random.seed", state, envir = .GlobalEnv)
  force(code)
}

.resolve.synthetic.n <- function(sampling, n) {
  fixed <- .synthetic.fixed.n(sampling)
  if (!is.null(fixed)) {
    if (is.null(n)) return(as.integer(fixed))
    n <- .synthetic.scalar.integer(n, "n", 1L)
    if (n != fixed) {
      stop("n disagrees with the sampling component's fixed size.",
           call. = FALSE)
    }
    return(n)
  }
  if (is.null(n)) stop("n is required for this sampling component.",
                       call. = FALSE)
  .synthetic.scalar.integer(n, "n", 1L)
}

.decorate.ssrhe.dataset <- function(object, spec) {
  metadata <- spec$metadata
  recipe.id <- spec$recipe.id %||% ""
  if (grepl("^S[0-9]{2}\\.V[1-3]$", recipe.id)) {
    object$domain <- "1D mixture"
    object$family <- metadata$shape.name
    object$variant <- metadata$variant.name
    object$dim <- 1L
    object$x <- object$predictors[, 1L]
    object$index.k <- NA_integer_
    object$curvature <- NA_real_
    object$coefficients <- NA_real_
  } else if (grepl("^ssrhe\\.(flat|quadform)\\.", recipe.id)) {
    object$domain <- if (metadata$family == "flat") {
      sprintf("%dD flat", metadata$intrinsic.dim)
    } else {
      sprintf("%dD quadform", metadata$intrinsic.dim)
    }
    object$family <- if (metadata$family == "flat") {
      "flat"
    } else {
      sprintf("quadform index %d", metadata$index.k)
    }
    object$variant <- if (metadata$family == "flat") {
      "uniform ambient sample"
    } else {
      sprintf("curvature %.2f", metadata$curvature)
    }
    object$dim <- as.integer(metadata$intrinsic.dim)
    object$x <- NA_real_
    object$index.k <- metadata$index.k %||% NA_integer_
    object$curvature <- metadata$curvature %||% NA_real_
    object$coefficients <- metadata$coefficients %||%
      rep(NA_real_, metadata$intrinsic.dim)
  }
  object
}

.materialize.synthetic.once <- function(
    spec, n, seed.plan, rng.policy, streams = NULL, attempt = 1L) {
  geometry <- spec$geometry.spec
  response <- spec$response.spec
  random.frame <- identical(
    geometry$parameters$frame, "random.orthonormal")
  if (rng.policy == "legacy") {
    frame.matrix <- if (random.frame) {
      .draw.synthetic.frame(geometry, seed.plan$frame.seed)
    } else {
      .draw.synthetic.frame(geometry)
    }
    set.seed(seed.plan$seed)
    sample.out <- .draw.synthetic.sampling(
      spec$sampling.spec, n, geometry)
  } else {
    sample.state <- .synthetic.state.for.attempt(
      streams$sampling, attempt)
    sample.out <- .with.synthetic.state(
      sample.state,
      .draw.synthetic.sampling(spec$sampling.spec, n, geometry))
    if (random.frame) {
      frame.state <- .synthetic.state.for.attempt(
        streams$geometry.frame, attempt)
      frame.matrix <- .with.synthetic.state(
        frame.state, .draw.synthetic.frame(geometry))
    } else {
      frame.matrix <- .draw.synthetic.frame(geometry)
    }
  }
  tag <- spec$registry.tag %||% ""
  if (tag == "G3d" && !is.null(sample.out$latent)) {
    colnames(sample.out$latent) <- c("u", "v")
  } else if (tag == "G4" && !is.null(sample.out$latent)) {
    colnames(sample.out$latent) <- c("tA", "")
  } else if (tag %in% c("G5", "G6") && !is.null(sample.out$latent)) {
    colnames(sample.out$latent) <- c("x1", "x2")
  }
  predictors <- if (is.null(sample.out$predictors)) {
    .embed.synthetic.geometry(
      geometry, sample.out$latent, frame.matrix = frame.matrix)
  } else {
    sample.out$predictors
  }
  if (tag %in% c("G3a", "G3b", "G3c", "G3d", "G4")) {
    colnames(predictors) <- c("x1", "x2", "x3")
  } else if (tag %in% c("G5", "G6")) {
    colnames(predictors) <- c("x1", "x2")
  } else if (tag == "G7") {
    colnames(predictors) <- paste0("p", seq_len(ncol(predictors)))
  }
  sample.out$predictors <- predictors
  truth.out <- .evaluate.synthetic.truth(
    spec$truth.spec, sample.out$latent, predictors)
  if (.synthetic.response.uses.rng(response)) {
    if (rng.policy == "legacy") {
      if (!is.null(seed.plan$response.seed)) {
        set.seed(seed.plan$response.seed)
      }
      response.out <- .draw.synthetic.response(
        response, truth.out$value, sample.out$region)
    } else {
      response.state <- .synthetic.state.for.attempt(
        streams$response, attempt)
      response.out <- .with.synthetic.state(
        response.state,
        .draw.synthetic.response(
          response, truth.out$value, sample.out$region))
    }
  } else {
    response.out <- .draw.synthetic.response(
      response, truth.out$value, sample.out$region)
  }
  list(sample = sample.out, frame = frame.matrix, truth = truth.out,
       response = response.out)
}

#' Materialize a synthetic-data specification
#'
#' @param spec A `synthetic_spec`.
#' @param n Sample size, or `NULL` for fixed-size sampling.
#' @param seed Scalar whole-number seed.
#' @param rng.policy Versioned RNG policy.
#' @param validate Whether to validate the returned dataset.
#' @return A `synthetic_dataset`.
#' @examples
#' spec <- synthetic.registry.spec("G1")
#' materialize.synthetic(spec, n = 20L, seed = 1L)
#' @export
materialize.synthetic <- function(
    spec, n = NULL, seed,
    rng.policy = c("named.stream.v1", "legacy"), validate = TRUE) {
  if (!inherits(spec, "synthetic_spec")) {
    stop("spec must be a synthetic_spec.", call. = FALSE)
  }
  rng.policy <- match.arg(rng.policy)
  n <- .resolve.synthetic.n(spec$sampling.spec, n)
  .validate.synthetic.support(spec$geometry.spec, spec$sampling.spec)
  random.frame <- identical(
    spec$geometry.spec$parameters$frame, "random.orthonormal")
  response.random <- .synthetic.response.uses.rng(spec$response.spec)
  frame.seed <- NULL
  seed.policy <- "named.stream.v1"
  if (!is.null(spec$compatibility)) {
    frame.seed <- spec$compatibility$frame.seed
    seed.policy <- spec$compatibility$seed.policy %||% seed.policy
  }
  # This resolves and validates every policy-derived seed before RNG state
  # changes or any draw is consumed.
  seed.plan <- .validate.synthetic.seed.plan(
    seed, frame.seed = frame.seed,
    need.frame = rng.policy == "legacy" && random.frame,
    need.response = rng.policy == "legacy" && response.random &&
      !identical(seed.policy, "ssrhe.1d.v1"))
  seed <- seed.plan$seed
  result <- .with.synthetic.rng.preserved({
    RNGversion(.synthetic.rng.version())
    if (rng.policy == "legacy") {
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    } else {
      RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    }
    streams <- if (rng.policy == "named.stream.v1") {
      .synthetic.named.streams(seed)
    } else {
      list(
        sampling = seed,
        geometry.frame = seed.plan$frame.seed,
        response = seed.plan$response.seed)
    }
    maximum.attempts <- if (spec$response.spec$family == "bernoulli") {
      spec$response.spec$parameters$parameters$maximum.attempts
    } else 1L
    minimum.positive <- if (spec$response.spec$family == "bernoulli") {
      spec$response.spec$parameters$parameters$minimum.positive
    } else 0L
    materialized <- NULL
    for (attempt in seq_len(maximum.attempts)) {
      materialized <- .materialize.synthetic.once(
        spec, n, seed.plan, rng.policy, streams, attempt)
      if (sum(materialized$response$value > 0) >= minimum.positive) break
    }
    if (sum(materialized$response$value > 0) < minimum.positive) {
      stop("Bernoulli minimum.positive was not reached within ",
           maximum.attempts, " attempts.", call. = FALSE)
    }
    list(materialized = materialized, streams = streams,
         attempts = attempt)
  })
  mat <- result$materialized
  object <- .new.synthetic.dataset(
    spec, n, seed, rng.policy, result$streams,
    mat$sample, mat$frame, mat$truth, mat$response,
    parameters = list(
      materialization.attempts = result$attempts,
      resolved.seeds = seed.plan),
    provenance = list(
      generator = "materialize.synthetic",
      geosmooth.version = tryCatch(
        as.character(utils::packageVersion("geosmooth")),
        error = function(e) NA_character_)))
  object <- .decorate.ssrhe.dataset(object, spec)
  if (validate) validate.synthetic.dataset(object)
  object
}
