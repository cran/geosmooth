# Serializable truth-estimand specifications and evaluators.

.new.synthetic.truth <- function(family, parameters, scope, estimand,
                                 evaluation.coordinates, version = 1L,
                                 subclass = "synthetic_truth_component") {
  scope <- match.arg(scope, c("population", "finite.design"))
  evaluation.coordinates <- match.arg(
    evaluation.coordinates,
    c("latent", "predictors", "graph.distance"))
  .new.synthetic.component(
    "truth", family,
    list(
      scope = scope,
      estimand = .synthetic.scalar.character(estimand, "estimand"),
      evaluation.coordinates = evaluation.coordinates,
      parameters = parameters
    ),
    version = version,
    subclass = c(paste0("synthetic_", gsub("\\.", "_", family), "_truth"),
                 subclass))
}

.synthetic.named.truth.registry <- function() {
  list(
    "intrinsic.smooth2.v1" = list(
      scope = "population", coordinates = "latent",
      estimand = "E[Y|U=u] = sin(pi*u1)*cos(pi*u2)"),
    "intrinsic.linear.first.v1" = list(
      scope = "population", coordinates = "latent",
      estimand = "E[Y|U=u] = u1"),
    "helix.sin.v1" = list(
      scope = "population", coordinates = "latent",
      estimand = "E[Y|T=t] = sin(t)"),
    "g4.ridge.v1" = list(
      scope = "population", coordinates = "predictors",
      estimand = "E[Y|X=x] = sin(pi*x1)"),
    "g5.smooth.v1" = list(
      scope = "population", coordinates = "predictors",
      estimand = "E[Y|X=x] = sin(pi*x1)*cos(pi*x2)"),
    "g6.sin.logit.v1" = list(
      scope = "finite.design", coordinates = "predictors",
      estimand = paste0(
        "E[Y_i|X_1:n] = clipped expit(alpha + amplitude*sin(pi*x_i1)), ",
        "with alpha solving the design prevalence")),
    "g7.smooth.v1" = list(
      scope = "population", coordinates = "predictors",
      estimand = "E[Y|X=x] = sin(pi*x1) + 0.5*x2"),
    "ssrhe.shared.smooth.v1" = list(
      scope = "finite.design", coordinates = "latent",
      estimand = paste0(
        "Shared SSRHE benchmark signal centered over the realized design"))
  )
}

#' Specify a versioned named truth evaluator
#' @param id Evaluator ID.
#' @param parameters Serializable evaluator parameters.
#' @return A synthetic truth component.
#' @examples
#' synthetic.truth.named("intrinsic.linear.first.v1")
#' @export
synthetic.truth.named <- function(id, parameters = list()) {
  id <- .synthetic.scalar.character(id, "id")
  registry <- .synthetic.named.truth.registry()
  entry <- registry[[id]]
  if (is.null(entry)) stop("Unknown named truth evaluator: ", id, call. = FALSE)
  if (.synthetic.prohibited(parameters)) {
    stop("Truth parameters must be serializable.", call. = FALSE)
  }
  .new.synthetic.truth(
    "named", list(evaluator.id = id, evaluator.parameters = parameters),
    entry$scope, entry$estimand, entry$coordinates)
}

#' Specify a polynomial truth
#' @param coefficients Named polynomial coefficients.
#' @param evaluation.coordinates Coordinates on which to evaluate.
#' @return A synthetic truth component.
#' @examples
#' synthetic.truth.polynomial(c(b0 = 1, b1 = 2, b11 = -0.5))
#' @export
synthetic.truth.polynomial <- function(
    coefficients, evaluation.coordinates = "latent") {
  coefficient.names <- names(coefficients)
  coefficients <- as.double(coefficients)
  names(coefficients) <- coefficient.names
  if (!length(coefficients) || is.null(names(coefficients)) ||
      any(!nzchar(names(coefficients))) || any(!is.finite(coefficients))) {
    stop("coefficients must be a finite named numeric vector.", call. = FALSE)
  }
  .new.synthetic.truth(
    "polynomial", list(coefficients = coefficients),
    "population", "Polynomial conditional mean",
    evaluation.coordinates)
}

#' Specify a Gaussian-mixture truth
#' @param centers Component centers, one per row.
#' @param scales Positive isotropic standard deviations or covariance matrices.
#' @param weights Nonnegative mixture weights.
#' @param normalize Normalization scope.
#' @param evaluation.coordinates Coordinates on which to evaluate.
#' @return A synthetic truth component.
#' @examples
#' synthetic.truth.gaussian.mixture(
#'   centers = rbind(-1, 1), scales = c(0.5, 0.5), weights = c(1, 1)
#' )
#' @export
synthetic.truth.gaussian.mixture <- function(
    centers, scales, weights,
    normalize = c("none", "sample.max"),
    evaluation.coordinates = "latent") {
  centers <- as.matrix(centers)
  storage.mode(centers) <- "double"
  weights <- as.double(weights)
  if (!nrow(centers) || any(!is.finite(centers)) ||
      length(weights) != nrow(centers) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("Invalid Gaussian-mixture centers or weights.", call. = FALSE)
  }
  weights <- weights / sum(weights)
  if (is.list(scales)) {
    if (length(scales) != nrow(centers)) {
      stop("Covariance list length must equal component count.",
           call. = FALSE)
    }
    scales <- lapply(scales, function(S) {
      S <- as.matrix(S)
      storage.mode(S) <- "double"
      if (!identical(dim(S), c(ncol(centers), ncol(centers))) ||
          any(!is.finite(S)) ||
          any(eigen(S, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        stop("Each covariance must be finite positive definite.",
             call. = FALSE)
      }
      dimnames(S) <- NULL
      S
    })
  } else {
    scales <- as.double(scales)
    if (length(scales) != nrow(centers) ||
        any(!is.finite(scales)) || any(scales <= 0)) {
      stop("Isotropic scales must be positive with one per component.",
           call. = FALSE)
    }
  }
  normalize <- match.arg(normalize)
  .new.synthetic.truth(
    "gaussian.mixture",
    list(centers = centers, scales = scales, weights = weights,
         normalize = normalize),
    if (normalize == "sample.max") "finite.design" else "population",
    if (normalize == "sample.max") {
      "Gaussian-mixture signal normalized by its realized design maximum"
    } else {
      "Unnormalized Gaussian-mixture signal"
    },
    evaluation.coordinates)
}

#' Specify a finite-design occupation-mixture truth
#'
#' The analytic Gaussian-mixture density is transformed to Bernoulli
#' occupation probabilities by
#' `probability.maximum * density^gamma / max(density^gamma)`.
#'
#' @param centers Component centers, one per row.
#' @param covariances Positive-definite covariance matrices.
#' @param weights Nonnegative mixture weights.
#' @param gamma Positive density-shape exponent.
#' @param probability.maximum Maximum realized occupation probability.
#' @param normalization Currently `"design.maximum"`.
#' @return A finite-design synthetic truth component.
#' @examples
#' synthetic.truth.occupation.mixture(
#'   centers = rbind(c(-1, 0), c(1, 0)),
#'   covariances = list(diag(2), diag(2)), weights = c(1, 1)
#' )
#' @export
synthetic.truth.occupation.mixture <- function(
    centers, covariances, weights, gamma = 1,
    probability.maximum = 0.65,
    normalization = "design.maximum") {
  mixture <- synthetic.truth.gaussian.mixture(
    centers = centers, scales = covariances, weights = weights,
    normalize = "none", evaluation.coordinates = "predictors")
  gamma <- .synthetic.scalar.double(gamma, "gamma", 0, strict = TRUE)
  probability.maximum <- .synthetic.scalar.double(
    probability.maximum, "probability.maximum", 0, strict = TRUE)
  if (probability.maximum > 1) {
    stop("probability.maximum must not exceed one.", call. = FALSE)
  }
  normalization <- .synthetic.scalar.character(
    normalization, "normalization")
  if (normalization != "design.maximum") {
    stop("Only design.maximum normalization is implemented.", call. = FALSE)
  }
  .new.synthetic.truth(
    "occupation.mixture",
    c(
      mixture$parameters$parameters,
      list(
        gamma = gamma,
        probability.maximum = probability.maximum,
        normalization = normalization)),
    "finite.design",
    paste0(
      "Bernoulli occupation probability obtained by design-maximum ",
      "normalization of a Gaussian-mixture density"),
    "predictors")
}

#' Specify a sinusoidal finite-design logit truth
#' @param amplitude Sinusoidal log-odds amplitude.
#' @param target.prevalence Target mean pre-clipping probability.
#' @param clip Probability bounds.
#' @return A synthetic truth component.
#' @examples
#' synthetic.truth.logit(amplitude = 1, target.prevalence = 0.4)
#' @export
synthetic.truth.logit <- function(
    amplitude = 1.5, target.prevalence = 0.5,
    clip = c(0.05, 0.95)) {
  amplitude <- .synthetic.scalar.double(amplitude, "amplitude")
  target.prevalence <- .synthetic.scalar.double(
    target.prevalence, "target.prevalence")
  clip <- as.double(clip)
  if (length(clip) != 2L || any(!is.finite(clip)) ||
      clip[1] >= clip[2] || target.prevalence <= clip[1] ||
      target.prevalence >= clip[2]) {
    stop("clip and target.prevalence are incompatible.", call. = FALSE)
  }
  synthetic.truth.named(
    "g6.sin.logit.v1",
    list(amplitude = amplitude, target.prevalence = target.prevalence,
         clip = clip))
}

.evaluate.synthetic.polynomial <- function(X, coefficients) {
  n <- nrow(X)
  value <- rep(unname(coefficients["b0"] %||% 0), n)
  for (j in seq_len(ncol(X))) {
    nm <- paste0("b", j)
    if (nm %in% names(coefficients)) value <- value + coefficients[nm] * X[, j]
  }
  for (j in seq_len(ncol(X))) {
    for (k in j:ncol(X)) {
      nm <- paste0("b", j, k)
      if (nm %in% names(coefficients)) {
        value <- value + coefficients[nm] * X[, j] * X[, k]
      }
    }
  }
  as.numeric(value)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L ||
      length(x) == 1L && is.na(x)) y else x
}

.evaluate.synthetic.named.truth <- function(id, pars, latent, predictors) {
  if (id == "intrinsic.smooth2.v1") {
    return(list(value = sin(pi * latent[, 1]) * cos(pi * latent[, 2]),
                parameters = list()))
  }
  if (id == "intrinsic.linear.first.v1") {
    return(list(value = latent[, 1], parameters = list()))
  }
  if (id == "helix.sin.v1") {
    return(list(value = sin(latent[, 1]), parameters = list()))
  }
  if (id == "g4.ridge.v1") {
    return(list(value = sin(pi * predictors[, 1]), parameters = list()))
  }
  if (id == "g5.smooth.v1") {
    return(list(
      value = sin(pi * predictors[, 1]) * cos(pi * predictors[, 2]),
      parameters = list()))
  }
  if (id == "g6.sin.logit.v1") {
    amplitude <- pars$amplitude %||% 1.5
    prevalence <- pars$target.prevalence %||% 0.5
    clip <- pars$clip %||% c(0.05, 0.95)
    eta <- amplitude * sin(pi * predictors[, 1])
    mean.p <- function(alpha) mean(stats::plogis(alpha + eta))
    alpha <- stats::uniroot(
      function(a) mean.p(a) - prevalence, c(-50, 50))$root
    value <- pmin(clip[2], pmax(clip[1], stats::plogis(alpha + eta)))
    return(list(
      value = value,
      parameters = list(alpha = alpha, realized.prevalence = mean(value))))
  }
  if (id == "g7.smooth.v1") {
    return(list(
      value = sin(pi * predictors[, 1]) + 0.5 * predictors[, 2],
      parameters = list()))
  }
  if (id == "ssrhe.shared.smooth.v1") {
    d <- ncol(latent)
    value <- sin(pi * latent[, 1])
    if (d >= 2L) {
      value <- value + 0.7 * cos(pi * latent[, 2]) +
        0.45 * latent[, 1]^2 -
        0.25 * latent[, 1] * latent[, 2]
    }
    if (d >= 3L) {
      value <- value + 0.35 * sin(pi * latent[, 3] / 2) +
        0.15 * latent[, 2] * latent[, 3]
    }
    if (d >= 4L) {
      value <- value + 0.20 * cos(pi * latent[, 4] / 2) -
        0.10 * latent[, 1] * latent[, 4]
    }
    value <- as.numeric(scale(value, center = TRUE, scale = FALSE))
    return(list(value = value, parameters = list()))
  }
  stop("No evaluator implementation for named truth: ", id, call. = FALSE)
}

.gaussian.density.matrix <- function(X, center, scale) {
  d <- ncol(X)
  if (d == 1L && length(scale) == 1L) {
    return(stats::dnorm(X[, 1], mean = center[1], sd = scale))
  }
  delta <- sweep(X, 2L, center, "-")
  if (length(scale) == 1L) {
    return(exp(-rowSums(delta^2) / (2 * scale^2)) /
             ((2 * pi)^(d / 2) * scale^d))
  }
  inv <- solve(scale)
  quad <- rowSums((delta %*% inv) * delta)
  determinant.value <- exp(as.numeric(determinant(scale)$modulus))
  exp(-quad / 2) / sqrt((2 * pi)^d * determinant.value)
}

.evaluate.synthetic.truth <- function(truth, latent, predictors) {
  p <- truth$parameters
  coordinates <- if (p$evaluation.coordinates == "latent") latent else predictors
  if (truth$family == "named") {
    out <- .evaluate.synthetic.named.truth(
      p$parameters$evaluator.id,
      p$parameters$evaluator.parameters,
      latent, predictors)
  } else if (truth$family == "polynomial") {
    out <- list(
      value = .evaluate.synthetic.polynomial(
        coordinates, p$parameters$coefficients),
      parameters = list())
  } else if (truth$family == "gaussian.mixture") {
    pars <- p$parameters
    dens <- vapply(seq_len(nrow(pars$centers)), function(k) {
      scale <- if (is.list(pars$scales)) pars$scales[[k]] else pars$scales[k]
      .gaussian.density.matrix(coordinates, pars$centers[k, ], scale)
    }, numeric(nrow(coordinates)))
    value <- as.numeric(dens %*% pars$weights)
    extra <- list()
    if (pars$normalize == "sample.max") {
      normalizer <- max(value)
      value <- value / normalizer
      extra$truth.normalizer <- normalizer
    }
    out <- list(value = value, parameters = extra)
  } else if (truth$family == "occupation.mixture") {
    pars <- p$parameters
    dens <- vapply(seq_len(nrow(pars$centers)), function(k) {
      .gaussian.density.matrix(
        coordinates, pars$centers[k, ], pars$scales[[k]])
    }, numeric(nrow(coordinates)))
    density.value <- as.numeric(dens %*% pars$weights)
    shaped <- density.value^pars$gamma
    normalizer <- max(shaped)
    if (!is.finite(normalizer) || normalizer <= 0) {
      stop("Occupation mixture has a nonpositive design maximum.",
           call. = FALSE)
    }
    value <- pars$probability.maximum * shaped / normalizer
    out <- list(
      value = value,
      parameters = list(
        analytic.density = density.value,
        design.maximum = normalizer,
        probability.maximum = pars$probability.maximum))
  } else {
    stop("Unsupported truth family: ", truth$family, call. = FALSE)
  }
  if (length(out$value) != nrow(predictors) || any(!is.finite(out$value))) {
    stop("Truth evaluator returned invalid values.", call. = FALSE)
  }
  out$value <- as.numeric(out$value)
  out
}
