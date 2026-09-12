# Response component specifications and draw algorithms.

.new.synthetic.response <- function(family, parameters, target,
                                    version = 1L) {
  .new.synthetic.component(
    "response", family,
    list(target = target, parameters = parameters),
    version = version,
    subclass = c(paste0("synthetic_", gsub("\\.", "_", family), "_response"),
                 "synthetic_response_component"))
}

#' Gaussian response specification
#' @param sd Response standard deviation.
#' @return A synthetic response component.
#' @examples
#' synthetic.response.gaussian(sd = 0.2)
#' @export
synthetic.response.gaussian <- function(sd) {
  sd <- .synthetic.scalar.double(sd, "sd", 0)
  .new.synthetic.response(
    "gaussian", list(sd = sd), "conditional.mean")
}

#' Heteroskedastic Gaussian response specification
#' @param base.sd Base standard deviation.
#' @param truth.multiplier Multiplier applied to truth.
#' @return A synthetic response component.
#' @examples
#' synthetic.response.heteroskedastic.gaussian(0.1, 0.2)
#' @export
synthetic.response.heteroskedastic.gaussian <- function(
    base.sd, truth.multiplier) {
  base.sd <- .synthetic.scalar.double(base.sd, "base.sd", 0)
  truth.multiplier <- .synthetic.scalar.double(
    truth.multiplier, "truth.multiplier")
  .new.synthetic.response(
    "heteroskedastic.gaussian",
    list(base.sd = base.sd, truth.multiplier = truth.multiplier),
    "conditional.mean")
}

#' Laplace response with Gaussian contamination
#' @param laplace.scale Laplace scale.
#' @param outlier.fraction Fraction contaminated.
#' @param outlier.sd Contamination standard deviation.
#' @param minimum.outliers Minimum contamination count.
#' @param count.rounding Count policy.
#' @param laplace.algorithm Versioned Laplace algorithm.
#' @return A synthetic response component.
#' @examples
#' synthetic.response.laplace.outlier(0.1, 0.05, 1)
#' @export
synthetic.response.laplace.outlier <- function(
    laplace.scale, outlier.fraction, outlier.sd,
    minimum.outliers = 0L, count.rounding = "round",
    laplace.algorithm = "uniform.inverse.v1") {
  laplace.scale <- .synthetic.scalar.double(
    laplace.scale, "laplace.scale", 0, strict = TRUE)
  outlier.fraction <- .synthetic.scalar.double(
    outlier.fraction, "outlier.fraction")
  if (outlier.fraction < 0 || outlier.fraction > 1) {
    stop("outlier.fraction must lie in [0, 1].", call. = FALSE)
  }
  outlier.sd <- .synthetic.scalar.double(outlier.sd, "outlier.sd", 0)
  minimum.outliers <- .synthetic.scalar.integer(
    minimum.outliers, "minimum.outliers", 0L)
  .new.synthetic.response(
    "laplace.outlier",
    list(laplace.scale = laplace.scale,
         outlier.fraction = outlier.fraction,
         outlier.sd = outlier.sd,
         minimum.outliers = minimum.outliers,
         count.rounding = count.rounding,
         laplace.algorithm = laplace.algorithm),
    "conditional.mean")
}

#' Clustered Gaussian response specification
#' @param residual.sd Residual standard deviation.
#' @param intraclass.correlation Intraclass correlation in `[0,1)`.
#' @return A synthetic response component.
#' @examples
#' synthetic.response.clustered.gaussian(0.2, 0.3)
#' @export
synthetic.response.clustered.gaussian <- function(
    residual.sd, intraclass.correlation) {
  residual.sd <- .synthetic.scalar.double(
    residual.sd, "residual.sd", 0, strict = TRUE)
  intraclass.correlation <- .synthetic.scalar.double(
    intraclass.correlation, "intraclass.correlation")
  if (intraclass.correlation < 0 || intraclass.correlation >= 1) {
    stop("intraclass.correlation must lie in [0, 1).", call. = FALSE)
  }
  .new.synthetic.response(
    "clustered.gaussian",
    list(residual.sd = residual.sd,
         intraclass.correlation = intraclass.correlation),
    "conditional.mean")
}

#' Bernoulli response specification
#' @param minimum.positive Minimum accepted positive responses.
#' @param maximum.attempts Maximum complete redraw attempts.
#' @return A synthetic response component.
#' @examples
#' synthetic.response.bernoulli(minimum.positive = 1L, maximum.attempts = 5L)
#' @export
synthetic.response.bernoulli <- function(
    minimum.positive = 0L, maximum.attempts = 1L) {
  minimum.positive <- .synthetic.scalar.integer(
    minimum.positive, "minimum.positive", 0L)
  maximum.attempts <- .synthetic.scalar.integer(
    maximum.attempts, "maximum.attempts", 1L)
  .new.synthetic.response(
    "bernoulli",
    list(minimum.positive = minimum.positive,
         maximum.attempts = maximum.attempts),
    "conditional.probability")
}

.synthetic.response.uses.rng <- function(response) {
  p <- response$parameters$parameters
  switch(
    response$family,
    gaussian = p$sd > 0,
    heteroskedastic.gaussian = TRUE,
    laplace.outlier = TRUE,
    clustered.gaussian = TRUE,
    bernoulli = TRUE,
    TRUE)
}

.draw.synthetic.response <- function(response, truth, region = NULL) {
  p <- response$parameters$parameters
  family <- response$family
  n <- length(truth)
  parameters <- list()
  if (family == "gaussian") {
    value <- if (p$sd == 0) truth else truth + stats::rnorm(n, sd = p$sd)
  } else if (family == "heteroskedastic.gaussian") {
    sd <- p$base.sd + p$truth.multiplier * truth
    if (any(!is.finite(sd)) || any(sd < 0)) {
      stop("Heteroskedastic response produced invalid standard deviations.",
           call. = FALSE)
    }
    value <- truth + stats::rnorm(n, sd = sd)
    parameters$response.sd <- sd
  } else if (family == "laplace.outlier") {
    u <- stats::runif(n) - 0.5
    error <- -p$laplace.scale * sign(u) * log1p(-2 * abs(u))
    count <- max(
      p$minimum.outliers,
      if (p$count.rounding == "round") round(p$outlier.fraction * n)
      else floor(p$outlier.fraction * n))
    index <- if (count) sample.int(n, count) else integer()
    if (length(index)) {
      error[index] <- error[index] + stats::rnorm(
        length(index), sd = p$outlier.sd)
    }
    value <- truth + error
    parameters$outlier.index <- index
  } else if (family == "clustered.gaussian") {
    if (is.null(region)) {
      stop("Clustered response requires cluster region labels.",
           call. = FALSE)
    }
    cluster <- as.integer(region)
    K <- max(cluster)
    tau <- sqrt(p$intraclass.correlation * p$residual.sd^2 /
                  (1 - p$intraclass.correlation))
    b <- stats::rnorm(K, sd = tau)
    eps <- stats::rnorm(n, sd = p$residual.sd)
    value <- truth + b[cluster] + eps
    parameters <- list(tau = tau, cluster.effect = b)
  } else if (family == "bernoulli") {
    if (any(truth < 0 | truth > 1)) {
      stop("Bernoulli truth must lie in [0, 1].", call. = FALSE)
    }
    value <- stats::rbinom(n, size = 1L, prob = truth)
  } else {
    stop("Unsupported response family: ", family, call. = FALSE)
  }
  list(value = as.numeric(value), parameters = parameters)
}
