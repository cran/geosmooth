# Retained legacy mixed-stratum sampling.

# Sampling component specifications and exact draw algorithms.

.new.synthetic.sampling <- function(family, parameters, version = 1L,
                                    subclass = "synthetic_sampling_component") {
  .new.synthetic.component(
    "sampling", family, parameters, version,
    subclass = c(paste0("synthetic_", gsub("\\.", "_", family), "_sampling"),
                 subclass))
}

#' Legacy G4 stratified sampling specification
#' @param fraction.a Fraction allocated to stratum A.
#' @param algorithm Versioned traversal algorithm.
#' @return A synthetic sampling component.
#' @examples
#' synthetic.sampling.stratified(fraction.a = 0.6)
#' @export
synthetic.sampling.stratified <- function(
    fraction.a = 0.5,
    algorithm = c("stratum.sequential.v1",
                  "legacy.g4.stratum.sequential.v1")) {
  fraction.a <- .synthetic.scalar.double(fraction.a, "fraction.a")
  if (fraction.a <= 0 || fraction.a >= 1) {
    stop("fraction.a must lie strictly between zero and one.", call. = FALSE)
  }
  .new.synthetic.sampling(
    "stratified.g4",
    list(fraction.a = fraction.a, algorithm = match.arg(algorithm)))
}

.draw.synthetic.g4 <- function(sampling, n, geometry) {
  p <- sampling$parameters
  gp <- geometry$parameters
  n <- .synthetic.scalar.integer(n, "n", 1L)
  out <- list(latent = NULL, predictors = NULL, latent.mask = NULL,
              region = NULL, parameters = list())
    nA <- max(1L, round(p$fraction.a * n))
    nB <- n - nA
    if (nA < 1L || nB < 1L) {
      stop("G4 requires at least one row in both declared strata.",
           call. = FALSE)
    }
    eta <- gp$eta
    tA <- stats::runif(nA, -1, 0)
    A <- cbind(tA, stats::rnorm(nA, sd = eta),
               stats::rnorm(nA, sd = eta))
    s1 <- stats::runif(nB, 0, 1)
    s2 <- stats::runif(nB, -0.5, 0.5)
    B <- cbind(s1, s2, 0)
    out$predictors <- rbind(A, B)
    out$latent <- rbind(cbind(tA, NA_real_), cbind(s1, s2))
    out$latent.mask <- !is.na(out$latent)
    out$region <- c(rep("A", nA), rep("B", nB))
    out$parameters$nA <- nA
    out$parameters$nB <- nB
  out
}

.synthetic.fixed.n <- function(sampling) sampling$parameters$fixed.n
.validate.synthetic.support <- function(geometry, sampling) {
  if (identical(geometry$family, "g4.segment.rectangle")) {
    if (!identical(sampling$family, "stratified.g4"))
      stop("G4 stratified geometry requires G4 stratified sampling.", call. = FALSE)
    return(invisible(TRUE))
  }
  dgraphs::validate.synthetic.sampling(sampling, geometry)
  invisible(TRUE)
}
