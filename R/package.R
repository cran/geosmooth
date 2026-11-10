#' geosmooth package
#'
#' Geometric smoothing and conditional expectation methods for coordinate data,
#' point clouds, and weighted graphs. The package provides local polynomial,
#' model-averaged, trend-filtering, graph low-pass, occupation-density, and
#' Hessian-energy methods.
#'
#' @section Start here:
#' Use \code{browseVignettes("geosmooth")} to open the installed HTML guides,
#' including \code{vignette("function-guide", package = "geosmooth")}.
#' The \href{https://pgajer.github.io/geosmooth/articles/function-guide.html}{online guide}
#' links to the full reference. The task links below work in installed help
#' without an internet connection.
#'
#' @section Finding a function:
#' Read \code{vignette("function-guide", package = "geosmooth")} for a
#' task-based guide to all exported functions, prediction limits, and examples.
#' \describe{
#'   \item{Local regression}{Start with \code{\link{fit.lps}} or
#'     \code{\link{fit.malps}}.}
#'   \item{Lifting trend filters}{Use \code{\link{fit.lpl.tf}} or
#'     \code{\link{fit.slpl.tf}} for residual penalties with optional
#'     prediction synchronization.}
#'   \item{Univariate trend filters}{Use \code{\link{fit.gaussian.trend.filtering}}
#'     for continuous responses or \code{\link{fit.logistic.trend.filtering}}
#'     for binary responses. Assess association and prediction with
#'     \code{\link{assess.gaussian.trend.filtering}} or
#'     \code{\link{assess.logistic.trend.filtering}}, respectively.}
#'   \item{Graph smoothing}{Use \code{\link{fit.metric.graph.lowpass}} for
#'     spectral smoothing or \code{\link{fit.graph.trend.filtering}} for
#'     graph-difference penalties.}
#'   \item{Hessian penalties}{Use \code{\link{fit.ssrhe.hessian.regression}}
#'     for quadratic energy or \code{\link{fit.ssrhe.hessian.l1.regression}}
#'     for an absolute-value penalty.}
#'   \item{Probability mass}{Use \code{\link{fit.density}} for support weights,
#'     \code{\link{fit.subject.od}} for repeated visits, and
#'     \code{\link{normalize.density}} to convert an existing fitted field.}
#'   \item{Synthetic examples}{Combine components with
#'     \code{\link{synthetic.spec}} and draw data with
#'     \code{\link{materialize.synthetic}}.}
#' }
#'
#' @seealso \code{vignette("function-guide", package = "geosmooth")},
#'   \code{vignette("synthetic-datasets", package = "geosmooth")}
#'
#' @examples
#' X <- matrix(seq(0, 1, length.out = 30), ncol = 1)
#' y <- sin(2 * pi * X[, 1])
#' fit <- fit.lps(X, y, support.grid = 8L, degree.grid = 1L,
#'                kernel.grid = "tricube", foldid = rep(1:3, 10))
#' head(predict(fit))
#'
#' @useDynLib geosmooth, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils head modifyList
"_PACKAGE"

.geosmooth.ge0.status <- function() {
    list(
        package = "geosmooth",
        native.support = "registered",
        native.stub = rcpp_geosmooth_native_stub(),
        vendored.ann = file.exists(system.file("licenses", "ANN-Copyright-Notice.txt",
                                               package = "geosmooth")),
        vendored.eigen = file.exists(system.file("include", "Eigen", "Core",
                                                 package = "geosmooth"))
    )
}
