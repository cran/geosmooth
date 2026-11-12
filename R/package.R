#' geosmooth: Geometric Methods for Nonparametric Regression and Density Estimation
#'
#' `geosmooth` estimates a smooth response, or a probability mass, from data
#' given as a matrix of coordinates (one row per observation) or as a weighted
#' graph. It provides local polynomial smoothing, model-averaged local
#' polynomial smoothing, several forms of trend filtering, graph low-pass
#' filtering, Hessian-energy regression, and estimators of probability mass
#' and occupation density.
#'
#' @section Start here:
#' \code{vignette("function-guide", package = "geosmooth")} is organized by
#' task. It lists every exported function, says which fitted models can
#' predict at new points, and gives examples.
#' \code{browseVignettes("geosmooth")} lists all installed guides. The
#' \href{https://pgajer.github.io/geosmooth/articles/function-guide.html}{online version of the guide}
#' links to the full reference.
#'
#' @section Finding a function:
#' The main functions, by task:
#' \describe{
#'   \item{Local regression}{\code{\link{fit.lps}} fits a local polynomial
#'     smoother and selects its settings by cross-validation.
#'     \code{\link{fit.malps}} averages the predictions of overlapping local
#'     polynomial fits.}
#'   \item{Lifting trend filters}{\code{\link{fit.lpl.tf}} penalizes the
#'     absolute residuals of local polynomial predictions.
#'     \code{\link{fit.slpl.tf}} adds a penalty that makes overlapping local
#'     predictions agree.}
#'   \item{One-dimensional trend filters}{\code{\link{fit.gaussian.trend.filtering}}
#'     fits a continuous response and \code{\link{fit.logistic.trend.filtering}}
#'     a binary response. \code{\link{assess.gaussian.trend.filtering}} and
#'     \code{\link{assess.logistic.trend.filtering}} test whether the response
#'     varies with the predictor and assess held-out prediction.}
#'   \item{Graph smoothing}{\code{\link{fit.metric.graph.lowpass}} smooths a
#'     response on the vertices of a graph by filtering its spectrum, and
#'     \code{\link{fit.graph.trend.filtering}} penalizes differences along
#'     the edges.}
#'   \item{Hessian penalties}{\code{\link{fit.ssrhe.hessian.regression}}
#'     penalizes estimated second derivatives with a quadratic penalty, and
#'     \code{\link{fit.ssrhe.hessian.l1.regression}} with an absolute-value
#'     penalty.}
#'   \item{Probability mass}{\code{\link{fit.density}} estimates probability
#'     mass over a fixed set of support points, \code{\link{fit.subject.od}}
#'     does so from the repeated visits of one subject, and
#'     \code{\link{normalize.density}} converts an existing fitted field to
#'     probability mass.}
#'   \item{Synthetic examples}{\code{\link{synthetic.spec}} describes a
#'     synthetic data set, and \code{\link{materialize.synthetic}} draws the
#'     data.}
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
#' @importFrom stats predict
#' @importFrom utils head modifyList tail
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
