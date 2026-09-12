#' geosmooth package
#'
#' Geometric smoothing and conditional expectation methods for coordinate data,
#' point clouds, and weighted graphs. The package provides local polynomial,
#' model-averaged, trend-filtering, graph low-pass, occupation-density, and
#' Hessian-energy methods.
#'
#' @useDynLib geosmooth, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils head modifyList
"_PACKAGE"

.geosmooth.ge0.status <- function() {
    list(
        package = "geosmooth",
        phase = "GE0",
        native.stub = rcpp_geosmooth_native_stub(),
        vendored.ann = file.exists(system.file("licenses", "ANN-Copyright-Notice.txt",
                                               package = "geosmooth")),
        vendored.eigen = file.exists(system.file("include", "Eigen", "Core",
                                                 package = "geosmooth"))
    )
}
