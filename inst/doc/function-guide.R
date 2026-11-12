## ----standalone-width, echo=FALSE, results='asis'-----------------------------
# Let pkgdown own the site layout; widen only the standalone installed vignette.
if (!identical(Sys.getenv("IN_PKGDOWN"), "true")) {
  cat("<style>body { max-width: 1120px; }</style>")
}

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
read.recipe <- getExportedValue("knitr", "read_chunk")
read.recipe(system.file("doc-tools/first-fit.R", package = "geosmooth"))
read.recipe(system.file("doc-tools/graph-workflow.R", package = "geosmooth"))

## ----first-fit-model----------------------------------------------------------
library(geosmooth)
set.seed(1)
X <- matrix(seq(0, 1, length.out = 60), ncol = 1)
truth <- sin(2 * pi * X[, 1])
y <- truth + rnorm(nrow(X), sd = 0.08)
fit <- fit.lps(
  X, y, foldid = rep(1:5, length.out = nrow(X)),
  support.grid = c(8L, 12L, 16L), degree.grid = 0:1,
  kernel.grid = "gaussian"
)
fit$selected

## ----first-fit-plot, echo=FALSE, fig.width=7.2, fig.height=4, out.width="100%", fig.cap="The points are noisy observations. The solid blue line is the selected fit; the dashed orange line is known only because these data are synthetic.", fig.alt="Noisy observations along a sine curve, with a fitted solid curve close to the dashed synthetic truth."----
plot(X[, 1], y, pch = 16, col = "#666666", xlab = "Coordinate",
     ylab = "Response", main = "A local polynomial fit to noisy observations")
lines(X[, 1], truth, col = "#D55E00", lty = 2, lwd = 2)
lines(X[, 1], fitted(fit), col = "#0072B2", lty = 1, lwd = 2)
legend("topright", c("Observations", "Fitted response", "Known synthetic truth"),
       pch = c(16, NA, NA), lty = c(NA, 1, 2),
       col = c("#666666", "#0072B2", "#D55E00"), bty = "n", cex = .8)

## ----first-predictions--------------------------------------------------------
head(fitted(fit))
new.X <- matrix(seq(0.1, 0.9, length.out = 5), ncol = 1)
predict(fit, newdata = new.X)

## ----shared-generics----------------------------------------------------------
malps.fit <- fit.malps(
  X, y, degree = 1L, support.type = "knn", support.size = 12L
)
new.y <- cos(2 * pi * X[, 1])
malps.refit <- refit(malps.fit, y = new.y)
S <- smoother.matrix(malps.fit)
max(abs(as.vector(S %*% new.y) - predict(malps.refit)))

## ----neighbor-probability, eval=FALSE-----------------------------------------
# fit <- fit.neighbor.probability(X, y, method = "case_corrected", k = 3,
#                                 id = record_id, group = participant_id,
#                                 normalization = "max", metric = "linf.simplex")
# predict(fit, newdata = query_X, query.group = query_participant_id,
#         type = "diagnostics")

## ----graph-workflow-----------------------------------------------------------
# Each row, vertex, and response refers to the same location.
locations <- matrix(seq(0, 1, length.out = 20), ncol = 1,
                    dimnames = list(sprintf("site-%02d", 1:20), "position"))
response <- sin(2 * pi * locations[, 1])
graph <- dgraphs::create.rknn.graph(locations, radius = .06,
                                    graph.detail = "minimal")
# dgraphs' new dgraph representation has public accessors. The second branch
# supports the pinned transitional dependency used by this development release.
if (inherits(graph, "dgraph")) {
  adjacency <- dgraphs::graph.adjacency(graph)
  edge.lengths <- dgraphs::graph.lengths(graph)
} else {
  adjacency <- graph$adj_list
  edge.lengths <- graph$weight_list
}
graph.fit <- fit.metric.graph.lowpass(
  adjacency, edge.lengths, response,
  n.eigenpairs = nrow(locations), eigen.solver = "dense",
  conductance.rule = "inverse.length.power", conductance.alpha = 1,
  conductance.epsilon = 1e-8, eta.grid = c(.001, .01, .1, 1)
)
head(data.frame(site = rownames(locations), response,
                fitted = fitted(graph.fit)))
graph.fit$gcv[c("eta.optimal", "gcv.optimal", "effective.df")]
# Reuse the fitted object's spectral basis, with the selected smoothing fixed.
second.fit <- refit(graph.fit, y = cos(2 * pi * locations[, 1]))
head(fitted(second.fit))

## ----gaussian-trend-filtering-assessment--------------------------------------
x.continuous <- seq(0, 1, length.out = 40)
y.continuous <- sin(2*pi*x.continuous) + 0.2*cos(seq_along(x.continuous))
gaussian.fit <- fit.gaussian.trend.filtering(
  x.continuous, y.continuous, lambda = c(Inf, 0.01), backend = "banded")
# Nine permutations keep this API example short; use more for inference.
gaussian.assessment <- assess.gaussian.trend.filtering(
  gaussian.fit, permutations = 9, inner.folds = 3, outer.folds = 3, seed = 42,
  additional.tests = c("absolute.deviation", "total.variation"))
gaussian.assessment$tests
gaussian.assessment$prediction$scores

## ----hessian-selection--------------------------------------------------------
hessian.fit <- fit.ssrhe.hessian.regression(
  X, y, k = 12L, tangent.dim = 1L,
  lambda.selection = "cv", lambda1.grid = c(0.01, 0.1),
  cv.control = list(foldid = rep(1:5, length.out = nrow(X)))
)
hessian.fit$lambda

## ----occupation---------------------------------------------------------------
support <- matrix(seq(0, 1, length.out = 6), ncol = 1)
occupation <- fit.subject.od(
  support, subject.index = c(1L, 3L, 3L, 6L), method = "empirical"
)
occupation$rho
sum(occupation$rho)

## ----help-navigation, eval=FALSE----------------------------------------------
# help(package = "geosmooth")
# help("geosmooth-package", package = "geosmooth")
# help("fit.lps", package = "geosmooth")
# help("predict.lps", package = "geosmooth")
# vignette("function-guide", package = "geosmooth")
# vignette("synthetic-datasets", package = "geosmooth")

