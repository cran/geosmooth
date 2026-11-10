## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----circle-------------------------------------------------------------------
library(geosmooth)

circle.spec <- synthetic.spec(
  geometry = dgraphs::synthetic.circle(radius = 1, ambient.dim = 2),
  sampling = dgraphs::synthetic.sampling.grid.interval(
    0, 2 * pi, endpoints = "exclude.lower"),
  truth = synthetic.truth.polynomial(c(b0 = 0)),
  response = synthetic.response.gaussian(sd = 0.05)
)
circle.data <- materialize.synthetic(circle.spec, n = 80, seed = 2026)
circle.data

## ----circle-views, fig.width=8, fig.height=4, out.width="100%", fig.cap="Left: the circle in its two observed coordinates. Right: responses versus latent angle, with the constant zero truth shown by a dashed line. These views describe the same 80 observations.", fig.alt="A circle of 80 points beside noisy responses scattered around a dashed zero line over angles from zero to two pi."----
old.par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(circle.data, view = "geometry", col = "#0072B2", asp = 1,
     main = "Observed geometry")
plot(circle.data, view = "response", col = "grey40", xlab = "Angle (radians)",
     main = "Response and truth")
par(old.par)

## ----quadform-----------------------------------------------------------------
surface.spec <- synthetic.spec(
  geometry = dgraphs::synthetic.quadform(
    intrinsic.dim = 2,
    ambient.dim = 5,
    forms = list(diag(c(1, -1))),
    frame = "random.orthonormal"
  ),
  sampling = dgraphs::synthetic.sampling.uniform.box(-1, 1),
  truth = synthetic.truth.polynomial(c(b0 = 0, b1 = 1, b22 = -0.5)),
  response = synthetic.response.gaussian(sd = 0.1)
)
surface.data <- materialize.synthetic(surface.spec, n = 100, seed = 17)
dim(surface.data$latent)
dim(surface.data$predictors)
head(as.data.frame(surface.data))

## ----surface-view, fig.width=6, fig.height=4.5, out.width="80%", fig.cap="Two-coordinate projection of a surface embedded in five dimensions. Point colors indicate the known, nonconstant conditional mean; the key uses its actual numerical range.", fig.alt="A projected surface point cloud colored by its known conditional mean, with a numeric color key."----
plot(surface.data, view = "geometry", color = "truth", main = "Projected surface")

## ----surface-fit--------------------------------------------------------------
surface.fit <- fit.lps(surface.data$predictors, surface.data$response,
  foldid = rep(1:5, length.out = surface.data$n),
  support.grid = c(20L, 30L), degree.grid = 1L,
  kernel.grid = "gaussian", coordinate.method = "local.pca", chart.dim = 2L)
surface.fit$selected
head(data.frame(response = surface.data$response, fitted = fitted(surface.fit),
                truth = surface.data$truth))

## ----registry-----------------------------------------------------------------
ids <- c("G1", "G2", "S01.V1")
do.call(rbind, lapply(ids, function(id) {
  spec <- synthetic.registry.spec(id)
  data.frame(id, geometry = spec$geometry.spec$family,
             truth = spec$truth.spec$parameters$estimand)
}))
registered <- synthetic.registry.spec("S01.V1")
registered.data <- materialize.synthetic(
  registered, n = 160,
  seed = synthetic.registry.seed("S01.V1", replicate = 1),
  rng.policy = "legacy"
)
synthetic.dataset.checksum(registered.data)

