## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----circle-------------------------------------------------------------------
library(geosmooth)

circle.spec <- synthetic.spec(
  geometry = synthetic.circle(radius = 1, ambient.dim = 2),
  sampling = synthetic.sampling.grid.interval(
    0, 2 * pi, endpoints = "exclude.lower"),
  truth = synthetic.truth.polynomial(c(b0 = 0)),
  response = synthetic.response.gaussian(sd = 0.05)
)
circle.data <- materialize.synthetic(circle.spec, n = 80, seed = 2026)
circle.data
plot(circle.data)

## ----quadform-----------------------------------------------------------------
surface.spec <- synthetic.spec(
  geometry = synthetic.quadform(
    intrinsic.dim = 2,
    ambient.dim = 5,
    forms = list(diag(c(1, -1))),
    frame = "random.orthonormal"
  ),
  sampling = synthetic.sampling.uniform.box(-1, 1),
  truth = synthetic.truth.polynomial(c(b0 = 0, b1 = 1, b22 = -0.5)),
  response = synthetic.response.gaussian(sd = 0.1)
)
surface.data <- materialize.synthetic(surface.spec, n = 100, seed = 17)
as.data.frame(surface.data)

## ----registry-----------------------------------------------------------------
head(synthetic.registry.ids())
registered <- synthetic.registry.spec("S01.V1")
registered.data <- materialize.synthetic(
  registered, n = 160,
  seed = synthetic.registry.seed("S01.V1", replicate = 1),
  rng.policy = "legacy"
)
synthetic.dataset.checksum(registered.data)

