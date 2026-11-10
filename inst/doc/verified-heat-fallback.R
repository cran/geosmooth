## -----------------------------------------------------------------------------
library(geosmooth)
adj <- list(2L, c(1L,3L), c(2L,4L), c(3L,5L), 4L)
w <- list(1, c(1,2), c(2,3), c(3,4), 4)
basis <- metric.graph.lowpass.basis(
  adj, w, weight.type="conductance", spectrum="full",
  numerical.policy="verified.fallback")
path <- apply.metric.graph.lowpass.path(
  basis, c(2,NA,4,NA,-1), eta.grid=c(0,.01,.1,1))
path$tuning
path$numerical$diagnostics[,c("eta","n.observed","fallback.used",
                              "candidate.residual","operator.error")]

