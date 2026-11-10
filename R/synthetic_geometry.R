# Retained statistical recipe and public geometry delegation.

.synthetic.g4.segment.rectangle <- function(eta = 0.02) {
  eta <- .synthetic.scalar.double(eta, "eta", 0)
  .new.synthetic.component(
    "geometry", "g4.segment.rectangle",
    list(
      ambient.dim = 3L,
      intrinsic.dim = NA_integer_,
      intrinsic.dim.by.region = c(A = 1L, B = 2L),
      codimension.by.region = c(A = 2L, B = 1L),
      eta = eta,
      declared.regions = c("A", "B")
    ),
    subclass = "synthetic_stratified_geometry")
}


.embed.synthetic.geometry <- function(geometry, latent, frame.matrix = NULL) {
  dgraphs::embed.synthetic.geometry(geometry, latent, frame.matrix)
}
.validate.synthetic.latent <- function(geometry, latent) {
  dgraphs::validate.synthetic.geometry(geometry, latent)
}
