# Convert a full-data summed-loss penalty to a training-fold penalty.
.tf.cv.penalty.factor <- function(training.mass, full.mass,
                                  convention = c("training_mass", "none")) {
    convention <- match.arg(convention)
    if (length(training.mass) != 1L || length(full.mass) != 1L ||
        !is.finite(training.mass) || !is.finite(full.mass) ||
        training.mass <= 0 || full.mass <= 0 || training.mass > full.mass)
        stop("CV requires positive training weight mass no greater than full weight mass.",
             call. = FALSE)
    if (convention == "training_mass") training.mass / full.mass else 1
}
