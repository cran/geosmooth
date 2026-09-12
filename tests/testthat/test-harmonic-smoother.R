test_that("perform.harmonic.smoothing returns native smoothing diagnostics", {
    adj.list <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
    weight.list <- list(1, c(1, 1), c(1, 1), 1)
    values <- c(0, 10, -10, 1)

    out <- suppressWarnings(perform.harmonic.smoothing(
        adj.list = adj.list,
        weight.list = weight.list,
        values = values,
        region.vertices = 1:3,
        max.iterations = 10,
        tolerance = 1e-10
    ))

    expect_named(out, c(
        "harmonic_predictions", "converged", "num_region", "num_boundary",
        "num_interior", "num_iterations", "max_change", "max_residual"
    ))
    expect_length(out$harmonic_predictions, length(values))
    expect_equal(out$num_region, 3L)
    expect_equal(out$num_iterations, 10L)
})

test_that("harmonic.smoother uses the geosmooth native backend", {
    adj.list <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
    weight.list <- list(1, c(1, 1), c(1, 1), 1)
    values <- c(0, 10, -10, 1)

    out <- harmonic.smoother(
        adj.list = adj.list,
        weight.list = weight.list,
        values = values,
        region.vertices = 1:4,
        max.iterations = 12,
        tolerance = 1e-10,
        record.frequency = 2,
        stability.window = 2,
        stability.threshold = 0.1
    )

    expect_s3_class(out, "harmonic_smoother")
    expect_named(out, c(
        "harmonic_predictions", "i_harmonic_predictions", "i_basins",
        "stable_iteration", "topology_differences", "basin_cx_differences",
        "converged", "num_region", "num_boundary", "num_interior"
    ))
    expect_length(out$harmonic_predictions, length(values))
    expect_equal(nrow(out$i_harmonic_predictions), length(values))
    expect_equal(length(out$i_basins), ncol(out$i_harmonic_predictions))
    expect_equal(out$topology_differences, out$basin_cx_differences)
    expect_true(out$stable_iteration >= 1L)

    smry <- summary(out)
    expect_s3_class(smry, "summary.harmonic_smoother")
    expect_equal(smry$stable_iteration, out$stable_iteration)
})
