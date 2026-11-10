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
        "harmonic_predictions", "converged", "status", "num_region", "num_boundary",
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
        "stable_iteration", "stability_detected", "num_iterations", "recorded_iterations",
        "max_change", "max_residual", "status", "topology_differences", "basin_cx_differences",
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

test_that("both harmonic solvers recover an independently calculated weighted solution", {
    adj <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
    weights <- list(1, c(1, 2), c(2, 4), c(4, 1), 1)
    values <- c(99, 1, -10, 4, -99)
    expected <- (1 / (2 + 1e-10) + 4 / (4 + 1e-10)) /
        (1 / (2 + 1e-10) + 1 / (4 + 1e-10))
    for (solver in list(perform.harmonic.smoothing, harmonic.smoother)) {
        fit <- solver(adj, weights, values, 2:4, tolerance = 1e-12)
        expect_true(fit$converged)
        expect_identical(fit$status, "converged")
        expect_equal(fit$harmonic_predictions[c(1, 2, 4, 5)], values[c(1, 2, 4, 5)])
        expect_equal(fit$harmonic_predictions[3], expected, tolerance = 1e-12)
        expect_lt(tail(fit$max_residual, 1), 1e-12)
    }
    expect_equal(get.region.boundary(adj, 2:4), c(2L, 4L))
    expect_equal(get.region.boundary(adj, 1:5), integer())
})

test_that("extrema stability is not numerical convergence on an alternating path", {
    adj <- list(2L, c(1L, 3L), c(2L, 4L), 3L)
    weights <- list(1, c(1, 1), c(1, 1), 1)
    values <- c(0, 2, -1, 1)
    basic <- perform.harmonic.smoothing(adj, weights, values, 1:4)
    expect_identical(basic$status, "no_boundary")
    expect_false(basic$converged)
    expect_identical(basic$num_iterations, 0L)
    expect_equal(basic$harmonic_predictions, values)
    every <- harmonic.smoother(adj, weights, values, 1:4,
                               max.iterations = 20, record.frequency = 1)
    alternate <- harmonic.smoother(adj, weights, values, 1:4,
                                   max.iterations = 20, record.frequency = 2)
    expect_equal(every$harmonic_predictions, alternate$harmonic_predictions)
    expect_false(every$converged)
    expect_false(alternate$converged)
    expect_false(every$stability_detected)
    expect_true(alternate$stability_detected)
    expect_gt(tail(alternate$max_residual, 1), 1)
    expect_match(paste(capture.output(print(alternate)), collapse = " "), "Numerical convergence: no")
    expect_match(paste(capture.output(summary(every)), collapse = " "), "not detected")
})

test_that("undetected stability is NA and final records use actual iteration numbers", {
    adj <- list(2L, c(1L, 3L), 2L)
    weights <- list(1, c(1, 1), 1)
    fit <- harmonic.smoother(adj, weights, c(0, 1, 0), 1:3,
                             max.iterations = 5, record.frequency = 2,
                             stability.window = 10)
    expect_false(fit$stability_detected)
    expect_identical(fit$stable_iteration, NA_integer_)
    expect_equal(fit$recorded_iterations, c(0L, 2L, 4L, 5L))
    expect_equal(fit$i_harmonic_predictions[, 4], fit$harmonic_predictions)
    expect_identical(fit$num_iterations, 5L)
    expect_identical(summary(fit)$converged, fit$converged)
    # Numerical convergence before a full extrema window does not invent detection.
    flat <- harmonic.smoother(adj, weights, rep(1, 3), 1:3)
    expect_true(flat$converged)
    expect_false(flat$stability_detected)
    expect_identical(flat$stable_iteration, NA_integer_)
})

test_that("no-interior regions and disconnected unanchored components have truthful status", {
    adj <- list(2L, c(1L, 3L), 2L, 5L, 4L)
    weights <- list(1, c(1, 1), 1, 1, 1)
    values <- c(1, 0, 2, 0, 1)
    for (solver in list(perform.harmonic.smoothing, harmonic.smoother)) {
        boundary.only <- solver(adj, weights, values, 2L)
        expect_true(boundary.only$converged)
        expect_identical(boundary.only$status, "no_interior")
        expect_identical(boundary.only$num_iterations, 0L)
        mixed <- solver(adj, weights, values, c(1L, 2L, 4L, 5L), max.iterations = 10)
        expect_false(mixed$converged)
        expect_identical(mixed$status, "iteration_limit")
        expect_equal(mixed$harmonic_predictions[2:3], values[2:3])
        expect_gt(tail(mixed$max_residual, 1), 0.5)
    }
})
