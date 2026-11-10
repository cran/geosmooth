test_that("density controls and terminal dots cannot be silently ignored", {
    X <- matrix(1:4, ncol = 1)
    weights <- c(1, 0, 2, 1)
    expect_error(fit.density(X, weights, made.up.control = 42), "unsupported.*made.up.control")
    expect_error(fit.density(X, weights, density.control = list(renomalize = FALSE)), "renomalize")
    expect_error(fit.density(X, weights, density.control = list(renormalize = TRUE, renormalize = FALSE)), "unique")
    expect_error(fit.density(X, weights, density.control = list(FALSE)), "nonempty names")
    for (bad in list(NA, 1, "TRUE", c(TRUE, FALSE))) {
        expect_error(fit.density(X, weights, density.control = list(renormalize = bad)), "TRUE or FALSE")
        expect_error(fit.density(X, weights, return.details = bad), "TRUE or FALSE")
    }
    expect_error(normalize.density(weights, renomalize = FALSE), "renomalize")
    expect_error(fit.density(X, weights, graph.control = list(walk.steps = 1)), "graph.control")
    expect_equal(fit.density(X, weights)$rho, weights / sum(weights))
    expect_equal(fit.density(X, weights, density.control = list(smoothness.auto.1d = FALSE))$rho,
                 weights / sum(weights))
})

test_that("random walk controls accept aliases but reject ambiguity and typos", {
    graph <- list(adj.list = list(2L, c(1L, 3L), 2L),
                  weight.list = list(1, c(1, 1), 1))
    run <- function(ctrl) fit.density(matrix(1:3), c(1, 0, 0),
        method = "graph_random_walk", graph = graph, graph.control = ctrl)
    expect_equal(run(list(walk.step = 1))$rho, c(0, 1, 0))
    expect_equal(run(list(walk_step = 1))$rho, c(0, 1, 0))
    expect_error(run(list(wlak.step = 2)), "wlak.step")
    expect_error(run(list(walk.step = 1, walk_step = 2)), "conflicting")
    expect_error(run(list(walk.step = 1, walk.steps = 2)), "conflicting")
    expect_error(run(list(normalize = "TRUE")), "TRUE or FALSE")
    expect_error(run(list(walk.steps = 2^32)), "integer")
    expect_error(run(list(walk.step.grid = 1:2)), "unsupported")
    # Zero lengths are supported by density affinities, unlike metric low-pass.
    graph$weight.list <- lapply(graph$weight.list, function(w) w * 0)
    expect_equal(run(list(affinity.method = "inverse_length", walk.step = 1))$rho, c(0, 1, 0))
})

test_that("public graph methods reject fractional and nonfinite indices before coercion", {
    adj <- list(2L, c(1L, 3L), 2L)
    weights <- list(1, c(1, 1), 1)
    for (bad in c(1.25, NA, NaN, Inf, 2^32, 0)) {
        malformed <- adj
        malformed[[1]] <- bad
        expect_error(perform.harmonic.smoothing(malformed, weights, 1:3, 1:2), "integer vertex")
        expect_error(harmonic.smoother(malformed, weights, 1:3, 1:2), "integer vertex")
        expect_error(get.region.boundary(malformed, 1:2), "integer vertex")
        expect_error(fit.metric.graph.lowpass(malformed, weights, 1:3), "integer vertex")
        expect_error(fit.density(matrix(1:3), c(1, 0, 0), method = "graph_random_walk",
            graph = list(adj.list = malformed, weight.list = weights)), "integer vertex")
        expect_error(normalize.density(1:3, adj.list = malformed), "integer vertex")
        expect_error(harmonic.smoother(adj, weights, 1:3, c(1, bad)), "integer vertex")
        expect_error(get.region.boundary(adj, c(1, bad)), "integer vertex")
    }
    expect_equal(get.region.boundary(lapply(adj, as.numeric), c(1, 2)), 2L)
    expect_error(harmonic.smoother(adj, weights, 1:3, 1:2, max.iterations = 1.5), "positive integer")
    expect_error(harmonic.smoother(adj, weights, 1:3, 1:2, record.frequency = Inf), "positive integer")
    expect_error(harmonic.smoother(adj, weights, c(1, NA, 2), 1:2), "finite")
})
