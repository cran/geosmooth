test_that("the graph tutorial preserves vertex alignment under joint permutation", {
    env <- new.env(parent = globalenv())
    recipe <- system.file("doc-tools/graph-workflow.R", package = "geosmooth")
    capture.output(sys.source(recipe, envir = env))
    p <- c(seq(2, 20, 2), seq(1, 19, 2))
    inverse <- order(p)
    adj <- lapply(env$adjacency[p], function(v) as.integer(inverse[v]))
    lengths <- env$edge.lengths[p]
    changed <- fit.metric.graph.lowpass(adj, lengths, env$response[p],
        n.eigenpairs = 20L, eigen.solver = "dense", eta.grid = c(.001, .01, .1, 1),
        conductance.rule = "inverse.length.power", conductance.alpha = 1,
        conductance.epsilon = 1e-8)
    expect_equal(fitted(changed)[inverse], fitted(env$graph.fit), tolerance = 1e-9)
    expect_equal(residuals(env$second.fit), unname(cos(2*pi*env$locations[, 1])) - fitted(env$second.fit))
    expect_equal(env$graph.fit$gcv$eta.optimal, changed$gcv$eta.optimal)
    expect_error(metric.graph.lowpass.operator(list(2L, integer()), list(1, numeric())), "reciprocal")
    op <- metric.graph.lowpass.operator(list(2L, 1L), list(2, 2),
        conductance.rule = "inverse.length.power", conductance.alpha = 1,
        conductance.epsilon = 1e-8)
    expect_equal(op$edge.table$conductance, 1/(2+1e-8))
    # Public dependency distance route must retain selected vertex order and
    # allow a shortest path to traverse a vertex outside the requested subset.
    distances <- .geosmooth.shortest.path(list(2L, c(1L,3L), 2L),
                                          list(2, c(2,3), 3), c(3L,1L))
    expect_equal(unname(distances), matrix(c(0,5,5,0), 2))
})
