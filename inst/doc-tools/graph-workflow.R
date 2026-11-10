## ---- graph-workflow
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
