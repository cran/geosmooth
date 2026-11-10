test_that("native symbols use the geosmooth package prefix", {
    root <- geosmooth.test.source.root()
    skip_if(is.null(root),
            "package source tree is unavailable in installed-package tests")
    rcpp <- readLines(file.path(root, "src", "RcppExports.cpp"), warn = FALSE)
    namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
    omp.compat <- readLines(file.path(root, "src", "omp_compat.h"), warn = FALSE)
    metric.lowpass <- readLines(
        file.path(root, "src", "metric_graph_lowpass.cpp"),
        warn = FALSE
    )

    expect_true(any(grepl("_geosmooth_", rcpp)))
    expect_false(any(grepl("_gflow_", rcpp)))
    expect_false(any(grepl("_gflow_", namespace)))
    expect_true(any(grepl("geosmooth_get_max_threads", omp.compat)))
    expect_true(any(grepl("geosmooth_get_max_threads", metric.lowpass)))
    expect_false(any(grepl("gflow_", omp.compat)))
    expect_false(any(grepl("gflow_", metric.lowpass)))
})

test_that("dgraphs exports the shortest-path contract required by geosmooth", {
        exports <- getNamespaceExports("dgraphs")
    expected <- if ("dgraph" %in% exports) {
        c("dgraph", "graph.geodesic.distances", "graph.adjacency", "graph.lengths")
    } else "shortest.path"
    expect_true(all(expected %in% exports))
})
