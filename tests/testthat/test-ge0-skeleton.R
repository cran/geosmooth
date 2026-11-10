test_that("GE0 vendored support assets are present", {
    root <- geosmooth.test.source.root()
    skip_if(is.null(root),
            "package source tree is unavailable in installed-package tests")

    expect_true(file.exists(file.path(root, "src", "ANN", "ANN.cpp")))
    expect_true(file.exists(file.path(root, "src", "ANN", "ANN.h")))
    expect_true(file.exists(file.path(root, "inst", "licenses",
                                      "ANN-Copyright-Notice.txt")))
    expect_true(file.exists(file.path(root, "inst", "licenses", "LGPL-2.1.txt")))
    expect_true(file.exists(file.path(root, "inst", "licenses", "MPL-2.0.txt")))
    expect_true(file.exists(file.path(root, "inst", "include", "Eigen", "Core")))
    expect_true(file.exists(file.path(root, "inst", "include", "geosmooth",
                                      "local_pca_charts.h")))
})

test_that("GE0 native scaffold is registered", {
    status <- geosmooth:::.geosmooth.ge0.status()

    expect_identical(status$native.stub, TRUE)
})
