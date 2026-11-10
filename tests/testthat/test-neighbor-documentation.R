test_that("the common interface and prediction method have executable help examples", {
    root <- geosmooth.test.source.root()
    skip_if(is.null(root), "No source tree in installed checks.")
    for (name in c("fit.neighbor.probability", "predict.neighbor.probability.fit")) {
        rd <- tools::parse_Rd(file.path(root, "man", paste0(name, ".Rd")))
        expect_true(any(vapply(rd, attr, "", which = "Rd_tag") == "\\examples"))
    }
    expect_true("fit.neighbor.probability" %in% getNamespaceExports("geosmooth"))
    expect_true(is.function(getS3method("predict", "neighbor.probability.fit")))
})
