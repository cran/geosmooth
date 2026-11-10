test_that("every exported function and registered S3 method has an Rd example section", {
    root <- geosmooth.test.source.root()
    skip_if(is.null(root),
            "package source tree is unavailable in installed-package tests")

    namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
    exports <- sub(
        "^export[(](.*)[)]$", "\\1",
        grep("^export[(]", namespace, value = TRUE)
    )
    s3.methods <- sub(
        "^S3method[(]([^,]+),([^,)]+).*$", "\\1.\\2",
        grep("^S3method[(]", namespace, value = TRUE)
    )
    documented.objects <- c(exports, s3.methods)

    rd.files <- list.files(
        file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE
    )
    documented <- unlist(lapply(rd.files, function(path) {
        rd <- tools::parse_Rd(path)
        tags <- vapply(rd, attr, "", which = "Rd_tag")
        if (!any(tags == "\\examples")) return(character())
        unlist(lapply(rd[tags == "\\alias"], as.character), use.names = FALSE)
    }), use.names = FALSE)

    missing <- setdiff(documented.objects, documented)
    expect_equal(
        missing, character(),
        info = paste(
            "Exports or S3 registrations without Rd examples:",
            paste(missing, collapse = ", ")
        )
    )
})
