geosmooth.test.source.root <- function() {
    starts <- c(
        tryCatch(
            testthat::test_path("..", ".."),
            error = function(e) NA_character_
        ),
        getwd()
    )
    starts <- unique(starts[!is.na(starts) & nzchar(starts)])

    for (start in starts) {
        root <- normalizePath(start, mustWork = FALSE)
        while (dir.exists(root)) {
            description <- file.path(root, "DESCRIPTION")
            if (file.exists(description)) {
                package <- tryCatch(
                    read.dcf(description, fields = "Package")[[1L]],
                    error = function(e) NA_character_
                )
                if (identical(package, "geosmooth")) return(root)
            }
            parent <- dirname(root)
            if (identical(parent, root)) break
            root <- parent
        }
    }

    NULL
}
