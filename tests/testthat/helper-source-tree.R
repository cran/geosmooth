geosmooth.test.source.root <- function(starts = NULL) {
    if (is.null(starts)) starts <- c(
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
                # Installed packages have DESCRIPTION too, but no canonical
                # Rcpp R source and their sources are replaced by lazy-load data.
                source.marker <- file.path(root, "R", "RcppExports.R")
                installed.marker <- file.path(root, "Meta", "package.rds")
                if (identical(package, "geosmooth") &&
                    file.exists(source.marker) && !file.exists(installed.marker)) {
                    # dirname() changes Windows separators while walking upward.
                    return(normalizePath(root, mustWork = TRUE))
                }
            }
            parent <- dirname(root)
            if (identical(parent, root)) break
            root <- parent
        }
    }

    NULL
}
