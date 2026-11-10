#!/usr/bin/env Rscript
# Run from the package source directory with geosmooth and its dependencies installed.
input <- normalizePath("vignettes/pca-pls-chart-comparison.Rmd", mustWork=TRUE)
output <- file.path("inst", "extdata", "vignette-results", "pca-pls-comparison.rds")
previous <- options(geosmooth.recompute.vignettes=TRUE)
work <- tempfile("pca-pls-render-")
dir.create(work)
environment <- new.env(parent=globalenv())
tryCatch({
    rmarkdown::render(input, output_dir=work, quiet=TRUE, envir=environment)
    provenance <- list(
        vignette.sha256=digest::digest(file=input, algo="sha256"),
        package.versions=vapply(c("geosmooth", "dgraphs", "Rcpp", "digest"),
            function(package) as.character(utils::packageVersion(package)), character(1)),
        R.version=R.version.string,
        platform=R.version$platform,
        seeds=as.vector(outer(271000+100*(1:5), 1:5, "+")),
        repetitions=5L, training=180L, tuning=60L, test=120L)
    payload <- list(results=environment$results,
                    illustrations=environment$illustrations,
                    provenance=provenance)
    stopifnot(nrow(payload$results)==200L, length(payload$illustrations)==5L)
    dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
    saveRDS(payload, output, version=3, compress="xz")
    message("Saved ", output)
}, finally={
    options(previous)
    unlink(work, recursive=TRUE)
})
