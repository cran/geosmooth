test_that("OD-CV5 all-method smoke runner produces uniform report artifacts", {
    out.dir <- file.path(tempdir(), paste0("odcv5-", Sys.getpid()))
    pkg.root <- normalizePath(file.path(getwd(), "..", ".."),
                              winslash = "/", mustWork = FALSE)
    if (!file.exists(file.path(pkg.root, "DESCRIPTION"))) {
        pkg.root <- normalizePath(".", winslash = "/", mustWork = TRUE)
    }
    script <- file.path(pkg.root, "scripts",
                        "run_od_cv5_all_method_smoke.R")
    skip_if_not(file.exists(script),
                "source-tree report runner is not included in the package tarball")
    res <- system2(
        file.path(R.home("bin"), "Rscript"),
        c(script, paste0("--out-dir=", out.dir)),
        stdout = TRUE,
        stderr = TRUE
    )
    expect_true(file.exists(file.path(out.dir,
                                      "od_cv5_all_method_smoke_report.html")),
                info = paste(res, collapse = "\n"))
    summary.path <- file.path(out.dir, "tables", "od_cv5_method_summary.csv")
    candidate.path <- file.path(out.dir, "tables", "od_cv5_candidate_table.csv")
    expect_true(file.exists(summary.path))
    expect_true(file.exists(candidate.path))
    summary <- utils::read.csv(summary.path, stringsAsFactors = FALSE)
    expect_setequal(
        summary$method,
        c("graph_random_walk", "chart_kernel", "local_likelihood_density",
          "local_likelihood_bernoulli", "lps_count",
          "lps_logistic_binary", "ps_lps_count")
    )
    expect_true(all(summary$status == "ok"))
    expect_true(all(is.finite(summary$visit.cv.neg.log.rho)))
    expect_true(all(summary$n.candidates > 0))
    html <- paste(readLines(file.path(out.dir,
                                      "od_cv5_all_method_smoke_report.html"),
                            warn = FALSE),
                  collapse = "\n")
    expect_true(grepl("Figure 1\\.", html))
    expect_true(grepl("contract smoke figure", html))
})
