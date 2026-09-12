test_that("grouped package test lanes stop on a failed expectation", {
  runner <- file.path(test_path("..", ".."), "scripts", "run_test_group.R")
  skip_if_not(file.exists(runner),
              "source-tree test runner is not included in the package tarball")
  source <- paste(readLines(runner, warn = FALSE), collapse = "\n")
  expect_match(
    source,
    "testthat::test_file\\([\\s\\S]*stop_on_failure = TRUE",
    perl = TRUE)
})
