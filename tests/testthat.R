library(testthat)
library(geosmooth)

validation.contexts <- gsub(
  "^test-|\\.R$", "",
  c(
    "test-lps-bandwidth-multiplier.R",
    "test-lps-binary-metric-consistency.R",
    "test-lps-binary-separation.R",
    "test-lps-binomial-na-consistency.R",
    "test-lps-nested-grouped-cv.R",
    "test-lps-ridge-alignment.R",
    "test-lps-tier0-correctness.R",
    "test-lps-tier0-correctness-extended.R",
    "test-lps-tier4-uncertainty.R",
    "test-odcv5-all-method-smoke-report.R"
  )
)

# These acceptance and scientific-validation studies remain available via
# `make test-validation` and `make test-all`, but are too CPU-intensive for
# package checks on shared local and external check farms.
test_check(
  "geosmooth",
  filter = paste(validation.contexts, collapse = "|"),
  invert = TRUE
)
