# Run from the repository root with:  Rscript tests/testthat.R
# Set SMAM_SLOW_TESTS=false to skip the slow fits (MRME exact, MRH, bootstraps).
library(testthat)
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
