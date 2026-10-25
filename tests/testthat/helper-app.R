`%||%` <- function(x, y) if (is.null(x)) y else x

# Load the app helpers in R/ the same way shiny does, without starting the app.
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(smam)
  library(ggplot2)
})
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

slow_tests <- function() {
  if (identical(tolower(Sys.getenv("SMAM_SLOW_TESTS", "true")), "false")) {
    testthat::skip("slow tests switched off")
  }
}

# A small simulated track in the prepared format: time in hours, then coordinates.
sim_track <- function(model = "MR", n = 40, Tend = 40, seed = 1) {
  set.seed(seed)
  tgrid <- seq(0, Tend, length.out = n)
  d <- switch(model,
    MR = smam::rMR(tgrid, 1, 0.5, 1, "m"),
    MRME = smam::rMRME(tgrid, 1, 0.5, 1, 0.05, "m"),
    MRH = smam::rMRH(tgrid, 4, 0.5, 0.1, 5, 0.8, "m"),
    MM = smam::rMM(tgrid, 1, 0.1, 1, 0.1, "m1"),
    BMME = smam::rBMME(tgrid, sigma = 1, delta = 0.5)
  )
  d <- d[, setdiff(names(d), "state")]
  names(d) <- c("time", paste0("X", seq_len(ncol(d) - 1)))
  d
}

fit_spec <- function(model, data, options = list(), start = NULL, use_segment = FALSE) {
  m <- smam_models[[model]]
  if (is.null(start)) {
    start <- vapply(m$params, `[[`, 0, "default")
    if (anyNA(start)) start <- NULL
  }
  list(model = model, start = start, options = options, use_segment = use_segment, data = data)
}

mrme_bounds <- function() {
  p <- smam_models$MRME$params
  list(lower = vapply(p, `[[`, 0, "lower"), upper = vapply(p, `[[`, 0, "upper"))
}
