test_that("the model registry is consistent", {
  for (id in names(smam_models)) {
    m <- smam_models[[id]]
    expect_equal(m$id, id)
    ids <- vapply(m$params, `[[`, "", "id")
    expect_false(anyDuplicated(ids) > 0)
    for (p in m$params) {
      expect_true(all(c("label", "latex", "unit", "meaning", "lower", "upper") %in% names(p)))
      expect_lt(p$lower, p$upper)
      if (!is.na(p$default)) expect_true(p$default >= p$lower && p$default <= p$upper, info = paste(id, p$id))
    }
    expect_true(id %in% names(model_math))
  }
  expect_setequal(unname(model_choices), names(smam_models))
})

test_that("derived quantities are sensible", {
  mr <- derived_quantities("MR", c(lamM = 2, lamR = 0.5, sigma = 1))
  expect_equal(mr$Value, c(0.5, 2, 0.2))
  mrh <- derived_quantities("MRH", c(lamM = 4, lamR = 0.5, lamH = 0.1, sigma = 5, p = 0.8))
  shares <- mrh$Value[grepl("share", mrh$Quantity)]
  expect_equal(sum(shares), 1, tolerance = 1e-3)
  mm <- derived_quantities("MM", c(lam1 = 1, lam2 = 0.1, sigma1 = 1, sigma2 = 0.1))
  expect_equal(nrow(mm), 3)
  bm <- derived_quantities("BMME", c(sigma = 2, delta = 0.1))
  expect_equal(bm$Value[2], signif(2 * sqrt(24), 4))
})

test_that("convergence codes are read the right way for each optimiser", {
  expect_true(convergence_ok("MR", "", 0))
  expect_false(convergence_ok("MR", "", 1))
  expect_true(convergence_ok("MRME", "exact", 4))
  expect_false(convergence_ok("MRME", "exact", -1))
  expect_true(convergence_ok("MRME", "approx", 0))
  expect_true(convergence_ok("MRH", "", 3))
  expect_match(convergence_note("MRME", "naive", 1), "success")
  expect_match(convergence_note("BMME", "", 52), "did not converge")
})

test_that("the generated R code parses for every model and option", {
  b <- mrme_bounds()
  cases <- list(
    list("MR", list(likelihood = "full")), list("MR", list(likelihood = "composite")),
    list("MRME", c(list(variant = "exact"), b)), list("MRME", c(list(variant = "naive"), b)),
    list("MRME", list(variant = "approx")), list("MRH", list(numThreads = 2)),
    list("MM", list()), list("BMME", list())
  )
  for (cs in cases) {
    for (seg in c(TRUE, FALSE)) {
      code <- fit_code(cs[[1]], c(1, 2, 3), cs[[2]], seg)
      expect_silent(parse(text = code))
      expect_match(code, smam_models[[cs[[1]]]]$fun %||% "", fixed = TRUE)
    }
  }
  expect_match(fit_code("BMME", NULL, list(), FALSE), "start = NULL")
})

test_that("the guide and the theme build", {
  expect_true(grepl("References", as.character(guide_panel())))
  for (id in names(smam_themes)) {
    withr_env <- Sys.getenv("SMAM_THEME")
    Sys.setenv(SMAM_THEME = id)
    css <- app_css()
    expect_false(grepl("\\{[a-z]+\\}", css), info = id)
    expect_s3_class(app_theme(), "bs_theme")
    expect_true(all(c("moving", "resting", "handling") %in% names(state_colours())), info = id)
    expect_s3_class(theme_smam(), "theme")
    Sys.setenv(SMAM_THEME = withr_env)
  }
})
