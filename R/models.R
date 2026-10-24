# Model registry: everything the UI needs to know about each smam model.
# Units assume the prepared data use hours for time and kilometres for
# coordinates (the app converts to these units on the Data tab).

param <- function(id, label, latex, default, unit, meaning,
                  lower = 1e-6, upper = 10) {
  list(id = id, label = label, latex = latex, default = default, unit = unit,
       meaning = meaning, lower = lower, upper = upper)
}

smam_models <- list(

  MR = list(
    id = "MR",
    name = "Moving-Resting (MR)",
    fun = "fitMR",
    summary = "Two states. The animal moves as a Brownian motion while moving and stays still while resting.",
    segment = TRUE,
    params = list(
      param("lamM", "&lambda;<sub>M</sub>", "\\lambda_M", 1, "1/hr",
            "rate of leaving the moving state"),
      param("lamR", "&lambda;<sub>R</sub>", "\\lambda_R", 1, "1/hr",
            "rate of leaving the resting state"),
      param("sigma", "&sigma;", "\\sigma", 1, "km/&radic;hr",
            "mobility (volatility) of the Brownian motion while moving")
    ),
    se = "hessian",
    states = c("moving", "resting")
  ),

  MRME = list(
    id = "MRME",
    name = "Moving-Resting with Measurement Error (MRME)",
    fun = "fitMRME",
    summary = "The MR process observed with Gaussian GPS measurement error. Recommended for real GPS data.",
    segment = TRUE,
    params = list(
      param("lamM", "&lambda;<sub>M</sub>", "\\lambda_M", 4, "1/hr",
            "rate of leaving the moving state"),
      param("lamR", "&lambda;<sub>R</sub>", "\\lambda_R", 0.4, "1/hr",
            "rate of leaving the resting state"),
      param("sigma", "&sigma;", "\\sigma", 1, "km/&radic;hr",
            "mobility (volatility) of the Brownian motion while moving"),
      param("sig_err", "&sigma;<sub>&epsilon;</sub>", "\\sigma_\\epsilon", 0.05, "km",
            "standard deviation of the measurement error")
    ),
    se = "bootstrap",
    states = c("moving", "resting")
  ),

  MRH = list(
    id = "MRH",
    name = "Moving-Resting-Handling (MRH)",
    fun = "fitMRH",
    summary = "Three states. After moving, the animal either rests or handles (for example a kill), both motionless. Computationally heavy: start with a short track.",
    segment = TRUE,
    params = list(
      param("lamM", "&lambda;<sub>M</sub>", "\\lambda_M", 4, "1/hr",
            "rate of leaving the moving state", 0.001, 10),
      param("lamR", "&lambda;<sub>R</sub>", "\\lambda_R", 0.5, "1/hr",
            "rate of leaving the resting state", 0.001, 10),
      param("lamH", "&lambda;<sub>H</sub>", "\\lambda_H", 0.1, "1/hr",
            "rate of leaving the handling state", 0.001, 10),
      param("sigma", "&sigma;", "\\sigma", 5, "km/&radic;hr",
            "mobility (volatility) of the Brownian motion while moving", 0.001, 10),
      param("p", "p", "p", 0.8, "",
            "probability that a moving bout is followed by resting (otherwise handling)",
            0.001, 0.999)
    ),
    se = "hessian",
    states = c("moving", "resting", "handling")
  ),

  MM = list(
    id = "MM",
    name = "Moving-Moving (MM)",
    fun = "fitMM",
    summary = "Two moving states with different mobility, for example travelling and foraging.",
    segment = FALSE,
    params = list(
      param("lam1", "&lambda;<sub>1</sub>", "\\lambda_1", 1, "1/hr",
            "rate of leaving moving state 1 (the faster one)"),
      param("lam2", "&lambda;<sub>2</sub>", "\\lambda_2", 0.1, "1/hr",
            "rate of leaving moving state 2 (the slower one)"),
      param("sigma1", "&sigma;<sub>1</sub>", "\\sigma_1", 1, "km/&radic;hr",
            "mobility in moving state 1 (larger)"),
      param("sigma2", "&sigma;<sub>2</sub>", "\\sigma_2", 0.1, "km/&radic;hr",
            "mobility in moving state 2 (smaller)")
    ),
    se = "bootstrap",
    states = NULL
  ),

  BMME = list(
    id = "BMME",
    name = "Brownian Motion with Measurement Error (BMME)",
    fun = "fitBMME",
    summary = "Constant movement as a Brownian motion, observed with Gaussian measurement error. Fast baseline model.",
    segment = TRUE,
    params = list(
      param("sigma", "&sigma;", "\\sigma", NA, "km/&radic;hr",
            "mobility (volatility) of the Brownian motion"),
      param("delta", "&delta;", "\\delta", NA, "km",
            "standard deviation of the measurement error")
    ),
    se = "hessian",
    states = NULL
  )
)

model_choices <- stats::setNames(names(smam_models),
                                 vapply(smam_models, `[[`, "", "name"))

# Interpretable quantities derived from the estimates.
derived_quantities <- function(model, est) {
  e <- as.list(est)
  rows <- switch(model,
    MR = , MRME = list(
      c("Mean duration of a moving bout", 1 / e$lamM, "hr"),
      c("Mean duration of a resting bout", 1 / e$lamR, "hr"),
      c("Long run share of time moving", e$lamR / (e$lamM + e$lamR), "")
    ),
    MRH = {
      cyc <- 1 / e$lamM + e$p / e$lamR + (1 - e$p) / e$lamH
      list(
        c("Mean duration of a moving bout", 1 / e$lamM, "hr"),
        c("Mean duration of a resting bout", 1 / e$lamR, "hr"),
        c("Mean duration of a handling bout", 1 / e$lamH, "hr"),
        c("Long run share of time moving", (1 / e$lamM) / cyc, ""),
        c("Long run share of time resting", (e$p / e$lamR) / cyc, ""),
        c("Long run share of time handling", ((1 - e$p) / e$lamH) / cyc, "")
      )
    },
    MM = list(
      c("Mean duration in moving state 1", 1 / e$lam1, "hr"),
      c("Mean duration in moving state 2", 1 / e$lam2, "hr"),
      c("Long run share of time in state 1", e$lam2 / (e$lam1 + e$lam2), "")
    ),
    BMME = list(
      c("Typical displacement in 1 hr (per coordinate)", e$sigma, "km"),
      c("Typical displacement in 24 hr (per coordinate)", e$sigma * sqrt(24), "km")
    )
  )
  data.frame(
    Quantity = vapply(rows, `[`, "", 1),
    Value = signif(as.numeric(vapply(rows, `[`, "", 2)), 4),
    Unit = vapply(rows, `[`, "", 3),
    stringsAsFactors = FALSE
  )
}

# How to read the convergence code returned by each fitting routine.
uses_nloptr <- function(model, variant) {
  model == "MRH" || (model == "MRME" && variant %in% c("exact", "naive"))
}

convergence_ok <- function(model, variant, code) {
  if (uses_nloptr(model, variant)) code > 0 else code == 0
}

convergence_note <- function(model, variant, code) {
  ok <- convergence_ok(model, variant, code)
  if (uses_nloptr(model, variant)) {
    paste0("nloptr status ", code, if (ok) " (success)" else " (did not converge, try other starting values)")
  } else {
    paste0("optim code ", code, if (ok) " (converged)" else " (did not converge, try other starting values)")
  }
}

# R code that reproduces a fit outside the app.
fit_code <- function(model, start, options, use_segment) {
  num <- function(x) paste0("c(", paste(signif(x, 6), collapse = ", "), ")")
  seg <- if (use_segment) ', segment = "segment"' else ""
  fun <- smam_models[[model]]$fun
  extra <- ""
  if (model == "MR") extra <- paste0(', likelihood = "', options$likelihood, '"')
  if (model == "MRME") {
    fun <- switch(options$variant, exact = "fitMRME", naive = "fitMRME_naive",
                  approx = "fitMRMEapprox")
    if (options$variant != "approx") {
      extra <- paste0(", lower = ", num(options$lower), ", upper = ", num(options$upper))
    }
  }
  if (model == "MRH") {
    extra <- paste0(", numThreads = ", options$numThreads)
  }
  start_txt <- if (is.null(start)) "NULL" else num(start)
  paste0(
    "library(smam)\n",
    "dat <- read.csv(\"prepared_data.csv\")  # downloaded from the Data tab\n",
    if (!use_segment) "dat$segment <- NULL\n" else "",
    "fit <- ", fun, "(dat, start = ", start_txt, seg, extra, ")\n",
    if (fun %in% c("fitMRME", "fitMR", "fitMRH", "fitMM", "fitBMME")) "estimate(fit)\n" else "fit$estimate\n"
  )
}
