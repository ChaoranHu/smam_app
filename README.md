# smam_app

A Shiny app that makes the R package [smam](https://github.com/ChaoranHu/smam)
(Statistical Modeling of Animal Movements) easy to use without writing code.

Live version: <https://huchaoran.shinyapps.io/smam_app/>

## What it does

The app has four steps, one tab each, plus a Guide tab with the model
definitions, notation and references.

1. **Data.** Load an example (mountain lion f109 from the package, or a simulated
   track) or upload a CSV file. Pick the time column and the coordinate columns
   and their units; the app converts to hours and kilometres. Optionally define
   segments, either from a column (0 means not used) or as a season window in
   every year from a date column (the same idea as `smam::seasonFilter`).
2. **Fit model.** Fit any of the five models in smam:

   | Model | smam function | Parameters |
   |---|---|---|
   | Moving-Resting (MR) | `fitMR` | $\lambda_M, \lambda_R, \sigma$ |
   | Moving-Resting with Measurement Error (MRME) | `fitMRME`, `fitMRME_naive`, `fitMRMEapprox` | $\lambda_M, \lambda_R, \sigma, \sigma_\epsilon$ |
   | Moving-Resting-Handling (MRH) | `fitMRH` | $\lambda_M, \lambda_R, \lambda_H, \sigma, p$ |
   | Moving-Moving (MM) | `fitMM` | $\lambda_1, \lambda_2, \sigma_1, \sigma_2$ |
   | Brownian Motion with Measurement Error (BMME) | `fitBMME` | $\sigma, \delta$ |

   Here $\lambda_M, \lambda_R, \lambda_H$ (and $\lambda_1, \lambda_2$) are the rates
   (per hour) of the exponential bout lengths in each state, so a bout lasts
   $1/\lambda$ hours on average; $\sigma$ (and $\sigma_1, \sigma_2$) is the Brownian
   motion mobility in $\text{km}/\sqrt{\text{hr}}$; $\sigma_\epsilon$ and $\delta$ are
   the standard deviations of the GPS measurement error in km; $p$ is the
   probability that a moving bout is followed by resting rather than handling.

   Results show each estimate with its meaning and unit, derived quantities
   (mean bout lengths, long run share of time in each state), standard errors
   and 95% intervals (Hessian, parametric bootstrap or Godambe information,
   depending on the model), and the equivalent R code.
3. **States.** Estimate the behavioural state at each fix for MR and MRH
   (`fitStateMR`, `fitViterbiMR`, `fitStateMRH`, `fitViterbiMRH`) and plot the
   track coloured by state.
4. **Simulate.** Simulate tracks with `rMR`, `rMRME`, `rMRH`, `rMM`, `rBMME`,
   download them, or send them to the Data tab to try the estimators.

Fitting MRME and MRH can take minutes. The app runs every fit in a background R
process (via `callr`), shows the optimizer log while it runs, keeps the
interface responsive and lets you cancel.

## Run locally

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "callr", "smam"))
shiny::runApp()
```

## Tests

```sh
Rscript tests/testthat.R                     # R tests; SMAM_SLOW_TESTS=false skips the slow fits
BASE_URL=http://127.0.0.1:8765 node tests/e2e/check.js both   # browser tests, needs playwright
```

The R tests cover data preparation, the model registry, every model, estimator,
standard error method and state method in the background worker, the job manager
and the figures. The browser tests click through every tab on a desktop and a phone
screen and check that the formulas render. GitHub Actions runs both on every pull
request, and runs the browser tests against the live site after each deploy, with
screenshots attached to each run.

## Themes

`R/theme.R` holds the colour themes (`teal`, `journal`, `telemetry`, `park`). Set the
environment variable `SMAM_THEME` to try one; the default is `park` (National park).

The paw print in the app icon (`www/smam-icon.svg`) is the
`paw` icon from [Ionicons](https://github.com/ionic-team/ionicons) (MIT licence).

## Files

* `app.R`: user interface and server logic.
* `R/models.R`: model registry (parameters, defaults, units, derived quantities).
* `R/worker.R`: the code that calls smam inside the background process.
* `R/jobs.R`: background job helpers.
* `R/data.R`: example data and data preparation.
* `R/plots.R`: figures.
* `R/guide.R`: Guide tab text and references.
* `R/theme.R`: colours, fonts and the matching ggplot theme.
* `tests/`: R tests (`testthat`) and browser tests (`e2e`).
