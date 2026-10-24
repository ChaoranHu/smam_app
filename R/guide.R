# Text for the Guide tab and the short model descriptions on the Fit tab.
# Math is written in LaTeX and rendered by MathJax.

model_math <- list(
  MR = shiny::tagList(
    shiny::p("The state \\(S(t) \\in \\{M, R\\}\\) alternates between moving \\(M\\) and resting \\(R\\).",
             "Each moving bout lasts an exponential time with rate \\(\\lambda_M\\) and each resting bout",
             "an exponential time with rate \\(\\lambda_R\\). The location is"),
    shiny::p("$$X(t) = X(0) + \\sigma\\, B\\big(M(t)\\big),$$"),
    shiny::p("$$M(t) = \\int_0^t \\mathbf{1}\\{S(u) = M\\}\\, du,$$"),
    shiny::p("where \\(B\\) is a standard Brownian motion, \\(M(t)\\) is the total time spent moving up to time \\(t\\),",
             "\\(\\mathbf{1}\\{\\cdot\\}\\) is the indicator function and \\(\\sigma\\) is the mobility parameter.",
             "The animal stays put while resting.")
  ),
  MRME = shiny::tagList(
    shiny::p("The true location \\(X(t)\\) follows the Moving-Resting (MR) process: the animal alternates between",
             "moving bouts (exponential with rate \\(\\lambda_M\\)), during which it moves as a Brownian motion with",
             "mobility \\(\\sigma\\), and resting bouts (exponential with rate \\(\\lambda_R\\)), during which it stays put.",
             "The GPS records the location with error:"),
    shiny::p("$$Y(t_i) = X(t_i) + \\epsilon_i,$$"),
    shiny::p("where the errors \\(\\epsilon_i \\sim N(0, \\sigma_\\epsilon^2 I)\\) are independent,",
             "\\(t_i\\) is the \\(i\\)th observation time, \\(Y(t_i)\\) the recorded location, \\(I\\) the identity matrix",
             "and \\(\\sigma_\\epsilon\\) the standard deviation of the measurement error.",
             "Estimation maximises a composite likelihood (Hu et al., 2021).")
  ),
  MRH = shiny::tagList(
    shiny::p("Three states: moving \\(M\\), resting \\(R\\) and handling \\(H\\). Bout lengths are exponential with rates",
             "\\(\\lambda_M, \\lambda_R, \\lambda_H\\). After a moving bout the animal rests with probability \\(p\\)",
             "or handles with probability \\(1 - p\\), and then moves again. It moves as \\(\\sigma B\\) while moving",
             "and stays put in \\(R\\) and \\(H\\).")
  ),
  MM = shiny::tagList(
    shiny::p("Two moving states \\(M_1\\) and \\(M_2\\) with bout rates \\(\\lambda_1, \\lambda_2\\). The location moves as",
             "a Brownian motion with mobility \\(\\sigma_1\\) in \\(M_1\\) and \\(\\sigma_2\\) in \\(M_2\\), where",
             "\\(\\sigma_1 > \\sigma_2\\):"),
    shiny::p("$$X(t) = X(0) + \\sigma_1 B_1\\big(M_1(t)\\big)$$"),
    shiny::p("$$\\qquad\\qquad + \\sigma_2 B_2\\big(M_2(t)\\big),$$"),
    shiny::p("where \\(B_1, B_2\\) are independent standard Brownian motions and \\(M_k(t)\\) is the time spent in state \\(M_k\\) up to \\(t\\).")
  ),
  BMME = shiny::tagList(
    shiny::p("The animal always moves:"),
    shiny::p("$$Y(t_i) = \\sigma B(t_i) + \\epsilon_i,$$"),
    shiny::p("where the errors \\(\\epsilon_i \\sim N(0, \\delta^2)\\) are independent, \\(\\sigma\\) is the mobility of the Brownian motion \\(B\\) and \\(\\delta\\) is the measurement error",
             "standard deviation. Leave the starting values empty to use a moment estimator.")
  )
)

references <- shiny::tags$ol(
  shiny::tags$li("Yan, J., Chen, Y., Lawrence-Apfel, K., Ortega, I. M., Pozdnyakov, V., Williams, S., and Meyer, T. (2014).",
                 "A moving-resting process with an embedded Brownian motion for animal movements.",
                 shiny::em("Population Ecology"), "56(2): 401 to 415.",
                 shiny::a("doi:10.1007/s10144-013-0428-8", href = "https://doi.org/10.1007/s10144-013-0428-8", target = "_blank")),
  shiny::tags$li("Pozdnyakov, V., Elbroch, L. M., Labarga, A., Meyer, T., and Yan, J. (2017).",
                 "Discretely observed Brownian motion governed by telegraph process: estimation.",
                 shiny::em("Methodology and Computing in Applied Probability."),
                 shiny::a("doi:10.1007/s11009-017-9547-6", href = "https://doi.org/10.1007/s11009-017-9547-6", target = "_blank")),
  shiny::tags$li("Pozdnyakov, V., Meyer, T., Wang, Y., and Yan, J. (2014).",
                 "On modeling animal movements using Brownian motion with measurement error.",
                 shiny::em("Ecology"), "95(2): 247 to 253.",
                 shiny::a("doi:10.1890/13-0532.1", href = "https://doi.org/10.1890/13-0532.1", target = "_blank")),
  shiny::tags$li("Pozdnyakov, V., Elbroch, L. M., Hu, C., Meyer, T., and Yan, J. (2020).",
                 "On estimation for Brownian motion governed by telegraph process with multiple off states.",
                 shiny::em("Methodology and Computing in Applied Probability"), "22: 1275 to 1291.",
                 shiny::a("doi:10.1007/s11009-020-09774-1", href = "https://doi.org/10.1007/s11009-020-09774-1", target = "_blank")),
  shiny::tags$li("Hu, C., Elbroch, L. M., Meyer, T., Pozdnyakov, V., and Yan, J. (2021).",
                 "Moving-resting process with measurement error in animal movement modeling.",
                 shiny::em("Methods in Ecology and Evolution."),
                 shiny::a("doi:10.1111/2041-210X.13694", href = "https://doi.org/10.1111/2041-210X.13694", target = "_blank"))
)

guide_panel <- function() {
  shiny::withMathJax(shiny::div(
    class = "container-md py-3", style = "max-width: 920px;",
    shiny::h3("What this app does"),
    shiny::p("This app is a point and click front end for the R package",
             shiny::a("smam", href = "https://github.com/ChaoranHu/smam", target = "_blank"),
             "(Statistical Modeling of Animal Movements). It fits continuous time movement models to",
             "GPS tracks, estimates the hidden behavioural state at each fix, and simulates tracks from the models."),

    shiny::h4("Workflow"),
    shiny::tags$ol(
      shiny::tags$li(shiny::strong("Data."), "Load an example or upload a CSV file. Choose the time column and the",
                     "coordinate columns and their units. The app converts time to hours and distance to kilometres,",
                     "which is the scale the package is tuned for. Optionally define segments (see below)."),
      shiny::tags$li(shiny::strong("Fit model."), "Pick a model, adjust starting values if needed and press Fit.",
                     "Fits run in the background, so the app stays responsive and you can cancel.",
                     "Standard errors are a separate, optional step because bootstrap methods can take a long time."),
      shiny::tags$li(shiny::strong("States."), "For the MR and MRH models, estimate whether the animal was moving,",
                     "resting or handling at each fix, using the fitted parameters."),
      shiny::tags$li(shiny::strong("Simulate."), "Generate tracks from any model to build intuition or to test",
                     "the estimators. A simulated track can be sent straight to the Data tab.")
    ),

    shiny::h4("Data format"),
    shiny::p("One row per GPS fix. One numeric column with the observation time and one or more numeric columns with",
             "planar coordinates (for example UTM easting and northing). Longitude and latitude should be projected",
             "to a planar system first. Times do not need to be equally spaced."),
    shiny::p("A ", shiny::strong("segment"), " column lets you fit only part of the data, or pool several separate",
             "periods. A value of 0 means the fix is not used, and each other integer marks one segment. The",
             "likelihoods of the segments are added, so for example all summers can be pooled while each summer",
             "is treated as its own track. With a date column the app can build such segments for a season",
             "window in every year, as the smam function seasonFilter does."),

    shiny::h4("Models and notation"),
    shiny::p("Throughout, \\(t\\) is time in hours, \\(X(t)\\) is the true location in kilometres, \\(B(t)\\) is a",
             "standard Brownian motion, and a rate \\(\\lambda\\) of an exponential bout length means the bout lasts",
             "\\(1/\\lambda\\) hours on average. A mobility \\(\\sigma\\) means that after moving for \\(h\\) hours the",
             "displacement along each coordinate has standard deviation \\(\\sigma\\sqrt{h}\\) km."),
    lapply(names(smam_models), function(id) {
      m <- smam_models[[id]]
      shiny::div(class = "mb-3",
                 shiny::h5(m$name),
                 shiny::p(class = "text-muted", m$summary),
                 model_math[[id]])
    }),
    shiny::p("For the MR model, the long run share of time spent moving is",
             "\\(\\lambda_R / (\\lambda_M + \\lambda_R)\\), because the mean moving bout is \\(1/\\lambda_M\\)",
             "and the mean resting bout is \\(1/\\lambda_R\\)."),

    shiny::h4("Standard errors"),
    shiny::tags$ul(
      shiny::tags$li("MR (full likelihood) and BMME: from the Hessian, returned with the fit."),
      shiny::tags$li("MRME: parametric bootstrap (recommended) or Godambe information. Both refit the model many",
                     "times, so start with a small number of bootstrap samples."),
      shiny::tags$li("MM: parametric bootstrap. MRH: numerical Hessian."),
      shiny::tags$li("Approximate 95% intervals are shown as estimate \\(\\pm 1.96 \\times\\) SE, and for bootstrap",
                     "methods also as the 2.5% and 97.5% bootstrap quantiles.")
    ),

    shiny::h4("Tips"),
    shiny::tags$ul(
      shiny::tags$li("MRME fits can take minutes on a few hundred fixes, and MRH fits and MRH state estimation",
                     "are slower still. Try a short segment first; any running computation can be cancelled."),
      shiny::tags$li("If the optimiser does not converge, change the starting values. The defaults for MRME are",
                     "those used for mountain lion f109 in the package vignette."),
      shiny::tags$li("Every fit shows the equivalent R code, so results can be reproduced in a script.")
    ),

    shiny::h4("References"),
    references,
    shiny::p(class = "text-muted small mt-3",
             "smam version ", as.character(utils::packageVersion("smam")),
             ". Package authors: Chaoran Hu, Vladimir Pozdnyakov, Jun Yan.")
  ))
}
