# smam web app: a point and click interface to the smam R package
# (Statistical Modeling of Animal Movements).
# Helpers live in R/ and are loaded automatically by shiny.

library(shiny)
library(bslib)
library(smam)
library(ggplot2)
library(DT)
library(callr)

options(shiny.maxRequestSize = 50 * 1024^2)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- UI helpers ------------------------------------------------------------

param_label <- function(p) {
  unit <- if (nzchar(p$unit)) paste0(" (", p$unit, ")") else ""
  HTML(paste0("<b>", p$label, "</b> ", p$meaning, unit))
}

start_inputs <- function(model) {
  m <- smam_models[[model]]
  lapply(m$params, function(p) {
    numericInput(paste0("start_", model, "_", p$id), param_label(p),
                 value = p$default, min = 0, step = 0.01)
  })
}

bound_inputs <- function(model) {
  m <- smam_models[[model]]
  lapply(m$params, function(p) {
    layout_columns(
      col_widths = c(4, 4, 4),
      div(class = "pt-2", HTML(p$label)),
      numericInput(paste0("lower_", model, "_", p$id), NULL, p$lower, min = 0),
      numericInput(paste0("upper_", model, "_", p$id), NULL, p$upper, min = 0)
    )
  })
}

model_panel <- function(model, ...) {
  conditionalPanel(sprintf("input.model == '%s'", model), ...)
}

sim_param_defaults <- list(
  MR = c(lamM = 1, lamR = 0.5, sigma = 1),
  MRME = c(lamM = 1, lamR = 0.5, sigma = 1, sig_err = 0.05),
  MRH = c(lamM = 4, lamR = 0.5, lamH = 0.1, sigma = 5, p = 0.8),
  MM = c(lam1 = 1, lam2 = 0.1, sigma1 = 1, sigma2 = 0.1),
  BMME = c(sigma = 1, delta = 0.5)
)

sim_inputs <- function(model) {
  m <- smam_models[[model]]
  conditionalPanel(
    sprintf("input.sim_model == '%s'", model),
    lapply(m$params, function(p) {
      numericInput(paste0("sim_", model, "_", p$id), param_label(p),
                   value = sim_param_defaults[[model]][[p$id]], min = 0, step = 0.01)
    }),
    if (model %in% c("MR", "MRME", "MRH", "MM")) {
      s0 <- switch(model,
        MR = , MRME = c(moving = "m", resting = "r"),
        MRH = c(moving = "m", resting = "r", handling = "h"),
        MM = c("moving state 1" = "m1", "moving state 2" = "m2"))
      selectInput(paste0("sim_", model, "_s0"), "State at time 0", s0)
    }
  )
}

# ---- UI --------------------------------------------------------------------

ui <- page_navbar(
  title = tags$span(span(class = "brand-mark", icon(active_theme()$brand_icon)), "smam",
                    span(class = "brand-sub d-none d-md-inline", "animal movement models")),
  window_title = "smam: animal movement models",
  id = "nav",
  theme = app_theme(),
  fillable = FALSE,
  header = tagList(
    withMathJax(),
    tags$head(tags$meta(name = "theme-color", content = active_theme()$pal$ink)),
    tags$style(HTML(app_css())),
    # On phones the results sit below the inputs, so jump to them after a button press.
    tags$script(HTML("
      $(document).on('click', '#fit_go, #state_go, #sim_go', function() {
        if (window.innerWidth >= 576) return;
        var main = $(this).closest('.bslib-sidebar-layout').children('.main')[0];
        if (main) setTimeout(function() { main.scrollIntoView({ behavior: 'smooth' }); }, 300);
      });
      // close the phone menu after picking a tab
      $(document).on('click', '.navbar-collapse.show .nav-link', function() {
        var el = $(this).closest('.navbar-collapse')[0];
        if (window.bootstrap && bootstrap.Collapse) bootstrap.Collapse.getOrCreateInstance(el).hide();
      });
    "))
  ),

  # ---- Data ----
  nav_panel(
    "1. Data", icon = icon("table"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always", max_height_mobile = "auto",
        width = 340,
        tags$h6("Data source", class = "mt-0"),
        radioButtons("src", NULL,
                     c("Example data" = "example", "Upload a CSV file" = "upload")),
        conditionalPanel("input.src == 'example'",
                         selectInput("example", "Example", example_choices)),
        conditionalPanel(
          "input.src == 'upload'",
          fileInput("file", NULL, accept = c(".csv", ".txt", ".tsv")),
          selectInput("sep", "Separator",
                      c(Comma = ",", Semicolon = ";", Tab = "\t", Whitespace = ""))
        ),
        accordion(
          open = FALSE,
          accordion_panel(
            "Columns and units",
            selectInput("time_col", "Time column", NULL),
            selectInput("time_unit", "Time unit in the file", names(time_factors)),
            selectizeInput("coord_cols", "Coordinate columns", NULL, multiple = TRUE),
            selectInput("coord_unit", "Coordinate unit in the file", names(coord_factors)),
            checkboxInput("shift_time", "Start time at 0", TRUE),
            checkboxInput("center", "Centre coordinates at the first fix", TRUE),
            helpText("smam works best with time in hours and distance in km.",
                     "The app converts to these units.")
          ),
          accordion_panel(
            "Segments",
            radioButtons("seg_mode", NULL,
                         c("Use all data as one track" = "none",
                           "Segment column in the file" = "column",
                           "Season window from a date column" = "season")),
            conditionalPanel("input.seg_mode == 'column'",
                             selectInput("seg_col", "Segment column", NULL),
                             helpText("0 = not used, other integers = segment id.")),
            conditionalPanel(
              "input.seg_mode == 'season'",
              selectInput("date_col", "Date column", NULL),
              textInput("date_format", "Date format", "%Y-%m-%d"),
              layout_columns(
                textInput("season_start", "From (MM-DD)", "06-01"),
                textInput("season_end", "To (MM-DD)", "08-31")
              ),
              helpText("Each year's window becomes one segment. A window may wrap",
                       "the new year, for example 12-01 to 02-28.")
            )
          )
        ),
        downloadButton("dl_prepared", "Download prepared data", class = "btn-outline-primary")
      ),
      tab_intro(1, "Load a track",
                "Pick an example or upload a CSV of GPS fixes. The app converts time to hours and distance to kilometres."),
      uiOutput("data_notes"),
      uiOutput("data_boxes"),
      navset_card_tab(
        nav_panel("Track", plotOutput("plot_track", height = "520px")),
        nav_panel("Coordinates over time", plotOutput("plot_coords", height = "420px")),
        nav_panel("Sampling intervals", plotOutput("plot_steps", height = "380px")),
        nav_panel("Prepared data", DTOutput("table_prepared")),
        nav_panel("Raw file", DTOutput("table_raw"))
      )
    )
  ),

  # ---- Fit ----
  nav_panel(
    "2. Fit model", icon = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always", max_height_mobile = "auto",
        width = 380,
        tags$h6("Model", class = "mt-0"),
        selectInput("model", NULL, model_choices, selected = "MRME"),
        uiOutput("model_summary"),
        tags$h6("Starting values"),
        lapply(names(smam_models), function(id) model_panel(id, start_inputs(id))),
        model_panel("MR", radioButtons("MR_likelihood", "Likelihood",
                                       c("Full" = "full", "Composite" = "composite"),
                                       inline = TRUE)),
        model_panel(
          "MRME",
          radioButtons("MRME_variant", "Estimator",
                       c("Composite likelihood (fitMRME, recommended)" = "exact",
                         "Naive composite likelihood (fitMRME_naive)" = "naive",
                         "Discrete error approximation (fitMRMEapprox, very slow)" = "approx")),
          conditionalPanel(
            "input.MRME_variant == 'approx'",
            div(class = "alert alert-warning py-2 small",
                "This estimator is very slow: one likelihood evaluation takes about 2 seconds for",
                "10 fixes, and a fit needs hundreds of them. Use it only on a short track.")
          )
        ),
        model_panel("MRH", numericInput("MRH_threads", "CPU threads", 2, min = 1, step = 1)),
        conditionalPanel(
          "input.model == 'MRH' || (input.model == 'MRME' && input.MRME_variant != 'approx')",
          accordion(
            open = FALSE,
            accordion_panel(
              "Optimizer bounds",
              layout_columns(col_widths = c(4, 4, 4), "", tags$small("lower"), tags$small("upper")),
              model_panel("MRME", bound_inputs("MRME")),
              model_panel("MRH", bound_inputs("MRH"))
            )
          )
        ),
        uiOutput("segment_toggle"),
        layout_columns(
          col_widths = c(8, 4),
          actionButton("fit_go", "Fit model", class = "btn-primary w-100", icon = icon("play")),
          actionButton("fit_cancel", "Cancel", class = "btn-outline-danger w-100")
        )
      ),
      tab_intro(2, "Fit a movement model",
                "Choose a model, set starting values and press Fit model. Fits run in the background and can be cancelled."),
      uiOutput("status_fit"),
      uiOutput("fit_mismatch"),
      card(card_header(textOutput("model_math_title", inline = TRUE)), uiOutput("model_math")),
      uiOutput("fit_results")
    )
  ),

  # ---- States ----
  nav_panel(
    "3. States", icon = icon("shoe-prints"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always", max_height_mobile = "auto",
        width = 360,
        tags$h6("Model", class = "mt-0"),
        radioButtons("state_model", NULL,
                     c("Moving-Resting (MR)" = "MR", "Moving-Resting-Handling (MRH)" = "MRH")),
        radioButtons("state_source", "Parameters",
                     c("From the latest fit" = "fit", "Enter manually" = "manual")),
        uiOutput("state_fit_note"),
        conditionalPanel(
          "input.state_source == 'manual'",
          conditionalPanel("input.state_model == 'MR'",
                           lapply(smam_models$MR$params, function(p)
                             numericInput(paste0("theta_MR_", p$id), param_label(p), p$default, min = 0))),
          conditionalPanel("input.state_model == 'MRH'",
                           lapply(smam_models$MRH$params, function(p)
                             numericInput(paste0("theta_MRH_", p$id), param_label(p), p$default, min = 0)))
        ),
        radioButtons("state_method", "Method",
                     c("Most likely state at each fix (forward backward)" = "fb",
                       "Most likely whole path (Viterbi)" = "viterbi")),
        conditionalPanel("input.state_model == 'MR'",
                         sliderInput("state_cutoff", "Call 'moving' when P(moving) exceeds",
                                     0.05, 0.95, 0.5, step = 0.05)),
        layout_columns(
          col_widths = c(8, 4),
          actionButton("state_go", "Estimate states", class = "btn-primary w-100", icon = icon("play")),
          actionButton("state_cancel", "Cancel", class = "btn-outline-danger w-100")
        ),
        helpText("Uses all fixes of the prepared data, in order.")
      ),
      tab_intro(3, "Estimate behavioural states",
                "Find out when the animal was moving, resting or handling, using the fitted parameters."),
      uiOutput("status_state"),
      uiOutput("state_results")
    )
  ),

  # ---- Simulate ----
  nav_panel(
    "4. Simulate", icon = icon("shuffle"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always", max_height_mobile = "auto",
        width = 360,
        tags$h6("Model", class = "mt-0"),
        selectInput("sim_model", NULL, model_choices, selected = "MR"),
        tags$h6("Parameters"),
        lapply(names(smam_models), sim_inputs),
        tags$h6("Observation times"),
        layout_columns(
          numericInput("sim_T", "Duration (hr)", 100, min = 1),
          numericInput("sim_n", "Number of fixes", 200, min = 3, step = 1)
        ),
        radioButtons("sim_spacing", NULL,
                     c("Regular spacing" = "regular", "Irregular (random) spacing" = "random"),
                     inline = TRUE),
        layout_columns(
          numericInput("sim_dim", "Dimensions", 2, min = 1, max = 3, step = 1),
          numericInput("sim_seed", "Random seed", 123, step = 1)
        ),
        actionButton("sim_go", "Simulate", class = "btn-primary w-100", icon = icon("dice")),
        uiOutput("sim_actions")
      ),
      tab_intro(4, "Simulate tracks",
                "Generate tracks from any model to build intuition, then send one to the Data tab to test the estimators."),
      uiOutput("sim_results")
    )
  ),

  nav_panel("Guide", icon = icon("book"), guide_panel()),
  nav_spacer(),
  nav_item(tags$a(icon("github"), "smam on GitHub",
                  href = "https://github.com/ChaoranHu/smam", target = "_blank"))
)

# ---- Server ----------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(job = NULL, fit = NULL, se = NULL, states = NULL,
                       sim = NULL, uploaded_sim = NULL)

  session$onSessionEnded(function() cancel_job(isolate(rv$job)))

  # ---- data loading ----

  raw_data <- reactive({
    if (input$src == "example") {
      if (!is.null(rv$uploaded_sim) && identical(input$example, "simulated")) return(rv$uploaded_sim)
      return(example_data(input$example))
    }
    req(input$file)
    sep <- input$sep
    tryCatch(
      utils::read.table(input$file$datapath, header = TRUE, sep = sep,
                        stringsAsFactors = FALSE, check.names = FALSE,
                        quote = "\"'", comment.char = ""),
      error = function(e) {
        showNotification(paste("Could not read the file:", conditionMessage(e)), type = "error")
        NULL
      }
    )
  })

  observeEvent(raw_data(), {
    d <- raw_data()
    req(d)
    nm <- names(d)
    map <- if (input$src == "example") example_mapping[[input$example]] else NULL
    if (is.null(map)) {
      tc <- guess_time_col(d)
      map <- list(time = tc, coords = guess_coord_cols(d, tc), time_unit = "hours",
                  coord_unit = "km", segment = "none")
    }
    updateSelectInput(session, "time_col", choices = nm, selected = map$time)
    updateSelectizeInput(session, "coord_cols", choices = nm, selected = map$coords)
    updateSelectInput(session, "time_unit", selected = map$time_unit)
    updateSelectInput(session, "coord_unit", selected = map$coord_unit)
    updateSelectInput(session, "seg_col", choices = nm,
                      selected = if (!is.null(map$segment_col)) map$segment_col else nm[length(nm)])
    date_like <- nm[vapply(d, function(x) is.character(x) || inherits(x, "Date"), TRUE)]
    updateSelectInput(session, "date_col", choices = nm,
                      selected = if (!is.null(map$date_col)) map$date_col else
                        if (length(date_like)) date_like[1] else nm[1])
    updateRadioButtons(session, "seg_mode", selected = map$segment)
  })

  segment_vector <- reactive({
    d <- raw_data()
    req(d)
    switch(input$seg_mode,
      none = NULL,
      column = {
        req(input$seg_col %in% names(d))
        v <- suppressWarnings(as.integer(as.numeric(d[[input$seg_col]])))
        validate(need(!all(is.na(v)), "The segment column must contain integers."))
        v[is.na(v)] <- 0L
        v
      },
      season = {
        req(input$date_col %in% names(d))
        dates <- as.Date(as.character(d[[input$date_col]]), format = input$date_format)
        validate(
          need(sum(!is.na(dates)) > 0, "No dates could be read. Check the date column and format."),
          need(grepl("^\\d{2}-\\d{2}$", input$season_start) && grepl("^\\d{2}-\\d{2}$", input$season_end),
               "Season dates must look like MM-DD.")
        )
        tryCatch(season_segments(dates, input$season_start, input$season_end),
                 error = function(e) validate(need(FALSE, conditionMessage(e))))
      }
    )
  })

  prepared <- reactive({
    d <- raw_data()
    req(d, input$time_col, input$coord_cols)
    req(input$time_col %in% names(d), all(input$coord_cols %in% names(d)))
    seg <- segment_vector()
    res <- tryCatch(
      prepare_data(d, input$time_col, input$coord_cols, input$time_unit,
                   input$coord_unit, input$shift_time, input$center, seg),
      error = function(e) e
    )
    validate(need(!inherits(res, "error"), if (inherits(res, "error")) conditionMessage(res)))
    if (!is.null(res$data$segment)) {
      validate(need(any(res$data$segment != 0), "No fixes fall in any segment."))
    }
    res
  })

  prepared_data <- reactive(prepared()$data)

  output$data_notes <- renderUI({
    notes <- prepared()$notes
    if (!length(notes)) return(NULL)
    div(class = "alert alert-info py-2", HTML(paste(notes, collapse = "<br>")))
  })

  output$data_boxes <- renderUI({
    d <- prepared_data()
    dt <- diff(d$time)
    nseg <- if (is.null(d$segment)) "none" else length(unique(d$segment[d$segment != 0]))
    used <- if (is.null(d$segment)) nrow(d) else sum(d$segment != 0)
    layout_columns(
      fill = FALSE,
      col_widths = breakpoints(xs = c(6, 6, 6, 6), lg = c(3, 3, 3, 3)),
      value_box(tagList(icon("location-dot"), "Fixes"), format(nrow(d), big.mark = ","),
                p(if (!is.null(d$segment)) paste(format(used, big.mark = ","), "in segments") else "GPS locations"),
                class = "vb"),
      value_box(tagList(icon("clock"), "Duration"), paste(signif(diff(range(d$time)), 4), "hr"),
                p(sprintf("%.1f days", diff(range(d$time)) / 24)), class = "vb"),
      value_box(tagList(icon("stopwatch"), "Median interval"), paste(signif(stats::median(dt), 3), "hr"),
                p(sprintf("range %.3g to %.3g hr", min(dt), max(dt))), class = "vb vb-accent"),
      value_box(tagList(icon("layer-group"), "Segments"), nseg, p(paste(length(coord_names(d)), "coordinate(s)")),
                class = "vb vb-accent")
    )
  })

  output$plot_track <- renderPlot({
    d <- prepared_data()
    plot_track(d, if (!is.null(d$segment)) "segment" else "time")
  }, res = 96)
  output$plot_coords <- renderPlot(plot_coords(prepared_data()), res = 96)
  output$plot_steps <- renderPlot(plot_steps(prepared_data()), res = 96)
  output$table_prepared <- renderDT({
    datatable(prepared_data(), rownames = FALSE, fillContainer = FALSE, class = "compact hover", options = list(pageLength = 10, scrollX = TRUE, dom = "tip", pagingType = "simple")) |>
      formatSignif(columns = c("time", coord_names(prepared_data())), digits = 6)
  })
  output$table_raw <- renderDT({
    req(raw_data())
    datatable(raw_data(), rownames = FALSE, fillContainer = FALSE, class = "compact hover", options = list(pageLength = 10, scrollX = TRUE, dom = "tip", pagingType = "simple"))
  })
  output$dl_prepared <- downloadHandler(
    filename = function() "prepared_data.csv",
    content = function(file) utils::write.csv(prepared_data(), file, row.names = FALSE)
  )

  # ---- background jobs ----

  launch <- function(task, spec, label) {
    if (job_running(rv$job)) {
      showNotification("Another computation is still running. Cancel it first.", type = "warning")
      return(FALSE)
    }
    rv$job <- start_job(task, spec, label)
    TRUE
  }

  observe({
    job <- rv$job
    req(job)
    invalidateLater(700)
    if (job_running(job)) return()
    res <- job_result(job)
    isolate({
      rv$job <- NULL
      if (!res$ok) {
        hint <- if (grepl("singular", res$error)) " Try more bootstrap samples or other starting values." else ""
        showNotification(paste0(job$label, " failed: ", res$error, hint), type = "error", duration = NULL)
        return()
      }
      v <- res$value
      if (job$task == "fit") {
        rv$fit <- c(job$spec[c("model", "start", "options", "use_segment", "data")], v)
        rv$se <- if (!is.null(v$vcov)) list(vcov = v$vcov, boot = NULL, method = "Hessian") else NULL
        showNotification(sprintf("%s finished in %s.", job$label, format_seconds(v$elapsed)),
                         type = "message")
      } else if (job$task == "se") {
        rv$se <- list(vcov = v$vcov, boot = v$boot, method = job$spec$method)
        showNotification(sprintf("Standard errors finished in %s.", format_seconds(v$elapsed)),
                         type = "message")
      } else if (job$task == "state") {
        rv$states <- c(job$spec[c("model", "method", "theta")], v)
      }
    })
  })

  status_ui <- function(tasks) {
    job <- rv$job
    if (is.null(job) || !job$task %in% tasks) return(NULL)
    invalidateLater(1000)
    card(
      class = "status-card",
      card_header(
        class = "d-flex justify-content-between align-items-center",
        span(span(class = "spinner-border spinner-border-sm text-primary me-2"), job$label),
        span(class = "text-muted small fw-normal", "running for ", format_seconds(job_elapsed(job)))
      ),
      if (nzchar(log <- job_log_tail(job))) tags$pre(class = "small-log mb-0", log)
      else card_body(class = "text-muted small", "Working. Press Cancel to stop.")
    )
  }
  output$status_fit <- renderUI(status_ui(c("fit", "se")))
  output$status_state <- renderUI(status_ui("state"))

  observeEvent(input$fit_cancel, {
    if (job_running(rv$job) && rv$job$task %in% c("fit", "se")) {
      cancel_job(rv$job)
      rv$job <- NULL
      showNotification("Computation cancelled.", type = "warning")
    }
  })
  observeEvent(input$state_cancel, {
    if (job_running(rv$job) && rv$job$task == "state") {
      cancel_job(rv$job)
      rv$job <- NULL
      showNotification("Computation cancelled.", type = "warning")
    }
  })

  # ---- fitting ----

  output$model_summary <- renderUI(p(class = "text-muted small", smam_models[[input$model]]$summary))
  output$model_math_title <- renderText(paste("Model:", smam_models[[input$model]]$name))
  output$model_math <- renderUI(withMathJax(model_math[[input$model]]))

  output$fit_mismatch <- renderUI({
    f <- rv$fit
    if (is.null(f) || identical(f$model, input$model)) return(NULL)
    div(class = "alert alert-warning py-2 small",
        sprintf("The estimates below are from the last fit (%s). Press Fit model to fit %s.",
                smam_models[[f$model]]$name, smam_models[[input$model]]$name))
  })

  has_segments <- reactive({
    d <- tryCatch(prepared_data(), error = function(e) NULL)
    !is.null(d) && !is.null(d$segment)
  })

  output$segment_toggle <- renderUI({
    if (!has_segments()) return(NULL)
    if (!smam_models[[input$model]]$segment) {
      return(div(class = "alert alert-warning py-2 small",
                 "This model has no segment option, so all fixes are used."))
    }
    checkboxInput("use_segment", "Fit the segments defined on the Data tab", TRUE)
  })

  collect_vals <- function(prefix, model) {
    ids <- vapply(smam_models[[model]]$params, `[[`, "", "id")
    vapply(ids, function(id) {
      v <- input[[paste0(prefix, model, "_", id)]]
      if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
    }, 0)
  }

  observeEvent(input$fit_go, {
    d <- tryCatch(prepared_data(), error = function(e) NULL)
    if (is.null(d)) {
      showNotification("Prepare the data on the Data tab first.", type = "error")
      return()
    }
    model <- input$model
    m <- smam_models[[model]]
    start <- collect_vals("start_", model)
    if (model == "BMME" && all(is.na(start))) {
      start <- NULL
    } else if (anyNA(start) || any(start <= 0)) {
      showNotification("All starting values must be positive numbers.", type = "error")
      return()
    }
    if (model == "MRH" && start[["p"]] >= 1) {
      showNotification("The switching probability p must be below 1.", type = "error")
      return()
    }
    use_seg <- m$segment && has_segments() && isTRUE(input$use_segment)
    options <- list()
    if (model == "MR") options$likelihood <- input$MR_likelihood
    if (model == "MRME") options$variant <- input$MRME_variant
    if (model %in% c("MRME", "MRH")) {
      options$lower <- unname(collect_vals("lower_", model))
      options$upper <- unname(collect_vals("upper_", model))
      if (anyNA(c(options$lower, options$upper)) || any(options$lower >= options$upper)) {
        showNotification("Each lower bound must be below its upper bound.", type = "error")
        return()
      }
      if (any(start < options$lower | start > options$upper)) {
        showNotification("Starting values must lie inside the optimizer bounds.", type = "error")
        return()
      }
    }
    if (model == "MRH") options$numThreads <- max(1, as.integer(input$MRH_threads))
    if (!use_seg) d$segment <- NULL
    spec <- list(model = model, start = if (is.null(start)) NULL else unname(start),
                 options = options, use_segment = use_seg, data = d)
    launch("fit", spec, paste("Fitting", m$name))
  })

  output$fit_results <- renderUI({
    f <- rv$fit
    if (is.null(f)) {
      return(empty_state("chart-line", "Choose a model and press Fit model. The estimates, what they mean and their standard errors appear here."))
    }
    m <- smam_models[[f$model]]
    variant <- if (is.null(f$options$variant)) "" else f$options$variant
    n_used <- if (isTRUE(f$use_segment)) sum(f$data$segment != 0) else nrow(f$data)
    tagList(
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center gap-2 flex-wrap",
          span(paste("Estimates:", m$name)),
          if (convergence_ok(f$model, variant, f$convergence))
            span(class = "badge rounded-pill text-bg-success", icon("check"), "converged")
          else span(class = "badge rounded-pill text-bg-warning", icon("triangle-exclamation"), "not converged")
        ),
        card_body(
          p(class = "text-muted small mb-1",
            sprintf("Fitted to %s fixes%s in %s. Log likelihood %s. Convergence: %s.",
                    format(n_used, big.mark = ","),
                    if (isTRUE(f$use_segment)) paste0(" in ", length(unique(f$data$segment[f$data$segment != 0])), " segments") else "",
                    format_seconds(f$elapsed), format(f$loglik, digits = 6),
                    convergence_note(f$model, variant, f$convergence))),
          div(class = "table-responsive param-wrap", tableOutput("est_table")),
          tags$h6("What the estimates mean", class = "section-label"),
          div(class = "table-responsive param-wrap", tableOutput("derived_table")),
          div(class = "d-flex gap-2",
              downloadButton("dl_estimates", "Download estimates", class = "btn-outline-primary btn-sm"))
        )
      ),
      card(card_header("Standard errors and confidence intervals"),
           card_body(uiOutput("se_have"), uiOutput("se_panel"))),
      card(card_header("Reproduce in R"), tags$pre(fit_code(f$model, f$start, f$options, f$use_segment)))
    )
  })

  est_df <- reactive({
    f <- rv$fit
    req(f)
    m <- smam_models[[f$model]]
    df <- data.frame(
      Parameter = vapply(m$params, function(p) p$label, ""),
      Meaning = vapply(m$params, function(p) p$meaning, ""),
      Unit = vapply(m$params, function(p) p$unit, ""),
      Estimate = signif(f$estimate, 5),
      stringsAsFactors = FALSE
    )
    se <- rv$se
    if (!is.null(se) && !is.null(se$vcov) && nrow(se$vcov) == length(f$estimate)) {
      s <- sqrt(pmax(diag(se$vcov), 0))
      df$SE <- signif(s, 4)
      df$`95% CI (normal)` <- sprintf("%.4g to %.4g", f$estimate - 1.96 * s, f$estimate + 1.96 * s)
      if (!is.null(se$boot)) {
        b <- as.matrix(se$boot)
        q <- apply(b, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
        df$`95% CI (bootstrap)` <- sprintf("%.4g to %.4g", q[1, ], q[2, ])
      }
    }
    df
  })

  # On screen the meaning sits under the symbol, which keeps the table narrow on phones.
  output$est_table <- renderTable({
    df <- est_df()
    df$Parameter <- paste0("<b>", df$Parameter, "</b><br><span class='small text-muted'>",
                           df$Meaning, "</span>")
    df$Meaning <- NULL
    df
  }, sanitize.text.function = identity,
                                  digits = 5, class = "param-table")

  output$derived_table <- renderTable({
    f <- rv$fit
    req(f)
    derived_quantities(f$model, f$estimate)
  }, digits = 4, class = "param-table")

  output$dl_estimates <- downloadHandler(
    filename = function() paste0("smam_", rv$fit$model, "_estimates.csv"),
    content = function(file) {
      df <- est_df()
      df$Parameter <- names(rv$fit$estimate)
      df$Meaning <- gsub("&[a-z]+;", "", df$Meaning)
      df$Unit <- gsub("&radic;", "sqrt ", df$Unit)
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  output$se_have <- renderUI({
    f <- rv$fit
    req(f)
    se <- rv$se
    if (!is.null(se)) {
      return(p(class = "text-success small",
               paste0("Shown in the table above (method: ", se$method, ").",
                      if (!is.null(se$boot)) " Bootstrap intervals use the 2.5% and 97.5% quantiles of the refitted estimates." else "")))
    }
    if (f$model == "BMME" || (f$model == "MR" && f$options$likelihood == "full")) {
      return(p(class = "text-muted", "The Hessian could not be inverted at this estimate, so no standard errors are available."))
    }
    NULL
  })

  output$se_panel <- renderUI({
    f <- rv$fit
    req(f)
    variant <- if (is.null(f$options$variant)) "" else f$options$variant
    switch(f$model,
      MR = if (f$options$likelihood != "full")
        p(class = "text-muted", "Standard errors come with the full likelihood fit. Refit with Likelihood = Full."),
      BMME = NULL,
      {
        if (f$model == "MRME" && variant == "approx") {
          return(p(class = "text-muted", "Not available for the approximation estimator. Use fitMRME for inference."))
        }
        if (isTRUE(f$use_segment)) {
          return(p(class = "text-muted", "The package computes these for a single track only. Refit without segments to get standard errors."))
        }
        methods <- switch(f$model,
          MRME = c("Parametric bootstrap (recommended)" = "pBootstrap", "Godambe information" = "Godambe"),
          MM = c("Parametric bootstrap" = "pBootstrap"),
          MRH = c("Numerical Hessian" = "Hessian"))
        tagList(
          layout_columns(
            col_widths = c(5, 3, 2, 2),
            selectInput("se_method", "Method", methods),
            conditionalPanel("input.se_method != 'Hessian'",
                             numericInput("se_nBS", "Bootstrap samples", 25, min = 2, step = 1)),
            numericInput("se_threads", "Threads", 2, min = 1, step = 1),
            div(class = "pt-4", actionButton("se_go", "Compute", class = "btn-primary"))
          ),
          helpText("Each bootstrap sample refits the model, so this takes about nBS times the fit time",
                   "divided by the number of threads.")
        )
      }
    )
  })

  observeEvent(input$se_go, {
    f <- rv$fit
    req(f)
    spec <- list(model = f$model, variant = f$options$variant, estimate = f$estimate,
                 data = f$data, use_segment = FALSE, method = input$se_method,
                 nBS = max(2, as.integer(input$se_nBS %||% 25)),
                 numThreads = max(1, as.integer(input$se_threads)))
    # Godambe estimates a 4 by 4 variability matrix from the bootstrap samples,
    # which is singular with only a few of them
    if (identical(spec$method, "Godambe") && spec$nBS < 10) {
      showNotification("Godambe information needs at least 10 bootstrap samples.", type = "error")
      return()
    }
    launch("se", spec, "Computing standard errors")
  })

  # ---- states ----

  state_theta <- reactive({
    model <- input$state_model
    if (input$state_source == "manual") return(collect_vals("theta_", model))
    f <- rv$fit
    if (is.null(f)) return(NULL)
    if (model == "MR" && f$model %in% c("MR", "MRME")) return(f$estimate[c("lamM", "lamR", "sigma")])
    if (model == "MRH" && f$model == "MRH") return(f$estimate)
    NULL
  })

  output$state_fit_note <- renderUI({
    if (input$state_source != "fit") return(NULL)
    th <- state_theta()
    if (is.null(th)) {
      return(div(class = "alert alert-warning py-2 small",
                 sprintf("No %s fit yet. Fit %s on the Fit tab or enter parameters manually.",
                         input$state_model,
                         if (input$state_model == "MR") "MR or MRME" else "MRH")))
    }
    note <- if (input$state_model == "MR" && rv$fit$model == "MRME")
      " (from the MRME fit; measurement error is ignored here)" else ""
    div(class = "alert alert-light py-2 small",
        paste0("Using ", paste(names(th), signif(th, 4), sep = " = ", collapse = ", "), note))
  })

  observeEvent(input$state_go, {
    d <- tryCatch(prepared_data(), error = function(e) NULL)
    if (is.null(d)) {
      showNotification("Prepare the data on the Data tab first.", type = "error")
      return()
    }
    th <- state_theta()
    if (is.null(th) || anyNA(th) || any(th <= 0)) {
      showNotification("Provide positive parameter values.", type = "error")
      return()
    }
    d$segment <- NULL
    spec <- list(model = input$state_model, method = input$state_method, theta = th,
                 cutoff = input$state_cutoff, data = d, use_segment = FALSE)
    launch("state", spec, paste("Estimating", input$state_model, "states"))
  })

  state_labels <- reactive({
    s <- rv$states
    req(s)
    res <- s$states
    if (s$model == "MR") ifelse(res$states == 1, "moving", "resting")
    else c("moving", "resting", "handling")[res$states + 1]
  })

  output$state_results <- renderUI({
    s <- rv$states
    if (is.null(s)) {
      return(empty_state("shoe-prints", "Press Estimate states to see the track coloured by behaviour and the probability of each state over time."))
    }
    lab <- state_labels()
    share <- round(100 * prop.table(table(factor(lab, levels = smam_models[[s$model]]$states))), 1)
    tagList(
      card(
        card_header(sprintf("%s states (%s)", s$model,
                            if (s$method == "fb") "forward backward" else "Viterbi")),
        p(class = "small text-muted",
          paste0("Share of fixes: ", paste(names(share), paste0(share, "%"), sep = " ", collapse = ", "),
                 ". Computed in ", format_seconds(s$elapsed), ".")),
        layout_columns(
          col_widths = c(6, 6),
          plotOutput("state_track", height = "460px"),
          div(
            if (s$method == "fb") plotOutput("state_prob", height = "260px"),
            plotOutput("state_path", height = "200px")
          )
        )
      ),
      card(card_header("State table"), DTOutput("state_table"),
           downloadButton("dl_states", "Download states", class = "btn-outline-primary btn-sm m-2"))
    )
  })

  output$state_track <- renderPlot({
    s <- rv$states
    req(s)
    cn <- setdiff(names(s$states), c("time", "segment", "p.m", "p.r", "p.h", "states"))
    plot_track(s$states[, c("time", cn)], "state", state_labels())
  }, res = 96)
  output$state_prob <- renderPlot({
    s <- rv$states
    req(s, s$method == "fb")
    plot_state_prob(s$states, smam_models[[s$model]]$states)
  }, res = 96)
  output$state_path <- renderPlot(plot_state_path(rv$states$states, state_labels()), res = 96)
  output$state_table <- renderDT({
    s <- rv$states
    req(s)
    df <- s$states
    df$state <- state_labels()
    datatable(df, rownames = FALSE, fillContainer = FALSE, class = "compact hover", options = list(pageLength = 10, scrollX = TRUE, dom = "tip", pagingType = "simple")) |>
      formatSignif(columns = setdiff(names(df), c("state", "states")), digits = 5)
  })
  output$dl_states <- downloadHandler(
    filename = function() paste0("smam_", rv$states$model, "_states.csv"),
    content = function(file) {
      df <- rv$states$states
      df$state <- state_labels()
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  # ---- simulation ----

  observeEvent(input$sim_go, {
    model <- input$sim_model
    par <- vapply(smam_models[[model]]$params, function(p) {
      v <- input[[paste0("sim_", model, "_", p$id)]]
      if (is.null(v) || is.na(v)) NA_real_ else v
    }, 0)
    names(par) <- vapply(smam_models[[model]]$params, `[[`, "", "id")
    n <- as.integer(input$sim_n)
    Tend <- input$sim_T
    dim <- as.integer(input$sim_dim)
    if (anyNA(par) || any(par <= 0) || is.na(n) || n < 3 || is.na(Tend) || Tend <= 0 ||
        is.na(dim) || dim < 1) {
      showNotification("Check the simulation settings: all values must be positive.", type = "error")
      return()
    }
    if (model == "MRH" && par[["p"]] >= 1) {
      showNotification("p must be below 1.", type = "error")
      return()
    }
    set.seed(input$sim_seed)
    tgrid <- if (input$sim_spacing == "regular") seq(0, Tend, length.out = n)
             else c(0, sort(stats::runif(n - 1, 0, Tend)))
    s0 <- input[[paste0("sim_", model, "_s0")]]
    res <- tryCatch(switch(model,
      MR = rMR(tgrid, par[["lamM"]], par[["lamR"]], par[["sigma"]], s0, dim = dim, state = TRUE),
      MRME = rMRME(tgrid, par[["lamM"]], par[["lamR"]], par[["sigma"]], par[["sig_err"]], s0,
                   dim = dim, state = TRUE),
      MRH = rMRH(tgrid, par[["lamM"]], par[["lamR"]], par[["lamH"]], par[["sigma"]], par[["p"]],
                 s0, dim = dim, state = TRUE),
      MM = rMM(tgrid, par[["lam1"]], par[["lam2"]], par[["sigma1"]], par[["sigma2"]], s0, dim = dim),
      BMME = rBMME(tgrid, dim = dim, sigma = par[["sigma"]], delta = par[["delta"]])
    ), error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Simulation failed:", conditionMessage(res)), type = "error")
      return()
    }
    st <- NULL
    if ("state" %in% names(res)) {
      st <- if (model == "MRH") c("moving", "resting", "handling")[res$state + 1]
            else ifelse(res$state == 1, "moving", "resting")
      res$state <- NULL
    }
    names(res)[-1] <- paste0("X", seq_len(ncol(res) - 1))
    res$time <- tgrid
    rv$sim <- list(model = model, par = par, data = res, states = st)
  })

  output$sim_actions <- renderUI({
    req(rv$sim)
    div(class = "d-grid gap-2 mt-2",
        actionButton("sim_use", "Use this track on the Data tab", class = "btn-success",
                     icon = icon("arrow-right")),
        downloadButton("dl_sim", "Download CSV", class = "btn-outline-primary"))
  })

  output$sim_results <- renderUI({
    if (is.null(rv$sim)) {
      return(empty_state("shuffle", "Set the model parameters and press Simulate. The simulated track and its data appear here."))
    }
    s <- rv$sim
    tagList(
      card(card_header(paste("Simulated", smam_models[[s$model]]$name, "track")),
           layout_columns(col_widths = c(6, 6),
                          plotOutput("sim_track", height = "440px"),
                          plotOutput("sim_coords", height = "440px"))),
      card(card_header("Simulated data"), DTOutput("sim_table"))
    )
  })

  sim_frame <- reactive({
    s <- rv$sim
    req(s)
    d <- s$data
    if (!is.null(s$states)) d$true_state <- s$states
    d
  })

  output$sim_track <- renderPlot({
    s <- rv$sim
    req(s)
    if (!is.null(s$states)) plot_track(s$data, "state", s$states, state_title = "True state") else plot_track(s$data, "time")
  }, res = 96)
  output$sim_coords <- renderPlot(plot_coords(rv$sim$data), res = 96)
  output$sim_table <- renderDT({
    d <- sim_frame()
    datatable(d, rownames = FALSE, fillContainer = FALSE, class = "compact hover", options = list(pageLength = 10, scrollX = TRUE, dom = "tip", pagingType = "simple")) |>
      formatSignif(columns = c("time", coord_names(rv$sim$data)), digits = 5)
  })
  output$dl_sim <- downloadHandler(
    filename = function() paste0("smam_simulated_", rv$sim$model, ".csv"),
    content = function(file) utils::write.csv(sim_frame(), file, row.names = FALSE)
  )

  observeEvent(input$sim_use, {
    req(rv$sim)
    rv$uploaded_sim <- sim_frame()
    choices <- c(example_choices, "Simulated track from the Simulate tab" = "simulated")
    updateRadioButtons(session, "src", selected = "example")
    updateSelectInput(session, "example", choices = choices, selected = "simulated")
    nav_select("nav", "1. Data")
    showNotification("The simulated track is now the active data set.", type = "message")
  })
}

shinyApp(ui, server)
