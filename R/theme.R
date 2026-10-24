# Look and feel: colours, fonts, component styles and the matching ggplot theme.
# One palette is used everywhere so the plots, tables and controls read as one design.

smam_pal <- list(
  ink      = "#17302B",  # navbar and headings
  text     = "#24312F",
  muted    = "#64716D",
  page     = "#F4F5F0",  # page background, a warm off white
  surface  = "#FFFFFF",  # cards and sidebar
  line     = "#E2E5DD",
  grid     = "#ECEEE8",
  primary  = "#1F6F5F",  # deep teal green
  tint     = "#E6F0EC",  # light primary for soft backgrounds
  accent   = "#D9822B",  # ochre
  blue     = "#2F6690",
  olive    = "#7A9A3A",
  danger   = "#B44A3C"
)

app_theme <- function() {
  inter <- font_link(
    "Inter",
    href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"
  )
  bs_theme(
    version = 5,
    bg = smam_pal$page, fg = smam_pal$text,
    primary = smam_pal$primary, secondary = smam_pal$muted,
    success = "#2E7D4F", info = smam_pal$blue, warning = smam_pal$accent,
    danger = smam_pal$danger,
    base_font = font_collection(inter, "system-ui", "-apple-system", "Segoe UI", "Roboto",
                                "Helvetica Neue", "Arial", "sans-serif"),
    heading_font = font_collection(inter, "system-ui", "-apple-system", "Segoe UI", "sans-serif"),
    code_font = font_collection("ui-monospace", "SFMono-Regular", "Menlo", "Consolas", "monospace"),
    "font-size-base" = "0.95rem",
    "headings-font-weight" = 600,
    "border-radius" = "0.6rem",
    "border-radius-sm" = "0.45rem",
    "border-radius-lg" = "0.9rem",
    "border-color" = smam_pal$line,
    "navbar-bg" = smam_pal$ink,
    "input-bg" = smam_pal$surface,
    "link-color" = smam_pal$primary
  )
}

app_css <- function() {
  p <- smam_pal
  css <- "
  body { background: {page}; letter-spacing: -0.005em; }

  /* navbar */
  .navbar.navbar-static-top, nav.navbar { background-color: {ink} !important; }
  .navbar .navbar-brand, .navbar .nav-link { color: rgba(255,255,255,.82) !important; }
  .navbar .navbar-brand, .navbar .nav-link.active, .navbar .nav-link:hover { color: #FFFFFF !important; }
  .navbar .navbar-toggler { border-color: rgba(255,255,255,.3); }
  .navbar .navbar-toggler-icon { filter: invert(1) brightness(2); }
  .navbar { box-shadow: 0 1px 0 rgba(0,0,0,.08), 0 6px 18px rgba(23,48,43,.12); }
  .navbar-brand { display: flex; align-items: center; gap: .55rem; font-weight: 700; letter-spacing: -0.01em; }
  .brand-mark { display: inline-flex; align-items: center; justify-content: center; width: 2rem; height: 2rem;
                border-radius: .55rem; background: {accent}; color: #fff; font-size: .95rem; }
  .brand-sub { font-weight: 400; opacity: .7; font-size: .9rem; }
  .navbar .nav-link { font-weight: 500; border-radius: .5rem; padding: .45rem .8rem !important; }
  .navbar .nav-link.active { background: rgba(255,255,255,.12); }

  /* layout: the sidebar runs full height from the navbar down */
  .container-fluid > .tab-content > .tab-pane:has(> .bslib-sidebar-layout) { margin-top: 0 !important; }
  .container-fluid > .tab-content > .tab-pane > .bslib-sidebar-layout { margin: 0 -12px !important;
     min-height: calc(100vh - 58px); }
  .bslib-sidebar-layout { --bslib-sidebar-bg: {surface}; --bslib-sidebar-border: 1px solid {line};
                          border: 0 !important; border-radius: 0 !important; }
  .bslib-sidebar-layout > .sidebar { box-shadow: 1px 0 0 {line}; }
  .bslib-sidebar-layout > .main { padding: 1.5rem 1.75rem; }
  /* an empty output (no status, no warning) should not leave a gap; it must stay displayed so shiny renders it */
  .bslib-sidebar-layout > .main > .shiny-html-output:empty { margin-top: -16px; }
  .sidebar .form-label, .sidebar label.control-label { font-weight: 500; font-size: .875rem; color: {text}; }
  .sidebar h6, .section-label { text-transform: uppercase; letter-spacing: .06em; font-size: .72rem;
                                font-weight: 600; color: {muted}; margin: 1rem 0 .5rem; }
  .sidebar .help-block, .sidebar .form-text { font-size: .8rem; color: {muted}; }
  .sidebar .accordion { --bs-accordion-bg: transparent; --bs-accordion-border-color: {line};
                        --bs-accordion-active-bg: {tint}; --bs-accordion-active-color: {ink}; }
  .sidebar .accordion-button { font-weight: 600; font-size: .9rem; }

  /* page headers inside each tab */
  .tab-intro { display: flex; gap: .9rem; align-items: flex-start; margin-bottom: 1.25rem; }
  .tab-intro .step { flex: 0 0 auto; width: 2.25rem; height: 2.25rem; border-radius: 50%;
                     background: {tint}; color: {primary}; font-weight: 700; display: flex;
                     align-items: center; justify-content: center; }
  .tab-intro h2 { font-size: 1.35rem; margin: 0 0 .15rem; color: {ink}; }
  .tab-intro p { margin: 0; color: {muted}; max-width: 46rem; }

  /* cards */
  .card { --bs-card-bg: {surface}; --bs-card-border-color: {line}; --bs-card-cap-bg: transparent;
          border-radius: .9rem; box-shadow: 0 1px 2px rgba(23,48,43,.04), 0 2px 8px rgba(23,48,43,.05);
          margin-bottom: 1.25rem; }
  .card-header { font-weight: 600; color: {ink}; padding: .85rem 1.1rem; border-bottom: 1px solid {line}; }
  .card-header .badge { font-weight: 600; }
  .card-body { padding: 1rem 1.1rem; }
  .nav-tabs.card-header-tabs .nav-link { color: {muted}; font-weight: 500; border: 0;
                                         border-bottom: 2px solid transparent; }
  .nav-tabs.card-header-tabs .nav-link.active { color: {primary}; background: transparent;
                                                border-bottom-color: {primary}; }

  /* value boxes */
  .bslib-value-box.vb { --bslib-color-bg: {surface}; --bslib-color-fg: {text};
                        background: {surface}; color: {text}; border: 1px solid {line}; }
  .bslib-value-box.vb { min-height: 0 !important; }
  .bslib-value-box.vb .value-box-area { padding: 1rem 1.15rem !important; }
  .bslib-value-box.vb .value-box-title { display: flex !important; flex-direction: row !important;
     align-items: center; gap: .1rem; margin-bottom: .2rem !important; }
  .bslib-value-box.vb .value-box-title { color: {muted}; font-size: .8rem; font-weight: 500;
                                         text-transform: uppercase; letter-spacing: .05em; }
  .bslib-value-box.vb .value-box-value { color: {ink}; font-weight: 700; font-size: 1.6rem; }
  .bslib-value-box.vb p { color: {muted}; font-size: .82rem; margin: 0; }
  .bslib-value-box.vb .value-box-title svg { color: {primary}; margin-right: .35rem; }
  .bslib-value-box.vb.vb-accent .value-box-title svg { color: {accent}; }

  /* buttons and inputs */
  .btn { font-weight: 500; }
  .btn-primary { box-shadow: 0 1px 2px rgba(23,48,43,.2); }
  .form-control, .form-select, .selectize-input { border-color: #D6DAD1; }
  .form-control:focus, .form-select:focus, .selectize-input.focus { border-color: {primary};
     box-shadow: 0 0 0 .2rem rgba(31,111,95,.15); }
  .irs--shiny .irs-bar, .irs--shiny .irs-single { background: {primary}; border-color: {primary}; }

  /* tables */
  .param-wrap table { width: 100% !important; margin-bottom: 1rem; --bs-table-bg: transparent; }
  .param-wrap table thead th { font-size: .72rem; text-transform: uppercase; letter-spacing: .05em;
                               color: {muted}; font-weight: 600; border-bottom: 1px solid {line} !important;
                               background: transparent; }
  .param-wrap table td { border-color: {grid} !important; vertical-align: middle; padding: .55rem .6rem; }
  .param-wrap table tbody tr:last-child td { border-bottom: 0; }
  #est_table td:first-child { min-width: 10rem; }
  #est_table td:nth-child(3) { font-weight: 600; font-variant-numeric: tabular-nums; color: {ink}; }
  .nav-tabs.card-header-tabs { flex-wrap: nowrap; overflow-x: auto; scrollbar-width: none; }
  .nav-tabs.card-header-tabs .nav-link { white-space: nowrap; }
  table.dataTable { --bs-table-bg: transparent; }
  .dataTables_scrollHead table.dataTable thead th { border-bottom: 1px solid {line} !important; }
  .dataTables_scrollBody table.dataTable, .dataTables_scrollBody table.dataTable thead th { border-top: 0 !important; border-bottom: 0 !important; }
  .dataTables_wrapper.no-footer .dataTables_scrollBody { border-bottom: 1px solid {line} !important; }
  table.dataTable > thead { border-bottom: 0 !important; }
  .dataTables_paginate .pagination { gap: .35rem; margin-top: .5rem; }
  .dataTables_paginate .page-link { border-radius: .5rem !important; border: 1px solid {line}; color: {primary};
     background: {surface}; font-weight: 500; padding: .3rem .8rem; }
  .dataTables_paginate .page-item.disabled .page-link { color: #B7BDB3; background: {surface}; }
  .card .datatables.html-widget { height: auto !important; flex: 0 0 auto !important; }
  .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_paginate { font-size: .8rem; color: {muted}; padding-top: .6rem; }
  table.dataTable thead th { font-size: .78rem; color: {muted}; font-weight: 600; }
  table.dataTable td { font-variant-numeric: tabular-nums; }
  .table-responsive table { min-width: 100%; }

  /* status, alerts, empty states, code */
  .status-card { border-left: 4px solid {primary} !important; }
  .small-log { font-size: .78rem; max-height: 260px; overflow-y: auto; margin: 0; background: #F7F8F5;
               padding: .75rem 1rem; border-radius: 0 0 .9rem .9rem; }
  .alert { border-radius: .7rem; border-width: 0; border-left: 4px solid; }
  .alert-warning { background: #FBF1E3; border-color: {accent}; color: #6B4310; }
  .alert-info { background: #EAF1F6; border-color: {blue}; color: #1E3F58; }
  .alert-light { background: {tint}; border-color: {primary}; color: {ink}; }
  .empty-state { text-align: center; padding: 2.5rem 1rem; color: {muted}; }
  .empty-state .icon { font-size: 1.6rem; color: {primary}; background: {tint}; width: 3.5rem; height: 3.5rem;
                       border-radius: 50%; display: inline-flex; align-items: center; justify-content: center;
                       margin-bottom: .75rem; }
  .empty-state p { margin: 0 auto; max-width: 28rem; }
  pre { white-space: pre-wrap; background: #F7F8F5; border: 1px solid {line}; border-radius: .6rem;
        padding: .8rem 1rem; font-size: .82rem; color: {ink}; }
  .card > pre { margin: 1rem 1.1rem; }

  .shiny-notification { background: {surface}; color: {text}; border: 1px solid {line}; border-left: 4px solid {primary};
                        border-radius: .7rem; box-shadow: 0 8px 24px rgba(23,48,43,.15); opacity: 1; }
  .shiny-notification-warning { border-left-color: {accent}; }
  .shiny-notification-error { border-left-color: {danger}; }

  /* guide */
  .guide { line-height: 1.6; }
  .guide h3 { color: {ink}; }
  .guide h4 { color: {ink}; margin-top: 2rem; padding-top: 1rem; border-top: 1px solid {line}; }
  .guide .model-block { background: {surface}; border: 1px solid {line}; border-radius: .9rem;
                        padding: 1rem 1.25rem; margin-bottom: 1rem; }
  .guide .model-block h5 { color: {primary}; }

  /* keep wide formulas and tables inside their card */
  .MathJax_Display, .MJXc-display, mjx-container[display] { overflow-x: auto; overflow-y: hidden; max-width: 100%; }

  /* phones: inputs first, then results */
  @media (max-width: 575.98px) {
    html body .bslib-sidebar-layout.bslib-sidebar-layout { display: flex !important; flex-direction: column !important; }
    html body .bslib-sidebar-layout.bslib-sidebar-layout > .sidebar { order: 1 !important; width: 100% !important;
       max-height: none !important; box-shadow: 0 1px 0 {line}; }
    html body .bslib-sidebar-layout.bslib-sidebar-layout > .main { order: 2 !important; padding: 1rem .75rem; }
    html body .bslib-sidebar-layout.bslib-sidebar-layout > .collapse-toggle { display: none !important; }
    .tab-intro { margin-bottom: .9rem; }
    #plot_track, #state_track, #sim_track, #sim_coords { height: 360px !important; }
    #plot_coords, #plot_steps { height: 300px !important; }
    #data_boxes .bslib-grid { grid-template-columns: repeat(2, minmax(0, 1fr)) !important; gap: .75rem !important; }
    #data_boxes .bslib-grid-item { grid-column: auto !important; }
    #data_boxes .card { margin-bottom: 0; }
    #shiny-notification-panel { left: .75rem; right: .75rem; bottom: .75rem; width: auto; max-width: none; }
    .tab-intro h2 { font-size: 1.15rem; }
    .tab-intro p { font-size: .88rem; }
    .bslib-value-box { min-height: 0 !important; }
    .bslib-value-box .value-box-area { padding: .6rem .8rem !important; }
    .bslib-value-box.vb .value-box-value { font-size: 1.3rem; }
    .bslib-value-box.vb .value-box-title { font-size: .68rem; }
    .param-wrap table td, .param-wrap table th { padding: .35rem .4rem; font-size: .85rem; }
    .navbar-brand { font-size: 1.05rem; }
    .card-header { padding: .7rem .85rem; }
    .card-body { padding: .8rem .85rem; }
  }
  "
  for (nm in names(p)) css <- gsub(paste0("{", nm, "}"), p[[nm]], css, fixed = TRUE)
  css
}

# Header shown at the top of each step.
tab_intro <- function(step, title, text) {
  div(class = "tab-intro",
      div(class = "step", step),
      div(h2(title), p(text)))
}

# Placeholder shown before there are results.
empty_state <- function(icon_name, text) {
  card(card_body(class = "empty-state",
                 div(class = "icon", icon(icon_name)),
                 p(text)))
}

# ggplot theme that matches the app.
theme_smam <- function(base_size = 13) {
  p <- smam_pal
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = p$surface, colour = NA),
      panel.background = ggplot2::element_rect(fill = p$surface, colour = NA),
      panel.grid.major = ggplot2::element_line(colour = p$grid, linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = p$muted, size = ggplot2::rel(0.85)),
      axis.title = ggplot2::element_text(colour = p$text, size = ggplot2::rel(0.9)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      legend.title = ggplot2::element_text(colour = p$text, size = ggplot2::rel(0.85)),
      legend.text = ggplot2::element_text(colour = p$muted, size = ggplot2::rel(0.85)),
      legend.key.size = ggplot2::unit(0.9, "lines"),
      plot.caption = ggplot2::element_text(colour = p$muted, size = ggplot2::rel(0.75), hjust = 0),
      plot.caption.position = "plot",
      plot.margin = ggplot2::margin(10, 14, 8, 8)
    )
}
