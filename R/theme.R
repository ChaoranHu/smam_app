# Look and feel: colours, fonts, component styles and the matching ggplot theme.
# One palette is used everywhere so the plots, tables and controls read as one design.

# Each theme is a palette, a pair of fonts and plot colours. The active one is
# chosen with the SMAM_THEME environment variable, falling back to the default.

topo_pattern <- function(colour) {
  # faint contour lines, like a topographic map, drawn as an inline SVG
  svg <- sprintf(paste0(
    "<svg xmlns='http://www.w3.org/2000/svg' width='360' height='120' viewBox='0 0 360 120'>",
    "<g fill='none' stroke='%s' stroke-width='1'>",
    "<path d='M-10 30 C 40 5, 90 55, 150 28 S 260 0, 370 35'/>",
    "<path d='M-10 52 C 50 25, 100 78, 160 50 S 270 22, 370 58'/>",
    "<path d='M-10 74 C 45 50, 110 100, 170 72 S 280 46, 370 80'/>",
    "<path d='M-10 96 C 55 74, 115 120, 180 94 S 285 70, 370 102'/>",
    "<ellipse cx='255' cy='30' rx='28' ry='10'/><ellipse cx='255' cy='30' rx='14' ry='5'/>",
    "</g></svg>"), colour)
  paste0("url(\"data:image/svg+xml,", utils::URLencode(svg, reserved = TRUE), "\")")
}

smam_themes <- list(

  # the design from October 2026: teal with an ochre accent
  teal = list(
    label = "Teal", brand_icon = "route",
    pal = list(ink = "#17302B", text = "#24312F", muted = "#64716D", page = "#F4F5F0",
               surface = "#FFFFFF", line = "#E2E5DD", grid = "#ECEEE8", primary = "#1F6F5F",
               tint = "#E6F0EC", accent = "#D9822B", blue = "#2F6690", olive = "#7A9A3A",
               danger = "#B44A3C", soft = "#F7F8F5"),
    base_font = "Inter", heading_font = "Inter", heading_weight = 600, heading_case = "none",
    fonts_href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap",
    radius = 0.9, pattern = FALSE,
    states = c(moving = "#D9822B", resting = "#2F6690", handling = "#7A9A3A"),
    series = c("#1F6F5F", "#D9822B", "#2F6690", "#8E5572", "#7A9A3A", "#B8860B", "#4F8FA8", "#C0574A"),
    gradient = c("#E3A04F", "#9DB36B", "#3E9078", "#1F6F5F", "#17302B")
  ),

  # A. a field notebook: moss green, parchment and rust, serif headings
  journal = list(
    label = "Field journal", brand_icon = "feather-pointed",
    pal = list(ink = "#26362A", text = "#2B2A26", muted = "#6F6A5D", page = "#F5F1E7",
               surface = "#FFFDF8", line = "#E3DCCB", grid = "#EEE8DA", primary = "#3F6B4A",
               tint = "#E8EEE2", accent = "#A4552B", blue = "#3D6A8A", olive = "#8A8F3C",
               danger = "#9E3B2E", soft = "#F8F4EA"),
    base_font = "Lato", heading_font = "Alegreya", heading_weight = 700, heading_case = "none",
    fonts_href = "https://fonts.googleapis.com/css2?family=Alegreya:wght@500;700&family=Lato:wght@400;700&display=swap",
    radius = 0.5, pattern = TRUE,
    states = c(moving = "#A4552B", resting = "#3D6A8A", handling = "#8A8F3C"),
    series = c("#3F6B4A", "#A4552B", "#3D6A8A", "#8A8F3C", "#7B4B6A", "#C08A3E", "#5E8C8C", "#6B5B45"),
    gradient = c("#D9A45B", "#B3A866", "#7C9A62", "#3F6B4A", "#22331F")
  ),

  # B. tracking software: navy, signal orange and blue, crisp sans serif
  telemetry = list(
    label = "Telemetry", brand_icon = "satellite-dish",
    pal = list(ink = "#0F1C2B", text = "#1C2733", muted = "#5D6A79", page = "#F1F3F6",
               surface = "#FFFFFF", line = "#DCE1E7", grid = "#EAEDF1", primary = "#1D6FA5",
               tint = "#E3EDF6", accent = "#E8743B", blue = "#1D6FA5", olive = "#2FA37F",
               danger = "#C8463D", soft = "#F5F7F9"),
    base_font = "Inter", heading_font = "Inter", heading_weight = 700, heading_case = "none",
    fonts_href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@500&display=swap",
    radius = 0.4, pattern = FALSE, mono_numbers = TRUE,
    states = c(moving = "#E8743B", resting = "#1D6FA5", handling = "#2FA37F"),
    series = c("#1D6FA5", "#E8743B", "#2FA37F", "#7A4FB0", "#D4A11E", "#C8463D", "#4BA3C7", "#5D6A79"),
    gradient = c("#F2C14E", "#E8743B", "#B8457B", "#5A3F9E", "#1B2F6B")
  ),

  # C. park service signage: brown, forest green and gold, condensed capitals
  park = list(
    label = "National park", brand_icon = "paw", brand_image = "smam-icon.svg",
    pal = list(ink = "#3A2F24", text = "#2E2A22", muted = "#6E6656", page = "#F3F0E6",
               surface = "#FFFFFF", line = "#E0D9C6", grid = "#ECE6D5", primary = "#4E6B2F",
               tint = "#E9EEDD", accent = "#C98A2E", blue = "#4A6D8C", olive = "#4E6B2F",
               danger = "#A8442F", soft = "#F7F4EC"),
    base_font = "Roboto", heading_font = "Oswald", heading_weight = 500, heading_case = "uppercase",
    fonts_href = "https://fonts.googleapis.com/css2?family=Oswald:wght@400;500;600&family=Roboto:wght@400;500;700&display=swap",
    radius = 0.3, pattern = TRUE,
    states = c(moving = "#B5542E", resting = "#4A6D8C", handling = "#4E6B2F"),
    series = c("#4E6B2F", "#B5542E", "#4A6D8C", "#C98A2E", "#6E4A6B", "#7A8B3A", "#3E7C7A", "#6E6656"),
    gradient = c("#E3B65A", "#C98A2E", "#8C7A3A", "#4E6B2F", "#2B3A1E")
  )
)

active_theme_id <- function() {
  id <- Sys.getenv("SMAM_THEME", "park")
  if (id %in% names(smam_themes)) id else "park"
}
active_theme <- function() smam_themes[[active_theme_id()]]

app_theme <- function() {
  th <- active_theme()
  pal <- th$pal
  fonts <- font_link(th$base_font, href = th$fonts_href)
  stack <- c("system-ui", "-apple-system", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif")
  heading_stack <- if (th$heading_font == "Alegreya") c("Georgia", "serif") else stack
  bs_theme(
    version = 5,
    bg = pal$page, fg = pal$text,
    primary = pal$primary, secondary = pal$muted,
    success = "#2E7D4F", info = pal$blue, warning = pal$accent,
    danger = pal$danger,
    base_font = do.call(font_collection, c(list(fonts, th$base_font), as.list(stack))),
    heading_font = do.call(font_collection, c(list(th$heading_font), as.list(heading_stack))),
    code_font = font_collection("JetBrains Mono", "ui-monospace", "SFMono-Regular", "Menlo", "Consolas", "monospace"),
    "font-size-base" = "0.95rem",
    "headings-font-weight" = th$heading_weight,
    "border-radius" = paste0(th$radius * 0.67, "rem"),
    "border-radius-sm" = paste0(th$radius * 0.5, "rem"),
    "border-radius-lg" = paste0(th$radius, "rem"),
    "border-color" = pal$line,
    "navbar-bg" = pal$ink,
    "input-bg" = pal$surface,
    "link-color" = pal$primary
  )
}

app_css <- function() {
  th <- active_theme()
  p <- c(th$pal, radius = paste0(th$radius, "rem"),
         heading = paste0("'", th$heading_font, "', ", if (th$heading_font == "Alegreya") "Georgia, serif" else "system-ui, sans-serif"),
         hcase = th$heading_case, hweight = th$heading_weight,
         pattern = if (isTRUE(th$pattern)) topo_pattern("rgba(255,255,255,0.07)") else "none")
  css <- "
  body { background: {page}; letter-spacing: -0.005em; }

  /* navbar: the selected tab is plain white text, the others are dimmed */
  .navbar.navbar-static-top, nav.navbar { background-color: {ink} !important; background-image: {pattern};
     background-size: 360px 120px; }
  .navbar .navbar-brand, .navbar .nav-link { color: rgba(255,255,255,.68) !important; }
  .navbar .navbar-brand, .navbar .nav-link.active, .navbar .nav-link:hover { color: #FFFFFF !important; }
  .navbar .nav-link.active { font-weight: 600; }
  .navbar .nav-link, .navbar .nav-link.active, .navbar .nav-link:focus { background: none !important;
     border: 0 !important; box-shadow: none !important; }
  .navbar .navbar-toggler { border-color: rgba(255,255,255,.3); }
  .navbar .navbar-toggler-icon { filter: invert(1) brightness(2); }
  .navbar { box-shadow: 0 1px 0 rgba(0,0,0,.08), 0 6px 18px rgba(0,0,0,.10); }
  .navbar-brand { display: flex; align-items: center; gap: .55rem; font-weight: 700; letter-spacing: -0.01em;
     font-family: {heading}; text-transform: {hcase}; }
  .brand-mark { display: inline-flex; align-items: center; justify-content: center; width: 2rem; height: 2rem;
                border-radius: calc({radius} * .6); background: {accent}; color: #fff; font-size: .95rem; overflow: hidden; }
  .brand-mark img { width: 100%; height: 100%; display: block; }
  .brand-sub { font-weight: 400; opacity: .7; font-size: .9rem; font-family: var(--bs-body-font-family); text-transform: none; }
  .navbar .nav-link { font-weight: 500; padding: .45rem .8rem !important; }
  h1, h2, h3, h4, h5, .card-header { font-family: {heading}; }
  .tab-intro h2, .guide h3, .guide h4 { text-transform: {hcase}; letter-spacing: 0.01em; }

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
          border-radius: {radius}; box-shadow: 0 1px 2px rgba(0,0,0,.04), 0 2px 8px rgba(0,0,0,.05);
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
  .btn-primary { box-shadow: 0 1px 2px rgba(0,0,0,.2); }
  .form-control, .form-select, .selectize-input { border-color: #D6DAD1; }
  .form-control:focus, .form-select:focus, .selectize-input.focus { border-color: {primary};
     box-shadow: 0 0 0 .2rem {tint}; }
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
  .small-log { font-size: .78rem; max-height: 260px; overflow-y: auto; margin: 0; background: {soft};
               padding: .75rem 1rem; border-radius: 0 0 {radius} {radius}; }
  .alert { border-radius: .7rem; border-width: 0; border-left: 4px solid; }
  .alert-warning { background: #FBF1E3; border-color: {accent}; color: #6B4310; }
  .alert-info { background: #EAF1F6; border-color: {blue}; color: #1E3F58; }
  .alert-light { background: {tint}; border-color: {primary}; color: {ink}; }
  .empty-state { text-align: center; padding: 2.5rem 1rem; color: {muted}; }
  .empty-state .icon { font-size: 1.6rem; color: {primary}; background: {tint}; width: 3.5rem; height: 3.5rem; align-self: center;
                       border-radius: 50%; display: inline-flex; align-items: center; justify-content: center;
                       margin-bottom: .75rem; }
  .empty-state p { margin: 0 auto; max-width: 28rem; }
  pre { white-space: pre-wrap; background: {soft}; border: 1px solid {line}; border-radius: .6rem;
        padding: .8rem 1rem; font-size: .82rem; color: {ink}; }
  .card > pre { margin: 1rem 1.1rem; }

  .shiny-notification { background: {surface}; color: {text}; border: 1px solid {line}; border-left: 4px solid {primary};
                        border-radius: .7rem; box-shadow: 0 8px 24px rgba(0,0,0,.12); opacity: 1; }
  .shiny-notification-warning { border-left-color: {accent}; }
  .shiny-notification-error { border-left-color: {danger}; }

  /* guide */
  .guide { line-height: 1.6; }
  .guide h3 { color: {ink}; }
  .guide h4 { color: {ink}; margin-top: 2rem; padding-top: 1rem; border-top: 1px solid {line}; }
  .guide .model-block { background: {surface}; border: 1px solid {line}; border-radius: {radius};
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
  if (isTRUE(th$mono_numbers)) css <- paste(css, "
  .bslib-value-box.vb .value-box-value, #est_table td:nth-child(3), table.dataTable td, .param-wrap td { font-family: 'JetBrains Mono', ui-monospace, monospace; font-size: .92em; }
  .bslib-value-box.vb .value-box-value { font-size: 1.45rem; }")
  for (nm in names(p)) css <- gsub(paste0("{", nm, "}"), p[[nm]], css, fixed = TRUE)
  css
}

# Mountain lion (cougar) tracks: two prints walking up to the right. Each
# print has four teardrop toes without claw marks over a heel pad with two
# lobes on top and three below. Drawn in the current text colour at icon
# size, like icon().
cougar_track_icon <- function() {
  print <- paste0(
    "<ellipse cx='21' cy='37' rx='8' ry='10.5' transform='rotate(-28 21 37)'/>",
    "<ellipse cx='39' cy='20' rx='8.5' ry='11' transform='rotate(-10 39 20)'/>",
    "<ellipse cx='61' cy='17' rx='8.5' ry='11' transform='rotate(8 61 17)'/>",
    "<ellipse cx='79' cy='32' rx='8' ry='10.5' transform='rotate(28 79 32)'/>",
    "<path d='M50 55C45 49 37 49 33 53C26 60 19 70 19 80C19 88 25 92 31 91",
    "C35 90 37 89 39 89C43 92 46 93 50 93C54 93 57 92 61 89C63 89 65 90 69 91",
    "C75 92 81 88 81 80C81 70 74 60 67 53C63 49 55 49 50 55Z'/>")
  HTML(paste0(
    "<svg class='cougar-track' viewBox='0 0 104 100' width='1.15em' height='1.1em' ",
    "fill='currentColor' aria-hidden='true' role='img' style='vertical-align:-0.18em'>",
    "<g transform='translate(2 42) scale(.56)'>", print, "</g>",
    "<g transform='translate(46 2) scale(.56)'>", print, "</g>",
    "</svg>"))
}

# Mark at the left of the navbar: the theme's own image (in www/) if it has
# one, otherwise a Font Awesome icon.
brand_mark <- function() {
  th <- active_theme()
  span(class = "brand-mark",
       if (!is.null(th$brand_image)) tags$img(src = th$brand_image, alt = "") else icon(th$brand_icon))
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
                 div(class = "icon", if (inherits(icon_name, "html")) icon_name else icon(icon_name)),
                 p(text)))
}

# ggplot theme that matches the app.
theme_smam <- function(base_size = 13) {
  p <- active_theme()$pal
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
