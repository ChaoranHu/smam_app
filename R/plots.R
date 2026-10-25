# ggplot2 figures shared by the Data, States and Simulate tabs.

# Plot colours come from the active theme (R/theme.R).
state_colours <- function() c(active_theme()$states, "not used" = "#C9CEC5")
series_colours <- function() active_theme()$series
time_gradient <- function() active_theme()$gradient

plot_track <- function(d, colour_by = c("time", "segment", "state"), states = NULL,
                       state_title = "Estimated state") {
  colour_by <- match.arg(colour_by)
  cn <- coord_names(d)
  if (length(cn) < 2) {
    return(plot_coords(d))
  }
  df <- data.frame(x = d[[cn[1]]], y = d[[cn[2]]], time = d$time)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_path(colour = "#C9CEC5", linewidth = 0.35)
  if (colour_by == "segment" && "segment" %in% names(d)) {
    df$segment <- ifelse(d$segment == 0, "not used", paste("segment", d$segment))
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = segment), size = 1.2) +
      ggplot2::scale_colour_manual(
        values = c("not used" = "#DADDD5",
                   stats::setNames(scales_hue(length(unique(d$segment[d$segment != 0]))),
                                   paste("segment", sort(unique(d$segment[d$segment != 0]))))),
        name = NULL)
  } else if (colour_by == "state" && !is.null(states)) {
    df$state <- states
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = state), size = 1.4) +
      ggplot2::scale_colour_manual(values = state_colours(), name = state_title)
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = time), size = 1.2) +
      ggplot2::scale_colour_gradientn(colours = time_gradient(), name = "time (hr)")
  }
  p + ggplot2::geom_point(data = df[1, ], shape = 24, size = 3.2, fill = "#FFFFFF",
                        colour = active_theme()$pal$ink, stroke = 0.9) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = paste(cn[1], "(km)"), y = paste(cn[2], "(km)"),
                  caption = "Triangle marks the first location") +
    theme_smam()
}

scales_hue <- function(n) {
  if (n <= 0) return(character())
  if (n <= length(series_colours())) return(series_colours()[seq_len(n)])
  grDevices::hcl(h = seq(15, 375, length.out = n + 1)[seq_len(n)], c = 70, l = 50)
}

plot_coords <- function(d) {
  cn <- coord_names(d)
  long <- do.call(rbind, lapply(cn, function(nm) {
    data.frame(time = d$time, value = d[[nm]], coordinate = nm)
  }))
  ggplot2::ggplot(long, ggplot2::aes(time, value, colour = coordinate)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = scales_hue(length(cn))) +
    ggplot2::labs(x = "time (hr)", y = "location (km)", colour = NULL) +
    theme_smam() +
    ggplot2::theme(legend.position = "top", legend.justification = "left")
}

plot_steps <- function(d) {
  dt <- diff(d$time)
  ggplot2::ggplot(data.frame(dt = dt), ggplot2::aes(dt)) +
    ggplot2::geom_histogram(bins = 40, fill = active_theme()$pal$primary, colour = "#FFFFFF", linewidth = 0.3) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "time between consecutive fixes (hr, log scale)", y = "count") +
    theme_smam()
}

plot_state_prob <- function(res, labels) {
  pcols <- intersect(c("p.m", "p.r", "p.h"), names(res))
  probs <- as.matrix(res[, pcols, drop = FALSE])
  probs <- probs / rowSums(probs)
  long <- do.call(rbind, lapply(seq_along(pcols), function(i) {
    data.frame(time = res$time, prob = probs[, i], state = labels[i])
  }))
  long$state <- factor(long$state, levels = labels)
  ggplot2::ggplot(long, ggplot2::aes(time, prob, fill = state)) +
    ggplot2::geom_area(position = "stack", alpha = 0.9, colour = "#FFFFFF", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = state_colours(), name = NULL) +
    ggplot2::labs(x = "time (hr)", y = "posterior probability of state") +
    theme_smam() +
    ggplot2::theme(legend.position = "top", legend.justification = "left")
}

plot_state_path <- function(res, states) {
  df <- data.frame(time = res$time, state = factor(states, levels = names(state_colours())))
  df <- droplevels(df)
  ggplot2::ggplot(df, ggplot2::aes(time, state, colour = state)) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = state_colours(), guide = "none") +
    ggplot2::labs(x = "time (hr)", y = NULL) +
    theme_smam()
}
