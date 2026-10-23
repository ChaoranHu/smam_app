# ggplot2 figures shared by the Data, States and Simulate tabs.

state_colours <- c(moving = "#D55E00", resting = "#0072B2", handling = "#009E73",
                   "not used" = "grey75")

plot_track <- function(d, colour_by = c("time", "segment", "state"), states = NULL) {
  colour_by <- match.arg(colour_by)
  cn <- coord_names(d)
  if (length(cn) < 2) {
    return(plot_coords(d))
  }
  df <- data.frame(x = d[[cn[1]]], y = d[[cn[2]]], time = d$time)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_path(colour = "grey70", linewidth = 0.3)
  if (colour_by == "segment" && "segment" %in% names(d)) {
    df$segment <- ifelse(d$segment == 0, "not used", paste("segment", d$segment))
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = segment), size = 1.2) +
      ggplot2::scale_colour_manual(
        values = c("not used" = "grey80",
                   stats::setNames(scales_hue(length(unique(d$segment[d$segment != 0]))),
                                   paste("segment", sort(unique(d$segment[d$segment != 0]))))),
        name = NULL)
  } else if (colour_by == "state" && !is.null(states)) {
    df$state <- states
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = state), size = 1.4) +
      ggplot2::scale_colour_manual(values = state_colours, name = "Estimated state")
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = time), size = 1.2) +
      ggplot2::scale_colour_viridis_c(name = "time (hr)")
  }
  p + ggplot2::geom_point(data = df[1, ], shape = 24, size = 3, fill = "white", colour = "black") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = paste(cn[1], "(km)"), y = paste(cn[2], "(km)"),
                  caption = "Triangle marks the first location") +
    ggplot2::theme_bw(base_size = 13)
}

scales_hue <- function(n) {
  if (n <= 0) return(character())
  grDevices::hcl(h = seq(15, 375, length.out = n + 1)[seq_len(n)], c = 100, l = 55)
}

plot_coords <- function(d) {
  cn <- coord_names(d)
  long <- do.call(rbind, lapply(cn, function(nm) {
    data.frame(time = d$time, value = d[[nm]], coordinate = nm)
  }))
  ggplot2::ggplot(long, ggplot2::aes(time, value, colour = coordinate)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::labs(x = "time (hr)", y = "location (km)", colour = NULL) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "top")
}

plot_steps <- function(d) {
  dt <- diff(d$time)
  ggplot2::ggplot(data.frame(dt = dt), ggplot2::aes(dt)) +
    ggplot2::geom_histogram(bins = 40, fill = "#0072B2", colour = "white") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "time between consecutive fixes (hr, log scale)", y = "count") +
    ggplot2::theme_bw(base_size = 13)
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
    ggplot2::geom_area(position = "stack", alpha = 0.85) +
    ggplot2::scale_fill_manual(values = state_colours, name = NULL) +
    ggplot2::labs(x = "time (hr)", y = "posterior probability of state") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "top")
}

plot_state_path <- function(res, states) {
  df <- data.frame(time = res$time, state = factor(states, levels = names(state_colours)))
  df <- droplevels(df)
  ggplot2::ggplot(df, ggplot2::aes(time, state, colour = state)) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = state_colours, guide = "none") +
    ggplot2::labs(x = "time (hr)", y = NULL) +
    ggplot2::theme_bw(base_size = 13)
}
