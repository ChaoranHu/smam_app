# Example data sets and helpers that turn a raw table into the format smam
# expects: first column time, then location coordinates, optionally followed
# by a segment column (0 = not used, other integers = segment id).

f109_prepared <- function() {
  d <- smam::f109raw
  d$dE <- d$dE / 1000
  d$dN <- d$dN / 1000
  tt <- d$t %% 1
  d$season <- ifelse(tt > 5 / 12 & tt <= 8 / 12, "summer",
                     ifelse(tt > 11 / 12 | tt <= 2 / 12, "winter", "other"))
  d$hr <- (d$t - 2009) * 365.25 * 24
  # year starting in March, so one winter is not split across two years
  d$year <- (d$t - 2 / 12) %/% 1
  d
}

example_data <- function(id) {
  switch(id,
    f109 = {
      d <- smam::f109
      d$date <- as.character(d$date)
      d
    },
    summer2012 = {
      d <- f109_prepared()
      d <- d[d$season == "summer" & d$year == 2012, c("hr", "dE", "dN")]
      rownames(d) <- NULL
      d
    },
    winter2011 = {
      d <- f109_prepared()
      d <- d[d$season == "winter" & d$year == 2011, c("hr", "dE", "dN")]
      rownames(d) <- NULL
      d
    },
    summers = {
      d <- f109_prepared()
      d$summer_year <- ifelse(d$season == "summer", d$year, 0)
      d <- d[, c("hr", "dE", "dN", "summer_year")]
      d
    },
    simMR = {
      set.seed(2024)
      tgrid <- sort(sample(seq(0, 100, length.out = 2000), 150))
      d <- smam::rMR(tgrid, lamM = 1, lamR = 0.5, sigma = 1, s0 = "m")
      names(d) <- c("time", "east", "north")
      d
    }
  )
}

example_choices <- c(
  "Mountain lion f109, summer 2012 (MRME paper example)" = "summer2012",
  "Mountain lion f109, winter 2011" = "winter2011",
  "Mountain lion f109, all summers pooled (segment = year)" = "summers",
  "Mountain lion f109, full record 2009 to 2012 (has a date column)" = "f109",
  "Simulated moving-resting track (150 points)" = "simMR"
)

# Suggested column mapping for each example, so the Data tab is ready to use.
example_mapping <- list(
  summer2012 = list(time = "hr", coords = c("dE", "dN"), time_unit = "hours",
                    coord_unit = "km", segment = "none"),
  winter2011 = list(time = "hr", coords = c("dE", "dN"), time_unit = "hours",
                    coord_unit = "km", segment = "none"),
  summers = list(time = "hr", coords = c("dE", "dN"), time_unit = "hours",
                 coord_unit = "km", segment = "column", segment_col = "summer_year"),
  f109 = list(time = "cumTime", coords = c("centerE", "centerN"), time_unit = "hours",
              coord_unit = "km", segment = "season", date_col = "date"),
  simMR = list(time = "time", coords = c("east", "north"), time_unit = "hours",
               coord_unit = "km", segment = "none")
)

time_factors <- c(hours = 1, minutes = 1 / 60, seconds = 1 / 3600, days = 24,
                  "decimal years" = 365.25 * 24)
coord_factors <- c(km = 1, m = 1 / 1000)

guess_time_col <- function(df) {
  nm <- names(df)
  hit <- grep("^(time|t|hr|hour|hours|cumtime|datetime|timestamp)$", nm, ignore.case = TRUE)
  if (length(hit)) return(nm[hit[1]])
  num <- nm[vapply(df, is.numeric, TRUE)]
  if (length(num)) num[1] else nm[1]
}

guess_coord_cols <- function(df, time_col) {
  num <- setdiff(names(df)[vapply(df, is.numeric, TRUE)], time_col)
  pref <- grep("east|north|^x|^y|^de$|^dn$|lon|lat|center", num, ignore.case = TRUE, value = TRUE)
  if (length(pref) >= 2) return(pref[1:2])
  utils::head(num, 2)
}

# Build segment ids from a date column: one segment per year for the dates
# between start (MM-DD) and end (MM-DD), using smam::seasonFilter.
season_segments <- function(dates, start, end) {
  ok <- !is.na(dates)
  tmp <- data.frame(date = dates[ok], row = which(ok))
  res <- smam::seasonFilter(tmp, start, end)
  seg <- integer(length(dates))
  seg[res$row] <- as.integer(res$BATCH)
  seg
}

# Returns list(data = prepared data.frame, notes = character messages).
prepare_data <- function(raw, time_col, coord_cols, time_unit = "hours",
                         coord_unit = "km", shift_time = TRUE, center = TRUE,
                         segment = NULL) {
  notes <- character()
  if (is.null(raw) || !nrow(raw)) stop("No data loaded.")
  if (!length(time_col) || !time_col %in% names(raw)) stop("Choose a time column.")
  if (!length(coord_cols)) stop("Choose at least one coordinate column.")
  if (time_col %in% coord_cols) stop("The time column cannot also be a coordinate column.")

  to_num <- function(x, nm) {
    v <- suppressWarnings(as.numeric(as.character(x)))
    if (all(is.na(v))) stop("Column '", nm, "' is not numeric.")
    v
  }
  time <- to_num(raw[[time_col]], time_col) * time_factors[[time_unit]]
  coords <- lapply(coord_cols, function(nm) to_num(raw[[nm]], nm) * coord_factors[[coord_unit]])
  names(coords) <- make.names(coord_cols, unique = TRUE)
  out <- data.frame(time = time, coords, check.names = FALSE)
  if (!is.null(segment)) out$segment <- as.integer(segment)

  keep <- stats::complete.cases(out)
  if (any(!keep)) notes <- c(notes, paste(sum(!keep), "rows with missing values were removed."))
  out <- out[keep, , drop = FALSE]
  if (is.unsorted(out$time)) {
    out <- out[order(out$time), , drop = FALSE]
    notes <- c(notes, "Rows were sorted by time.")
  }
  dup <- duplicated(out$time)
  if (any(dup)) {
    out <- out[!dup, , drop = FALSE]
    notes <- c(notes, paste(sum(dup), "rows with a repeated time were removed."))
  }
  if (nrow(out) < 3) stop("At least 3 valid observations are needed.")
  if (shift_time) out$time <- out$time - out$time[1]
  if (center) {
    for (nm in names(coords)) out[[nm]] <- out[[nm]] - out[[nm]][1]
  }
  rownames(out) <- NULL
  list(data = out, notes = notes)
}

coord_names <- function(d) setdiff(names(d), c("time", "segment"))
