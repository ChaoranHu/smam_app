test_that("every example loads and prepares with its suggested mapping", {
  for (id in example_choices) {
    raw <- example_data(id)
    expect_s3_class(raw, "data.frame")
    expect_gt(nrow(raw), 10)
    map <- example_mapping[[id]]
    expect_true(all(c(map$time, map$coords) %in% names(raw)), info = id)
    seg <- switch(map$segment,
      none = NULL,
      column = as.integer(raw[[map$segment_col]]),
      season = season_segments(as.Date(raw[[map$date_col]]), "06-01", "08-31")
    )
    res <- prepare_data(raw, map$time, map$coords, map$time_unit, map$coord_unit,
                        TRUE, TRUE, seg)
    d <- res$data
    expect_equal(names(d)[1], "time", info = id)
    expect_false(is.unsorted(d$time), info = id)
    expect_false(anyDuplicated(d$time) > 0, info = id)
    expect_equal(d$time[1], 0, info = id)
    expect_equal(unname(unlist(d[1, coord_names(d)])), rep(0, length(map$coords)), info = id)
    if (!is.null(seg)) expect_true(any(d$segment != 0), info = id)
  }
})

test_that("the f109 summer 2012 example matches the package data", {
  d <- prepare_data(example_data("summer2012"), "hr", c("dE", "dN"))$data
  # about three months of fixes, in km
  expect_gt(diff(range(d$time)), 24 * 60)
  expect_lt(max(abs(d$dE)), 100)
})

test_that("units are converted to hours and kilometres", {
  raw <- data.frame(t = c(0, 60, 120, 240), x = c(0, 1000, 2000, 500), y = c(0, 0, 500, 500))
  d <- prepare_data(raw, "t", c("x", "y"), time_unit = "minutes", coord_unit = "m")$data
  expect_equal(d$time, c(0, 1, 2, 4))
  expect_equal(d$x, c(0, 1, 2, 0.5))
  d2 <- prepare_data(raw, "t", c("x", "y"), time_unit = "days", shift_time = FALSE, center = FALSE)$data
  expect_equal(d2$time, raw$t * 24)
})

test_that("messy input is cleaned with notes", {
  raw <- data.frame(t = c(3, 1, 2, 2, NA, 5), x = c(1, 2, 3, 4, 5, NA), y = 1:6)
  res <- prepare_data(raw, "t", c("x", "y"))
  expect_equal(res$data$time, c(0, 1, 2))
  expect_true(any(grepl("missing", res$notes)))
  expect_true(any(grepl("sorted", res$notes)))
  expect_true(any(grepl("repeated", res$notes)))
})

test_that("bad input gives a clear error", {
  raw <- data.frame(t = 1:5, x = 1:5, y = letters[1:5])
  expect_error(prepare_data(raw, "t", character()), "coordinate")
  expect_error(prepare_data(raw, "t", c("t", "x")), "cannot also")
  expect_error(prepare_data(raw, "t", c("x", "y")), "not numeric")
  expect_error(prepare_data(raw[1:2, ], "t", "x"), "At least 3")
  expect_error(prepare_data(NULL, "t", "x"), "No data")
})

test_that("the CSV in the repository reads and maps automatically", {
  raw <- utils::read.table(file.path(root, "test.csv"), header = TRUE, sep = ",")
  tc <- guess_time_col(raw)
  expect_equal(tc, "time")
  expect_equal(guess_coord_cols(raw, tc), c("east", "north"))
  expect_equal(nrow(prepare_data(raw, tc, c("east", "north"))$data), nrow(raw))
})

test_that("season windows become one segment per year, also across the new year", {
  dates <- seq(as.Date("2010-01-01"), as.Date("2012-12-31"), by = "week")
  s <- season_segments(dates, "06-01", "08-31")
  expect_equal(length(unique(s[s != 0])), 3)
  expect_true(all(format(dates[s != 0], "%m") %in% c("06", "07", "08")))
  w <- season_segments(dates, "12-01", "02-28")
  expect_true(all(format(dates[w != 0], "%m") %in% c("12", "01", "02")))
  expect_gte(length(unique(w[w != 0])), 3)
})
