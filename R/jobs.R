# A small background job manager built on callr. Each job runs smam_worker()
# in its own R process, writes the optimizer output to a log file that the UI
# shows while it runs, and can be cancelled.

start_job <- function(task, spec, label) {
  log_file <- tempfile(fileext = ".log")
  proc <- callr::r_bg(
    smam_worker,
    args = list(task = task, spec = spec),
    stdout = log_file, stderr = "2>&1",
    supervise = TRUE
  )
  list(proc = proc, log = log_file, task = task, label = label,
       started = Sys.time(), spec = spec)
}

job_running <- function(job) {
  !is.null(job) && job$proc$is_alive()
}

job_elapsed <- function(job) {
  as.numeric(difftime(Sys.time(), job$started, units = "secs"))
}

# Last few lines of the job log; long numeric vectors printed by nloptr are
# trimmed so the panel stays readable.
job_log_tail <- function(job, n = 12) {
  if (is.null(job) || !file.exists(job$log)) return("")
  lines <- tryCatch(readLines(job$log, warn = FALSE), error = function(e) character())
  lines <- lines[nzchar(trimws(lines))]
  lines <- substr(lines, 1, 160)
  paste(utils::tail(lines, n), collapse = "\n")
}

# Returns list(ok = TRUE, value = ...) or list(ok = FALSE, error = "message").
job_result <- function(job) {
  tryCatch(
    list(ok = TRUE, value = job$proc$get_result()),
    error = function(e) {
      msg <- conditionMessage(e)
      parent <- e$parent
      if (!is.null(parent)) msg <- conditionMessage(parent)
      list(ok = FALSE, error = msg)
    }
  )
}

cancel_job <- function(job) {
  if (!is.null(job) && job$proc$is_alive()) job$proc$kill()
  invisible(NULL)
}

format_seconds <- function(s) {
  s <- round(s)
  if (s < 60) return(paste0(s, " s"))
  sprintf("%d min %02d s", s %/% 60, s %% 60)
}
