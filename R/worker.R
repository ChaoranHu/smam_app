# Code that runs in a background R process (via callr), so long fits do not
# freeze the app. callr sends only the function body to the child process,
# so smam_worker() must be self contained: no calls to other helpers in R/.

smam_worker <- function(task, spec) {
  suppressPackageStartupMessages(library(smam))
  t0 <- Sys.time()

  dat <- spec$data
  seg <- if (isTRUE(spec$use_segment)) "segment" else NULL
  # columns the model sees: time, coordinates and (optionally) the segment
  if (is.null(seg)) dat <- dat[, setdiff(names(dat), "segment"), drop = FALSE]

  pnames <- list(
    MR   = c("lamM", "lamR", "sigma"),
    MRME = c("lamM", "lamR", "sigma", "sig_err"),
    MRH  = c("lamM", "lamR", "lamH", "sigma", "p"),
    MM   = c("lam1", "lam2", "sigma1", "sigma2"),
    BMME = c("sigma", "delta")
  )[[spec$model]]

  as_cov <- function(v) {
    # bootstrap based estimators return list(cov, BS_detail) when detailBS = TRUE
    if (is.list(v) && !is.null(v$cov)) {
      list(vcov = v$cov, boot = v$BS_detail)
    } else {
      list(vcov = as.matrix(v), boot = NULL)
    }
  }

  out <- switch(task,

    fit = {
      o <- spec$options
      fit <- switch(spec$model,
        MR = fitMR(dat, start = spec$start, segment = seg,
                   likelihood = o$likelihood),
        MRME = switch(o$variant,
          exact = fitMRME(dat, start = spec$start, segment = seg,
                          lower = o$lower, upper = o$upper, print_level = 1),
          naive = fitMRME_naive(dat, start = spec$start, segment = seg,
                                lower = o$lower, upper = o$upper),
          approx = fitMRMEapprox(dat, start = spec$start, segment = seg)
        ),
        MRH = fitMRH(dat, start = spec$start, segment = seg,
                     numThreads = o$numThreads,
                     lower = o$lower, upper = o$upper, print_level = 1),
        MM = fitMM(dat, start = spec$start),
        BMME = fitBMME(dat, start = spec$start, segment = seg)
      )
      est <- as.numeric(fit$estimate)
      names(est) <- pnames
      v <- NULL
      if (spec$model == "MR" && !all(is.na(fit$varest))) v <- fit$varest
      if (spec$model == "BMME") {
        # the segmented fit names this element varest, the plain fit var.est
        vb <- if (!is.null(fit$var.est)) fit$var.est else fit$varest
        if (!is.null(vb) && !all(is.na(vb))) v <- vb
      }
      list(estimate = est, loglik = fit$loglik, convergence = fit$convergence,
           vcov = v, boot = NULL)
    },

    se = {
      est <- unname(spec$estimate)
      threads <- spec$numThreads
      v <- switch(spec$model,
        MRME = if (spec$variant == "naive") {
          if (spec$method == "Godambe") {
            estVarMRMEnaive_Godambe(est, dat, nBS = spec$nBS, numThreads = threads)
          } else {
            estVarMRMEnaive_pBootstrap(est, dat, nBS = spec$nBS, detailBS = TRUE,
                                       numThreads = threads)
          }
        } else {
          if (spec$method == "Godambe") {
            estVarMRME_Godambe(est, dat, nBS = spec$nBS, numThreads = threads)
          } else {
            estVarMRME_pBootstrap(est, dat, nBS = spec$nBS, detailBS = TRUE,
                                  numThreads = threads)
          }
        },
        MM = estVarMM(est, dat, nBS = spec$nBS, detailBS = TRUE,
                      numThreads = threads),
        MRH = vcov(structure(list(estimate = est, data = dat), class = "smam_mrh"),
                   numThreads = threads)
      )
      as_cov(v)
    },

    state = {
      theta <- unname(spec$theta)
      res <- switch(paste(spec$model, spec$method),
        "MR fb"       = fitStateMR(dat, theta, cutoff = spec$cutoff),
        "MR viterbi"  = fitViterbiMR(dat, theta, cutoff = spec$cutoff),
        "MRH fb"      = fitStateMRH(dat, theta),
        "MRH viterbi" = fitViterbiMRH(dat, theta)
      )
      list(states = res)
    }
  )

  out$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}
