# These call smam_worker() directly, which is what the background process runs.

expect_fit <- function(res, model) {
  k <- length(smam_models[[model]]$params)
  expect_length(res$estimate, k)
  expect_true(all(is.finite(res$estimate)), info = model)
  expect_true(all(res$estimate > 0), info = model)
  expect_true(is.finite(res$loglik), info = model)
  expect_false(is.null(res$convergence), info = model)
  expect_gte(res$elapsed, 0)
  expect_s3_class(derived_quantities(model, res$estimate), "data.frame")
}

test_that("MR fits with full and composite likelihood, with and without segments", {
  d <- sim_track("MR", n = 60, Tend = 60)
  full <- smam_worker("fit", fit_spec("MR", d, list(likelihood = "full")))
  expect_fit(full, "MR")
  expect_true(convergence_ok("MR", "", full$convergence))
  expect_equal(dim(full$vcov), c(3, 3))
  comp <- smam_worker("fit", fit_spec("MR", d, list(likelihood = "composite")))
  expect_fit(comp, "MR")
  d$segment <- rep(c(1, 0, 2), each = 20)
  seg <- smam_worker("fit", fit_spec("MR", d, list(likelihood = "full"), use_segment = TRUE))
  expect_fit(seg, "MR")
  expect_false(isTRUE(all.equal(seg$estimate, full$estimate)))
})

test_that("BMME fits with a moment start, plain and segmented, with a Hessian", {
  d <- sim_track("BMME", n = 80, Tend = 80)
  res <- smam_worker("fit", fit_spec("BMME", d))
  expect_fit(res, "BMME")
  expect_equal(dim(res$vcov), c(2, 2))
  d$segment <- rep(1:2, each = 40)
  seg <- smam_worker("fit", fit_spec("BMME", d, use_segment = TRUE))
  expect_fit(seg, "BMME")
  expect_equal(dim(seg$vcov), c(2, 2))
})

test_that("MRME approximation fits", {
  # fitMRMEapprox takes about 2 s per likelihood evaluation on 10 fixes, so a full
  # fit runs for many minutes; run it only on request
  if (!identical(Sys.getenv("SMAM_APPROX_TEST"), "true")) skip("set SMAM_APPROX_TEST=true to run")
  d <- sim_track("MRME", n = 6, Tend = 6)
  res <- smam_worker("fit", fit_spec("MRME", d, list(variant = "approx")))
  expect_fit(res, "MRME")
})

test_that("MRME exact and naive composite likelihood fit", {
  slow_tests()
  d <- sim_track("MRME", n = 30, Tend = 30)
  for (v in c("exact", "naive")) {
    res <- smam_worker("fit", fit_spec("MRME", d, c(list(variant = v), mrme_bounds())))
    expect_fit(res, "MRME")
  }
})

test_that("MM fits", {
  slow_tests()
  d <- sim_track("MM", n = 60, Tend = 60)
  res <- smam_worker("fit", fit_spec("MM", d))
  expect_fit(res, "MM")
})

test_that("MRH fits and gives a Hessian", {
  # one MRH likelihood takes seconds even on 8 fixes, so a fit takes many minutes
  if (!identical(Sys.getenv("SMAM_MRH_TEST"), "true")) skip("set SMAM_MRH_TEST=true to run")
  d <- sim_track("MRH", n = 5, Tend = 5)
  p <- smam_models$MRH$params
  res <- smam_worker("fit", fit_spec("MRH", d, list(
    numThreads = 2, lower = vapply(p, `[[`, 0, "lower"), upper = vapply(p, `[[`, 0, "upper"))))
  expect_fit(res, "MRH")
  se <- smam_worker("se", list(model = "MRH", numThreads = 2, estimate = res$estimate, data = d))
  expect_equal(dim(se$vcov), c(5, 5))
})

test_that("standard errors: MRME bootstrap and Godambe, MM bootstrap", {
  slow_tests()
  d <- sim_track("MRME", n = 25, Tend = 25)
  est <- c(lamM = 1, lamR = 0.5, sigma = 1, sig_err = 0.05)
  for (method in c("pBootstrap", "Godambe")) {
    for (variant in c("exact", "naive")) {
      se <- smam_worker("se", list(model = "MRME", variant = variant, method = method,
                                   estimate = est, nBS = 2, numThreads = 2, data = d))
      expect_equal(dim(se$vcov), c(4, 4), info = paste(method, variant))
      if (method == "pBootstrap") expect_equal(ncol(as.matrix(se$boot)), 4)
    }
  }
  dm <- sim_track("MM", n = 30, Tend = 30)
  se <- smam_worker("se", list(model = "MM", method = "pBootstrap", nBS = 2, numThreads = 2,
                               estimate = c(1, 0.1, 1, 0.1), data = dm))
  expect_equal(dim(se$vcov), c(4, 4))
})

test_that("states: MR forward backward and Viterbi", {
  d <- sim_track("MR", n = 50, Tend = 50)
  for (method in c("fb", "viterbi")) {
    res <- smam_worker("state", list(model = "MR", method = method, theta = c(1, 0.5, 1),
                                     cutoff = 0.5, data = d))$states
    expect_equal(nrow(res), nrow(d), info = method)
    expect_true(all(res$states %in% c(0, 1)), info = method)
    if (method == "fb") {
      expect_true(all(abs(res$p.m + res$p.r - 1) < 1e-6))
      expect_s3_class(plot_state_prob(res, c("moving", "resting")), "ggplot")
    }
  }
})

test_that("states: MRH forward backward and Viterbi", {
  slow_tests()
  d <- sim_track("MRH", n = 6, Tend = 6)
  for (method in c("fb", "viterbi")) {
    res <- smam_worker("state", list(model = "MRH", method = method,
                                     theta = c(4, 0.5, 0.1, 5, 0.8), data = d))$states
    expect_equal(nrow(res), nrow(d), info = method)
    expect_true(all(res$states %in% 0:2), info = method)
  }
})
