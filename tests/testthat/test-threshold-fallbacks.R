test_that("one calibration pair is the default minimum", {
  expect_equal(formals(estimate_adduct_clustering_thresholds)$min_pairs, 1)
  expect_equal(formals(estimate_isotope_clustering_thresholds)$min_pairs, 1)
})

test_that("adduct and isotope defaults are returned when calibration fails", {
  fail_estimation <- function(...) stop("not enough usable pairs")

  adduct <- msmica_threshold_estimate_with_fallback(
    fail_estimation, list(), "adduct", 10, 0.4
  )
  expect_equal(adduct$adduct_correlation_time_threshold, 10)
  expect_equal(adduct$adduct_correlation_r_threshold, 0.4)
  expect_equal(adduct$summary$threshold_method, "fixed_default_fallback")
  expect_match(adduct$summary$fallback_reason, "not enough usable pairs")

  isotope <- msmica_threshold_estimate_with_fallback(
    fail_estimation, list(), "isotope", 5, 0.7
  )
  expect_equal(isotope$isotopic_correlation_time_threshold, 5)
  expect_equal(isotope$isotopic_correlation_r_threshold, 0.7)
  expect_equal(isotope$summary$threshold_method, "fixed_default_fallback")
  expect_match(isotope$summary$fallback_reason, "not enough usable pairs")
})

test_that("finite empirical thresholds remain unchanged", {
  estimate <- function(...) list(
    adduct_correlation_time_threshold = 3,
    adduct_correlation_r_threshold = 0.63,
    summary = data.frame(threshold_method = "adduct_empirical")
  )
  result <- msmica_threshold_estimate_with_fallback(
    estimate, list(), "adduct", 10, 0.4
  )
  expect_equal(result$adduct_correlation_time_threshold, 3)
  expect_equal(result$adduct_correlation_r_threshold, 0.63)
  expect_equal(result$summary$threshold_method, "adduct_empirical")
})

test_that("fallback summaries retain observed calibration-pair counts", {
  fail_estimation <- function(...) stop("only 3 calibration pairs were found; at least 20 are required")
  result <- msmica_threshold_estimate_with_fallback(
    fail_estimation, list(), "isotope", 5, 0.7
  )
  expect_equal(result$summary$n_pairs, 3)
  expect_equal(result$summary$n_pairs_for_threshold, 3)
})
