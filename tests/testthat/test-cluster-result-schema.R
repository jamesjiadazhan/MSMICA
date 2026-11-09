test_that("empty cluster results retain the downstream schema", {
  candidates <- data.frame(
    Mono_mass = numeric(),
    Adduct_annotated = character(),
    mz_annotated = numeric(),
    time_annotated = numeric(),
    mz_time_annotated = character()
  )

  result <- ensure_msmica_cluster_result_schema(data.frame(), candidates)

  expect_equal(nrow(result), 0L)
  expect_true(all(c(
    "Mono_mass", "Adduct_annotated", "mz_annotated", "time_annotated",
    "mz_time_annotated", "correlation", "adduct_corr_cluster",
    "mz_isotope", "time_isotope", "MSMICA_identification", "Probability"
  ) %in% names(result)))
})

test_that("nonempty rows without evidence columns retain rows without fabricated evidence", {
  clusters <- data.frame(correlation = 0.8, Probability = 99)
  candidates <- data.frame(Mono_mass = numeric())

  result <- ensure_msmica_cluster_result_schema(clusters, candidates)
  expect_equal(result$correlation, clusters$correlation)
  expect_equal(result$Probability, clusters$Probability)
  expect_true(is.na(result$mz_isotope))
  expect_true(is.na(result$adduct_corr_cluster))
})

test_that("partially bound empty cluster results gain only missing columns", {
  partial <- data.frame(Mono_mass = numeric(), mz_time_annotated = character())
  result <- ensure_msmica_cluster_result_schema(
    partial,
    data.frame(Adduct_annotated = character())
  )

  expect_true(all(c(
    "Mono_mass", "mz_time_annotated", "correlation", "mz_isotope",
    "time_isotope", "Probability"
  ) %in% names(result)))
  expect_equal(nrow(result), 0L)
})

test_that("cluster rows without isotope correlations keep rows and gain NA evidence", {
  partial <- data.frame(Mono_mass = 100, mz_time_annotated = "100_10")
  result <- ensure_msmica_cluster_result_schema(
    partial,
    data.frame(Adduct_annotated = character())
  )

  expect_equal(nrow(result), 1L)
  expect_true(is.na(result$correlation))
  expect_equal(result$Mono_mass, 100)
})
