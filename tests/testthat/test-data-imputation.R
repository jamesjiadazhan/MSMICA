# Exercise the preprocessing expressions from the actual algorithm body,
# stopping before normalization and excluding unrelated database setup.
run_algorithm_imputation <- function(input, method) {
  expressions <- as.list(body(MSMICA_algorithm))[-1]
  text <- vapply(expressions, function(x) paste(deparse(x), collapse = " "), character(1))
  start <- which(text == "met_raw_wide_original = met_raw_wide")
  end <- which(grepl("^met_raw_wide = quantile_normalization", text))
  stopifnot(length(start) == 1L, length(end) == 1L)
  env <- new.env(parent = environment(MSMICA_algorithm))
  env$met_raw_wide <- input
  env$imputation_method <- method
  for (i in seq.int(start + 1L, end - 1L)) eval(expressions[[i]], env)
  env$met_raw_wide
}

test_that("the algorithm imputes half minimum before log2 transformation", {
  raw <- data.frame(mz = c(100, 200, 300), time = c(10, 20, 30),
                    a = c(1024, 0.25, 1), b = c(2048, 0.5, 2),
                    c = c(NA_real_, NA_real_, NA_real_))
  result <- run_algorithm_imputation(raw, "half_min")
  expect_equal(result$c, log2(c(512, 0.125, 0.5)))
  expect_equal(as.matrix(result[, 1:2]), as.matrix(raw[, 1:2]))
  expect_equal(as.matrix(result[, 3:4]), log2(as.matrix(raw[, 3:4])))
  expect_identical(names(result), names(raw))
})

test_that("the algorithm still passes log2 intensities to QRILC", {
  skip_if_not_installed("imputeLCMD")
  input <- data.frame(mz = seq_len(20) + 100, time = seq_len(20),
                      a = seq_len(20) + 10, b = seq_len(20) + 20,
                      c = seq_len(20) + 30)
  input$a[c(2, 5)] <- NA_real_
  input$b[c(3, 8)] <- NA_real_
  input$c[c(4, 9)] <- NA_real_
  logged <- data_transformation(input)
  set.seed(17)
  expected <- imputeLCMD::impute.QRILC(as.data.frame(logged[, 3:5]))[[1]]
  set.seed(17)
  result <- run_algorithm_imputation(input, "QRILC")
  expect_equal(as.matrix(result[, 3:5]), as.matrix(expected))
})
