test_that("half minimum agrees with halving raw intensities before logging", {
  raw <- data.frame(mz = c(100, 200, 300), time = c(10, 20, 30),
                    a = c(1024, 0.25, 1), b = c(2048, 0.5, 2),
                    c = c(NA_real_, NA_real_, NA_real_))
  result <- data_imputation(data_transformation(raw), "half_min")
  expect_equal(result$c, log2(c(512, 0.125, 0.5)))
  expect_equal(result[, 1:2], raw[, 1:2])
  expect_equal(as.matrix(result[, 3:4]), log2(as.matrix(raw[, 3:4])))
  expect_identical(names(result), names(raw))
})

test_that("entirely missing features remain missing without infinite imputations", {
  input <- data.frame(mz = c(100, 200), time = c(10, 20),
                      a = c(NA_real_, 10), b = c(NA_real_, 11))
  result <- data_imputation(input, "half_min")
  expect_true(all(is.na(result[1, 3:4])))
  expect_equal(unlist(result[2, 3:4], use.names = FALSE), c(10, 11))
})

test_that("QRILC continues to receive log2 intensities", {
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
  result <- data_imputation(logged, "QRILC")
  expect_equal(as.matrix(result[, 3:5]), as.matrix(expected))
})

test_that("NA disables imputation as documented", {
  input <- data.frame(mz = 100, time = 10, a = 5, b = NA_real_)
  expect_identical(data_imputation(input, NA), input)
})
