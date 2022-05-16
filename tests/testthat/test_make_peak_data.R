context("test: make_peak_data")

test_that("Testing peak data making function", {

  # Test make_peak_data input checks--------------------------------------------

  # MZ should be numeric
  expect_error(
    make_peak_data(MZ = c("22")),
    "MZ must be numeric."
  )

  # MZ should be non-zero and positive
  expect_error(
    make_peak_data(MZ = c(0)),
    "MZ should contain non-zero and positive values."
  )

  # Intensity should be numeric
  expect_error(
    make_peak_data(MZ = c(22, 23), Intensity = c("100", 200)),
    "Intensity must be numeric."
  )

  # MZ should be non-zero and positive
  expect_error(
    make_peak_data(MZ = c(22, 23), Intensity = c(100, -200)),
    "Intensity should contain non-zero and positive values."
  )

  # MZ and Intensity should be the same length
  expect_error(
    make_peak_data(MZ = c(22, 23), Intensity = c(100, 200, 300)),
    "Both MZ and Intensity should be the same length."
  )

  # Make a peak_data object-----------------------------------------------------

  PeakData <- make_peak_data(c(22,23,24), c(100,300,100))
  expect_true(inherits(PeakData, "peak_data"))

})
