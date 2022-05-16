context("test: get_ms1ft")

test_that("Testing get ms1ft function", {

  # Download required test files------------------------------------------------

  source(system.file("tests/testthat/file_check.R", package = "pspecterlib"))
  downfolder <- file_check("downfolder")

  # Test get_ms1ft input checks-------------------------------------------------

  # The MS1FT file path must exist
  expect_error(
    get_ms1ft("/bad/path"),
    "MS1FT file path must exist."
  )

  # The MS1FT file must be an ms1ft
  expect_error(
    get_ms1ft(file.path(downfolder, "BottomUp.mzid")),
    "MS1FT file must have the .ms1ft extension."
  )

  # If provided, the targets file must exist
  expect_error(
    get_ms1ft(
      MS1FTPath = file.path(downfolder, "TopDown.ms1ft"),
      TargetsPath = "/bad/path"
    ),
    "Targets file path must exist."
  )

  # The targets file must be a tsv
  expect_error(
    get_ms1ft(
      MS1FTPath = file.path(downfolder, "TopDown.ms1ft"),
      TargetsPath = file.path(downfolder, "BottomUp.mzML")
    ),
    "TargetsPath must be a tsv."
  )

  # Trim columns should be a true or false
  expect_error(
    get_ms1ft(file.path(downfolder, "TopDown.ms1ft"), TrimColumns = NA),
    "TrimColumns must be a TRUE or FALSE."
  )

  # Create a ms1ft object-------------------------------------------------------

  # Make an example with both files
  MS1FT <- get_ms1ft(
    MS1FTPath = file.path(downfolder, "TopDown.ms1ft"),
    TargetsPath = file.path(downfolder, "TopDown_IcTarget.tsv")
  )
  expect_true(inherits(MS1FT, "ms1ft"))

  # Make an example with no targets
  MS1FT_2 <- get_ms1ft(
    MS1FTPath = file.path(downfolder, "TopDown.ms1ft")
  )
  expect_true(inherits(MS1FT_2, "ms1ft"))
  expect_true(length(unique(MS1FT_2$ProteinName)) == 1)

})
