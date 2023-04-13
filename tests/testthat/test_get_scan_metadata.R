context("test: get_scan_metadata")

test_that("Testing get scan metadata function", {

  # Download required test files------------------------------------------------

  source(system.file("tests/testthat/file_check.R", package = "pspecterlib"))
  downfolder <- file_check("downfolder")
  BU_ScanMetadata <- file_check("BU")
  RAW_ScanMetadata <- file_check("RAW")

  # Test get_scan_metadata input checks-----------------------------------------

  # Path to mzML file must exist
  expect_error(
    get_scan_metadata("this/is/bad/path"),
    "MS file path must exist."
  )

  # Path must be to an mzML file
  expect_error(
    get_scan_metadata(file.path(downfolder, "BottomUp.mzid")),
    "MS file must be a mzML, mzXML, or raw"
  )

  # Path to ID file must exist
  expect_error(
    get_scan_metadata(MSPath = file.path(downfolder, "BottomUp.mzML"), IDPath = "bad/id/path"),
    "ID file path must exist."
  )

  # Path must be to an mzid file
  expect_error(
    get_scan_metadata(MSPath = file.path(downfolder, "BottomUp.mzML"), IDPath = file.path(downfolder, "BottomUp.mzML")),
    "ID file must be an mzID file."
  )

  # Create a scan_metadata object-----------------------------------------------

  # Create BU_ScanMetadata using the files
  BU_ScanMetadata <- get_scan_metadata(file.path(downfolder, "BottomUp.mzML"), file.path(downfolder, "BottomUp.mzid"))

  # The bottom up example should have 5525 rows and 17 columns.
  expect_equal(nrow(BU_ScanMetadata), 5525)
  expect_equal(ncol(BU_ScanMetadata), 17)

  # The attributes should have the MSPath, the IDPath, the file type, a mzRpwiz object,
  # an MSnbase object, and a list of MS1 Scans
  expect_equal(names(attributes(BU_ScanMetadata)$pspecter), c("MSPath", "IDPath", "MSFileType", "mzRpwiz", "MSnbase", "MS1Scans" ))

  # MSPath and IDPath should be the path to MS and ID file
  expect_equal(attributes(BU_ScanMetadata)$pspecter$MSPath, file.path(downfolder, "BottomUp.mzML"))
  expect_equal(attributes(BU_ScanMetadata)$pspecter$IDPath, file.path(downfolder, "BottomUp.mzid"))

  # The MSFileType should be XML
  expect_equal(attributes(BU_ScanMetadata)$pspecter$MSFileType, "XML")

  # mzRpwiz should not be NULL
  expect_true(!is.null(attributes(BU_ScanMetadata)$pspecter$mzRpwiz))

  # MSnbase should not be NULL as well
  expect_true(!is.null(attributes(BU_ScanMetadata)$pspecter$MSnbase))

  # Use: raw file example
  library(rawrr)
  rawfile <- file.path(path.package(package = 'rawrr'), 'extdata', 'sample.raw')
  RAW_ScanMetadata <- get_scan_metadata(MSPath = rawfile)

  # The MS file type should be RAW
  expect_equal(attributes(RAW_ScanMetadata)$pspecter$MSFileType, "RAW")

  # mzRpwiz should be NULL
  expect_true(is.null(attributes(RAW_ScanMetadata)$pspecter$mzRpwiz))

  # MSnbase should be NULL as well
  expect_true(is.null(attributes(RAW_ScanMetadata)$pspecter$MSnbase))

})
