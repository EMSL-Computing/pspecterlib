context("test: get_protein_table")

test_that("Testing get protein table function", {

  # Download required test files------------------------------------------------

  source(system.file("tests/testthat/file_check.R", package = "pspecterlib"))
  downfolder <- file_check("downfolder")
  BU_ScanMetadata <- file_check("BU")
  FASTA_PATH <- file.path(downfolder, "QC_Shew.fasta")

  # Test get_protein_table input checks-----------------------------------------

  # Scan Metadata object should be a scan_metadata object
  expect_error(
    get_protein_table(data.frame(test = 1)),
    "ScanMetadata must be a scan_metadata object from get_scan_metadata."
  )

  # FASTA Path must exist
  expect_error(
    get_protein_table(
      ScanMetadata = BU_ScanMetadata,
      FASTAPath = "/bad/path"
    ),
    "FASTA file path must exist."
  )

  # FASTA Path must point to a fasta file
  expect_error(
    get_protein_table(
      ScanMetadata = BU_ScanMetadata,
      FASTAPath = file.path(downfolder, "BottomUp.mzML")
    ),
    "FASTA file must have the .FASTA or .FA extension."
  )

  # ID data must be included in Scan Metadata object
  Test_ScanMetadata <- get_scan_metadata(file.path(downfolder, "BottomUp.mzML"))
  expect_error(
    get_protein_table(
      ScanMetadata = Test_ScanMetadata,
      FASTAPath = file.path(downfolder, "QC_Shew.fasta")
    ),
    "No ID data included. Please add identification data to the ScanMetadata object."
  )

  # Q Value maximum must be between 0 and 1
  expect_error(
    get_protein_table(
      ScanMetadata = BU_ScanMetadata,
      FASTAPath = file.path(downfolder, "QC_Shew.fasta"),
      QValueMaximum = 1.2
    ),
    "QValueMaximum must be a single value between 0-1."
  )

  # Score maximum must be between 0 and 1
  expect_error(
    get_protein_table(
      ScanMetadata = BU_ScanMetadata,
      FASTAPath = file.path(downfolder, "QC_Shew.fasta"),
      ScoreMaximum = -2
    ),
    "ScoreMaximum must be a single value between 0-1."
  )

  # Remove contaminants must be TRUE or FALSE
  expect_error(
    get_protein_table(
      ScanMetadata = BU_ScanMetadata,
      FASTAPath = FASTA_PATH,
      RemoveContaminants = "FALSE"
    ),
    "RemoveContaminants must be a TRUE or FALSE."
  )

  # Create a protein_table object-----------------------------------------------

  # Our example BU_ScanMetadata object has NA in Q-Value. Let's create examples of none and both.
  All_SM <- No_Q <- No_Score <- BU_ScanMetadata
  All_SM$`Q Value` <- runif(length(All_SM$`Q Value`))
  No_Score$Score <- NA
  No_Score$`Q Value` <- All_SM$`Q Value`

  # No Q value should trigger a message
  expect_message(
    suppressWarnings({test <- get_protein_table(
      ScanMetadata = No_Q,
      FASTAPath = FASTA_PATH,
      QValueMaximum = 0.5,
      ScoreMaximum = 0.001
    )}),
    "NA values found in Q Values. Filter Ignored."
  )

  # No score should trigger a message
  expect_message(
    suppressWarnings({test <- get_protein_table(
      ScanMetadata = No_Score,
      FASTAPath = FASTA_PATH,
      QValueMaximum = 0.5,
      ScoreMaximum = 0.001,
      RemoveContaminants = FALSE
    )}),
    "NA values found in Score. Filter Ignored."
  )

  # Now, create an example with all
  suppressWarnings({ProteinTable <- get_protein_table(
    ScanMetadata = All_SM,
    FASTAPath = FASTA_PATH
  )})
  expect_true(inherits(ProteinTable, "protein_table"))

})
