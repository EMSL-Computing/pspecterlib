context("test: multiple modifications")

test_that("Testing multiple modifications", {
  
  # Check inputs----------------------------------------------------------------
  
  # Sequence should be a string
  expect_error(
    multiple_modifications(TRUE, "Methyl,(3)[1]"),
    "Sequence must be a string."
  )
  
  # Sequence should be an acceptable sequence
  expect_error(
    multiple_modifications("TRUE", "Methyl,(3)[1]"),
    "Sequence is not an acceptable peptide/protein sequence. See ?pspecterlib::is_sequence for more details.",
    fixed = T
  )
  
  # Modification should be a string
  expect_error(
    multiple_modifications("TRICITIES", TRUE),
    "Modification must be a string."
  )
  
  # Write an incorrect modification format
  expect_error(
    multiple_modifications("TRICITIES", "Modify[3rd]Position"),
    "The proper Modification format is PTM,Residue(Positions)Number of Modifications, separated by semicolon.",
    fixed = T
  )
  
  # Unmodified string should be a logical
  expect_error(
    multiple_modifications("TRICITIES", "Methyl,(1^,2,3,4,5,6,7^,8,9)[3];1.00727,(2,4,9)[1]", "FALSE"),
    "ReturnUnmodified must be TRUE or FALSE."
  )
  
  # Modification should be in the glossary
  expect_error(
    multiple_modifications("TRICITIES", "Banana,(1^,2,3,4,5,6,7^,8,9)[3];1.00727,(2,4,9)[1]"),
    "Modification Banana is not in the backend glossary."
  )
  
  # Run test function-----------------------------------------------------------
  
  # Try to test everything at once - should run without error
  Glossary <- data.table::fread(system.file("extdata", "Unimod_v20220602.csv", package = "pspecterlib"))
  Complex_Query <- multiple_modifications("TRICITIES", "Methyl,(1^,2,3,4,5,6,7^,8,9)[2,3];1.00727,(2,4,9)[1]", TRUE,
                         AlternativeGlossary = Glossary[Glossary$Modification == "Methyl",])
  Complex_Query2 <- multiple_modifications("TRICITIES", "Methyl,(1^,2,3,4,5,6,7^,8,9)[2,3];1.00727,(2,4,9)[1]", TRUE)
  
  

})