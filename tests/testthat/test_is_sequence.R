context("test: is_sequence")

test_that("Testing sequence checks", {

  # Test is_sequence returns TRUE for only a correct amino acid sequence--------

  # Nothing should return a FALSE
  expect_false(is_sequence(NULL))

  # A non-string should return a FALSE
  expect_false(is_sequence(1))

  # Spaces should return FALSE
  expect_false(is_sequence("TEST TEST"))

  # Non-letters should return FALSE
  expect_false(is_sequence("TEST.THAT"))

  # Sequences cannot contain B, J, O, U, or X
  expect_false(is_sequence("TESTB"))
  
  # Sequences with brackets must be in glossary if not numeric
  expect_false(is_sequence("TEST[55.32.]"))
  expect_false(is_sequence("TOEST[Acetyl]"))

  # Otherwise, the sequence should be TRUE
  expect_true(is_sequence("TESTTHAT"))
  
  # Sequence can have brackets
  expect_true(is_sequence("TEST[Acetyl]"))
  
  # Ultimate test
  expect_true(is_sequence("T[55]E[Methyl]S[Acetyl]T[23.2345]"))

})
