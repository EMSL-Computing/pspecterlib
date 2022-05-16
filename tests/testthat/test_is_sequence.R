context("test: is_sequence")

test_that("Testing sequence checks", {

  # Test is_sequence returns TRUE for only a correct amino acid sequence--------

  # Nothing should return a FALSE
  expect_false(is_sequence(NULL))

  # A non-string should return an error
  expect_error(
    is_sequence(1),
    "The provided Sequence must be a single string. Example: TEST"
  )

  # A vector longer than 1 should return an error
  expect_error(
    is_sequence(c("TEST", "TEST")),
    "The provided Sequence must be a single string. Example: TEST"
  )

  # Spaces should return FALSE
  expect_false(is_sequence("TEST TEST"))

  # Sequence should not have brackets and should return FALSE
  expect_false(is_sequence("TEST[ACETYL]"))

  # Non-letters should return FALSE
  expect_false(is_sequence("TEST.THAT"))

  # A sequence cannot be a single residue
  expect_false(is_sequence("T"))

  # Sequences cannot contain B, J, O, U, or X
  expect_false(is_sequence("TESTB"))

  # Otherwise, the sequence should be TRUE
  expect_true(is_sequence("TESTTHAT"))

})
