context("test: make_ptm")

test_that("Testing make ptm function", {

  # Test make_ptm input checks--------------------------------------------------

  # Name must be a string
  expect_error(
    make_ptm(1),
    "Name must be a vector of strings describing each modification."
  )

  # AMU Change must be a vector of masses
  expect_error(
    make_ptm("Test", "32"),
    "AMU_Change must be a vector of masses written as a numeric."
  )

  # N Position should be a vector of integers
  expect_error(
    make_ptm("Test", 32, "2"),
    "N_Position should be a vector of integers."
  )

  # The length of each vector should be the same
  expect_error(
    make_ptm(
      Name = c("Acetyl", "Methyl", "Test"),
      AMU_Change = c(42.010565, 14.015650, 5),
      N_Position = c(2, 3)
    ),
    "The length of Names, AMU_Change, and N_Position must be the same."
  )

  expect_error(
    make_ptm(
      Name = c("Acetyl", "Methyl"),
      AMU_Change = c(42.010565, 14.015650),
      N_Position = c(2, 3),
      Molecular_Formula = list(list("H" = 2, "C" = 2, "O" = 1))
    ),
    "Molecular_Formula must be the same length as Names, AMU_Change, and N_Position."
  )

  # Make a peak_data object-----------------------------------------------------

  TestPTM <- make_ptm(
    Name = "Acetyl",
    AMU_Change = 43.045,
    N_Position = 2,
    Molecular_Formula = list(list("C" = 2, "H" = 3, "O" = 1))
  )
  expect_true(inherits(TestPTM, "modifications_pspecter"))

})

