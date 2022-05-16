context("test: make_molecule")

test_that("Confirm that make_molecule is working correctly and producing expected output", {

  # Test make_molecule input checks---------------------------------------------

  # A list must be passed to atom list
  expect_error(
    make_molecule(6),
    "AtomList should be a list where the names are the element, and the values are the counts of that element."
  )

  # FormulaOnly must be TRUE or FALSE
  expect_error(
    make_molecule(list("C" = 3), NA),
    "FormulaOnly must be a TRUE or FALSE."
  )

  # Atom List cannot be length 0
  expect_error(
    make_molecule(list()),
    "AtomList cannot be length 0."
  )

  # Create a molecule-----------------------------------------------------------

  # Create a molecule object
  Molecule1 <- make_molecule(list("C" = 2, "H" = 4, "O" = -1))

  # The colnames should be the input elements and the formula
  expect_equal(colnames(Molecule1), c("C", "H", "O", "Formula"))

  # There should be only one row
  expect_equal(nrow(Molecule1), 1)

  # The class should include molecule_pspecter
  expect_true(inherits(Molecule1, "molecule_pspecter"))

  # Formula only
  expect_equal(make_molecule(list("C" = 2, "H" = 4), FormulaOnly = TRUE), "C2H4")

})
