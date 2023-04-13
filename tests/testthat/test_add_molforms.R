context("test: add_molforms")

test_that("Testing add molecule functionality", {

  # Test add_molecules input checks---------------------------------------------

  # Make an example molecule
  Molecule1 <- as.molform("CHO3N-1")

  # All objects passed to add_molecules must be molecule_pspecter objects
  expect_error(
    add_molforms(Molecule1, list("C" = 3, "H" = 2)),
    "All input objects must be molform objects from as.molform."
  )

  # CapNegative must be a true or false
  expect_error(
    add_molforms(Molecule1, CapNegatives = NA),
    "CapNegatives must be a TRUE or FALSE."
  )

  # Run the function------------------------------------------------------------

  # Add molecules with and without negatives capped
  AddMol1 <- add_molforms(Molecule1, Molecule1, CapNegatives = TRUE) %>% collapse_molform()
  AddMol2 <- add_molforms(Molecule1, Molecule1, CapNegatives = FALSE) %>% collapse_molform()

  # Double check formulas
  expect_equal(AddMol1, "C2H2O6")
  expect_equal(AddMol2, "C2H2N-2O6")

})
