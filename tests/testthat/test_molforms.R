context("test: molform functions")

test_that("Testing molform functions", {

  # Test as.molform input checks------------------------------------------------
  
  # MolForm must be a string
  expect_error(
    as.molform(TRUE),
    "MolForm must be a single string."
  )
  
  # MolForm must be real elements
  expect_error(
    as.molform("An1Dr2E3W3"),
    "The elements An, Dr, E are not supported at this time."
  )
  
  # Test collapse_molform-------------------------------------------------------
  
  # Molform should be the proper object
  expect_error(
    collapse_molform("C2"),
    "molform should be an object from the molform class."
  )
  
  # 1's should be removed
  expect_equal(collapse_molform(as.molform("C2H2OONa")), "C2H2NaO2")
  
  # Test get_mw-----------------------------------------------------------------
  
  # Molform should be the proper object
  expect_error(
    get_mw("C2"),
    "molform should be an object from the molform class."
  )
  
  # Calculate mw of glucose
  expect_equal(round(get_mw(as.molform("C6H12O6")), 4), 180.1561)
  
  # Calculate mw of methyl peptide modification
  expect_equal(round(get_mw(as.molform("CH2")), 4), 14.0266)
  
  # Test get_monoisotopic-------------------------------------------------------
  
  # Molform should be the proper object
  expect_error(
    get_monoisotopic("C2"),
    "molform should be an object from the molform class."
  )
  
  # 2 carbons should be a simple 24
  expect_equal(get_monoisotopic(as.molform("C2")), 24)
  
  # Test add_molforms input checks----------------------------------------------

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

  # Add molecules with and without negatives capped
  AddMol1 <- add_molforms(Molecule1, Molecule1, CapNegatives = TRUE) %>% collapse_molform()
  AddMol2 <- add_molforms(Molecule1, Molecule1, CapNegatives = FALSE) %>% collapse_molform()

  # Double check formulas
  expect_equal(AddMol1, "C2H2O6")
  expect_equal(AddMol2, "C2H2N-2O6")

  # Test multiply_molforms input checks-----------------------------------------
  
  # Molform should be the proper object
  expect_error(
    multiply_molforms("CH2O", 6),
    "molform should be a molform object from as.molform."
  )
  
  # Scalar should be a single numeric
  expect_error(
    multiply_molforms(as.molform("CH2O"), "6"),
    "scalar should be a single numeric."
  )
  
  # If scalar is 0, none is returned
  expect_null(multiply_molforms(as.molform("CH2O"), 0))
  
  # Multiplying by 6 should give C6H12O6
  expect_equal(collapse_molform(multiply_molforms(as.molform("CH2O"), 6)), "C6H12O6")
  
  # Test get_aa_molform---------------------------------------------------------
  
  # Molform should be the proper object
  expect_error(
    get_aa_molform("BADSEQUENCE"),
    "BADSEQUENCE is not an acceptable amino acid sequence. See ?is_sequence.", fixed = T
  )
  
  # 
  
  
  
  
})
