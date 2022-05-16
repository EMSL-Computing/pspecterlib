context("test: make_mass_modified_ion")

test_that("Testing mass modified ion object creator function", {

  # Test make_mass_modified_ion input checks------------------------------------

  # Ion must be an a, b, c, x, y, and or z
  expect_error(
    make_mass_modified_ion(Ion = c("a", "b", "c", "d", "x", "y", "z")),
    "Ion must be any combination of a, b, c, x, y, or z."
  )

  # Symbol must be "+", "++", "-", "--", "^", and or "^^"
  expect_error(
    make_mass_modified_ion(
      Ion = c("a", "a", "a"),
      Symbol = c("+", "++", "+++")
    )
  )

  # AMU Change must be numeric
  expect_error(make_mass_modified_ion("a", "+", "32"), "AMU_Change must be non-zero numerics.")

  # All three input vectors must be of the same length
  expect_error(
    make_mass_modified_ion(
      Ion = c("a", "b", "c"),
      Symbol = c("+", "-", "^"),
      AMU_Change = c(1, 2)
    ),
    "Ion, Symbol, and AMU_Change must all be the same length."
  )

  # Make a mass modified ion----------------------------------------------------

  # Run an example that shouldn't error out
  MassModifiedIon <- make_mass_modified_ion("a", "+", 1.007)
  expect_true(inherits(MassModifiedIon, "modified_ion"))

  # Error out if duplicate ions are generated
  expect_error(
    make_mass_modified_ion(c("a", "a"), c("+", "+"), c(2, 3)),
    "Duplicate ions detected: a+"
  )

})
