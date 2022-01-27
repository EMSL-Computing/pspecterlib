#' Creates a modified_ion object for get_matched_peaks
#'
#' @description The modified ion allows for additional ion groups (an a-c or x-z
#'      ion with a mass shift) to be tested.
#'
#'
#' @param Ion A vector of ions to modify. Must be a, b, c, x, y, z. Required.
#' @param Symbol An additional symbol to distinguish this ion from the original
#'    a-c or x-z. Acceptable inputs are +, ++, -, --, ^, or ^^. Must be the same
#'    length as the Ion vector. Required.
#' @param AMU_Change A vector of the AMU mass changes to apply. Must be the same
#'    length as the Ion vector. Required.
#'
#' @details Depending on the analysis, users may want to investigate ion groups
#'    (a-c or x-z) with an additional proton, water loss, etc. This allows for that
#'    functionality, as well as a user-defined symbol for that change.
#'
#' @examples
#' \dontrun{
#'
#' make_mass_modified_ion(Ion = c("b", "z"),
#'                        Symbol = c("+", "+"),
#'                        AMU_Change = c(1.00727647, 1.00727647))
#'
#' }
#' @export
make_mass_modified_ion <- function(Ion,
                                   Symbol,
                                   AMU_Change) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Assert that ion is a-c or x-z
  if (!is.character(Ion) | all(unique(Ion) %in% c("a", "b", "c", "x", "y", "z")) == FALSE) {
    stop("Ion must be any combination of a, b, c, x, y, or z.")
  }

  # Assert that symbol is an acceptable input
  if (!is.character(Symbol) | all(unique(Symbol) %in% c("+", "++", "-", "--", "^", "^^")) == FALSE) {
    stop("Symbol must be any combination of +, ++, -, --, ^, or ^^.")
  }

  # Assert that AMU Change is a nonzero numeric
  if (!is.numeric(AMU_Change) | any(AMU_Change == 0)) {
    stop("AMU_Change must be non-zero numerics.")
  }

  # Assert that all three vectors are of the same length
  if (length(Ion) != length(Symbol) | length(Ion) != length(AMU_Change)) {
    stop("Ion, Symbol, and AMU_Change must all be the same length.")
  }

  #####################
  ## GENERATE OBJECT ##
  #####################

  # Generate modified ion object
  ModifiedIon <- data.table::data.table(Ion = Ion, Modified_Ion = paste0(Ion, Symbol), AMU_Change)

  # Test for duplicates
  IonCounts <- table(ModifiedIon$Modified_Ion, dnn = "Ion") %>% data.frame()
  if (any(IonCounts$Freq > 1)) {
    stop(paste("Duplicate ions detected:", paste(IonCounts[IonCounts$Freq > 1, "Ion"], collapse = ", ")))
  }

  # Otherwise, return object
  class(ModifiedIon) <- c(class(ModifiedIon), "modified_ion")

  return(ModifiedIon)

}





