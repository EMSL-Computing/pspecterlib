#' Creates a modifications_pspecter object for get_matched_peaks
#'
#' @description The modifications objects allows for the generation of a stable
#'     object for testing post-translational methods (PTMs) as more are discovered.
#'
#' @param Name A vector of the names of modifications to check. Required.
#' @param AMU_Change A vector of the AMU mass changes to apply. Required.
#' @param N_Position A vector of vectors of N positions to apply a modification. Required.
#' @param Molecular_Formula A list of lists, with atoms written as the name, and the
#'     number of atoms written as the element list(list("C" = 2, "H" = 3, "O" = 1)). Required.
#'
#' @details The modification_object informs PSpecteR where modifications change
#'    mass values for fragment calculations, as well as in the calculation of expected
#'    abundances (thus why the molecular formula is needed). Users can choose to not
#'    include the molecular formula, though isotopic percentages may be slightly off.
#'
#' @examples
#' \dontrun{
#'
#' # Make an acetyl modification on the 2nd residue
#' PTM <- make_ptm(Name = "Acetyl", AMU_Change = 43.045, N_Position = 2,
#'                 Molecular_Formula = list(list("C" = 2, "H" = 3, "O" = 1)))
#'
#' }
#' @export
make_ptm <- function(Name,
                     AMU_Change,
                     N_Position,
                     Molecular_Formula = NULL) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Assert that Names is a vector of characters
  if (is.null(Name) || !is.character(Name)) {
    stop("Names must be a vector of characters describing each modification.")
  }

  # Assert that AMU_Change is a vector of masses
  if (is.null(AMU_Change) || !is.numeric(AMU_Change)) {
    stop("AMU_Change must be a vector of masses written as a numeric.")
  }

  # Assert that N_Position is a vector integers
  if (is.null(N_Position) || !is.numeric(N_Position)) {
    stop("N_Position should be a vector of integers.")
  }

  # Make N_Position a vector of integers
  N_Position <- N_Position %>% abs() %>% round()

  # Assert that the three vectors are of the same length
  if (length(Name) != length(AMU_Change & length(Name) != length(N_Position))) {
    stop("The length of Names, AMU_Change, and N_Position must be the same.")
  }

  # Run additional checks if Molecular_Formula is not NULL
  if (is.null(Molecular_Formula) == FALSE) {

    # Assert that the length is the same
    if (length(Name) != length(Molecular_Formula)) {
      stop("Molecular_Formula must be the same length as Names, AMU_Change, and N_Position.")
    }

  }

  ######################
  ## BUILD THE OBJECT ##
  ######################

  # Build modification data table
  PTM <- data.table::data.table(
          "Name" = Name,
          "AMU Change" = AMU_Change,
          "N Position" = N_Position
        )


  # If Molecular_Formula is not NULL
  if (is.null(Molecular_Formula) == FALSE) {

    # Get molecular formula objects
    MolForms <- lapply(1:length(Molecular_Formula), function(el) {
      make_molecule(Molecular_Formula[[el]])
    })

    # Add names
    names(MolForms) <- Name

    # Add attribute
    attr(PTM, "pspecter")$MolForms <- MolForms

    # Add molecular formula to the dataframe
    PTM$`Molecular Formula` <- lapply(1:length(MolForms), function(el) {unlist(MolForms[[el]][1, "Formula"])}) %>% unlist()

  }

  # Add class
  class(PTM) <- c(class(PTM), "modifications_pspecter")

  # Return results
  return(PTM)

}
