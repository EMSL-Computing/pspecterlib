#' Makes a Molecule Dataframe, Helper Function to get_peptide_coverage
#'
#' @description This simple function generates molecule data frames from lists for easy addition
#'
#' @param AtomList A list of atom counts
#'
#' @details A data frame is derived from the list, with each atom count and the formula
#'
#' @examples
#' \dontrun{
#'
#' # Generate a molecule from a list of atoms
#' make_molecule(AtomList = list("C" = 6, "H" = 12, "O" = 6, "N" = 0))
#'
#' }
#'
#' @export
make_molecule <- function(AtomList, FormulaOnly = FALSE) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # AtomList should be a list
  if (!is.list(AtomList)) {
    stop('AtomList should be a list where the names are the element, and the values are the counts of that element. Example: AtomList = list("C" = 6, "H" = 12, "O" = 6, "N" = 0)')
  }

  # If Atom List is length 0, return a carbon with 0
  if (length(AtomList) == 0) {
    res <- list("Molecule" = data.frame("C" = 0, "Formula" = "C0"))
    class(res) <- "molecule_pspecter"
    return(res)
  }


  #############################
  ## MAKE MOLECULE DATAFRAME ##
  #############################

  # Remove 0's
  AtomList <- AtomList[AtomList != 0]

  # Define a simple collapse function
  collapseFUN <- function(x) {paste0(colnames(x), x[1,], collapse = "")}

  # Convert to a dataframe
  Molecule <- AtomList %>% as.data.frame()

  # Make duplicate row names the same
  colnames(Molecule) <- colnames(Molecule) %>%
    gsub(pattern = ".", replacement = "&", fixed = T) %>%
    gsub(pattern = "&.+", replacement = "")

  # Collapse down names
  Molecule <- rowsum(t(Molecule), group = colnames(Molecule), na.rm = TRUE) %>% t() %>% as.data.frame()

  # Make the formula
  Molecule$Formula <- Molecule %>% collapseFUN()

  # If formula only, return only the formula
  if (FormulaOnly) {return(Molecule$Formula)}

  ###################
  ## RETURN OBJECT ##
  ###################

  class(Molecule) <- c(class(Molecule), "molecule_pspecter")

  # Return the molecule data frame
  return(Molecule)

}
