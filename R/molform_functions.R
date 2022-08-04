#' Make a molecular formula object
#' 
#' @description This function generates a molecular formula object. 
#' 
#' @param MolForm A molecular formula written as a string like "C6H12O6". Only 
#'    elements up to uranium are supported. 
#'     
#' @details A list of length 92 (Hydrogen to Uranium) is stored in the attributes 
#'    and used for quick mathematics. 
#'    
#' @examples
#' \dontrun{
#' 
#' # Generate a molecular formula object
#' as.molform("C6H12O6")
#' 
#' # This should fail, as 'Ot' and 'Lp' are not elements 
#' as.molform("C6H12Ot5Lp2")
#' 
#' }
#' 
#' @export
as.molform <- function(MolForm) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that the molecular formula is a string
  if (!is.character(MolForm) | length(MolForm) != 1) {
    stop("MolForm must be a single string.")
  }
  
  # Split out elements 
  Elements <- strsplit(MolForm, "[0-9]") %>% unlist() %>% .[. != ""]
  
  # List acceptable elements
  AcceptableElements <- c(
   "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", 
   "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
   "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
   "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
   "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", 
   "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
   "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
   "Fr", "Ra", "Ac", "Th", "Pa", "U"
  )
  
  # Check for non-established elements
  NonEstablished <- Elements[Elements %in% AcceptableElements == FALSE]
  
  if (length(NonEstablished) != 0) {
    stop(paste("The elements", paste0(NonEstablished, sep = ", "), 
               "are not supported at this time."))
  }
  
  #########################
  ## GENERATE ATTRIBUTES ##
  #########################
  
  
  
}    