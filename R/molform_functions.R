#' Make a molecular formula object
#' 
#' @description This function generates a molecular formula object. 
#' 
#' @param MolForm A molecular formula written as a string like "C6H12O6". Only 
#'    elements up to uranium are supported. 
#'     
#' @details A list of length 92 (Hydrogen to Uranium) is stored in the attributes 
#'    and used for quick mathematics. All numerics are rounded to the nearest whole
#'    number.
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
  Elements <- strsplit(MolForm, "-?[0-9]") %>% unlist() %>% .[. != ""]
  
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
    stop(paste("The elements", paste0(NonEstablished, collapse = ", "), 
               "are not supported at this time."))
  }
  
  #########################
  ## GENERATE ATTRIBUTES ##
  #########################
  
  # Pull elemental counts
  Counts <- strsplit(MolForm, "[[:alpha:]]") %>% unlist() %>% .[. != ""] %>% as.numeric()
  Coutns <- round(Counts)
  
  # Fill an atomic vector
  AtomList <- rep(0, 92)
  names(AtomList) <- AcceptableElements
  
  # Add counts
  AtomList[Elements] <- Counts
  
  # Add the attribute
  attr(MolForm, "AtomList") <- AtomList
  
  # Add a class
  class(MolForm) <- c(class(MolForm), "molform")
  
  # Return the object
  return(MolForm)
  
}    

#' Add molecular formula objects
#' 
#' @description Add molform objects together
#' 
#' @param ... Any number of molecules from the as.molform class
#' @param CapNegatives A TRUE/FALSE to indicate whether negative elements should be
#'     capped at 0. Useful for modifications where elements are lost. Default is TRUE.
#' 
#' @examples
#' \dontrun{
#' 
#' add_molforms(as.molform("C6H12O6"), as.molform("H2Fr2Sr2Ar5"), as.molform("H-1C-2"))
#' 
#' }
#' @export
add_molforms <- function(..., CapNegatives = TRUE) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check for extra molecules
  molecule_objects <- list(...)
  
  # Iterate through all combinations and determine that it is an object of the
  # molecule_pspecter class
  ObjectTest <- lapply(molecule_objects, function(molecule) {
    inherits(molecule, "molform")
  }) %>% unlist() %>% all()
  
  # Check that all objects are molecule_pspecter objects
  if (ObjectTest == FALSE) {
    stop("All inputs must be molform objects from as.molform")
  }
  
  # Check that CapNegatives is a TRUE or FALSE
  if (!is.logical(CapNegatives) | is.na(CapNegatives)) {
    stop("CapNegatives must be a TRUE or FALSE.")
  }
  
  ################
  ## ADD VALUES ##
  ################
  
  # Pull and add all AtomLists
  Total <- rep(0, 92)
  names(Total) <- names(attr(molecule_objects[[1]], "AtomList"))
  for (item in molecule_objects) {
    Total <- Total + attr(item, "AtomList")
  }
  
  # Remove negatives if enabled
  if (CapNegatives) {Total <- Total[Total>0]} else {Total <- Total[Total != 0]}
  
  # Return NULL if nothing is left
  if (length(Total) == 0) {return(NULL)}
  
  # Build a molform object
  Total <- lapply(1:length(Total), function(x) {paste(names(Total)[x], Total[x], sep = "")}) %>% 
    unlist() %>% 
    paste0(collapse = "")
  
  return(as.molform(Total))

}

#' Multiply molecular formula objects 
#' 
#' @description Multiply molecular formulas by a constant
#' 
#' @param molform An object of the as.molform class
#' @param scalar A positive integer determining the number of times to multiply 
#'    a molecular formula by. 
#' 
#' @examples
#' \dontrun{
#' 
#' multiply_molforms(as.molform("C6H12O6"), 3)
#' 
#' }
#' @export
multiply_molforms <- function(molform, scalar) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that molform is a molecular formula object
  if (!inherits(molform, "molform")) {
    stop("molform should be an object of the as.molform class.")
  }
  
  # Check that the scalar is a single numeric
  if (!is.numeric(scalar) | length(scalar) != 1) {
    stop("scalar should be a single numeric.")
  }
  scalar <- abs(round(scalar))
  if (scalar == 0) {return(NULL)}
  
  #####################
  ## MULTIPLY VALUES ##
  #####################
  
  # Multiply the object
  Total <- attr(molform, "AtomList") * scalar
  Total <- Total[Total != 0]
  
  # Build a molform object
  Total <- lapply(1:length(Total), function(x) {paste(names(Total)[x], Total[x], sep = "")}) %>% 
    unlist() %>% 
    paste0(collapse = "")
  
  return(as.molform(Total))
  
}

#' Calculate a molecular formula for an amino acid sequence
#' 
#' @param sequence The amino acid sequence. 
#' 
#' @examples 
#' \dontrun{
#' 
#' getAAMolForm("TESTTESTTESTTEST")
#' 
#' }
#' @export
getAAMolForm <- function(sequence) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  if (!is_sequence(sequence)) {
    stop(paste0(sequence, " is not an acceptable amino acid sequence. See ?is_sequence."))
  }
  
  #####################################
  ## CALCULATE THE MOLECULAR FORMULA ##
  #####################################
  
  # Initialize a dictionary of molecular formulas 
  MolFormDict <- function(x) {
    switch(x, 
           "I" = as.molform("C6H13N1O2"), 
           "L" = as.molform("C6H13N1O2"), 
           "K" = as.molform("C6H14N2O2"), 
           "M" = as.molform("C5H11N1O2S1"), 
           "F" = as.molform("C9H11N1O2"), 
           "T" = as.molform("C4H9N1O3"), 
           "W" = as.molform("C11H12N2O2"), 
           "V" = as.molform("C5H11N1O2"), 
           "R" = as.molform("C6H14N4O2"), 
           "H" = as.molform("C6H9N3O2S"), 
           "A" = as.molform("C3H7N1O2"), 
           "N" = as.molform("C4H8N2O3"), 
           "D" = as.molform("C4H7N1O4"), 
           "C" = as.molform("C3H7N1O2S1"), 
           "E" = as.molform("C5H9N1O4"), 
           "Q" = as.molform("C5H10N2O3"), 
           "G" = as.molform("C2H5N1O2"), 
           "P" = as.molform("C5H9N1O2"), 
           "S" = as.molform("C3H7N1O3"), 
           "Y" = as.molform("C9H11N1O3")
    )}
  
  # Get counts of amino acids 
  AACounts <- strsplit(sequence, "") %>% 
    unlist() %>% 
    table(dnn = "AA") %>% 
    data.frame() %>%
    dplyr::mutate(AA = as.character(AA))
  
  # Generate a molform object
  Total <- as.molform("")
  for (row in 1:nrow(AACounts)) {
    Total <- add_molforms(Total, multiply_molforms(MolFormDict(AACounts[row, "AA"]), AACounts[row, "Freq"]))
  }
  
  # Remove waters 
  Total <- add_molforms(Total, multiply_molforms(as.molform("H-2O-1"), nchar(sequence) - 1))

  # Remove zeroes
  Total <- Total[Total != 0]
  
  # Build a molform object
  Total <- lapply(1:length(Total), function(x) {paste(names(Total)[x], Total[x], sep = "")}) %>% 
    unlist() %>% 
    paste0(collapse = "")
  
  return(as.molform(Total))
  
}




