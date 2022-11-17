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
  Counts <- round(Counts)
  
  # Fill an atomic vector
  AtomList <- rep(0, length(AcceptableElements))
  names(AtomList) <- AcceptableElements
  AtomList[Elements] <- Counts
  
  # Add a class
  class(AtomList) <- c(class(AtomList), "molform")
  
  # Return the object
  return(AtomList)
  
}    

#' Collapse a molform object to a simple molecular formula
#' 
#' @description A simple printing option for an object of the as.molform class. 
#' 
#' @examples
#' \dontrun{
#' print(as.molform("C6H12O6"))
#' }
#' @export
collapse_molform <- function(molform) {
  
  ###################
  ## CHECK INPUTS ##
  ##################
  
  # Molform should be an object of the molform class
  if (!inherits(molform, "molform")) {
    stop("molform should be an object from the molform class.")
  }
  
  ######################
  ## COLLAPSE MOLFORM ##
  ######################
  
  collapsed <- molform[molform != 0]
  return(paste0(names(collapsed), collapsed, sep = "", collapse = ""))
  
}

#' Get average mass from the molecular formula 
#' 
#' @description Calculate the average mass for a molform object.
#'
#' @param molform An object of the as.molform class
#' 
#' @examples
#' \dontrun{
#' 
#' get_mw(as.molform("C6H12O6"))
#' 
#' }
#' @export
get_mw <- function(molform) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Molform should be an object of the molform class
  if (!inherits(molform, "molform")) {
    stop("molform should be an object from the molform class.")
  }
  
  ######################
  ## CALCULATE RESULT ##
  ######################
  
  AverageMasses <- c(1.007940721855, 4.00260164839362, 6.9400370562, 9.012182, 10.811027768, 
     12.0107358985, 14.0067430888, 15.99940530471, 18.998403, 20.1800463039, 
     22.98977, 24.3050518651, 26.981538, 28.085413282599, 30.973762, 
     32.0660849916, 35.45253851, 39.098301439102, 39.947676594063, 
     40.07802265518, 44.95591, 47.8667497092, 50.9414719975, 51.99613763897, 
     54.93805, 55.84515013384, 58.693356464998, 58.9332, 63.5456439019, 
     65.3955668951, 69.72307154608, 72.6127589578, 74.921596, 78.9593889728, 
     79.9035286243, 83.7993250788, 85.4676637502, 87.6166459844, 88.905848, 
     91.2236473858, 92.906378, 95.9312915794, 97.907216, 101.0649451127, 
     102.905504, 106.4153280285, 107.86815069743, 112.4115526683, 
     114.8180858507, 118.7101106395, 121.7597883042, 126.904468, 127.6031253845, 
     131.2924806487, 132.905447, 137.32688569073, 138.9054486831, 
     140.11572154786, 140.907648, 144.236126978, 144.912744, 150.3663440037, 
     151.964366222, 157.2521192511, 158.925343, 162.4970300387, 164.930319, 
     167.2563010726, 168.934211, 173.0376918021, 174.9667175726, 178.4849709358, 
     180.9478759364, 183.8417786794, 186.20670567, 190.2248610803, 
     192.216053791, 195.07779131504, 196.966552, 200.5991493631, 204.38331701508, 
     207.21689158, 208.980383, 208.982416, 209.987131, 222.01757, 
     223.019731, 226.025403, 227.027747, 231.035879, 232.03805, 238.028913066965
  )
  
  return((molform * AverageMasses) %>% sum())
  
}

#' Get monoisotopic mass from the molecular formula 
#' 
#' @description Calculate the monoisotopic mass for a molform object. No charge is assumed.
#'
#' @param molform An object of the as.molform class
#' 
#' @examples
#' \dontrun{
#' 
#' get_monoisotopic(as.molform("C6H12O6"))
#' 
#' }
#' @export
get_monoisotopic <- function(molform) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Molform should be an object of the molform class
  if (!inherits(molform, "molform")) {
    stop("molform should be an object from the molform class.")
  }
  
  ######################
  ## CALCULATE RESULT ##
  ######################
  
  MonoisotopicMasses <- c(1.007825, 3.016029, 6.015122, 9.012182, 
                      10.012937, 12, 14.003074, 15.994915, 18.998403, 19.99244, 22.98977, 
                      23.985042, 26.981538, 27.976927, 30.973762, 31.972071, 34.968853, 
                      35.967546, 38.963707, 39.962591, 44.95591, 45.952629, 49.947163, 
                      49.94605, 54.93805, 53.939615, 58.9332, 57.935348, 62.929601, 
                      63.929147, 68.925581, 69.92425, 74.921596, 73.922477, 78.918338, 
                      77.920386, 84.911789, 83.913425, 88.905848, 89.904704, 92.906378, 
                      91.90681, 97.907216, 95.907598, 102.905504, 101.905608, 106.905093, 
                      105.906458, 112.904061, 111.904821, 120.903818, 119.90402, 126.904468, 
                      123.905896, 132.905447, 129.90631, 137.907107, 135.907144, 140.907648, 
                      141.907719, 144.912744, 143.911995, 150.919846, 151.919788, 158.925343, 
                      155.924278, 164.930319, 161.928775, 168.934211, 167.933894, 174.940768, 
                      173.94004, 179.947466, 179.946706, 184.952956, 183.952491, 190.960591, 
                      189.95993, 196.966552, 195.965815, 202.972329, 203.973029, 208.980383, 
                      208.982416, 209.987131, 222.01757, 223.019731, 226.025403, 227.027747, 
                      232.03805, 231.035879, 234.040946)
  
  return((molform * MonoisotopicMasses) %>% sum())
  
}


# Initialize a list of molecular formulas 
MolFormDict <- list(
  "I" = as.molform("C6H13N1O2"), 
  "L" = as.molform("C6H13N1O2"), 
  "K" = as.molform("C6H14N2O2"), 
  "M" = as.molform("C5H11N1O2S1"), 
  "F" = as.molform("C9H11N1O2"), 
  "T" = as.molform("C4H9N1O3"), 
  "W" = as.molform("C11H12N2O2"), 
  "V" = as.molform("C5H11N1O2"), 
  "R" = as.molform("C6H14N4O2"), 
  "H" = as.molform("C6H9N3O2"), 
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
)

#' Add molecular formula objects
#' 
#' @description Add molform objects together
#' 
#' @param ... Any number of molform objects from as.molform. 
#' @param CapNegatives A TRUE/FALSE to indicate whether negative elements should be
#'     capped at 0. Useful for modifications where elements are lost. Default is TRUE.
#' 
#' @examples
#' \dontrun{
#' 
#' add_molforms(as.molform("C6H12O6"), as.molform("H2K2F2Na5"), as.molform("H-1C-2"))
#' 
#' }
#' @export
add_molforms <- function(..., CapNegatives = TRUE) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that all inputs are molecular formula objects
  InputTest <- lapply(list(...), function(x) {inherits(x, "molform")}) %>% unlist() %>% all()
  if (!InputTest) {
    stop("All input objects must be molform objects from as.molform.")
  }
  
  # Check that CapNegatives is a TRUE or FALSE
  if (!is.logical(CapNegatives) | is.na(CapNegatives)) {
    stop("CapNegatives must be a TRUE or FALSE.")
  }
  
  ################
  ## ADD VALUES ##
  ################
  
  # Reduce results 
  Added <- Reduce(`+`, list(...))
  
  # Pull and add all AtomLists
  if (!CapNegatives) {return(Added)} else{
    Added[Added < 0] <- 0
    return(Added)
  }

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
  
  # Check that molform is an appropriate object
  if (!inherits(molform, "molform")) {
    stop("molform should be a molform object from as.molform.")
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
  return(molform * scalar)
  
}

#' Calculate a molecular formula for an amino acid sequence
#' 
#' @param sequence The amino acid sequence. 
#' 
#' @examples 
#' \dontrun{
#' 
#' get_aa_molform("TESTTESTTESTTEST")
#' 
#' }
#' @export
get_aa_molform <- function(sequence) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  if (!is_sequence(sequence)) {
    stop(paste0(sequence, " is not an acceptable amino acid sequence. See ?is_sequence."))
  }
  
  #####################################
  ## CALCULATE THE MOLECULAR FORMULA ##
  #####################################
  
  # Get counts of amino acids 
  AA <- strsplit(sequence, "") %>% unlist()
  
  # Generate a water
  Water <- as.molform("H2O1")
  
  # Add all counts together
  return(Reduce(`+`, MolFormDict[AA]) - Water * (nchar(sequence) - 1))
  
}

# Stress test 
# stress <- rep("ACDEFGHIKLNPQRSTVWYACDEFGHIKLNPQRSTVWY", 100000)
# tictoc::tic(); test <- lapply(stress, BRAIN::getAtomsFromSeq); tictoc::toc() 11.783 sec
# tictoc::tic(); test <- lapply(stress, getAAMolForm); tictoc::toc() 14.865 sec

RelativeAbundances <- data.frame(
  element = c("H", "H", "He", "He", "Li", "Li", "Be", "B", "B", "C", "C", "N", "N", "O", "O", 
              "O", "F", "Ne", "Ne", "Ne", "Na", "Mg", "Mg", "Mg", "Al", "Si", "Si", "Si", "P", 
              "S", "S", "S", "S", "Cl", "Cl", "Ar", "Ar", "Ar", "K", "K", "K", "Ca", "Ca", "Ca", 
              "Ca", "Ca", "Ca", "Sc", "Ti", "Ti", "Ti", "Ti", "Ti", "V", "V", "Cr", "Cr", "Cr", 
              "Cr", "Mn", "Fe", "Fe", "Fe", "Fe", "Co", "Ni", "Ni", "Ni", "Ni", "Ni", "Cu", 
              "Cu", "Zn", "Zn", "Zn", "Zn", "Zn", "Ga", "Ga", "Ge", "Ge", "Ge", "Ge", "Ge", 
              "As", "Se", "Se", "Se", "Se", "Se", "Se", "Br", "Br", "Kr", "Kr", "Kr", "Kr", 
              "Kr", "Kr", "Rb", "Rb", "Sr", "Sr", "Sr", "Sr", "Y", "Zr", "Zr", "Zr", "Zr", "Zr", 
              "Nb", "Mo", "Mo", "Mo", "Mo", "Mo", "Mo", "Mo", "Tc", "Ru", "Ru", "Ru", "Ru", 
              "Ru", "Ru", "Ru", "Rh", "Pd", "Pd", "Pd", "Pd", "Pd", "Pd", "Ag", "Ag", "Cd", 
              "Cd", "Cd", "Cd", "Cd", "Cd", "Cd", "Cd", "In", "In", "Sn", "Sn", "Sn", "Sn", 
              "Sn", "Sn", "Sn", "Sn", "Sn", "Sn", "Sb", "Sb", "Te", "Te", "Te", "Te", "Te", 
              "Te", "Te", "Te", "I", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", "Xe", 
              "Cs", "Ba", "Ba", "Ba", "Ba", "Ba", "Ba", "Ba", "La", "La", "Ce", "Ce", "Ce", 
              "Ce", "Pr", "Nd", "Nd", "Nd", "Nd", "Nd", "Nd", "Nd", "Pm", "Sm", "Sm", "Sm", 
              "Sm", "Sm", "Sm", "Sm", "Eu", "Eu", "Gd", "Gd", "Gd", "Gd", "Gd", "Gd", "Gd", 
              "Tb", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy", "Ho", "Er", "Er", "Er", "Er", 
              "Er", "Er", "Tm", "Yb", "Yb", "Yb", "Yb", "Yb", "Yb", "Yb", "Lu", "Lu", "Hf", 
              "Hf", "Hf", "Hf", "Hf", "Hf", "Ta", "Ta", "W", "W", "W", "W", "W", "Re", "Re", 
              "Os", "Os", "Os", "Os", "Os", "Os", "Os", "Ir", "Ir", "Pt", "Pt", "Pt", "Pt", 
              "Pt", "Pt", "Au", "Hg", "Hg", "Hg", "Hg", "Hg", "Hg", "Hg", "Tl", "Tl", "Pb", 
              "Pb", "Pb", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "U", "U"),
  isotope = c("H1", "H2", "He3", "He4", "Li6", "Li7", "Be9", "B10", "B11", 
              "C12", "C13", "N14", "N15", "O16", "O17", "O18", "F19", "Ne20", 
              "Ne21", "Ne22", "Na23", "Mg24", "Mg25", "Mg26", "Al27", "Si28", 
              "Si29", "Si30", "P31", "S32", "S33", "S34", "S36", "Cl35", "Cl37", 
              "Ar36", "Ar38", "Ar40", "K39", "K40", "K41", "Ca40", "Ca42", 
              "Ca43", "Ca44", "Ca46", "Ca48", "Sc45", "Ti46", "Ti47", "Ti48", 
              "Ti49", "Ti50", "V50", "V51", "Cr50", "Cr52", "Cr53", "Cr54", 
              "Mn55", "Fe54", "Fe56", "Fe57", "Fe58", "Co59", "Ni58", "Ni60", 
              "Ni61", "Ni62", "Ni64", "Cu63", "Cu65", "Zn64", "Zn66", "Zn67", 
              "Zn68", "Zn70", "Ga69", "Ga71", "Ge70", "Ge72", "Ge73", "Ge74", 
              "Ge76", "As75", "Se74", "Se76", "Se77", "Se78", "Se80", "Se82", 
              "Br79", "Br81", "Kr78", "Kr80", "Kr82", "Kr83", "Kr84", "Kr86", 
              "Rb85", "Rb87", "Sr84", "Sr86", "Sr87", "Sr88", "Y89", "Zr90", 
              "Zr91", "Zr92", "Zr94", "Zr96", "Nb93", "Mo92", "Mo94", "Mo95", 
              "Mo96", "Mo97", "Mo98", "Mo100", "Tc98", "Ru96", "Ru98", "Ru99", 
              "Ru100", "Ru101", "Ru102", "Ru104", "Rh103", "Pd102", "Pd104", 
              "Pd105", "Pd106", "Pd108", "Pd110", "Ag107", "Ag109", "Cd106", 
              "Cd108", "Cd110", "Cd111", "Cd112", "Cd113", "Cd114", "Cd116", 
              "In113", "In115", "Sn112", "Sn114", "Sn115", "Sn116", "Sn117", 
              "Sn118", "Sn119", "Sn120", "Sn122", "Sn124", "Sb121", "Sb123", 
              "Te120", "Te122", "Te123", "Te124", "Te125", "Te126", "Te128", 
              "Te130", "I127", "Xe124", "Xe126", "Xe128", "Xe129", "Xe130", 
              "Xe131", "Xe132", "Xe134", "Xe136", "Cs133", "Ba130", "Ba132", 
              "Ba134", "Ba135", "Ba136", "Ba137", "Ba138", "La138", "La139", 
              "Ce136", "Ce138", "Ce140", "Ce142", "Pr141", "Nd142", "Nd143", 
              "Nd144", "Nd145", "Nd146", "Nd148", "Nd150", "Pm145", "Sm144", 
              "Sm147", "Sm148", "Sm149", "Sm150", "Sm152", "Sm154", "Eu151", 
              "Eu153", "Gd152", "Gd154", "Gd155", "Gd156", "Gd157", "Gd158", 
              "Gd160", "Tb159", "Dy156", "Dy158", "Dy160", "Dy161", "Dy162", 
              "Dy163", "Dy164", "Ho165", "Er162", "Er164", "Er166", "Er167", 
              "Er168", "Er170", "Tm169", "Yb168", "Yb170", "Yb171", "Yb172", 
              "Yb173", "Yb174", "Yb176", "Lu175", "Lu176", "Hf174", "Hf176", 
              "Hf177", "Hf178", "Hf179", "Hf180", "Ta180", "Ta181", "W180", 
              "W182", "W183", "W184", "W186", "Re185", "Re187", "Os184", "Os186", 
              "Os187", "Os188", "Os189", "Os190", "Os192", "Ir191", "Ir193", 
              "Pt190", "Pt192", "Pt194", "Pt195", "Pt196", "Pt198", "Au197", 
              "Hg196", "Hg198", "Hg199", "Hg200", "Hg201", "Hg202", "Hg204", 
              "Tl203", "Tl205", "Pb204", "Pb206", "Pb207", "Pb208", "Bi209", 
              "Po209", "At210", "Rn222", "Fr223", "Ra226", "Ac227", "Th232", 
              "Pa231", "U234", "U235", "U238"),
  mass = c(1.007825, 2.014102, 3.016029, 4.002603, 6.015122, 7.016004, 9.012182, 10.012937, 
             11.009305, 12, 13.003355, 14.003074, 15.000109, 15.994915, 16.999132, 17.99916, 
             18.998403, 19.99244, 20.993847, 21.991386, 22.98977, 23.985042, 24.985837, 
             25.982593, 26.981538, 27.976927, 28.976495, 29.97377, 30.973762, 31.972071, 
             32.971458, 33.967867, 35.967081, 34.968853, 36.965903, 35.967546, 37.962732, 
             39.962383, 38.963707, 39.963999, 40.961826, 39.962591, 41.958618, 42.958767, 
             43.955481, 45.953693, 47.952534, 44.95591, 45.952629, 46.951764, 47.947947, 
             48.947871, 49.944792, 49.947163, 50.943964, 49.94605, 51.940512, 52.940654, 
             53.938885, 54.93805, 53.939615, 55.934942, 56.935399, 57.93328, 58.9332, 
             57.935348, 59.930791, 60.93106, 61.928349, 63.92797, 62.929601, 64.927794, 
             63.929147, 65.926037, 66.927131, 67.924848, 69.925325, 68.925581, 70.924705, 
             69.92425, 71.922076, 72.923459, 73.921178, 75.921403, 74.921596, 73.922477, 
             75.919214, 76.919915, 77.91731, 79.916522, 81.9167, 78.918338, 80.916291, 
             77.920386, 79.916378, 81.913485, 82.914136, 83.911507, 85.91061, 84.911789, 
             86.909183, 83.913425, 85.909262, 86.908879, 87.905614, 88.905848, 89.904704, 
             90.905645, 91.90504, 93.906316, 95.908276, 92.906378, 91.90681, 93.905088, 
             94.905841, 95.904679, 96.906021, 97.905408, 99.907477, 97.907216, 95.907598, 
             97.905287, 98.905939, 99.90422, 100.905582, 101.90435, 103.90543, 102.905504, 
             101.905608, 103.904035, 104.905084, 105.903483, 107.903894, 109.905152, 
             106.905093, 108.904756, 105.906458, 107.904183, 109.903006, 110.904182, 
             111.902757, 112.904401, 113.903358, 115.904755, 112.904061, 114.903878, 
             111.904821, 113.902782, 114.903346, 115.901744, 116.902954, 117.901606, 
             118.903309, 119.902197, 121.90344, 123.905275, 120.903818, 122.904216, 
             119.90402, 121.903047, 122.904273, 123.902819, 124.904425, 125.903306, 
             127.904461, 129.906223, 126.904468, 123.905896, 125.904269, 127.90353, 
             128.904779, 129.903508, 130.905082, 131.904154, 133.905395, 135.90722, 
             132.905447, 129.90631, 131.905056, 133.904503, 134.905683, 135.90457, 
             136.905821, 137.905241, 137.907107, 138.906348, 135.907144, 137.905986, 
             139.905434, 141.90924, 140.907648, 141.907719, 142.90981, 143.910083, 
             144.912569, 145.913112, 147.916889, 149.920887, 144.912744, 143.911995, 
             146.914893, 147.914818, 148.91718, 149.917271, 151.919728, 153.922205, 
             150.919846, 152.921226, 151.919788, 153.920862, 154.922619, 155.92212, 
             156.923957, 157.924101, 159.927051, 158.925343, 155.924278, 157.924405, 
             159.925194, 160.92693, 161.926795, 162.928728, 163.929171, 164.930319, 
             161.928775, 163.929197, 165.93029, 166.932045, 167.932368, 169.93546, 
             168.934211, 167.933894, 169.934759, 170.936322, 171.936378, 172.938207, 
             173.938858, 175.942568, 174.940768, 175.942682, 173.94004, 175.941402, 
             176.94322, 177.943698, 178.945815, 179.946549, 179.947466, 180.947996, 
             179.946706, 181.948206, 182.950224, 183.950933, 185.954362, 184.952956, 
             186.955751, 183.952491, 185.953838, 186.955748, 187.955836, 188.958145, 
             189.958445, 191.961479, 190.960591, 192.962924, 189.95993, 191.961035, 
             193.962664, 194.964774, 195.964935, 197.967876, 196.966552, 195.965815, 
             197.966752, 198.968262, 199.968309, 200.970285, 201.970626, 203.973476, 
             202.972329, 204.974412, 203.973029, 205.974449, 206.975881, 207.976636, 
             208.980383, 208.982416, 209.987131, 222.01757, 223.019731, 226.025403, 
             227.027747, 232.03805, 231.035879, 234.040946, 235.043923, 238.050783),
  abundance = c(99.9885, 0.0115, 0.000137, 99.999863, 7.59, 92.41, 100, 19.9, 80.1, 
                98.93, 1.07, 99.632, 0.368, 99.757, 0.038, 0.205, 100, 90.48, 0.27, 
                9.25, 100, 78.99, 10, 11.01, 100, 92.2297, 4.6832, 3.0872, 100, 
                94.93, 0.76, 4.29, 0.02, 75.78, 24.22, 0.3365, 0.0632, 99.6003, 
                93.2581, 0.0117, 6.7302, 96.941, 0.647, 0.135, 2.086, 0.004, 0.187, 
                100, 8.25, 7.44, 73.72, 5.41, 5.18, 0.25, 99.75, 4.345, 83.789, 
                9.501, 2.365, 100, 5.845, 91.754, 2.119, 0.282, 100, 68.0769, 
                26.2231, 1.1399, 3.6345, 0.9256, 69.17, 30.83, 48.63, 27.9, 4.1, 
                18.75, 0.62, 60.108, 39.892, 20.84, 27.54, 7.73, 36.28, 7.61, 100, 
                0.89, 9.37, 7.63, 23.77, 49.61, 8.73, 50.69, 49.31, 0.35, 2.28, 
                11.58, 11.49, 57, 17.3, 72.17, 27.83, 0.56, 9.86, 7, 82.58, 100, 
                51.45, 11.22, 17.15, 17.38, 2.8, 100, 14.84, 9.25, 15.92, 16.68, 
                9.55, 24.13, 9.63, 100, 5.54, 1.87, 12.76, 12.6, 17.06, 31.55, 
                18.62, 100, 1.02, 11.14, 22.33, 27.33, 26.46, 11.72, 51.839, 
                48.161, 1.25, 0.89, 12.49, 12.8, 24.13, 12.22, 28.73, 7.49, 
                4.29, 95.71, 0.97, 0.66, 0.34, 14.54, 7.68, 24.22, 8.59, 
                32.58, 4.63, 5.79, 57.21, 42.79, 0.09, 2.55, 0.89, 4.74, 7.07, 
                18.84, 31.74, 34.08, 100, 0.09, 0.09, 1.92, 26.44, 4.08, 21.18, 
                26.89, 10.44, 8.87, 100, 0.106, 0.101, 2.417, 6.592, 7.854, 
                11.232, 71.698, 0.09, 99.91, 0.185, 0.251, 88.45, 11.114, 100, 
                27.2, 12.2, 23.8, 8.3, 17.2, 5.7, 5.6, 100, 3.07, 14.99, 11.24, 
                13.82, 7.38, 26.75, 22.75, 47.81, 52.19, 0.2, 2.18, 14.8, 20.47, 
                15.65, 24.84, 21.86, 100, 0.06, 0.1, 2.34, 18.91, 25.51, 24.9, 
                28.18, 100, 0.14, 1.61, 33.61, 22.93, 26.78, 14.93, 100, 0.13, 
                3.04, 14.28, 21.83, 16.13, 31.83, 12.76, 97.41, 2.59, 0.16, 5.26, 
                18.6, 27.28, 13.62, 35.08, 0.012, 99.988, 0.12, 26.5, 14.31, 30.64, 
                28.43, 37.4, 62.6, 0.02, 1.59, 1.96, 13.24, 16.15, 26.26, 40.78, 
                37.3, 62.7, 0.014, 0.782, 32.967, 33.832, 25.242, 7.163, 100, 
                0.15, 9.97, 16.87, 23.1, 13.18, 29.86, 6.87, 29.524, 70.476, 
                1.4, 24.1, 22.1, 52.4, 100, 100, 100, 100, 100, 100, 100, 100, 
                100, 0.0055, 0.72, 99.2745) / 100
)

#' Calculate an isotope profile using isopat
#' 
#' @description Generates an isotope profile using [isopat](https://github.com/cran/isopat)
#' 
#' @param molform An object of the as.molform class
#' @param min_abundance Minimum abundance for calculating isotopes. Default is 1. 
#' @param limit See ?isopat::isopattern for more details. 0.1 appears to calculate enough isotopes. 
#' 
#' @examples
#' \dontrun{
#' 
#' calculate_iso_profile(molform = as.molform("C6H12O6"), min_abundance = 1)
#' 
#' }
#' @export
calculate_iso_profile <- function(molform, min_abundance = 1, limit = 0.1) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Molform should be an object of the molform class
  if (!inherits(molform, "molform")) {
    stop("molform should be an object from the molform class.")
  }
  
  # Min abundance should be a numeric between 0 and 100 
  if (!is.numeric(min_abundance) | length(min_abundance) != 1) {
    stop("min_abundance should be a single numeric.")
  }
  if (min_abundance < 0 | min_abundance > 100) {
    stop("The range of min_abundance is between 0 and 100.")
  }
  
  ###############################
  ## CALCULATE ISOTOPE PROFILE ##
  ###############################
  
  # Collapse the molform
  MolForm <- collapse_molform(molform)
  
  # Calculate isotope profile 
  IsoProfile <- isopat::isopattern(RelativeAbundances, MolForm, limit) %>%
    data.table::data.table() %>%
    dplyr::select(mass, abundance) %>%
    dplyr::mutate(
      abundance = abundance / max(abundance) * 100,
      massbin = round(mass)
    ) %>%
    dplyr::group_by(massbin) %>%
    dplyr::filter(abundance == max(abundance) & abundance > min_abundance) %>%
    dplyr::ungroup() %>%
    dplyr::select(-massbin) %>%
    dplyr::mutate(
      isotope = 0:(nrow(.)-1),
      isolabel = paste("M+", isotope, sep = ""),
      isolabel = ifelse(isolabel == "M+0", "M", isolabel)
    )
  
  return(IsoProfile)

}

#' Ensures a string is an acceptable amino acid sequence
#'
#' @description A simple test function to ensure the provided sequence is an
#'    acceptable amino acid sequence
#'
#' @param Sequence An amino acid sequence which should be a single string. Required.
#'
#' @details The output will either be a "TRUE" acceptable sequence, or "FALSE"
#'    unacceptable sequence. An acceptable sequence cannot be NULL,
#'    will not have any spaces or non-letter characters, will be longer than 1
#'    amino acid, and will not contain the letters B, J, O, U, X, or Z. Proforma
#'    notation (i.e. "TES[Acetyl]T") is not accepted, though you can convert the 
#'    proforma string to a sequence with convert_proforma. Modifications should be of the
#'    pspecter_ptm class as defined by the "make_ptm" class and can be applied to
#'    the get_matched_peaks algorithm.
#'
#' @examples
#' \dontrun{
#'
#' # An acceptable sequence will return TRUE
#' is_sequence("TEST")
#'
#' # Unacceptable sequences will return FALSE
#' is_sequence("TESTSEQUENCE")
#' is_sequence("TES[Acetyl]T")
#' is_sequence("TIRED $CIENTIST WORKING L8")
#' is_sequence("T")
#'
#' }
#'
#' @export
is_sequence <- function(Sequence) {
  
  ###########################
  ## RUN THROUGH THE TESTS ##
  ###########################
  
  # If sequence is NULL, return FALSE
  if (is.null(Sequence)) {return(FALSE)}
  
  # If the sequence is not a string or longer than 1 string, throw error
  if (is.character(Sequence) == FALSE || length(Sequence) > 1) {
    stop("The provided Sequence must be a single string. Example: TEST")
  }
  
  # If the sequence has spaces, return FALSE.
  if (grepl("[[:space:]]", Sequence)) {return(FALSE)}
  
  # If the sequence has brackets, return FALSE.
  if (grepl("\\[", Sequence) | grepl("\\]", Sequence)) {return(FALSE)}
  
  # If the sequence has anything besides letters, return FALSE.
  if (grepl("[^aA-zZ]", Sequence)) {return(FALSE)}
  
  # If the number of characters in the sequence is less than 2, return FALSE.
  if (nchar(Sequence) < 2) {return(FALSE)}
  
  # If the sequence contains B, J, O, U, or X, return FALSE.
  if (grepl("[BbJjOoUuXxZz]", Sequence)) {return(FALSE)}
  
  # Unless something else comes up that needs to expand this function, then at
  # this step, it is an acceptable sequence.
  return(TRUE)
  
}

#' Convert a proforma string to a PTM object
#' 
#' @description Converts proforma strings to a PTM object that get_matched_peaks uses
#'     to calculate fragmentation patterns  
#'     
#' @author Degnan, David. Flores, Javier.
#' 
#' @param proforma A string written in the format "M.S[Methyl]S[22]S[23].V"
#' 
#' @examples
#' \dontrun{
#' 
#' convert_proforma("M.(S)[Acetyl]ATNNIAQARKLVEQLRIEAGIERIKVSKAASDLMSYCEQHARNDPLLVGVPASENPFKDK(KPCIIL)[-52.9879].")
#'
#' convert_proforma("M.SS[Methyl]S.V")
#' 
#' convert_proforma("M.S[Methyl]S[22]S[23].V")
#' 
#' convert_proforma("TESTTESTTEST")
#' 
#' }
#' @export
convert_proforma <- function(proforma) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Proteoform should be a string
  if (is.null(proforma) || !is.character(proforma)) {
    stop("proforma must be a vector of characters.")
  }
  
  ###########################
  ## CHECK SIMPLE SEQUENCE ##
  ###########################
  
  # If no brackets, test the sequence is acceptable and then return that sequence 
  if (!grepl("\\[|\\]", proforma)) {
    if (is_sequence(proforma)) {return(proforma)} else {stop("The detected sequence is not acceptable.")}
  }
  
  ###################
  ## LOAD GLOSSARY ##
  ###################
  
  # Load backend glossary
  Glossary <- data.table::fread(
    system.file("extdata", "Unimod_v20220602.csv", package = "ProteoMatch")
  )
  
  #####################################
  ## OTHERWISE, BUILD THE PTM OBJECT ##
  #####################################
  
  ## SEPARATE THE MODIFICATIONS ##
  
  # Grab the data within each of the hard brackets []
  Bracketed_Data <- proforma %>%
    gsub(pattern = "[", replacement = "^[", fixed = T) %>%
    gsub(pattern= "]", replacement = "]^", fixed = T) %>%
    strsplit("^", fixed = T) %>%
    unlist() %>%
    lapply(function(x) {if (grepl("[", x, fixed = T)) {x}}) %>%
    unlist()
  
  # Separate out the modifications from the mass changes
  Modifications <- NULL
  MassChanges <- NULL
  
  # Separate modifications and mass changes
  for (Mod in Bracketed_Data) {
    Mod <- gsub("\\[|\\]", "", Mod)
    StringTest <- suppressWarnings(is.na(as.numeric(Mod)))
    if (StringTest) {Modifications <- c(Modifications, Mod)} else {MassChanges <- c(MassChanges, as.numeric(Mod))}
  }
  
  # Check if modification is in the glossary
  if (!is.null(Modifications)) {
    
    # Ensure all modifications are in the library
    if (all(Modifications %in% Glossary$Modification) == FALSE) {
      
      stop(paste("Modification", Modifications[Modifications %in% Glossary$Modification == FALSE],
                 "is not in our library. If no input error is found, submit an issue",
                 "request on github.")
      )
      
    }
  }
  
  ## CLEAN THE SEQUENCE ##
  
  # First, use the proforma annotation
  Sequence <- proforma
  
  # Then, remove each of the bracketed pieces of information
  for (Mod in Bracketed_Data) {
    Sequence <- gsub(Mod, "", Sequence, fixed = TRUE)
  }
  
  # Then, remove any of the remaining parenthesis
  Sequence <- gsub("\\(|\\)", "", Sequence)
  
  # Get the number of periods in the sequence
  NumPeriods <- stringr::str_count(Sequence, "\\.")
  
  # There shouldn't be more than 2 periods at this point in the pipeline
  if (NumPeriods > 2) {
    stop(paste0("There are too many periods in the input proteoform", Proteoform))
  }
  
  # Split by the periods and take the largest chunk
  SeqSplit <- strsplit(Sequence, "\\.") %>% unlist()
  SeqPosition <- lapply(SeqSplit, nchar) %>% unlist() %>% which.max()
  Sequence <- SeqSplit[SeqPosition]
  
  ## CONVERT SEQUENCE TO POSITIONS ##
  
  SeqWithMods <- proforma
  
  # First, convert the mass changes to a simple "MassChange" 
  for (mod in MassChanges) {
    SeqWithMods <- gsub(paste0("[", mod, "]"), "[MassChange]", SeqWithMods, fixed = T)
  }
  
  # Get sequence cleaned with modifications 
  SeqWithMods <- SeqWithMods %>% 
    strsplit("\\.") %>% 
    unlist() %>% 
    .[SeqPosition] %>% 
    gsub(pattern = "\\(|\\)", replacement = "") 
  
  # Get modification starts and finishes 
  PTM_Object <- data.frame(
    Start = stringr::str_locate_all(SeqWithMods, "\\[")[[1]][,1],
    Finish = stringr::str_locate_all(SeqWithMods, "\\]")[[1]][,1]
  ) %>%
    dplyr::mutate(
      `Mod Number` = 1:nrow(.),
      Span = Finish - Start,
      `Rolling Span` = cumsum(Span),
      `Rolling Span` = ifelse(is.na(`Rolling Span`) & `Mod Number` == 1, Span, `Rolling Span`),
      `N Position` = Finish - `Mod Number` - `Rolling Span`,
      Name = gsub("\\[|\\]", "", Bracketed_Data),
      `AMU Change` = suppressWarnings(as.numeric(Name)),
    ) %>%
    dplyr::select(Name, `AMU Change`, `N Position`) %>%
    dplyr::mutate(
      `AMU Change` = purrr::map2(`AMU Change`, Name, function(`AMU Change`, Name) {
        ifelse(is.na(`AMU Change`), unlist(Glossary[Glossary$Modification == Name, "Mass Change"]), `AMU Change`)
      }) %>% unlist()
    )
  
  rownames(PTM_Object) <- 1:nrow(PTM_Object)
  
  ######################
  ## BUILD THE OBJECT ##
  ######################
  
  # Add proforma string input 
  attr(PTM_Object, "pspecter")$proforma <- proforma
  attr(PTM_Object, "pspecter")$cleaned_sequence <- Sequence
  attr(PTM_Object, "pspecter")$modifications <- Bracketed_Data
  attr(PTM_Object, "pspecter")$mass_changes <- MassChanges
  attr(PTM_Object, "pspecter")$PTMs <- Modifications
  
  # Add class
  class(PTM_Object) <- c(class(PTM_Object), "modifications_pspecter")
  
  # Return results
  return(PTM_Object)
  
}

