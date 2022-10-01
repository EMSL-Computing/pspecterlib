#' Matches Calculated Fragments to the Experimental Spectrum and creates the "Matched Peak" object
#'
#' @description Returns the "matched_peaks" object with experimental spectrum
#'    annotated with matches to calculated fragments.
#'
#' @param ScanMetadata  Object of the scan_metadata class from get_scan_metadata. Required,
#'    unless an alternative spectrum, sequence, and charge is provided (not recommended for beginners).
#' @param PeakData Object of the peak_data class from get_peak_data. Required,
#'    unless an alternative spectrum, sequence, and charge is provided (not recommended for beginners).
#' @param PPMThreshold The ppm error threshold between calculated fragments. Default is 10. Required.
#' @param IonGroups Determine which ion types to calculate. a, b, c, x, y, z are supported. Default
#'    is c("a", "b", "c", "x", "y", "z"). Required.
#' @param CalculateIsotopes A logical which indicates whether isotopes should be calculated.
#'    FALSE = Faster Calculations. Default is TRUE. Required.
#' @param MinAbundance Minimum abundance for calculating isotopes. Default is 0.1.
#' @param CorrelationScore A minimum correlation score to filter isotopes by. Range is 0 to 1.
#'    Default is 0. There is a 3 peak minimum to calculate a correlation score. Required.
#' @param AlternativeIonGroups A "modified_ion" object from "make_mass_modified ions." Default is NULL.
#' @param PTMs A modifications_pspecter object to test modified sequences. Default is NULL.
#'
#' @details
#' The data.table outputted by this function contains 17 columns.
#' \tabular{ll}{
#' M/Z \tab The calculated M/Z value of the fragment \cr
#' \tab \cr
#' Ion \tab The ion's type (a, b, c, x, y, or z) with the ion's position, oriented by terminus: N-terminus (a-c) or C-terminus (x-z) \cr
#' \tab \cr
#' Type \tab The ion's type (a, b, c, x, y, or z) with modified ion annotations (z*) if applicable \cr
#' \tab \cr
#' Position \tab The ion's position, oriented by terminus \cr
#' \tab \cr
#' Z \tab The charge of the fragment \cr
#' \tab \cr
#' Sequence \tab The peptide sequence of the fragment \cr
#' \tab \cr
#' N Position \tab The ion's position, oriented by only the N-terminus \cr
#' \tab \cr
#' General Type \tab The ion's type, removing modified ion annotation (z* would be z) \cr
#' \tab \cr
#' Modifications \tab Any PTMs assigned to this fragment. If none, the string will be "" \cr
#' \tab \cr
#' M/Z Tolerance \tab Based on the inputted PPM Tolerance, this value indicates how far off the calculated M/Z and Experimental M/Z can be. \cr
#' \tab \cr
#' M/Z Experimental \tab The experimental M/Z value that was matched to the calculated M/Z value for that fragment. \cr
#' \tab \cr
#' Intensity Experimental \tab The experimental intensity for the experimental M/Z value \cr
#' \tab \cr
#' PPM Error \tab A calculated value to indicate how far off the experimental and calculated value are from each other, in parts per million (PPM) \cr
#' \tab \cr
#' Molecular Formula \tab The formula of the sequence at that fragment. Used to determine isotopic percentages \cr
#' \tab \cr
#' Isotope \tab Annotation of isotopes in the M+n format, where M+0 is the non-isotope peak, and each successive isotope is M+1, M+2, etc. \cr
#' \tab \cr
#' Isotopic Percentage \tab The calculated intensity, proportional to M+0. For example, if M+1 has an isotopic percentage of 0.5, it is half the size of M+0. \cr
#' \tab \cr
#' Correlation Score \tab If at least two isotopes are recorded, a cosine correlation score of calculated and experimental intensities is determined for these 3+ fragments. \cr
#' \tab \cr
#' Residue \tab Name of the C-terminal residue. Used in downstream functions. \cr
#' \tab \cr
#' }
#'
#'
#' Objects of the class "matched_peaks" contain attributes that are referenced by downstream functions.
#'     The attributes are the input parameters as well as coverage (number of residues with an ion annotation,
#'     over the length of the peptide minus 1), and the median ppm error for all the fragments in
#'     the spectra.
#'
#' @examples
#' \dontrun{
#'
#' # Test bottom up data
#' BU_Peak <- get_peak_data(ScanMetadata = BU_ScanMetadata, ScanNumber = 31728)
#' BU_Match <- get_matched_peaks(ScanMetadata = BU_ScanMetadata, PeakData = BU_Peak)
#'
#' # Test bottom up data with a mass modified ion and a PTM
#' BU_Match2 <- get_matched_peaks(
#'   ScanMetadata = BU_ScanMetadata,
#'   PeakData = BU_Peak,
#'   PTMs = make_ptm(Name = c("Acetyl", "Methyl"),
#'                   AMU_Change = c(42.010565, 14.015650),
#'                   N_Position = c(2, 16),
#'                   Molecular_Formula = list(list("C" = 2, "H" = 2, "O" = 1),
#'                                            list("C" = 1, "H" = 2))),
#'   AlternativeIonGroups = make_mass_modified_ion(Ion = c("b", "z"),
#'                                                 Symbol = c("+", "+"),
#'                                                 AMU_Change = c(1.00727647, 1.00727647))
#' )
#'
#' # Test with top down data
#' TD_Peak <- get_peak_data(ScanMetadata = TD_ScanMetadata, ScanNumber = 5709)
#' TD_Match <- get_matched_peaks(ScanMetadata = TD_ScanMetadata, PeakData = TD_Peak)
#'
#' }
#'
#' @export
get_matched_peaks <- function(ScanMetadata = NULL,
                              PeakData = NULL,
                              PPMThreshold = 10,
                              IonGroups = c("a", "b", "c", "x", "y", "z"),
                              CalculateIsotopes = TRUE,
                              MinimumAbundance = 0.1,
                              CorrelationScore = 0,
                              AlternativeIonGroups = NULL,
                              PTMs = NULL,
                              ...) {

  .get_matched_peaks(
    ScanMetadata = ScanMetadata,
    PeakData = PeakData,
    PPMThreshold = PPMThreshold,
    IonGroups = IonGroups,
    CalculateIsotopes = CalculateIsotopes,
    MinimumAbundance = MinimumAbundance,
    CorrelationScore = CorrelationScore,
    AlternativeIonGroups = AlternativeIonGroups,
    PTMs = PTMs,
    ...
  )

}

.get_matched_peaks <- function(ScanMetadata,
                               PeakData,
                               PPMThreshold,
                               IonGroups,
                               CalculateIsotopes,
                               MinimumAbundance,
                               CorrelationScore,
                               AlternativeIonGroups,
                               PTMs,
                               AlternativeSequence = NULL,
                               AlternativeSpectrum = NULL,
                               AlternativeCharge = NULL,
                               CorrelationScore_FilterNA = FALSE,
                               ChargeThresh = 5,
                               ChargeThresh2 = 10,
                               DebugMode = TRUE) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # ScanMetadata and PeakData may be NULL if Alternative Sequence and Alternative Spectrum
  # and an Alternative Charge are provided
  if (is.null(AlternativeSequence) & is.null(AlternativeSpectrum) & is.null(AlternativeCharge)) {

    # Assert that ScanMetadata is a ScanMetadata object.
    if ("scan_metadata" %in% class(ScanMetadata) == FALSE) {
      stop("ScanMetadata must be a scan_metadata object generated by get_scan_metadata.")
    }

    # Assert that PeakData is a PeakData object
    if ("peak_data" %in% class(PeakData) == FALSE) {
      stop("PeakData must be a peak_data object generated by get_peak_data.")
    }

  }

  # Assert that PPM Threshold is a number
  if (is.numeric(PPMThreshold) == FALSE || length(PPMThreshold) > 1) {
    stop("PPMThreshold must be a single number.")
  }

  # PPM Threshold must be positive to make sense
  PPMThreshold <- abs(PPMThreshold)

  # Assert that ion groups is a vector no longer than length 6
  if (length(IonGroups) > 6) {
    stop("IonGroups cannot be longer than length 6.")
  }

  # Assert that the IonGroups vector contains the letters a-c and x-z only
  if (match(IonGroups, c("a", "b", "c", "x", "y", "z")) %>% is.na() %>% any()) {
    stop("IonGroups must only contain the letters a, b, c, x, y, and z.")
  }

  # Assert that the CalculateIsotopes parameter is a single logical value
  if (is.na(CalculateIsotopes) || is.logical(CalculateIsotopes) == FALSE || length(CalculateIsotopes) > 1) {
    stop("CalculateIsotopes must be a single logical value TRUE or FALSE.")
  }

  # Assert that Minimum Abundance is a single number
  if (is.numeric(MinimumAbundance) == FALSE || length(MinimumAbundance) > 1) {
    stop("MinimumAbundance must be a single numeric value. For example, 0.1.")
  }

  # Convert MinimumAbundance to a positive number
  MinimumAbundance <- abs(MinimumAbundance)

  # Assert that MinimumAbundance is in the range of 0 and 100
  if (MinimumAbundance > 100) {
    stop("MinimumAbundance must be between 0 and 100.")
  }

  # Assert that the Correlation Score is a single number
  if (is.numeric(CorrelationScore) == FALSE || length(CorrelationScore) > 1) {
    stop("CorrelationScore must be a single numeric value. For example, 0.1.")
  }

  # Convert Correlation Score to a positive number
  CorrelationScore <- abs(CorrelationScore)

  # Assert that Correlation Score is in the range of 0 and 1
  if (CorrelationScore > 1) {
    stop("CorrelationScore must be between 0 and 1.")
  }

  # Assert that the alternative sequence is a real sequence
  if (is.null(AlternativeSequence) == FALSE) {

    if (is_sequence(AlternativeSequence) == FALSE) {
      stop("The provided AlternativeSequence is not acceptable. See is_sequence for more details.")
    }

  }

  # Assert that the alternative spectrum is a real spectrum
  if (is.null(AlternativeSpectrum) == FALSE) {

    if ("peak_data" %in% class(AlternativeSpectrum) == FALSE) {
      stop("AlternativeSpectrum must be made with make_peak_data.")
    } else {PeakData <- AlternativeSpectrum}

  }

  # Assert that alternative charge is a single number
  if (is.null(AlternativeCharge) == FALSE) {

    if (length(AlternativeCharge) > 1 | is.numeric(AlternativeCharge) == FALSE) {
      stop("AlternativeCharge must be a single number.")
    }

    # Round to the nearest positive integer
    AlternativeCharge <- AlternativeCharge %>% abs() %>% round()

  }

  # PTMS must be of a defined class
  if (is.null(PTMs) == FALSE) {
    if ("modifications_pspecter" %in% class(PTMs) == FALSE) {
      stop("PTMs must be of the class 'modifications_pspecter' from make_ptm.")
    }
  }

  # Modified ions must be of a defined class
  if (is.null(AlternativeIonGroups) == FALSE) {
    if ("modified_ion" %in% class(AlternativeIonGroups) == FALSE) {
      stop("AlernativeIonGroups must be of the class 'modified_ion' from make_mass_modified_ion.")
    }
  }


  ###################################################################
  ## 0. DEFINE FUNCTION TO REMOVE EXTRANEOUS PEAK MATCHING OPTIONS ##
  ###################################################################
  
  # First, take the minimum PPM spacing in peak data
  PeakSpacing <- (PeakData$`M/Z` - dplyr::lag(PeakData$`M/Z`, default = dplyr::first(PeakData$`M/Z`))) /
    (dplyr::lag(PeakData$`M/Z`, default = dplyr::first(PeakData$`M/Z`))) * 1e6
  PPMBinSize <- min(PeakSpacing[PeakSpacing != 0])

  # Create a function to subset down the matching dataframe
  cleanCalculatedFragments <- function(Fragments) {

    # First, remove all peaks less than the min peak MZ and more than the max peak MZ
    Fragments <- Fragments %>% dplyr::filter(`M/Z` > min(PeakData$`M/Z`) & `M/Z` < max(PeakData$`M/Z`))

    #  Second, remove charge 1 less than the first n positions (ChargeThresh), 
    #  and charge 2 less than the second n positions (ChargeThresh2)
    toRm <- c(which(Fragments$Z > 1 & Fragments$Position <= ChargeThresh),
              which(Fragments$Z > 2 & Fragments$Position <= ChargeThresh2)) %>%
      unique() %>%
      sort()
    if (length(toRm) > 0) {Fragments <- Fragments[-toRm,]}

    
    # First, remove peaks that would never match 
    Fragments <- Fragments %>%
      dplyr::mutate(
        `PPM High` = Fragments$`M/Z` + (PPMThreshold/1e6 * Fragments$`M/Z`),
        `PPM Low` =  Fragments$`M/Z` - (PPMThreshold/1e6 * Fragments$`M/Z`),
        Within = purrr::map2(`PPM Low`, `PPM High`, function(Low, High) {
          nrow(PeakData[PeakData$`M/Z` >= Low & PeakData$`M/Z` <= High,]) != 0
        }) %>% unlist()
      ) %>%
      dplyr::filter(Within == TRUE) %>%
      dplyr::select(-c(`PPM Low`, `PPM High`, Within))
    
    # Second take the minimum charge peak within each ppm bin to prioritize smaller charges. 
    BinVal <- 0 # This is to count bins
    Fragments <- Fragments %>%
      dplyr::arrange(`M/Z`) %>%
      dplyr::mutate(
        PPM = (`M/Z` - dplyr::lag(`M/Z`, default = dplyr::first(`M/Z`))) / (dplyr::lag(`M/Z`, default = dplyr::first(`M/Z`))) * 1e6,
        Flag = PPM < PPMBinSize,
        Bin = lapply(Flag, function(Flag) {
          if (Flag == FALSE) {BinVal <<- BinVal + 1}
          BinVal
        }) %>% unlist()
      ) %>%
      dplyr::group_by(Bin) %>%
      dplyr::slice(which.min(Z)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c(PPM, Flag, Bin))
    
    return(Fragments)
    
  }

  ##########################
  ## 1. DEFINE ALL INPUTS ##
  ##########################

  # Get Scan Number
  ScanNumber <- attr(PeakData, "pspecter")$ScanNumber

  # Get the sequence
  if (is.null(AlternativeSequence)) {
    Sequence <- ScanMetadata[ScanMetadata$`Scan Number` == ScanNumber, "Sequence"] %>% unlist()
  } else {Sequence <- AlternativeSequence}

  # Get the precursor charge
  if (is.null(AlternativeCharge)) {
    PrecursorCharge <- ScanMetadata[ScanMetadata$`Scan Number` == ScanNumber, "Precursor Charge"] %>% unlist()
  } else {PrecursorCharge <- AlternativeCharge}

  #################################
  ## 2. CALCULATE BASE FRAGMENTS ##
  #################################

  if (DebugMode) {message("......Calculating Base Fragments")}

  # Calculate fragment data with MSnbase
  Fragments <- MSnbase::calculateFragments(sequence = Sequence, type = IonGroups,
                                           z = 1:PrecursorCharge) %>% data.table::data.table()

  # Rename Fragments
  colnames(Fragments) <- c("M/Z", "Ion", "Type", "Position", "Z", "Sequence")

  # Exclude N-deamidated and C-dehydrated specific modifications
  Fragments <- dplyr::filter(Fragments, !grepl("[.*_]", Fragments$Type))

  # Label the N-position. Remember that x,y,z fragments are determined from the C-terminus
  Fragments$`N Position` <- ifelse(Fragments$Type %in% c("a", "b", "c"),
                                   Fragments$Position, (nchar(Sequence) + 1) - Fragments$Position)

  ####################################################
  ## 3. ADD IONS WITH MASS MODIFICATIONS (OPTIONAL) ##
  ####################################################

  if (DebugMode) {message("......Adding Ions with Mass Modifications")}

  # Add Mass Modified Ions if they exist
  if (!(is.null(AlternativeIonGroups))) {

    # Iterate through Ion Choices
    NewFragments <- do.call(rbind, lapply(1:nrow(AlternativeIonGroups), function(row) {

      # Get the ions
      getIon <- AlternativeIonGroups$Ion[row]

      # Subset fragments
      subFrag <- Fragments[Fragments$Type == getIon,]

      # Proceed only if there's any fragments
      if (nrow(subFrag) > 0) {

        # Change type
        subFrag$Type <- AlternativeIonGroups$Modified_Ion[row]

        # Change ion
        subFrag$Ion <- paste0(subFrag$Type, subFrag$Position)

        # Change mass
        subFrag$`M/Z` <- subFrag$`M/Z` + AlternativeIonGroups$AMU_Change[row]

        return(subFrag)

      }

    }))

    # Add new fragments
    Fragments <- rbind(Fragments, NewFragments)

  }

  # Add a general type for faster subsetting
  Fragments$`General Type` <- Fragments$Type %>% substr(1, 1)

  ###############################################################
  ## 4. ADD POST-TRANSLATIONAL MODIFICATION WEIGHTS (OPTIONAL) ##
  ###############################################################

  if (DebugMode) {message("......Adding PTMs")}

  # Add blank column to track modifications
  Fragments$Modifications <- ""

  # Apply modifications
  if (is.null(PTMs) == FALSE) {

    # Split fragments by a,b,c ions and x,y,z ions
    FragmentsABC <- Fragments %>% subset(subset = grepl("a|b|c", `General Type`))
    FragmentsXYZ <- Fragments %>% subset(subset = grepl("x|y|z", `General Type`))

    # Apply each PTM
    ApplyPTM <- lapply(1:nrow(PTMs), function(row) {

      # Pull N Position Modification
      NPosition <- PTMs[row, "N Position"] %>% unlist()

      # Apply modifications to ABC MZ
      MZabc <- FragmentsABC[FragmentsABC$`N Position` %in% NPosition:max(FragmentsABC$Position), "M/Z"] %>% unlist()
      Zabc <- FragmentsABC[FragmentsABC$`N Position` %in% NPosition:max(FragmentsABC$Position), "Z"] %>% unlist()
      FragmentsABC[FragmentsABC$`N Position` %in% NPosition:max(FragmentsABC$Position), "M/Z"] <<-
        MZabc + (PTMs[row, "AMU Change"] %>% unlist() / Zabc)
      FragmentsABC[FragmentsABC$`N Position` %in% NPosition:max(FragmentsABC$Position), "Modifications"] <<-
        lapply(FragmentsABC[FragmentsABC$`N Position` %in% NPosition:max(FragmentsABC$Position), "Modifications"], function(Mod) {
          paste(Mod, paste0(PTMs[row, "Name"], "=", PTMs[row, "AMU Change"], "@", substr(Sequence, PTMs[row, "N Position"], PTMs[row, "N Position"]),
                            PTMs[row, "N Position"]))
        }) %>% unlist()

      # Apply modifications to XYZ MZ
      MZxyz <- FragmentsXYZ[FragmentsXYZ$`N Position` %in% 1:NPosition, "M/Z"] %>% unlist()
      Zxyz <- FragmentsXYZ[FragmentsXYZ$`N Position` %in% 1:NPosition, "Z"] %>% unlist()
      FragmentsXYZ[FragmentsXYZ$`N Position` %in% 1:NPosition, "M/Z"] <<-
        MZxyz + (PTMs[row, "AMU Change"] %>% unlist() / Zxyz)
      FragmentsXYZ[FragmentsXYZ$`N Position` %in% 1:NPosition, "Modifications"] <<-
        lapply(FragmentsXYZ[FragmentsXYZ$`N Position` %in% 1:NPosition, "Modifications"], function(Mod) {
          paste(Mod, paste0(PTMs[row, "Name"], "=", PTMs[row, "AMU Change"], "@", substr(Sequence, PTMs[row, "N Position"], PTMs[row, "N Position"]),
                            PTMs[row, "N Position"]))
        }) %>% unlist()

    })
    rm(ApplyPTM)

    # Combine fragments
    Fragments <- rbind(FragmentsABC, FragmentsXYZ)

    # Clean up leading and trailing whitespace
    Fragments$Modifications <- trimws(Fragments$Modifications)

  }

  # Trim down potential fragments to match
  Fragments <- cleanCalculatedFragments(Fragments)

  ###############################
  ## 5. ADD MOLECULAR FORMULAS ##
  ###############################

  if (DebugMode) {message("......Generating Molecular Formulas")}

  # Extract all unique Sequence and Modifications
  MolFormDF <- Fragments %>%
    dplyr::select(Sequence, Modifications) %>%
    unique()
  
  # Remove sequences with a single amino acid
  MolFormDF <- MolFormDF %>% 
    dplyr::mutate(Count = nchar(Sequence) > 1) %>% 
    dplyr::filter(Count) %>% 
    dplyr::select(-Count)

  # Iterate through, getting sequences and modifications and combining them
  MolFormDF$`Molecular Formula` <- lapply(1:nrow(MolFormDF), function(row) {

    # Step one: get sequence and modifications
    Seq <- MolFormDF$Sequence[row]
    Mod <- MolFormDF$Modifications[row]

    # Step two: convert sequence to molecule object
    Atoms <- get_aa_molform(Seq)

    # Step three: add mod if it exists
    if (Mod != "") {

      # Split out modifications
      ModNames <- Mod %>% strsplit("=") %>% lapply(function(x) {head(x, 1)}) %>% unlist()

      # Add to mass
      for (ModName in ModNames) {
        # Add modifications as we do in proteomatch 
        browser()
        Atoms <- add_molforms(Atoms, molform)
      }
    }

    # Step four: add slight changes depending on ion. This was removed to improve speed,
    # since slight changes in Molecular Formula won't drastically change isotopic distributions.
    #AccountForIon <- switch(IonType,
    #                        "a" = make_molecule(list("C" = -1, "H" = -1, "O" = -2)),
    #                        "b" = make_molecule(list("H" = -1, "O" = -1)),
    #                        "c" = make_molecule(list("H" = 2, "N" = 1, "O" = -1)),
    #                        "x" = make_molecule(list("C" = 1, "H" = -1)),
    #                        "y" = make_molecule(list("H" = 1)),
    #                        "z" = make_molecule(list("H" = -2, "N" = -1))
    #)
    #Atoms <- add_molecules(Atoms, AccountForIon)

    return(collapse_molform(Atoms))

  }) %>% unlist()
  
  # Add Molecular Formula
  Fragments <- merge(
    Fragments,
    MolFormDF %>% dplyr::select(Sequence, `Molecular Formula`),
    by = "Sequence"
  )
  
  ###########

  #######################################
  ## 6. CALCULATE ISOTOPES (OPTIONAL)  ##
  #######################################

  if (DebugMode) {message("......Calculating Isotopes")}

  # Add Isotope and Isotopic Percentage
  Fragments$Isotope <- "M"
  Fragments$`Isotopic Percentage` <- NA

  # Calculate Isotopes
  if (CalculateIsotopes) {

    # Get all unique Molecular Formulas
    MolForms <- Fragments$`Molecular Formula` %>% unique()
    
    # Get isotopic distributions
    IsotopeList <- do.call(dplyr::bind_rows, lapply(MolForms, function(MolForm) {
      
      # Get Isotope Relative Abundances
      IsotopeResults <- calculate_iso_profile(as.molform(MolForm), min_abundance = MinimumAbundance)
      
    }))
    
    # Rename the rows of the isotope list
    colnames(IsotopeList) <- c("M/Z", "")

    browser()
    
    # Filter out all monoiosotopic peaks
    IsotopeListFilter <- IsotopeList %>% dplyr::filter(Isotope != "M")

    # If there are Isotopes, apply them
    if (nrow(IsotopeListFilter) > 0) {

      # Get isotope fragments
      IsoFragments <- do.call(rbind, lapply(unique(IsotopeListFilter$`Molecular Formula`), function(Formula) {

        # Subset down to relevant fragments
        subFrag <- Fragments[Fragments$`Molecular Formula` == Formula,]

        # Determine number of isotopes to add
        Isos <- IsotopeListFilter[IsotopeListFilter$`Molecular Formula` == Formula, "Isotope"] %>% unlist()
        IsoNum <- length(Isos)

        do.call(rbind, lapply(1:IsoNum, function(num) {
          subFrag$`M/Z` <- subFrag$`M/Z` + (1.008665 * num / subFrag$Z)
          subFrag$Isotope <- Isos[num]
          return(subFrag)
        }))

      }))

      # Bind to fragments
      Fragments <- rbind(Fragments, IsoFragments)

      # Remove isotopic percentage
      Fragments <- Fragments %>% dplyr::select(-`Isotopic Percentage`)

      # Add actual isotopic percentage values
      Fragments <- merge(Fragments, IsotopeList, by = c("Molecular Formula", "Isotope")) %>%
        dplyr::arrange(`M/Z`)

      # Clean up fragments
      Fragments <- cleanCalculatedFragments(Fragments)
    }

  }

  ########################
  ## 7. MATCH FRAGMENTS ##
  ########################

  if (DebugMode) {message("......Matching Fragments")}

  # Determine the theoertical mz tolerance
  Fragments$`M/Z Tolerance` <- Fragments$`M/Z` * (PPMThreshold / 1e6)

  # For each theoretical peak, find the closest index in ms, where ms = theoretical
  LeftIndex <- findInterval(Fragments$`M/Z`, PeakData$`M/Z`, rightmost.closed = FALSE, all.inside = TRUE)

  # Compute mz differences (absolute) to closest element to each side, smaller to the left and next greater to the right:
  Fragments$`Left Difference` <- abs(PeakData$`M/Z`[LeftIndex] - Fragments$`M/Z`)
  Fragments$`Right Difference` <- abs(PeakData$`M/Z`[LeftIndex + 1] - Fragments$`M/Z`)
  Fragments$`Closest Index` <- LeftIndex

  # Set closest index as right side one, if difference is smaller:
  RightIndexBest <- which(Fragments$`Right Difference` < Fragments$`Left Difference`)
  Fragments$`Closest Index`[RightIndexBest] <- Fragments$`Closest Index`[RightIndexBest] + 1
  Fragments$`M/Z Difference` <- abs(PeakData$`M/Z`[Fragments$`Closest Index`] - Fragments$`M/Z`)

  # Keep only matches within the tolerance
  Fragments <- Fragments[which(Fragments$`M/Z Difference` < Fragments$`M/Z Tolerance`), ]
  Fragments$`M/Z Experimental` <- PeakData$`M/Z`[Fragments$`Closest Index`]
  Fragments$`Intensity Experimental` <- PeakData$Intensity[Fragments$`Closest Index`]

  # Remove non-necessary rows moving forward
  Fragments <- Fragments %>% dplyr::select(-c(`Left Difference`, `Right Difference`, `Closest Index`, `M/Z Difference`))

  # Calculate PPM Error
  Fragments$`PPM Error` <- ((Fragments$`M/Z Experimental` - Fragments$`M/Z`) / Fragments$`M/Z`) * 1e6

  ###############################################
  ## 8. CALCULATE COSINE SIMILARITY (OPTIONAL) ##
  ###############################################

  if (DebugMode) {message("......Calculating Cosine Similarity")}

  # Assign a blank score to all Fragments
  Fragments$`Correlation Score` <- NA

  # Calculate Cosine Similarity if Isotopes Enabled
  if (CalculateIsotopes) {

    # Create an identifier for each family of isotopes
    Fragments <- Fragments %>% dplyr::mutate(ID = paste(Ion, Z))

    # Get correlation score
    CS <- Fragments %>%
      dplyr::group_by(ID) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        `Correlation Score` = purrr::map(data, function(x) {
          if (nrow(x) < 3) {return(NA)} else {
            lsa::cosine(
              x = x$`Isotopic Percentage` %>% unlist(),
              y = x$`Intensity Experimental` %>% unlist()
            )[1,1]
          }
        }) %>% unlist()
      ) %>%
      dplyr::select(-data)

    # Merge to main dataframe
    Fragments <- Fragments %>% dplyr::select(-`Correlation Score`)
    Fragments <- merge(Fragments, CS, by = "ID")

    # Separate no score from correaltion scores
    NAFragments <- Fragments[is.na(Fragments$`Correlation Score`),]
    CorrFragments <- Fragments %>% dplyr::filter(`Correlation Score` >= CorrelationScore)

    # Filtering NA correlation score values
    if (CorrelationScore_FilterNA) {
      Fragments <- CorrFragments
    } else {
      Fragments <- rbind(NAFragments, CorrFragments)
    }

    # Remove ID
    Fragments <- Fragments %>% dplyr::select(-ID)

  }

  # If all the fragments were removed, warn user
  if (nrow(Fragments) == 0) {
    message("All fragments removed with the Correlation Score filter.")
    return(NULL)
  }

  ################################
  ## 8. ADD PEPTIDE POSITIONING ##
  ################################

  # Determine the peptide residue and position for each ion
  Fragments$Residue <- lapply(1:nrow(Fragments), function(row) {
    Seq <- Fragments[row, "Sequence"] %>% unlist()
    if (Fragments[row, "General Type"] %>% unlist() %in% c("a", "b", "c")) {
      substr(Seq, start = nchar(Seq), stop = nchar(Seq))
    } else {substr(Seq, start = 1, stop = 1)}
  }) %>%
    unlist() %>%
    paste0(Fragments$`N Position`)

  ##################
  ## BUILD OBJECT ##
  ##################

  # Reorder Fragments
  Fragments <- Fragments[,c("PPM Error", "Ion", "Z", "Isotope", "M/Z", "M/Z Experimental", "M/Z Tolerance",
               "Isotopic Percentage", "Intensity Experimental", "Correlation Score",
               "Type", "General Type", "Modifications", "Molecular Formula",
               "Position", "N Position", "Residue", "Sequence")] %>%
    dplyr::arrange(`M/Z`)

  attr(Fragments, "pspecter")$Coverage <- ((Fragments$`N Position` %>% unique() %>% .[.!= 1] %>% length()) / (nchar(Sequence) - 1) * 100) %>% round(2) %>% paste0("%")
  attr(Fragments, "pspecter")$MedianPPMError <- median(Fragments$`PPM Error`)
  attr(Fragments, "pspecter")$Sequence <- Sequence

  # If ScanMetadata is not NULL, add attributes
  if (is.null(ScanMetadata) == FALSE) {
    attr(Fragments, "pspecter")$MSPath <- attr(ScanMetadata, "pspecter")$MSPath
    attr(Fragments, "pspecter")$IDPath <- attr(ScanMetadata, "pspecter")$IDPath
  }

  attr(Fragments, "pspecter")$PPMThreshold <- PPMThreshold
  attr(Fragments, "pspecter")$IonGroups <- IonGroups
  attr(Fragments, "pspecter")$IsotopesIncluded <- CalculateIsotopes
  attr(Fragments, "pspecter")$IsotopicPercentageFilter <- IsotopicPercentage
  attr(Fragments, "pspecter")$CorrelationScoreFilter <- CorrelationScore

  # Add mass modified ions and PTMs if they are not NULL
  if (is.null(AlternativeIonGroups) == FALSE) {
    attr(Fragments, "pspecter")$AlternativeIonGroups <- AlternativeIonGroups
  }

  if (is.null(PTMs) == FALSE) {
    attr(Fragments, "pspecter")$PTMs = PTMs
  }

  # Add the matched peaks class
  class(Fragments) <- c(class(Fragments), "matched_peaks")

  return(Fragments)

}

