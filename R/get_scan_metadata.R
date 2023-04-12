#' Creates the "Scan Metadata" object from an MS file
#'
#' @description Aligns both MS and ID metadata per scan number. MS data is required. ID data is optional.
#'
#' @param MSPath Data path to the MS file. Currently, mzML, mzXML, and raw are supported.
#' @param IDPath Data path to the ID file. Currently, mzID is supported.
#'
#' @details
#' The data.table outputted by this function contains 16 main columns of data from both the MS and ID file.
#' \tabular{ll}{
#' Scan Number \tab The number of the MS scan where scan 1 was taken at the first retention time. \cr
#' \tab \cr
#' MS Level \tab Indiciates whether the scan is a MS1 (precursor) or MS2. \cr
#' \tab \cr
#' Retention Time \tab The time, in minutes, that the fragment was eluted. \cr
#' \tab \cr
#' Precursor M/Z \tab The M/Z value selected from the MS1 (precursor) scan to fragment for the MS2 scan. Will be NA for MS1 scans. \cr
#' \tab \cr
#' Precursor Charge \tab The Z value, charge, from the selected MS1 (precursor) scan. Will be NA for MS1 scans. \cr
#' \tab \cr
#' Precursor Scan \tab The MS1 scan number that the MS2 scan was fractioned from. Will be NA for MS1 scans. \cr
#' \tab \cr
#' Activation Method \tab A string to indicate which activation (dissociation) method was used to fragment the proteins. \cr
#' \tab \cr
#' Sequence \tab The amino acid sequence  (fragment) assigned to the MS2 scan. Only included if IDPath is supplied. \cr
#' \tab \cr
#' Protein ID \tab The name of the protein that the amino acid sequence is from. Only included if the IDPath is supplied. \cr
#' \tab \cr
#' Calculated Mass \tab The calculated, or expected, mass of the peptide fragment. Only included if the IDPath is supplied. \cr
#' \tab \cr
#' Experimental Mass \tab The experimental, or actual, mass that was matched in the MS2 file. Only included if the IDPath is supplied. \cr
#' \tab \cr
#' Score \tab A measure of how well the calculated and experimental spectra match. A lower score is better. Only included if the IDPath is supplied. \cr
#' \tab \cr
#' Q Value \tab An adjusted p-value that also measures how well the calculated and experimental spectra match. A lower score is preferred. Only included if the IDPath is supplied. \cr
#' \tab \cr
#' Decoy \tab A logical (TRUE/FALSE/NA) indicating whether the identified protein is a protein decoy or not. \cr
#' \tab \cr
#' Modifications \tab A string indicating whether a modification was identified or not. The format is "Name@@Position=MW & Name2@@Position2=MW2" etc \cr
#' \tab \cr
#' }
#'
#'
#' Objects of the class "scan_metadata" contain attributes that are referenced by downstream functions.
#'    The attributes are as follows:
#' \tabular{ll}{
#' MSPath \tab The path the user supplied in the MSPath parameter. \cr
#' \tab \cr
#' IDPath \tab The path the user supplied in the IDPath parameter. Default is NULL. \cr
#' \tab \cr
#' MSFileType \tab A flag to indicate whether the inputted MS file is XML-based (will be read with the mzR package) or RAW (will be read with the rawrr package). \cr
#' \tab \cr
#' MSnbase \tab An object from the MSnbase package used to pull XIC data when the MSFileType is XML. If not, then this attribute will be NULL. \cr
#' \tab \cr
#' MS1Scans \tab An integer vector listing all the MS1 scans in order. Used by the MS1Plots function. \cr
#' \tab \cr
#' }
#'
#' @examples
#' \dontrun{
#'
#' # Create a temporary directory and copy example data there
#' tmpdir <- tempdir()
#'
#' files <- c(
#'  "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzML",
#'  "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzid",
#'  "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown.mzML",
#'  "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown.mzid"
#' )
#'
#' for (file in files) {
#'   download.file(file, file.path(tmpdir, tail(unlist(strsplit(file, "/")), 1)))
#' }
#'
#' # Test with bottom up data
#' BU_ScanMetadata <- get_scan_metadata(MSPath = file.path(tmpdir, "BottomUp.mzML"),
#'                                      IDPath = file.path(tmpdir, "BottomUp.mzid"))
#'
#' # Test with top down data
#' TD_ScanMetadata <- get_scan_metadata(MSPath = file.path(tmpdir, "TopDown.mzML"),
#'                                      IDPath = file.path(tmpdir, "TopDown.mzid"))
#'
#' # Test with raw top down data
#' library(rawrr)
#' rawfile <- file.path(path.package(package = 'rawrr'), 'extdata', 'sample.raw')
#' RAW_ScanMetadata <- get_scan_metadata(MSPath = rawfile)
#'
#' }
#'
#' @export
get_scan_metadata <- function(MSPath,
                              IDPath = NULL) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Assert that the MS file exists
  if (file.exists(MSPath) == FALSE) {stop("MS file path must exist.")}

  # Assert that the MS file extension is mzML, mzXML, or raw
  if (grepl(".mzML|.mzXML|.mzml|.mzxml|.raw", MSPath) == FALSE) {stop("MS file must be a mzML, mzXML, or raw")}

  # If ID Path is not NULL, check it
  if (is.null(IDPath) == FALSE) {

    # Assert that the ID file exists
    if (file.exists(IDPath) == FALSE) {stop("ID file path must exist.")}

    # The only acceptable ID file extension is mzID
    if (grepl(".mzid|.mzID", IDPath) == FALSE) {stop("ID file must be an mzID file.")}

  }

  ############################
  ## DEFINE LOCAL FUNCTIONS ##
  ############################

  # Pull the activation method from the filter string
  getActMet <- function(filterStringVector) {
    filterStringVector %>%
      lapply(function(filter) {
        ifelse(grepl("@", filter), toupper(gsub("[^A-z]|\\[|\\]", "", unlist(strsplit(filter, "@"))[2])), NA)
      }) %>%
      unlist()
  }

  ##################
  ## PULL MS DATA ##
  ##################

  # Set MS File Type to NULL
  MSFileType <- NULL

  # Read mzML or mzXML file
  if (grepl(".mzML|.mzXML|.mzml|.mzxml", MSPath)) {

    # Set the MS File Type to XML
    MSFileType <- "XML"

    # Use mzR to read the file
    mzRpwiz <- mzR::openMSfile(MSPath, backend = "pwiz")

    # Read scan data
    ScanMetadata <- mzRpwiz %>%
      mzR::header() %>%
      data.table::data.table()

    # Scale retention time to minutes
    ScanMetadata$retentionTime <- ScanMetadata$retentionTime / 60

    # If there are MS1 scans in the file, change the precursor charge to 1
    if (1 %in% ScanMetadata$msLevel %>% unique()) {
      ScanMetadata[ScanMetadata$msLevel == 1, "precursorCharge"] <- 1
    }

    # Pull the activation method
    ScanMetadata$`Activation Method` <- getActMet(ScanMetadata$filterString)

    # Rename columns
    ScanMetadata <- data.table::data.table(
      "Scan Number" = ScanMetadata$acquisitionNum,
      "MS Level" = ScanMetadata$msLevel,
      "Retention Time" = round(ScanMetadata$retentionTime, 4),
      "Precursor M/Z" = round(ScanMetadata$precursorMZ, 4),
      "Precursor Charge" = ScanMetadata$precursorCharge,
      "Precursor Scan" = ScanMetadata$precursorScanNum,
      "Activation Method" = ScanMetadata$`Activation Method`
    )

    # Read raw file
  } else if (grepl(".raw", MSPath)) {

    # Set the MS File Type to RAW
    MSFileType <- "RAW"

    # Use rawDiag to read the file
    ScanMetadata <- rawrr::readIndex(MSPath) %>% data.table::data.table()

    # Change MS levels to numeric
    ScanMetadata$MSOrder <- data.table::fifelse(ScanMetadata$MSOrder == "Ms", 1, 2)

    # Change Precursor Charge to 1
    ScanMetadata$charge <- data.table::fifelse(ScanMetadata$charge == 0, 1, ScanMetadata$charge)

    # Change Precursor Mass to NA when MS Level is 1
    if (1 %in% ScanMetadata$MSOrder %>% unique()) {
      ScanMetadata[ScanMetadata$MSOrder == 1, "precursorMass"] <- NA
    }

    # Pull the activation method
    ScanMetadata$`Activation Method` <- getActMet(ScanMetadata$scanType)

    # Rename columns
    ScanMetadata <- data.table::data.table(
      "Scan Number" = ScanMetadata$scan,
      "MS Level" = ScanMetadata$MSOrder,
      "Retention Time" = round(ScanMetadata$rtinseconds / 60, 4),
      "Precursor M/Z" = round(ScanMetadata$precursorMass, 4),
      "Precursor Charge" = ScanMetadata$charge,
      "Precursor Scan" = ScanMetadata$masterScan,
      "Activation Method" = ScanMetadata$`Activation Method`
    )

  }

  ##################
  ## PULL ID DATA ##
  ##################

  # Fill in NAs for the rest of the dataframe if no ID file is provided or if the file
  # does not exist or is not an mzid.
  if (is.null(IDPath) || file.exists(IDPath) == FALSE || grepl(".mzid", IDPath) == F) {

    # Add empty ID columns to ScanMetadata
    ScanMetadata$Sequence <- ScanMetadata$`Protein ID` <- ScanMetadata$`Calculated Mass` <- NA
    ScanMetadata$`Experimental Mass` <- ScanMetadata$Score <- NA
    ScanMetadata$`Q Value` <- ScanMetadata$Decoy <- ScanMetadata$Description <- NA
    ScanMetadata$`Peptide Start Position` <- NA
    ScanMetadata$Modifications <- NA

    # Order ScanMetadata with MS2 first, followed by MS1
    ScanMetadata <- ScanMetadata[order(-ScanMetadata$`MS Level`),]
    ScanMetadata$Order <- 1:nrow(ScanMetadata)

  } else {

    # Pull ID data for both Protein Spectra Matches (PSMS), and score data
    ID <- mzR::openIDfile(IDPath)
    PSMS <- ID %>% mzR::psms() %>% data.table::data.table()
    Score <- ID %>% mzR::score() %>% data.table::data.table()
    MergedGroup <- dplyr::left_join(PSMS, Score, by = "spectrumID") %>% unique()

    # If modifications exist, convert amino acid sequences to ProForma strings
    if (nrow(mzR::modifications(ID)) > 0) {
      
      # Pull and format modifications
      MOD <- mzR::modifications(ID) %>%
        dplyr::mutate(
          spectrumID = spectrumID, 
          Proteoform = gsub("Pep_|[0-9]", "", peptideRef), # Remove extra text and numerics from peptide references
          Name = lapply(1:length(spectrumID), function(pos) {
            theName <- .$name[pos]
            ifelse(is.null(theName) || is.na(theName) || theName != "", theName, .$mass[pos])
          }) %>% unlist()
        ) %>%
        dplyr::select(spectrumID, Proteoform, Name) %>%
        unique() %>%
        dplyr::mutate(
          Proteoform = lapply(1:length(spectrumID), function(pos) {
            fix <- gsub("+", paste0("[", .$Name[pos], "]"), .$Proteoform[pos], fixed = T)
            fix <- gsub("-", paste0("[", .$Name[pos], "]"), fix, fixed = T)
            return(fix)
          }) %>% unlist()
        )  %>%
        dplyr::select(-Name)
      
      # Merged MODs and replace sequences with non NA Proteoforms
      MergedGroup <- dplyr::left_join(MergedGroup, MOD, by = "spectrumID") %>% 
        unique() %>%
        dplyr::mutate(sequence = ifelse(is.na(Proteoform), sequence, Proteoform))
      
    } 
    
    
    
    # Try to pull Score Data, and if it doesn't work, extract the text. Note: will add back
    # in when I have a tangible, shareable example.
    #tryCatch({
    
    #},

    #error = function(e) {

      # Read mzid file as a text file
    #  suppressWarnings({textID <- readLines(IDPath)})

      # Grab Spec E Value and Q Value
    # EV <- QV <- c()
    # for (line in textID) {
    #   if (grepl("SpecEValue", line)) {
    #     splitLine <- unlist(strsplit(line, " |=|\""))
    #     EV <- c(EV, as.numeric(splitLine[grep("value", splitLine) + 2]))
    #   } else if (grepl(":QValue", line)) {
    #     splitLine <- unlist(strsplit(line, " |=|\""))
    #     QV <- c(QV, as.numeric(splitLine[grep("value", splitLine) + 2]))
    #   }
    # }

    # # Make Score DataTable
    # Score <<- data.table::data.table(
    #   "MS.GF.SpecEValue" = EV %>% as.numeric(),
    #   "MS.GF.QValue" = QV
    # )

    #})

    # Q Value should be converted to NA if NULL
    if (is.null(MergedGroup$MS.GF.QValue)) {MergedGroup$MS.GF.QValue <- NA}

    # Rename columns
    IDData <- data.table::data.table(
      "Scan Number" = MergedGroup$acquisitionNum,
      "Sequence" = MergedGroup$sequence,
      "Protein ID" = MergedGroup$DatabaseAccess,
      "Calculated Mass" = round(MergedGroup$experimentalMassToCharge, 4),
      "Experimental Mass" = round(MergedGroup$calculatedMassToCharge, 4),
      "Score" = MergedGroup$MS.GF.SpecEValue,
      "Q Value" = round(as.numeric(MergedGroup$MS.GF.QValue), 2),
      "Decoy" = MergedGroup$isDecoy,
      "Description" = MergedGroup$DatabaseDescription,
      "Peptide Start Position" = MergedGroup$start
    )

    # Merge ScanMetdata and IDData
    ScanMetadata <- merge(ScanMetadata, IDData, by = "Scan Number", all = TRUE) %>% unique()

    # Order by Score
    ScanMetadata <- ScanMetadata[order(ScanMetadata$Score),]
    ScanMetadata$Order <- 1:nrow(ScanMetadata)

  }

  ##################
  ## BUILD OBJECT ##
  ##################

  # Add hidden attributes that the users don't need to know about: MSPath, IDPath,
  # MSFileType, mzRpwiz (object from mzR), MSnbase (object from MSnbase), MS1Scans.
  attr(ScanMetadata, "pspecter") = list(
    "MSPath" = MSPath,
    "IDPath" = IDPath,
    "MSFileType" = MSFileType,
    "mzRpwiz" = NULL,
    "MSnbase" = NULL,
    "MS1Scans" = ScanMetadata[ScanMetadata$`MS Level` == 1, "Scan Number"]
  )

  # If file type is XML, then store the mzRpwiz and MSnbase objects
  if (MSFileType == "XML") {
    attr(ScanMetadata, "pspecter")$mzRpwiz = mzRpwiz
    attr(ScanMetadata, "pspecter")$MSnbase = MSnbase::readMSData(attr(ScanMetadata, "pspecter")$MSPath, mode = "onDisk")
  }

  # Store the class information of this object
  class(ScanMetadata) <- c(class(ScanMetadata), "scan_metadata")

  return(ScanMetadata)
}
