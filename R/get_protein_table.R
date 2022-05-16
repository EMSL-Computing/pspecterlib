#' Counts the Proteins in ScanMetadata and generate "Protein Table" object
#'
#' @description A simple function which returns a table of protein counts from ScanMetadata
#'
#' @param ScanMetadata Object of the scan_metadata class from get_scan_metadata. Required.
#' @param FASTAPath File path to the database file, in FASTA format. Required.
#' @param QValueMaximum A maximum q-value to filter proteins on. The values range
#'    from 0-1, and NULL means no filter is applied. Default is NULL.
#' @param ScoreMaximum A maximum score to filter proteins on. The values range
#'    from 0-1. Remember most scores are small (i.e. 1e-20 or smaller). NULL means no filter is applied.
#'    Default is NULL.
#' @param RemoveContaminants A True/False to indicate whether contaminants will be removed
#'    or not. Default is TRUE.
#'
#' @details
#' The data.table outputted by this function contains 5 columns of data from the ID file.
#' \tabular{ll}{
#' Protein \tab The name of the protein (protein ID) from the ID file \cr
#' \tab \cr
#' Number of Peptides \tab The count of peptides identified for that protein. \cr
#' \tab \cr
#' Median Q Value \tab The median Q Value for the peptides identified for that protein. \cr
#' \tab \cr
#' Median Score \tab The median score for the peptides identified for that protein. \cr
#' \tab \cr
#' Description \tab The protein's description as provided in the ID file. \cr
#' \tab \cr
#' }
#'
#' Objects of the class "protein_table" contain attributes of the input variables: QValueMaximum, ScoreMaximum, and RemoveContaminants
#'
#' @examples
#' \dontrun{
#'
#' # Test bottom up data. Note that for these example files, the same FASTA file
#' # has been used for both bottom-up and top-down
#' tmpdir <- tempdir()
#' download.file("https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/QC_Shew.fasta",
#'               file.path(tmpdir, "QC_Shew.fasta"))
#'
#' ProteinTable1 <- get_protein_table(BU_ScanMetadata, FASTAPath = file.path(tmpdir, "QC_Shew.fasta"))
#' ProteinTable2 <- get_protein_table(BU_ScanMetadata, FASTAPath = file.path(tmpdir, "QC_Shew.fasta"), QValueMaximum = 0.1, ScoreMaximum = 0.00003)
#'
#' }
#'
#'
#' @export
get_protein_table <- function(ScanMetadata,
                              FASTAPath,
                              QValueMaximum = NULL,
                              ScoreMaximum = NULL,
                              RemoveContaminants = TRUE) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Check object
  if ("scan_metadata" %in% class(ScanMetadata) == FALSE) {
    stop("ScanMetadata must be a scan_metadata object from get_scan_metadata.")
  }

  # Assert that the FASTA file exists
  if (file.exists(FASTAPath) == FALSE) {stop("FASTA file path must exist.")}

  # Assert that the MS file extension is mzML, mzXML, or raw
  if (grepl(".FASTA|.FA|.fasta|.fa", FASTAPath) == FALSE) {stop("FASTA file must have the .FASTA or .FA extension.")}

  # Check for ID Data
  if (is.null(attr(ScanMetadata, "pspecter")$IDPath)) {
    stop("No ID data included. Please add identification data to the ScanMetadata object.")
  }

  # Assert that QValueMaximum is a single value between 0-1
  if (is.null(QValueMaximum) == FALSE) {
    if (!is.numeric(QValueMaximum) || length(QValueMaximum) > 1 || QValueMaximum < 0 || QValueMaximum > 1) {
      stop("QValueMaximum must be a single value between 0-1.")
    }
  }

  # Assert that ScoreMaximum is a single value between 0-1
  if (is.null(ScoreMaximum) == FALSE) {
    if (!is.numeric(ScoreMaximum) || length(ScoreMaximum) > 1 || ScoreMaximum < 0 || ScoreMaximum > 1) {
      stop("ScoreMaximum must be a single value between 0-1.")
    }
  }

  # Remove Contaminants must be a true or false
  if (is.na(RemoveContaminants) || is.logical(RemoveContaminants) == FALSE || length(RemoveContaminants) > 1) {
    stop("RemoveContaminants must be a TRUE or FALSE.")
  }

  #####################
  ## READ FASTA FILE ##
  #####################

  # Read input FASTA file
  FASTA <- seqinr::read.fasta(FASTAPath, seqtype = "AA", as.string = T)

  #####################
  ## MAKE DATA TABLE ##
  #####################

  # Extract relevant information
  ProteinData <- ScanMetadata %>% as.data.frame() %>% dplyr::select(c(`Protein ID`, Description, `Q Value`, Score))

  # Remove Protein Data
  ProteinData <- ProteinData[is.na(ProteinData$`Protein ID`) == FALSE,]

  # Remove contaminants if TRUE
  if (RemoveContaminants) {
    ProteinData <- ProteinData %>% dplyr::filter(grepl("Contaminant", `Protein ID`) == FALSE)
  }

  # Determine whether filters can be applied given a Q Value and a Score
  if (is.na(ProteinData$`Q Value`) %>% any()) {
    message("NA values found in Q Values. Filter Ignored.")
    QValueMaximum <- NULL
  }

  if (is.na(ProteinData$Score) %>% any()) {
    message("NA values found in Score. Filter Ignored.")
    ScoreMaximum <- NULL
  }

  # Apply filters that are not NULL
  if (is.null(QValueMaximum) == FALSE) {
    ProteinData <- ProteinData %>% dplyr::filter(`Q Value` <= QValueMaximum)
  }

  if (is.null(ScoreMaximum) == FALSE) {
    ProteinData <- ProteinData %>% dplyr::filter(Score <= ScoreMaximum)
  }

  # Rename the first column
  colnames(ProteinData)[1] <- "Protein"

  # Get counts and average score per protein
  ProteinTable <- ProteinData %>%
    dplyr::group_by(Protein) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      "Number of Peptides" = purrr::map(data, function(x) {nrow(x)}) %>% unlist(),
      "Median Q Value" = purrr::map(data, function(x) {x$`Q Value` %>% median(na.rm = T)}) %>% unlist(),
      "Median Score" = purrr::map(data, function(x) {x$Score %>% median(na.rm = T)}) %>% unlist(),
      "Description" = purrr::map(data, function(x) {x$Description[1]}) %>% unlist(),
      "Literature Sequence" = purrr::map(Protein, function(x) {FASTA[[Protein]][[1]]})
    ) %>%
    dplyr::select(-data) %>%
    dplyr::arrange(desc(`Number of Peptides`)) %>%
    data.table::data.table()

  ###################
  ## STORE RESULTS ##
  ###################

  # Create the final table
  FinalTable <- ProteinTable %>% dplyr::select(-`Literature Sequence`)

  # Add hidden attributes of the user inputs
  attr(FinalTable, "pspecter")$QValueMaximum <- ifelse(is.null(QValueMaximum), 1, QValueMaximum)
  attr(FinalTable, "pspecter")$ScoreMaximum <- ifelse(is.null(ScoreMaximum), 1, ScoreMaximum)
  attr(FinalTable, "pspecter")$RemoveContaminants <- RemoveContaminants
  attr(FinalTable, "pspecter")$LitSeq <- ProteinTable %>% dplyr::select(c(Protein, `Literature Sequence`))

  # Store the class information of this object
  class(FinalTable) <- c(class(FinalTable), "protein_table")

  # Return object
  return(FinalTable)

}
