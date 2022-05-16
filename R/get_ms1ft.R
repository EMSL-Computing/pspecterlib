#' Generates the "ms1ft" object, which is Top-Down MS1 Feature data
#'
#' @description Reads in ms1ft from the ProMex algorithm in InformedProteomics, and
#'    associates features to protein data from the Ic_Targets file from MSPathFinderT.
#'    All of this data is used to generate the promex feature plot.
#'
#' @param MS1FTPath Path to the MS1FT txt file, which is a tsv. Required.
#' @param TargetsPath Optional path to an IC targets file which links MS1FT features to proteins. Default is NULL.
#' @param TrimColumns A logical to indicate whether all columns or only columns relevant to
#'     promex_feature_plot should be returned. Default is TRUE.
#'
#' @details
#' The data.table outputted by this function contains many columns. Listed below are the main ones.
#' \tabular{ll}{
#' FeatureID \tab An identified mass across elution (retention) times is given a feature ID, typically a single number. \cr
#' \tab \cr
#' MinElutionTime \tab The smallest retention time of the feature \cr
#' \tab \cr
#' MaxElutionTime \tab The largest retention time of the feature \cr
#' \tab \cr
#' MonoMass \tab The monoisotopic mass of the feature. \cr
#' \tab \cr
#' Abundance \tab The relative abundance of the feature. \cr
#' \tab \cr
#' ProteinName \tab If the Targets file is provided, the name of the protein associated to that feature. \cr
#' \tab \cr
#' }
#'
#' @examples
#' \dontrun{
#'
#' # Test top down data
#' tmpdir <- tempdir()
#' download.file("https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown.ms1ft",
#'               file.path(tmpdir, "TopDown.ms1ft"))
#' download.file("https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown_IcTarget.tsv",
#'               file.path(tmpdir, "TopDown_IcTarget.tsv"))
#'
#' # Generate object
#' MS1FT <- get_ms1ft(MS1FTPath = file.path(tmpdir, "TopDown.ms1ft"),
#'                    TargetsPath = file.path(tmpdir, "TopDown_IcTarget.tsv"))
#'
#' }
#'
#' @export
get_ms1ft <- function(MS1FTPath,
                      TargetsPath = NULL,
                      TrimColumns = TRUE) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Assert that the MS1FT File Exists
  if (file.exists(MS1FTPath) == FALSE) {stop("MS1FT file path must exist.")}

  # Assert that the MS1FT file extension is ms1ft
  if (grepl(".ms1ft", MS1FTPath) == FALSE) {stop("MS1FT file must have the .ms1ft extension.")}

  # If Targets Path is not NULL, check it
  if (is.null(TargetsPath) == FALSE) {

    # Assert that the Targets file exists
    if (file.exists(TargetsPath) == FALSE) {stop("Targets file path must exist.")}

    # The only acceptable Targets file extension is tsv
    if (grepl(".tsv", TargetsPath) == FALSE) {stop("TargetsPath must be a tsv.")}

  }

  # Trim columns must be a true or false
  if (is.na(TrimColumns) || is.logical(TrimColumns) == FALSE || length(TrimColumns) > 1) {
    stop("TrimColumns must be a TRUE or FALSE.")
  }

  ###############################
  ## GENERATE MS1FT DATA TABLE ##
  ###############################

  # Read MS1FT File
  MS1FT <- data.table::fread(MS1FTPath, sep = "\t", header = TRUE)

  # Add protein data if TargetsPath is not NULL
  if (is.null(TargetsPath) == FALSE) {

    # Read targets file, extract protein data, and merge with MS1FT, keeping only MS1FT the same length
    Targets <- data.table::fread(TargetsPath, sep = "\t", header = TRUE)

    # Make FeatureID columns match up for merging
    colnames(Targets)[colnames(Targets) == "Ms1Features"] <- "FeatureID"
    MS1FT <- merge(MS1FT, Targets, by = "FeatureID", all.x = T)

    # Rename is NA proteins as "Protein Not Found"
    MS1FT$ProteinName <- ifelse(is.na(MS1FT$ProteinName), "Protein Not Found", MS1FT$ProteinName)

  }

  # Trim Columns down if "TrimColumns" is true
  if (TrimColumns) {

    # Trim table down
    if (is.null(TargetsPath)) {

      MS1FT <- MS1FT %>% dplyr::select(c(FeatureID, MinElutionTime, MaxElutionTime, MonoMass, Abundance))
      MS1FT$ProteinName <- "Protein Not Found"

    } else {

      MS1FT <- MS1FT %>% dplyr::select(c(FeatureID, MinElutionTime, MaxElutionTime, MonoMass, Abundance, ProteinName))

    }

  }

  ###################
  ## RETURN OBJECT ##
  ###################

  # Add attributes
  attr(MS1FT, "pspecter")$MS1FTPath = MS1FTPath
  attr(MS1FT, "pspecter")$TargetsPath = TargetsPath

  # Add class information
  class(MS1FT) <- c(class(MS1FT), "ms1ft")

  return(MS1FT)

}
