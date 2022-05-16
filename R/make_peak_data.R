#' Create a peak_data object from an M/Z and Intensity vector
#'
#' @description Generate a peak_data object. The MZ and Intensity vectors must be
#'     of the same length.
#'
#' @param MZ A vector of MZ values. Positive non-zero values only. Required.
#' @param Intensity A vector Intensity values. Positive non-zero values only. Required.
#'
#' @details
#' The data.table outputted by this function contains both the M/Z and Intensity vectors of the spectra.
#'
#' @examples
#' \dontrun{
#'
#' make_peak_data("MZ" = c(300, 301, 302), "Intensity" = c(50, 100, 50))
#'
#' }
#' @export
make_peak_data <- function(MZ,
                           Intensity) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Make a function to check MZ and Intensity
  .check_input <- function(vector, name) {

    # Both have to be numeric
    if (!is.numeric(vector)) {
      stop(paste(name, "must be numeric."))
    }

    # Both should be positive and nonzero
    if (any(vector < 0) | any(0 %in% vector)) {
      stop(paste(name, "should contain non-zero and positive values."))
    }

  }

  # Check inputs
  .check_input(MZ, "MZ")
  .check_input(Intensity, "Intensity")

  # Both vectors should be the same length
  if (length(MZ) != length(Intensity)) {
    stop("Both MZ and Intensity should be the same length.")
  }

  ##################
  ## BUILD OBJECT ##
  ##################

  # Generate peaks data.table
  Peaks <- data.table::data.table(
    "M/Z" = MZ,
    "Intensity" = Intensity,
    "Abundance" =  round(Intensity / max(Intensity) * 100, 4)
  ) %>%
    dplyr::arrange(`M/Z`)

  # Keep track of the inputs in the attributes
  attr(Peaks, "pspecter") <- list(
    "MSPath" = NULL,
    "ScanNumber" = NULL,
    "MinimumAbundance" = NULL,
    "TotalNumberPeaks" = NULL,
    "NumberPeaksPostFilter" = NULL,
    "PercentagePeaksRemain" = NULL
  )

  # Assign the peak data class
  class(Peaks) <- c(class(Peaks), "peak_data")

  # Return peaks
  return(Peaks)

}
