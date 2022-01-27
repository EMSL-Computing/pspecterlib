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
#'    amino acid, and will not contain the letters B, J, O, U, X, or Z. Modification
#'    notation (i.e. "TES[Acetyl]T") is not accepted. Modifications should be of the
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


