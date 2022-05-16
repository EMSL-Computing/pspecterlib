# Checks for and downloads required files

file_check <- function(test) {

  # Get download folder name
  if (Sys.info()["sysname"] == "Windows") {
    downfolder <- dirname("~")
  } else {
    downfolder <- path.expand("~")
    downfolder <- file.path(downfolder, "Downloads")
    downfolder <- paste0(downfolder, .Platform$file.sep)
  }

  # Add pspecter test
  downfolder <- file.path(downfolder, "PSpecterTest")

  # Return path if test is downfolder
  if (test == "downfolder") {
    return(downfolder)
  }

  # If the folder does not exist, download test data
  if (!dir.exists(downfolder)) {

    dir.create(downfolder)

    files <- c(
      "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzML",
      "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzid",
      "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown.ms1ft",
      "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown_IcTarget.tsv",
      "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/QC_Shew.fasta"
    )

    for (file in files) {
      download.file(file, file.path(downfolder, tail(unlist(strsplit(file, "/")), 1)))
    }

  }

  # Now create and return the appropriate scan metadata object
  if (test == "BU") {
    testdata <- get_scan_metadata(file.path(downfolder, "BottomUp.mzML"), file.path(downfolder, "BottomUp.mzid"))
  } else if (test == "RAW") {
    library(rawrr)
    rawfile <- file.path(path.package(package = 'rawrr'), 'extdata', 'sample.raw')
    testdata <- get_scan_metadata(MSPath = rawfile)
  } else {
    stop("Test name not recognized ")
  }

  return(testdata)

}
