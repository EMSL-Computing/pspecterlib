## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(pspecterlib)
library(dplyr)
library(DT)
library(ggplot2)

## ----copy data----------------------------------------------------------------
# Create a temporary directory and copy example data there
tmpdir <- tempdir()

# Pull example bottom up data filepath 
files <- c(
 "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzML",
 "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzid",
 "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown.mzML",
 "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown.mzid",
 "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/QC_Shew.fasta"
)

# Download files to temporary directory
for (file in files) {
  download.file(file, file.path(tmpdir, tail(unlist(strsplit(file, "/")), 1)))
}

# Test with raw top down data
library(rawrr)
rawfile <- file.path(path.package(package = 'rawrr'), 'extdata', 'sample.raw')
RAW_ScanMetadata <- get_scan_metadata(MSPath = rawfile)

## ----load data----------------------------------------------------------------
# Create the initial scan_metadata object to pass to downstream modules 
BU_ScanMetadata <- get_scan_metadata(MSPath = file.path(tmpdir, "BottomUp.mzML"),
                                     IDPath = file.path(tmpdir, "BottomUp.mzid"))

# Display first 6 entries in ScanMetadata table
BU_ScanMetadata %>%
  head() %>% 
  datatable(options = list(scrollX = TRUE))

## ----visualize scan metadata, fig.width = 6-----------------------------------
scan_metadata_plot(BU_ScanMetadata, XVar = "Precursor M/Z", YVar = "Retention Time",
  LabVar = "Score", Interactive = TRUE, MSFilter = 2, ScanNumFilter = c(32000, 34500))

## ----peak data----------------------------------------------------------------
# Use the scan number of the lowest e-score
BU_Peak <- get_peak_data(BU_ScanMetadata, 31728, MinAbundance = 1)

# View the first few peaks
BU_Peak %>%
  head() %>% 
  datatable()

