  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/EMSL-Computing/pspecterlib/branch/master/graph/badge.svg)](https://app.codecov.io/gh/EMSL-Computing/pspecterlib?branch=master)
  <!-- badges: end -->

# The PSpecteR Library (pspecterlib)
Allows users to use the backend functionality of the PSpecteR proteomics visualization application independent of the app. 

**Highlights of PSpecteR functionality include:**

1. Calculating and visualizing peptide/protein fragmentation patterns on experimental spectra
2. Generating extracted ion chromatograms (XICs)
3. Mapping identified peptides to protein sequences
4. Testing alternative peptides and mass modified ions 
5. Visualizing output from MSPathFinderT

...for both top-down and bottom-up proteomics. 

*Inputs:* MS file (mzML or ThermoFisher raw), ID file (optional mzid), FASTA file (optional) 

# How to install 

`devtools::install_github("EMSL-Computing/pspecterlib")`
