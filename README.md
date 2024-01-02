  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/EMSL-Computing/pspecterlib/branch/master/graph/badge.svg)](https://app.codecov.io/gh/EMSL-Computing/pspecterlib?branch=master)
  <!-- badges: end -->

# The PSpecteR Library (pspecterlib)

<img src="PSpecteR_Logo.png" alt="PSpecteR_Logo" style="width:320px;margin-right: 40px;" align="left"/>

Allows for the use of backend functionality of the PSpecteR proteomics visualization 
application independent of the app. 

**Highlights of PSpecteR functionality include:**

* Calculating and visualizing peptide/protein fragmentation patterns on experimental spectra 
* Generating extracted ion chromatograms (XICs)
* Mapping identified peptides to protein sequences
* Testing alternative peptides and mass modified ions  
* Visualizing output from MSPathFinderT

...for both top-down and bottom-up proteomics. 

*Inputs:* MS file (mzML or ThermoFisher raw), ID file (optional mzid), FASTA file (optional) 

# How to install 

`devtools::install_github("EMSL-Computing/pspecterlib")`

Read our vignette [here](https://emsl-computing.github.io/pspecterlib/)

# What does "PSpecteR" mean? 

Excellent question! We get the name from **P**roteomics Mass **Spect**rum Insp**e**cto**r**!
The last R is capitalized for the language it's written in. Now we know what you're 
thinking. Why isn't it PSpectoR, then? Well, our logo is a ghost inspector named
Inspector Specter. So we chose "specter". Now you know, and are empowered with 
this knowledge.





