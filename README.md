# README #

### What is this repository for? ###

This repository contains a workflow for the analysis fo single-cell RNASeq data using R/bioconductor. 
The main steps included are:

* Quality control
* Data normalization
* Feature selection
* Clustering
* Differential expression analysis
* Detection of rare cell subtypes with CellSIUS 
* Visualization

[![DOI](https://zenodo.org/badge/189391647.svg)](https://zenodo.org/badge/latestdoi/189391647)

### How do I get set up? ###

Check the vignette (vignettes/workflow_vignette.html) for a detailed description. All code can be run directly form the R-Markdown document (vignettes/workflow_vignette.Rmd).

IMPORTANT: This workflow was built and tested using Bioconductor 3.5. Because some of the packages it uses (especially scater) changed between Bioconductor 3.5 and 3.6, currently, the only supported version is R3.4.1 with Bioconductor 3.5.

### Need help? ###

Please contact Rebekka Wegmann [wegmann@imsb.biol.ethz.ch](mailto:wegmann@imsb.biol.ethz.ch) or Marilisa Neri [marilisa.neri@novartis.com](mailto:marilisa.neri@novartis.com)
