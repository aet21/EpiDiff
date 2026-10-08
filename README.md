---
title: "EpiDiff: DNA methylation based prediction of differentiation"
author: "Andrew E. Teschendorff, Huige Tong"
date: "`r Sys.Date()`"
package: "`r pkg_ver('EpiDiff')`"
output:
  BiocStyle::html_document
bibliography: EpiDiff.bib
vignette: >
  %\VignetteIndexEntry{Epigenetic Predictor of Differentiation State}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# System Requirements
The **EpiDiff** R package only imports that randomForest R-package is installable on all UNIX, MAC-OSX and Windows platforms. **EpiDiff** was developed and tested on Ubuntu 24.04.5 LTS and R version 4.6.1 .

# Installation guide
You can install with
```{r load, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
library(devtools);
devtools::install_github("aet21/EpiDiff");
```
or download the .tar.gz file from **aet21/EpiDiff** and install from command line using:
*R CMD INSTALL EpiDiff_X.Y.Z.tar.gz*
replacing X, Y and Z with current version number. Installing it on a normal desktop computer only takes seconds.

# Demonstration
**EpiDiff** is an R-package to estimate differentiation state of a DNA methylation sample (@EpiDiffPaper). Please refer to the vignette within the EpiDiff package, which provides a detailed demo of how EpiDiff works on a small real dataset. Briefly, the input to the main function is a normalized DNA methylation data matrix, output will include the estimated differentiation index (DI) for each sample (column) in the data matrix. Estimation of DI for a 1000-sample data matrix without needing to impute any CpGs takes 1 second.

# Instructions for use
Currently, only **EpiDiff** can be applied to WGBS and snmC-Seq data, the package only supports application to Illumina based DNA methylation beadarray platforms. Instructions on how to use EpiDiff are provided in the vignette of the EpiDiff R-package, where we provide a detailed demo of how to apply it to a small real Illumina DNAm dataset.


# Sessioninfo

```{r sessionInfo, echo=FALSE}
sessionInfo()
```

# References



