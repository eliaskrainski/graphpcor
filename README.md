
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Models for correlation matrices based on graphs

<img src="graphpcor_logo.png" style="width:3in" />

Graphs are used to represent dependency and to define models for
correlation matrices. The Penalized Complexity - PC prior for the graph
models are used to define priors that penalizes the contraction of the
correlation structure from a simpler one. Models are implemented using
the ‘cgeneric’ interface in the ‘INLA’ package
(<https://www.r-inla.org>) so that ‘INLA’ can be used to build and fit
complex data models using these as building blocks.

## Installing dependencies

The ‘INLA’ package is a suggested one, but you will need it for actually
fitting a model. You can install it with

    install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE) 

There are two suggested packages, ‘graph’ and ‘Rgraphviz’, that are on
Bioconductor, and you can install those with:

    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE) 
