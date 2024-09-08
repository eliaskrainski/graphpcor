# corGraphs

Models for correlation matrices based on graphs.

The 'INLA' package is a suggested one, but you will need it 
for actually fitting a model. You can install it with
```
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE) 
```

There are two suggested packages, ‘graph’ and ‘Rgraphviz’, 
that are on Bioconductor, and you can install those with:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE) 
```
