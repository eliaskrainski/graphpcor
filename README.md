# corGraphs

Models for correlation matrices based on graphs.

There are two suggested packages, ‘graph’ and ‘Rgraphviz’, 
that are on Bioconductor, and you can install those with:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE) 
```
