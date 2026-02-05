rm -rf ~/temp/graphpcor*
mkdir ~/temp/graphpcor/
cp -r DESCRIPTION NAMESPACE R/ man/ src/ demo/ ~/temp/graphpcor/
cd vignettes/
mkdir ~/temp/graphpcor/vignettes/
cp preamble.tex references.bib treepcor.Rmd ~/temp/graphpcor/vignettes/
cd ~/temp/
R CMD build graphpcor
R CMD check graphpcor_*.tar.gz --as-cran
