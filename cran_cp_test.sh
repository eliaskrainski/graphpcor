rm -rf ~/temp/graphpcor*
mkdir ~/temp/graphpcor/
git checkout main 
cp -r DESCRIPTION NAMESPACE R/ man/ demo/ ~/temp/graphpcor/
cd vignettes/
mkdir ~/temp/graphpcor/vignettes/
cp preamble.tex references.bib treepcor.Rmd ~/temp/graphpcor/vignettes/
git checkout devel
cd ~/temp/
R CMD build graphpcor
R CMD check graphpcor_*.tar.gz --as-cran
