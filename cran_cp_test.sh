./main_branch_copy_to.sh
rm -rf ~/temp/graphpcor*
mkdir ~/temp/graphpcor/
cd ../main_branch_graphpcor/
cp -r DESCRIPTION NAMESPACE data/ demo/ man/ R/ src/ ~/temp/graphpcor/
cd vignettes/
mkdir ~/temp/graphpcor/vignettes/
cp preamble.tex references.bib treepcor.Rmd ~/temp/graphpcor/vignettes/
cd ~/temp/
R CMD build graphpcor
R CMD check graphpcor_*.tar.gz --as-cran
