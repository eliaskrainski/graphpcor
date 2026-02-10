cp NAM* DESCR* READ*md   ../main_branch_graphpcor/
cp R/*R                  ../main_branch_graphpcor/R/
cp data/*                ../main_branch_graphpcor/data/
cp demo/*                ../main_branch_graphpcor/demo/
cp examples/*R           ../main_branch_graphpcor/examples/
cp man/*R                ../main_branch_graphpcor/man/
cp vignettes/*Rmd        ../main_branch_graphpcor/vignettes/
cd ../main_branch_graphpcor/
head -n -1 NAMESPACE > .namespacetemp; mv .namespacetemp NAMESPACE
sed -i '/useDynLib/d' R/basepcor_utils.R
git add DESCRIPTION NAMESPACE README.Rmd README.md 
git add R/* data/* demo/* examples/* man/* vignettes/*Rmd
git commit -m 'update main'
git push
cd ../graphpcor/
