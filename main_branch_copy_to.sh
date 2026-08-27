cp NAM* DESCR* READ*md   ../main_branch_graphpcor/
cp R/*R                  ../main_branch_graphpcor/R/
cp data/*                ../main_branch_graphpcor/data/
cp demo/*                ../main_branch_graphpcor/demo/
cd examples
cp Germany4.R PBC2.R PBC3.R SAheart_inla.R SAheart_stan.R swiss.R ../../main_branch_graphpcor/examples/
cd ..
cd src; cp *c *h ../../main_branch_graphpcor/src/; cd ..
cp man/*Rd               ../main_branch_graphpcor/man/
cp vignettes/*.Rmd        ../main_branch_graphpcor/vignettes/
cd ../main_branch_graphpcor/
## head -n -1 NAMESPACE > .namespacetemp; mv .namespacetemp NAMESPACE
## sed -i 's/useDynLib, .registration = TRUE//' R/basepcor_utils.R

git add DESCRIPTION NAMESPACE README.Rmd README.md 
git add R/* data/* demo/* examples/* man/* vignettes/*Rmd
git commit -m 'update main'
git push

cd ../graphpcor/
