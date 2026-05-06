cd src/
cp FUNCTIONS cgeneric_graphpcor.c  cgeneric.h  cgeneric_LKJ.c  cgeneric_pc_correl.c cgeneric_treepcor.c  cgeneric_Wishart.c  graphpcor.h  graphpcor_utils.c  graphpcor_utils.h ../../inla_branch_graphpcor/src/
cd ../../inla_branch_graphpcor/
git add src/* 
git commit -m 'update inla branch'
git push
cd ../graphpcor/
