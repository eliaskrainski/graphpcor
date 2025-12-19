library(graphpcor)

bc <- basecor(c(-1,-1,0.5), p=3)

R <- bc$base
R

bc
basecor(R)

all.equal(bc$theta, basecor(R)$theta, tol = 1e-4)
