library(graphpcor)
bc <- basepcor(c(-1,-1), p=3, itheta = c(2,3))
all.equal(bc, basepcor(bc$base, itheta =c(2,3)))
