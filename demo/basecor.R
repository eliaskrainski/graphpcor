library(graphpcor)

## A correlation matrix
c0 <- matrix(c(1,.9,.5, .9,1,.2, .5,.2,1), 3)

## build the 'basecor'
pc.c0 <- basecor(c0) ## base as matrix
pc.c0

## elements
pc.c0$base
pc.c0$theta
pc.c0$I0

## from 'theta' 
th0 <- pc.c0$theta
pc.th0 <- basecor(th0) ## base as vector
pc.th0

## numerically the same
all.equal(c0, pc.th0$base)

## from a numeric vector (theta)
th2 <- c(-1, -0.5)
bth2 <- basecor(th2, p = 3, itheta = c(2,3))
bth2

## from the correlation matrix
b2 <- basecor(bc$base, itheta = c(2,3))

all.equal(th2, b2$theta, tol = 1e-4)

