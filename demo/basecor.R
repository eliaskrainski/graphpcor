
## A correlation matrix
c0 <- matrix(c( 1.0,  0.8, -0.5,
                0.8,  1.0, -0.4,
               -0.5, -0.4,  1.0), 3)

## build the 'basecor'
b3 <- basecor(c0) ## base as matrix
b3

## elements
b3$base
b3$theta

## from 'theta' 
th3 <- b3$theta
bth3 <- basecor(th3, p = 3) ## base as vector
bth3

## numerically the same
all.equal(c0, bth3$base, tol = 1e-4)

## from a numeric vector (theta)
th2 <- c(-1, -0.5)
b2 <- basecor(th2, p = 3, iLtheta = c(2,3))
b2

## from the correlation matrix
b2c <- basecor(b2$base, iLtheta = c(2,3))

all.equal(th2, b2c$theta, tol = 1e-4)

## Hessian around the base 
hessian(b3)
hessian(b3, ifixed = 3)
hessian(b2)
