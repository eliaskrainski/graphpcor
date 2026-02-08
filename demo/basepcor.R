## A correlation matrix
c0 <- matrix(c( 1.0,  0.8, -0.5,
                0.8,  1.0, -0.4,
               -0.5, -0.4,  1.0), 3)

## p = 3, m = 3
b1 <- basepcor(c0)
b1

## p = 3, m = 2
b2 <- basepcor(c0, iLtheta = c(2,3))
b2

all.equal(b1$base, b2$base)

## Hessian
hessian(b2)
hessian(b1, ifixed = 3)
hessian(b1)

## p = 4, m = 4
th4 <- c(0.5,-1,0.5,-0.3)
ith4 <- c(2,3,8,12)
b44 <- basepcor(th4, p = 4, iLtheta = ith4)
b44

Sparse(round(solve(b44$base), 4), zeros.rm = TRUE)

## p = 4, m = 3 (with some common theta)
th3 <- c(0.5, -1, -0.3)
ip3 <- c(1, 2, 1, 3) ## 1st == 3rd
b43 <- basepcor(th3, p = 4, iLtheta = ith4, iparams = ip3)

all.equal(b44$base, b43$base) ## TRUE

## parameter dimension is now reduced
hessian(b44)
hessian(b43)

## If a subset of the parameters are known (fixed), then the
## Hessian is only computed with respect to the unknown ones
hessian(basepcor(th4, p=4, iLtheta = ith4), ifixed = 2:3)
hessian(basepcor(th4, p=4, iLtheta = ith4), ifixed = 1)
hessian(basepcor(th4, p=4, iLtheta = ith4), ifixed = 3)

hessian(basepcor(th3, p=4, iLtheta = ith4,
                 iparams = ip3), ifixed = 3)
hessian(basepcor(th3, p=4, iLtheta = ith4,
                 iparams = ip3), ifixed = 2:3)
hessian(basepcor(th3, p=4, iLtheta = ith4,
                 iparams = ip3), ifixed = 1:2)

## check
hessian(basepcor(th3, p=4, iLtheta = ith4,
                 iparams = ip3), ifixed = NULL)

