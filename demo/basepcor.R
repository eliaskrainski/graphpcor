
## p = 3, m = 2
bc <- basepcor(c(-1,-1), p = 3, itheta = c(2,3))
bc

round(solve(bc$base), 4)

all.equal(bc,
          basepcor(bc$base, itheta = c(2,3)))

## p = 4, m = 4
th2 <- c(0.5,-1,0.5,-0.3)
ith2 <- c(2,3,8,12)
b2 <- basepcor(th2, p = 4, itheta = ith2)
b2

Sparse(solve(b2$base), zeros.rm = TRUE)

all.equal(th2, basepcor(b2$base, itheta = ith2)$theta)

## Hessian around the base 
hessian(b2)

## p = 4, m = 3 with some common theta
th3 <- c(0.5, -1, -0.3)
ip <- c(1, 2, 1, 3) ## 1st == 3rd
b3 <- basepcor(th3, p = 4, itheta = ith2, iparams = ip)

all.equal(b2$base, b3$base) ## TRUE

## but the parameter dimension is now reduced
hessian(b3)

## If a subset of the parameters are known (fixed), then the
## Hessian is only computed with respect to the unknown ones
hessian(basepcor(th2, p=4, itheta = ith2, iunknown = 1))
hessian(basepcor(th2, p=4, itheta = ith2, iunknown = c(2,3)))
hessian(basepcor(th2, p=4, itheta = ith2, iunknown = c(2,4)))

hessian(basepcor(th3, p=4, itheta = ith2,
                 iparams = ip, iunknown = 1))
hessian(basepcor(th3, p=4, itheta = ith2,
                 iparams = ip, iunknown = 1:2))
hessian(basepcor(th3, p=4, itheta = ith2,
                 iparams = ip, iunknown = 3))

## check
all.equal(b2, basepcor(th2, p = 4, itheta = ith2,
                       iunknown = 1:4))
all.equal(b3, basepcor(th3, p = 4, itheta = ith2,
                       iparams = ip, iunknown = 1:3))
