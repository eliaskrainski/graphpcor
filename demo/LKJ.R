## a cgeneric model for the LKJ prior
## for theta from the CPC parametrization
cglkj <- cgeneric(
    model = "LKJ", n = 3, eta = 2,
    useINLAprecomp = FALSE)##!is.na(packageCheck("INLA", "25.12.16")))

## correlation matrix, p = 3
cc <- matrix(c(1,.8,-.625, 0.8,1,-.5, -0.625,-.5,1), 3)

## CPC parametrization: C(theta)
(bb <- basecor(cc))
th <- bb$theta

## p(theta | eta)
cgeneric_prior(cglkj, theta = th)

## precision inverse
solve(cgeneric_Q(cglkj, theta = th))
solve(cgeneric_Q(cglkj, theta = c(0,0,0)))
solve(cgeneric_Q(cglkj, theta = c(1,1,1)))


