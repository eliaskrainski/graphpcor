library(graphpcor)

## graph in Example 2.6 of the GMRF book
g <- graphpcor(x1~x2+x3, x2~x4, x3~x4)
(ne <- dim(g))

## Laplacian
(G <- Laplacian(g))

## base model (theta for L)
theta0l <- rep(-0.3, ne[2])

## define initial L
theta.low <- rep(-1, ne[2])
idx.low <- (!is.zero(G) & lower.tri(G))
L1a <- diag(ne[1]); L1a[idx.low] <- theta.low
L1a

## the C code to fill.in
idx.fill <- which(is.zero(G) & (!is.zero(t(chol(G + diag(ne[1]))))))
nfi <- length(idx.fill)
nfi

if(nfi) {
    ij.fill <- cbind(
        row(G)[idx.fill],
        col(G)[idx.fill])
    print(ij.fill)
    Lf <- .C("fillL",
             as.integer(ne[1]),
             as.integer(nfi),
             as.integer(ij.fill[, 1]-1),
             as.integer(ij.fill[, 2]-1),
             l=L1a)$l
    print(Lf)
    Q0 <- crossprod(Lf)
} else {
    Q0 <- crossprod(L1a)
}
Q0

V0 <- chol2inv(t(Lf))
V0

variance(g, theta = c(log(diag(V0))/2, theta.low))

cov2cor(V0)
##           [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.9486833 0.9128709 0.7745967
## [2,] 0.9486833 1.0000000 0.8660254 0.8164966
## [3,] 0.9128709 0.8660254 1.0000000 0.7071068
## [4,] 0.7745967 0.8164966 0.7071068 1.0000000


## define the sigma parameters
logsigmas <- log(c(0.3, 0.7, 1.5, 0.9))

sr <- exp(logsigmas)/sqrt(diag(V0))
sr
##[1] 0.07745967 0.28577380 1.06066017 0.90000000

V1 <- diag(sr) %*% V0 %*% diag(sr)
V1
##           [,1]      [,2]      [,3]      [,4]
## [1,] 0.0900000 0.1992235 0.4107919 0.2091411
## [2,] 0.1992235 0.4900000 0.9093267 0.5143928
## [3,] 0.4107919 0.9093267 2.2500000 0.9545942
## [4,] 0.2091411 0.5143928 0.9545942 0.8100000

stopifnot(all.equal(sqrt(diag(V1)), exp(logsigmas)))

round(chol(V1), 4)
##      [,1]   [,2]   [,3]   [,4]
## [1,]  0.3 0.6641 1.3693 0.6971
## [2,]  0.0 0.2214 0.0000 0.2324
## [3,]  0.0 0.0000 0.6124 0.0000
## [4,]  0.0 0.0000 0.0000 0.5196

## build the cgeneric model
## Note: here 'model' is a 'graphpcor'
cmodel <- cgeneric(
    model = g, ## use the graphpcor
    lambda = 1,
    base = theta0l,
    sigma.prior.reference = rep(1, ne[1]),
    sigma.prior.probability = rep(0.5, ne[1]),
    debug = 1e9)

prior(cmodel, theta = theta0l)

mu(cmodel)

initial(cmodel)

graph(cmodel) ## the structure of the precision

cmodel$f$cgeneric$data$doubles$thetabasescaled
## [1] -0.2190589 -0.2256071 -0.2690858 -0.2829668
drop(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]) %*% theta0l)
## [1] -0.2190589 -0.2256071 -0.2690858 -0.2829668

round(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]), 3)
##        [,1]   [,2]   [,3]   [,4]
## [1,]  0.928 -0.112 -0.078 -0.008
## [2,] -0.112  0.945 -0.035 -0.046
## [3,] -0.078 -0.035  1.014 -0.004
## [4,] -0.008 -0.046 -0.004  1.001

drop(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]) %*% theta.low)
## [1] -0.7301962 -0.7520235 -0.8969525 -0.9432227

drop(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]) %*% theta.low) -
    cmodel$f$cgeneric$data$doubles$thetabase
## [1] -0.5111374 -0.5264165 -0.6278668 -0.6602559

theta1 <- c(d = logsigmas,
            l = theta.low)

Q1c <- prec(cmodel, theta = theta1)

round(Q1c, 2)
round(Q <- prec(g, theta = theta1), 2)

##          [,1]     [,2]     [,3]    [,4]
## [1,] 166.6667 -45.1754 -12.1716  0.0000
## [2,] -45.1754  24.4898   0.0000 -3.8881
## [3,] -12.1716   0.0000   2.6667  0.0000
## [4,]   0.0000  -3.8881   0.0000  3.7037

prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = theta1)

dataf <- list(
    i = 1:ne[1],
    y = rep(NA, ne[1]))

library(INLA)

fit1 <- inla(
    formula = y ~ 0 + f(i, model = cmodel),
    family = 'poisson',
    data = dataf,
    control.mode = list(theta = theta1, fixed = TRUE)
)

all.equal(Q1c, prec(fit1))

