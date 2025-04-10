library(graphpcor)

g <- graphpcor(x ~ y, y ~ v, v ~ z, z ~ x)

plot(g)

class(g)

g

ne <- dim(g)
ne

summary(g)

## sometimes we need it
G <- Laplacian(g)
G

## alternatively
all.equal(G,
          Laplacian(graphpcor(x~y+z,y~v,v~z)))

## compact, but different ordering
Laplacian(graphpcor(x~y+z,v~y+z))

graphpcor(x1~x2+x3, x4~x2+x3)
Laplacian(graphpcor(x1~x2+x3, x4~x2+x3))

g <- graphpcor(x1~x2+x3, x2~x4, x3~x4) ## compact ordered
(G <- Laplacian(g)) ## the graph in Example 2.6 of the GMRF book

graphpcor(G) ## dag from a matrix

all.equal(graphpcor(G), g) ## TRUE if compact ordered

## base model (theta for L)
theta0l <- rep(-0.5, ne[2])

## build a Cholesky (for precision)
L0 <- chol(g, theta = theta0l)
L0

round(tcrossprod(L0), 4)

## the precision for a correlation matrix
Q0 <- prec(g, theta = theta0l)
Q0

C0 <- solve(Q0)
C0

## the Hese matrix around a base model
## using numDeriv::hessian
hessian(function(x) graphpcor:::KLD10(vcov(g, theta=x), C0), x = theta0l)

## using hessian for a `corgtraph` returns more stuff:
hessian(g, base = theta0l)

## using different ways to specify base model and different decomposition
all.equal(
    hessian(g, base = theta0l, decomposition = 'svd'),
    hessian(g, base = C0, decomposition = 'svd')
)

## variance method for graphpcor computes the correlation
##  if only theta for lower of L is provided
all.equal(
    C0,
    vcov(g, theta = theta0l)
)

vcov(g, theta = theta0l)
vcov(g, theta = c(log(c(2,1,4,0.5)), theta0l))

## the 'iid' case would be
vcov(g, theta = rep(0, ne[2]))
vcov(g, theta = rep(0, sum(ne)))

## marginal variance specified throught standard errors
sigmas <- c(0.3, 0.7, 1.2, 0.5)

vcov(g, theta = c(log(sigmas), theta0l))
vcov(g, theta = c(log(sigmas), rep(0, ne[2])))
vcov(g, theta = c(log(sigmas), rep(-1, ne[2])))
all.equal(vcov(g, theta =  rep(-1, ne[2])),
          cov2cor(vcov(g, theta = c(log(sigmas), rep(-1, ne[2])))))

vcov(g, theta = rep(1, ne[2])) ## no edge 2~3 but high covariance!!!
vcov(g, theta = c(log(sigmas), rep(1, ne[2]))) ## no edge 2~3 but high covariance!!!

vcov(g, theta = c(-.5,-.5,-.5,.5))
vcov(g, theta = c(-5,-5,-5,5))

## build the cgeneric model
## Note: here 'model' is a 'graphpcor'
library(INLA)
cmodel <- cgeneric(
    model = g, ## a `graphpcor` in model argument
    lambda = 1,
    base = theta0l, 
    sigma.prior.reference = rep(1, ne[1]),
    sigma.prior.probability = rep(0.5, ne[1]))

## Note: another way: using the Laplacian
all.equal(
    cmodel,
    cgeneric(
        model = "graphpcor", ## model now is acharacter
        graph = G, ## using G as a graph
        lambda = 1,
        base = theta0l,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## yet another way: a binary matrix
all.equal(
    cmodel,
    cgeneric(
        model = "graphpcor",
        graph = G!=0, ## any binary matrix works
        lambda = 1,
        base = theta0l,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## specify the base model from C0 (SAME generated from theta0l)
all.equal(
    cmodel,
    cgeneric(
        model = g,
        lambda = 1,
        base = C0,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## compatible base model
Qc <- Q0; Qc[1,3] <- Qc[3,1] <- 0
Qc

## all non-zero in Qc are also non-zero in Q0
all(which(!is.zero(Qc)) %in%
    which(!is.zero(Q0)))

### also chol
all(which(!is.zero(chol(Qc))) %in%
    which(!is.zero(chol(Q0))))

Cc <- chol2inv(Qc)
Cc

## the base model elemens (lconst, thetabase, hHeg)are now different
all.equal(cmodel,
          cgeneric(
              model = g,
              lambda = 1,
              base = Cc,
              sigma.prior.reference = rep(1, ne[1]),
              sigma.prior.probability = rep(0.5, ne[1])))

## incompatible base model
Qi <- Q0; Qi[4,1] <- Qi[1, 4] <- -0.5
Qi

Vi <- chol2inv(chol(Qi))
Ci <- cov2cor(Vi)
Ci

if(FALSE) { ## this give an error because Ci is not compatible
    cgeneric(
        model = g,
        lambda = 1,
        base = Ci,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
}

graph(cmodel)

## define some model
theta1 <- c(
    d = log(sigmas),
    l = rep(-1, ne[2])
)

prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = theta1)

V1 <- vcov(g, theta = theta1)
V1

C1 <- cov2cor(vcov(g, theta = theta1))
C1

round(solve(V1), 2)

Q1c <- prec(cmodel, theta = theta1)
Q1c

round(Q1c, 2)

all.equal(solve(as.matrix(Q1c)), vcov(g, theta = theta1))
all.equal(as.matrix(Q1c), prec(g, theta = theta1))

dataf <- list(
    i = 1:ne[1],
    y = rep(NA, ne[1]))

library(INLA)

fit0 <- inla(
    formula = y ~ 0 + f(i, model = cmodel),
    family = 'poisson',
    data = dataf,
    control.mode = list(theta = theta1, fixed = TRUE)
)

all.equal(Q1c, prec(fit0))

