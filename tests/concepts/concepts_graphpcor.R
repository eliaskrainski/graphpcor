library(graphpcor)

par(mfrow = c(2, 3), mar = c(0,0,0,0))
plot(graphpcor(x~y+v, z~y+v))
plot(graphpcor(x~y,x~v,z~y,z~v))
plot(graphpcor(x~y, v~x, y~z, z~v))
plot(graphpcor(y~x, v~x, z~y, z~v))
plot(graphpcor(y~x+z, v~z+x))
plot(graphpcor(y~x+z, v~x, z~v))

## the graph in Example 2.6 of the GMRF book
g <- graphpcor(x ~ y + v, z ~ y + v)

class(g)
g

par(mfrow=c(1,1))
plot(g)

summary(g) ## the graph: nodes and edges (nodes ordered as given)

ne <- dim(g)
ne

## sometimes we need it
G <- Laplacian(g)
G 

## alternatively
all.equal(G,
          Laplacian(graphpcor(x~y, v~x, y~z, z~v)))
all.equal(G,
          Laplacian(graphpcor(x~y, x~v, z~y, v~z)))

plot(graphpcor(G)) ## from a matrix

## base model (theta for lower triangle Cholesky)
theta0l <- rep(-0.5, ne[2])

## vcov() method for graphpcor computes the correlation
##  if only theta for lower of L is provided
C0 <- vcov(g, theta = theta0l)
C0

## the precision for a correlation matrix
Q0 <- cgeneric_Q(g, theta = theta0l)
Q0

all.equal(C0, as.matrix(solve(Q0)))

## the Hessian matrix around a base model
H0 <- hessian(g, x = theta0l)
H0

## a base model can also be a matrix
## however it shall give a precision with
## same sparse pattern as the graph
all.equal(H0, hessian(g, x = C0))

## the 'iid' case would be
vcov(g, theta = rep(0, ne[2]))
vcov(g, theta = rep(0, sum(ne)))

## marginal variance specified throught standard errors
sigmas <- c(0.3, 0.7, 1.2, 0.5)
## the covariance
vcov(g, theta = c(log(sigmas), rep(0, ne[2]))) ## IID
vcov(g, theta = c(log(sigmas), theta0l))

vcov(g, theta = rep(-3, ne[2])) ## no edge 2~3 but high correlation!!!

## build the cgeneric model
## Note: here 'model' is a 'graphpcor'
cmodel <- cgeneric(
    model = g, ## a `graphpcor` in model argument
    lambda = 1,
    base = theta0l, 
    sigma.prior.reference = rep(1, ne[1]),
    sigma.prior.probability = rep(0.5, ne[1]))

## Another way: using C0 matrix for base model
all.equal(
    cmodel,
    cgeneric(
        model = g, 
        lambda = 1,
        base = C0,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)


## Note: another way: using the Laplacian
all.equal(
    cmodel,
    cgeneric(
        model = G, ## using G as a graph
        lambda = 1,
        base = theta0l,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## yet another way: a binary matrix
all.equal(
    cmodel,
    cgeneric(
        model = G!=0, ## any binary matrix works
        lambda = 1,
        base = theta0l,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## compatible base model
Qc <- Q0; Qc[1,3] <- Qc[3,1] <- 0
Q0
Qc

## all non-zero in Qc are also non-zero in Q0
all(which(!is.zero(Qc)) %in%
    which(!is.zero(Q0)))

### also chol
all(which(!is.zero(chol(Qc))) %in%
    which(!is.zero(chol(Q0))))

Cc <- chol2inv(Qc)
Cc

## the base model elemens (thetabasescaled, thetab, hHneg and H) are now different
all.equal(
    cmodel,
    cgeneric(
        model = g,
        lambda = 1,
        base = Cc,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

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

cgeneric_graph(cmodel)

## define some model parameters
theta1 <- c(
    d = log(sigmas),
    l = rep(-1, ne[2])
)

cgeneric_prior(cmodel, theta = rnorm(sum(ne)))
cgeneric_prior(cmodel, theta = rnorm(sum(ne)))
cgeneric_prior(cmodel, theta = theta1)

V1 <- vcov(g, theta = theta1)
V1

C1 <- cov2cor(vcov(g, theta = theta1))
C1

round(solve(V1), 2)

Q1c <- cgeneric_Q(cmodel, theta = theta1)
Q1c

round(Q1c, 2)

all.equal(solve(as.matrix(Q1c)),
          vcov(g, theta = theta1))

all.equal(Q1c, 
          cgeneric_Q(g, theta = theta1))

dataf <- list(
    i = 1:ne[1],
    y = rep(NA, ne[1]))

library(INLA)

fit0 <- inla(
    formula = y ~ 0 + f(i, model = cmodel),
    family = 'poisson',
    data = dataf,
    control.compute = list(config = TRUE),
    control.mode = list(theta = theta1, fixed = TRUE)
)

all.equal(Q1c, cgeneric_Q(fit0))

