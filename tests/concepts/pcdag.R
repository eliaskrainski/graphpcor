library(graphpcor)

g <- dag(x ~ y, y ~ v, v ~ z, z ~ x)

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
          Laplacian(dag(x~y+z,y~v,v~z)))

## compact, but different ordering 
Laplacian(dag(x~y+z,v~y+z))

g <- dag(x1~x2+x3, x2~x4, x3~x4) ## compact ordered
(G <- Laplacian(g)) ## the graph in Example 2.6 of the GMRF book

dag(G) ## dag from a matrix

all.equal(dag(G), g) ## TRUE if compact ordered

## base model (theta for L)
theta0l <- rep(-0.5, ne[2])

## build a Cholesky (for precision)
L0 <- chol(g, theta = theta0l)
L0

round(tcrossprod(L0), 4)

## the precision for a correlation matrix
Q0 <- precision(g, theta = theta0l)
Q0

C0 <- solve(Q0)
C0

all.equal(C0, 
          variance(g, theta = theta0l)) ## variance method for dag

## the 'iid' case would be
variance(g, theta = rep(0, ne[2]))

## define full variance
sigmas <- c(0.3, 0.7, 1.2, 0.5)

variance(g, theta = c(log(sigmas), rep(0, ne[2])))
variance(g, theta = c(log(sigmas), rep(-1, ne[2])))
variance(g, theta = c(log(sigmas), rep(1, ne[2]))) ## no edge 2~3 but high covariance!!!

## build the cgeneric model
## Note: here 'model' is a 'dag' 
cmodel <- cgeneric(
    model = g, ## use the dag
    lambda = 1,
    base = theta0l,
    sigma.prior.reference = rep(1, ne[1]),
    sigma.prior.probability = rep(0.5, ne[1]))

## Note: another way: using the Laplacian
all.equal(
    cmodel,
    cgeneric(
        model = "pcdag", ## as character
        graph = G, ## using G as the 'corgraph' model
        lambda = 1,
        base = theta0l,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## yet another way: a binary matrix
all.equal(
    cmodel,
    cgeneric(
        model = "pcdag",
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

V1 <- variance(g, theta = theta1)
V1

C1 <- cov2cor(variance(g, theta = theta1))
C1

round(solve(V1), 2)

Q1c <- precision(cmodel, theta = theta1)
Q1c

round(Q1c, 2)

all.equal(solve(as.matrix(Q1c)), variance(g, theta = theta1))
all.equal(as.matrix(Q1c), precision(g, theta = theta1))

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

all.equal(Q1c, precision(fit0))

