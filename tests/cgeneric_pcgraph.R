library(INLA)
library(graphpcor)

### for concepts: test/concepts/pcgraph.R
### for details: tests/check/check_pcgraph.R
### for details: tests/detailed/detailed_pcgraph.R

## the graph in Example 2.6 of the GMRF book
g <- graph(x1 ~ x2+x3, x2~x4, x3~x4)
class(g)
g
summary(g)
Laplacian(g)

(ne <- dim(g))

## define the cgeneric model (see test/concepts/pcgraph.R for other options)
theta.base <- rep(-1, ne[2])
cmodel <- cgeneric(
    model = g,
    lambda = 1,
    base = theta.base,
    sigma.prior.reference = rep(1, ne[1]),
    sigma.prior.probability = rep(0.5, ne[1]))

sigmas <- c(5, 0.5, 1, 0.1)
thetaL <- c(-5, 1, 2, -0.1)
theta1 <- c(log(sigmas), thetaL)

Vg <- variance(g, theta = theta1)
Vg

cov2cor(variance(g, theta = theta.base)) ## base model correlation
cov2cor(Vg) ## correlation to be used to sample data

precision(cmodel, theta = theta1)

## some data
nrep <- 100
nd <- nrep * ne[1]

xx <- matrix(rnorm(nd), nrep) %*% chol(Vg)
cov(xx)

theta.y <- log(5)
datar <- data.frame(
    r = rep(1:nrep, ne[1]),
    i = rep(1:ne[1], each = nrep),
    y = rnorm(nd, 1 + xx, exp(-2*theta.y))
)

m1 <- y ~ f(i, model = cmodel, replicate = r)
fit <- inla(
    formula = m1,
    data = datar
)
fit$cpu.used

grep("function evaluations =", fit$logfile, value = TRUE)
grep("fn-calls=", fit$logfile, value = TRUE)

rbind(true = c(theta.y, theta1),
      fit = fit$mode$theta) 
fit$mode$theta - c(theta.y, theta1)

precision(cmodel, theta = fit$mode$theta[-1])

round(Vfit <- variance(g, theta = fit$mode$theta[-1]), 2)
round(Vg, 2)
round(cov(xx), 2)

cor(xx)
cov2cor(Vfit)
