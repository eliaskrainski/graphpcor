## for details see
## https://link.springer.com/article/10.1007/s10260-025-00788-y

library(INLA)
library(graphpcor)

## define the tree graph
g <- treepcor(
    p1 ~ p2 + c1 - c2,
    p2 ~ c3 + c4)
class(g)
g
summary(g)

(np <- dim(g))

sigmas <- c(4, 2, 1, 0.5)
thetal <- c(0, 1)
theta1 <- c(log(sigmas), thetal)

Vg <- vcov(g, theta = theta1)
Vg

cov2cor(Vg)

(Qg <- prec(cmodel, theta = theta1))
all.equal(Vg, as.matrix(solve(Qg)))

## some data
nrep <- 100
nd <- nrep * np[1]

xx <- matrix(rnorm(nd), nrep) %*% chol(Vg)
cov(xx)

theta.y <- log(10)
datar <- data.frame(
    r = rep(1:nrep, np[1]),
    i = rep(1:np[1], each = nrep),
    y = rnorm(nd, 1 + xx, exp(-2*theta.y))
)

## define the cgeneric model
cmodel <- cgeneric(
    model = g,
    lambda = 10,
    sigma.prior.reference = rep(1, np[1]),
    sigma.prior.probability = rep(0.05, np[1]),
    debug = TRUE)

graph(cmodel)
initial(cmodel)
prior(cmodel, theta = rep(0, sum(np)))
prior(cmodel, theta = rep(1, sum(np)))

np
prec(cmodel, theta = rep(0, sum(np)))
(Qc <- prec(cmodel, theta = theta1))
all.equal(Vg, as.matrix(solve(Qc)))

m1 <- y ~ f(i, model = cmodel, replicate = r)
pfix <- list(prec = list(initial = 10, fixed = TRUE))
fit <- inla(
    formula = m1,
    data = datar,
    control.family = list(hyper = pfix)
)
fit$cpu.used

grep("function evaluations =", fit$logfile, value = TRUE)
grep("fn-calls=", fit$logfile, value = TRUE)

rbind(true = c(theta1),
      fit = fit$mode$theta)

## true covariance
round(Vg, 2)

## observed covariance
round(cov(xx), 2)

## fitted covariance (from posterior mode)
round(Vfit <- vcov(g, theta = fit$mode$theta), 2)

## fitted correlation (from posterior mode)
round(Cfit <- vcov(g, theta = fit$mode$theta[np[1]+1:np[2]]), 2)
round(cov2cor(Vfit), 2)

## true correlation
round(cov2cor(Vg), 2)

## observed correlation
round(cor(xx), 2)

