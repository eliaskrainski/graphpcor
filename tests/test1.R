
library(corGraphs)
library(INLA)

dag1 <- list(p1 ~ c1 + c2 + c3)
np <- length(dag1)

(theta.p <- 0.5)

mcorr <- cov2cor(dag_covariance(dag1, theta.p))

(nc <- nrow(mcorr))
(theta.c <- seq(0.5, -0.5, length = nc))

ss <- diag(exp(-1.0 * theta.c))
mcov <- ss %*% mcorr %*% ss

round(100 * mcorr)
round(mcov, 1)

n <- 3000
m <- nrow(mcov)

ll <- chol(mcov)
xx <- matrix(rnorm(n * m), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(w1 = runif(n), w2 = runif(n))
dataf$y <- xx[,1] + xx[,2] * dataf$w1 + xx[,3] * dataf$w2

dataf$idx1 <- rep(1L, n)
dataf$idx2 <- rep(2L, n)
dataf$idx3 <- rep(3L, n)
dataf$repl <- 1:n

ff <- y ~ 0 +
    f(idx1, model = rGmodel, replicate = repl) +
    f(idx2, w1, copy = "idx1", replicate = repl) +
    f(idx3, w2, copy = "idx1", replicate = repl)

d1plot <- GraphPlot(dag1, base=0)

par(mar = c(1, 1, 1, 1))
plot(d1plot$gr, nodeAttrs = d1plot$nAttrs)

hfix <- list(prec = list(initial = 10, fixed = TRUE))

cGmodel <- cgeneric_dag_model(
    dag = dag1,
    lambda = 5,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.05, nc)
)

ff <- y ~ 0 +
    f(idx1, model = cGmodel, replicate = repl) +
    f(idx2, w1, copy = "idx1", replicate = repl) +
    f(idx3, w2, copy = "idx1", replicate = repl)

fit <- inla(
    formula = ff, 
    data = dataf,
    control.family = list(list(hyper = hfix)),
    verbose = !TRUE)

fit$cpu

rbind(true = c(theta.c, theta.p),
      cg = fit$mode$theta)

plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)

mcorr.fit <- cov2cor(dag_covariance(dag1, fit$mode$theta[nc+1:np]))

round(100 * cor(xx))
round(100 * mcorr.fit)

ss.fit <- diag(exp(-1.0 * fit$mode$theta[1:nc]))
mcov.fit <- ss %*% mcorr.fit %*% ss

round(cov(xx), 1)
round(mcov.fit, 1)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
