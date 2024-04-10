
library(corGraphs)
library(INLA)

inla.setOption(safe = FALSE,
               num.threads = 6)

dcg1 <- list(p1 ~ c1 + c2 + c3)
np <- length(dcg1)

(theta.p <- 0.5)

mcov0 <- dcg_covariance(dcg1, theta.p)
mcov0
mcorr <- cov2cor(mcov0)
round(mcorr, 2)

(nc <- nrow(mcorr))
(theta.c <- seq(0.5, -0.5, length = nc))

ss <- diag(exp(theta.c))
mcov <- ss %*% mcorr %*% ss

round(100 * mcorr)
round(mcov, 1)

n <- 5000
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

d1plot <- GraphPlot(dcg1, base=0)

par(mar = c(1, 1, 1, 1))
plot(d1plot$gr, nodeAttrs = d1plot$nAttrs)

hfix <- list(prec = list(initial = 10, fixed = TRUE))

gmodel <- dcg_model(
    dcg = dcg1,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.05, nc),
    lambda = 2,
    iprior = 3,
    useINLAprecomp = FALSE,
    debug =  0 ## bigger debug prints llllooootttttsss of details
)

ff <- y ~ 0 +
    f(idx1, model = gmodel, replicate = repl) +
    f(idx2, w1, copy = "idx1", replicate = repl) +
    f(idx3, w2, copy = "idx1", replicate = repl)

fit <- inla(
    formula = ff, 
    data = dataf,
    control.family = list(list(hyper = hfix)),
##    control.mode = list(theta = rep(0, 4), restart = !FALSE, fixed = TRUE),
    verbose = !TRUE)

fit$cpu.used

rbind(true = c(theta.c, theta.p),
      cg = fit$mode$theta)

plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)

mcorr.fit <- cov2cor(dcg_covariance(dcg1, fit$mode$theta[nc+1:np]))

round(100 * cor(xx))
round(100 * mcorr.fit)

ss.fit <- diag(exp(fit$mode$theta[1:nc]))
mcov.fit <- ss %*% mcorr.fit %*% ss

round(cov(xx), 1)
round(mcov.fit, 1)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)

