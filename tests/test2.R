
library(INLA)
library(corGraphs)

inla.setOption(
    safe = FALSE,
    num.threads = 6
)

g <- dtg(
    p1 ~ p2 + c1 + c2,
    p2 ~ c3 + c4)
g
d <- dim(g)
d

plot(g)

(theta.p <- rep(0, d[2]))

mcorr <- cov2cor(variance(g, theta = theta.p))
round(mcorr * 100)

(theta.c <- (0.5:d[1] - d[1]/2)/2)

dd <- diag(exp(theta.c))
mcov <- dd %*% mcorr %*%dd

round(mcov, 1)

n <- 3000

ll <- chol(mcov)
xx <- matrix(rnorm(n * d[1]), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(
    i = rep(1:d[1], each = n),
    r = rep(1:n, d[1]),
    y = rpois(n * d[1], exp(1 + xx))
)

gmodel <- cgeneric(
    model = g,
    lambda = 2,
    sigma.prior.reference = rep(1, d[1]),
    sigma.prior.probability = rep(0.1, d[1]),
    iprior = 3,
    useINLAprecomp = FALSE,
    debug = 1### if debug>999 and inla(..., verbose = TRUE) prints looooooottttssss of details    
)

ff <- y ~ 0 + factor(i) +
    f(i, model = gmodel, replicate = r, vb.correct = FALSE)

fit <- inla(
    formula = ff,
    family = "poisson",
    data = dataf,
    control.inla = list(int.strategy = "eb"),
##    control.mode = list(theta = rep(0, 6),  restart = FALSE, fixed = !TRUE),
    verbose = !TRUE) 

fit$cpu.used

fit$misc$nf

rbind(true = c(theta.c, theta.p), 
      cg = fit$mode$theta)

plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)

tail(fit$logfile, 30)

cc.fit <- cov2cor(variance(g, theta = fit$mode$theta[d[1]+1:d[2]]))

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(fit$mode$theta[1:d[1]]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
