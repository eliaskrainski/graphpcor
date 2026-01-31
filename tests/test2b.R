
library(graphpcor)
library(INLA)

inla.setOption(
    safe = FALSE,
    num.threads = 6
)

g <- treepcor(
    p1 ~ p2 + c1 + c2,
    p2 ~ c3 + c4 + p3,
    p3 ~ c5)
g
d <- dim(g)

nc <- d[1]
np <- d[2]

plot(g)

(theta.p <- seq(1/2, -1/2, length = np))

mcorr <- cov2cor(vcov(g, theta = theta.p))
round(mcorr * 100)

(theta.ch <- seq(-1/2, 1/2, length = nc))

dd <- diag(exp(theta.ch))
mcov <- dd %*% mcorr %*%dd

round(mcov, 1)

n <- 3000

ll <- chol(mcov)
xx <- matrix(rnorm(n * nc), n) %*% ll

cov(xx)
cor(xx)

c(n, nc)

dataf <- data.frame(
    i = rep(1:nc, each = n),
    r = rep(1:n, nc),
    y = as.vector(xx)##rpois(n * nc, exp(1 + xx))
)
head(dataf, 3)

gmodel <- cgeneric(
    model = g,
    lambda = 1,
    sigma.prior.reference = rep(5, nc),
    sigma.prior.probability = rep(0.2, nc),
    debug = 0 ### if debug>999 and inla(..., verbose = true) prints looooooottttssss of details
    )

ff <- y ~ 0 + factor(i) +
    f(i, model = gmodel, replicate = r, vb.correct = FALSE)

fit <- inla(
    formula = ff,
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
    data = dataf
)

fit$cpu.used

rbind(true = c(theta.ch, theta.p),
      cg = fit$mode$theta)


plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)

tail(fit$logfile, 30)

mcorr.fit <- cov2cor(vcov(
    g,
    theta = fit$mode$theta[nc+1:np]))
mcorr.fit

q.fit <- prec(gmodel, theta = fit$mode$theta)
mcov.fit <- solve(q.fit)

mcov
round(cov(xx), 2)
round(mcov.fit, 2)

cc.fit <- cov2cor(as.matrix(mcov.fit))

round(mcorr*100)
round(cor(xx)*100)
round(cc.fit*100)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
