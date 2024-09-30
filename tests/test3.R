
library(graphpcor)
library(INLA)

g <- dtg(
    p1 ~ p2 + p3 + c1 + c2,
    p2 ~ -c3 + c4,
    p3 ~ c5 -c6)

g

plot(g)

d <- dim(g)
d

(theta.p <- c(0, -0.33, 0.33))

mcov0 <- variance(g, theta = theta.p)
mcorr <- cov2cor(mcov0)
round(mcorr,2)

nc <- d[1]
(theta.c <- seq(-1, 1, length = nc))

dd <- diag(exp(theta.c))
mcov <- dd %*% mcorr %*% dd

dd^2
round(mcov, 1)

n <- 500

ll <- chol(mcov)
xx <- matrix(rnorm(n * nc), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(
    i = rep(1:nc, eac = n),
    r = rep(1:n, nc),
    y = as.vector(xx)##rpois(n*nc, exp(1 + xx))
)

gmodel <- cgeneric(
    model = g,
    lambda = 5,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc)
)

initial(gmodel)

cff <- y ~ 0 + factor(i) +
    f(i, model = gmodel, replicate = r, vb.correct = FALSE)

cfit <- inla(
    formula = cff,
    control.family=list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
    data = dataf,
    control.mode = list(
        theta = rep(0, d[1] + d[2])
    ), 
    control.inla = list(int.strategy = "eb")
    )

cfit$cpu.used

rbind(true = c(theta.c, theta.p), 
      cg = cfit$mode$theta)

plot(cfit, F, F, F, F, F, F, plot.opt.trace = TRUE)

cc.fit <- cov2cor(variance(g, theta = cfit$mode$theta[nc+1:d[2]]))

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(cfit$mode$theta[1:nc]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
