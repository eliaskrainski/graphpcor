
library(corGraphs)
library(INLA)

dcgp <- list(
    p1 ~ p2 + c1 + c2,
    p2 ~ p3 + c3 + c4,
    p3 ~ c5 + c6)

d2plot <- GraphPlot(dcgp, base=0)

par(mar = c(1, 1, 1, 1))
plot(d2plot$gr, nodeAttrs = d2plot$nAttrs)

dcg <- list(
    p1 ~ p2 + p3 + c1 + c2,
    p2 ~ -c3 + c4,
    p3 ~ c5 - c6)

(np <- length(dcg))
(theta.p <- c(0, -0.33, 0.33))

mcov <- corGraphs:::dcg_covariance(dcg, theta.p)
mcorr <- cov2cor(mcov)

nc <- nrow(mcov)
(theta.c <- rep(1, nc))

round(mcov, 1)
round(mcorr * 100)

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

gmodel <- dcg_model(
    dcg = dcg,
    lambda = 5,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc)
)

cff <- y ~ 0 + factor(i) +
    f(i, model = gmodel, replicate = r, vb.correct = FALSE)

cfit <- inla(
    formula = cff,
    ##    family = "poisson",
    control.family=list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
    data = dataf,
    control.mode = list(
        theta = rep(0, nc + np), 
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE) ### if true prints looooooottttssss of details

cfit$cpu

rbind(true = c(theta.c, theta.p), 
      cg = cfit$mode$theta)

plot(cfit, F, F, F, F, F, F, plot.opt.trace = TRUE)

cc.fit <- cov2cor(dcg_covariance(dcg, cfit$mode$theta[nc+1:np]))

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(cfit$mode$theta[1:nc]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
