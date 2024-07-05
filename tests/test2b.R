
library(corGraphs)
library(INLA)

inla.setOption(
    safe = FALSE,
    num.threads = 6
)

dcg <- list(
    p1 ~ p2 + c1 + c2,
    p2 ~ c3 + c4 + p3,
    p3 ~ c5)
(np <- length(dcg))
(theta.p <- seq(1/2, -1/2, length = np))

mcorr <- cov2cor(dcg_covariance(dcg, theta.p))
round(mcorr * 100)

(nc <- nrow(mcorr))

dgplot <- GraphPlot(dcg, base=0)

par(mar = c(1, 1, 1, 1))
plot(dgplot$gr, nodeAttrs = dgplot$nAttrs)

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

gmodel <- dcg_model(
    dcg = dcg,
    lambda = 1,
    sigma.prior.reference = rep(5, nc),
    sigma.prior.probability = rep(0.2, nc),
    iprior = 3,
    debug = 0 ### if debug>999 and inla(..., verbose = true) prints looooooottttssss of details
    )

ff <- y ~ 0 + factor(i) +
    f(i, model = gmodel, replicate = r, vb.correct = FALSE)

fit <- inla(
    formula = ff,
##    family = "poisson",
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
    data = dataf,
  #  control.mode = list(
#        theta = c(theta.ch, theta.p)*0.1,
 #       restart = TRUE, fixed = !TRUE),
##        restart = !TRUE, fixed = TRUE),
##    verbose = TRUE,
##    inla.call = "remote",
    control.inla = list(int.strategy = "eb")
) 

fit$cpu.used

rbind(true = c(theta.ch, theta.p),
      cg = fit$mode$theta)


plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)

tail(fit$logfile, 30)

mcorr.fit <- cov2cor(dcg_covariance(
    dcg,
    fit$mode$theta[nc+1:np]))
mcorr.fit

q.fit <- cgeneric_get(gmodel, "Q", fit$mode$theta, FALSE)
mcov.fit <- solve(q.fit)

mcov
round(cov(xx), 2)
round(mcov.fit, 2)

cc.fit <- cov2cor(mcov.fit)

round(mcorr*100)
round(cor(xx)*100)
round(cc.fit*100)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
