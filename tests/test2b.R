
library(corGraphs)
library(INLA)

inla.setOption(safe = FALSE,
               num.threads = 6)

dcg <- list(
    p1 ~ p2 + c1 + c2,
    p2 ~ c3 + c4)
np <- length(dcg)
nc <- 4

dgplot <- GraphPlot(dcg, base=0)

par(mar = c(1, 1, 1, 1))
plot(dgplot$gr, nodeAttrs = dgplot$nAttrs)

dc.el <- corGraphs:::dcg_e2covariance(corGraphs:::dcg_elements(dcg))
dc.el

(np <- length(dcg))
(theta.p <- c(0, 1))
(theta.ch <- seq(-nc/2, nc/2, length = nc))

v2 <- exp(theta.p * 2.0)
q2 <- theta.ch^2
vp <- sapply(1:np, function(j)
    sum(v2[dc.el$iv[[j]]]))

mcov <- matrix(
    theta.ch[row(diag(nc))] *
    vp[dc.el$itop] *
    theta.ch[col(diag(nc))],
    nc) + diag(1/q2)
mcov

mcorr <- cov2cor(mcov)
round(mcorr * 100)

round(mcov, 1)

n <- 300

ll <- chol(mcov)
xx <- matrix(rnorm(n * nc), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(
    i = rep(1:nc, eac = n),
    r = rep(1:n, nc),
    y = as.vector(xx)##rpois(n * nc, exp(1 + xx))
)

gmodel <- dcg_model(
    dcg = dcg,
    lambda = 5,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc),
    iprior = 3,
    sfixed = FALSE,
    useINLAprecomp = FALSE,
    debug = 0 ### if debug>999 and inla(..., verbose = true) prints looooooottttssss of details
    )

str(gmodel)

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
    verbose = TRUE,
##    inla.call = "remote",
    control.inla = list(int.strategy = "eb")
) 

fit$cpu

rbind(true = c(theta.ch, theta.p),
      cg = fit$mode$theta)

plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)

tail(fit$logfile, 30)

v2fit <- exp(fit$mode$theta[nc + 1:np] * 2.0) 
q2fit <- fit$mode$theta[1:nc]^2 
vpfit <- sapply(1:np, function(j)
    sum(v2fit[dc.el$iv[[j]]]))

mcov.fit <- matrix(
    fit$mode$theta[row(diag(nc))] *
    vpfit[dc.el$itop] *
    fit$mode$theta[col(diag(nc))],
    nc) + diag(1/q2fit)

mcov
round(cov(xx), 2)
round(mcov.fit, 2)

cc.fit <- cov2cor(mcov.fit)

round(mcorr*100)
round(cor(xx)*100)
round(cc.fit*100)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
