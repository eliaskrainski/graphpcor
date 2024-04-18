
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

(np <- length(dcg))
(theta.p <- c(0, 0))

mcorr <- cov2cor(dcg_covariance(dcg, theta.p))
round(mcorr * 100)

(nc <- nrow(mcorr))
(theta.c <- (0.5:nc - nc/2)/2)

dd <- diag(exp(theta.c))
mcov <- dd %*% mcorr %*%dd

round(mcov, 1)

n <- 3000

ll <- chol(mcov)
xx <- matrix(rnorm(n * nc), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(
    i = rep(1:nc, eac = n),
    r = rep(1:n, nc),
    y = rpois(n * nc, exp(1 + xx))
)

ff <- y ~ 0 + factor(i) +
    f(i, model = rGmodel, replicate = r, vb.correct = FALSE)

gmodel <- dcg_model(
    dcg = dcg,
    lambda = 2,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc),
    iprior = 3,
    useINLAprecomp = FALSE,
    debug = 0### if debug>999 and inla(..., verbose = TRUE) prints looooooottttssss of details    
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

cc.fit <- cov2cor(dcg_covariance(dcg, fit$mode$theta[nc+1:np]))

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(fit$mode$theta[1:nc]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
