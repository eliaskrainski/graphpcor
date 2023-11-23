
library(corGraphs)
library(INLA)

inla.setOption(safe = FALSE,
               num.threads = 6)

dag <- list(
    p1 ~ p2 + c1 + c2,
    p2 ~ c3 + c4)
np <- length(dag)
nc <- 4

dgplot <- GraphPlot(dag, base=0)

par(mar = c(1, 1, 1, 1))
plot(dgplot$gr, nodeAttrs = dgplot$nAttrs)

(np <- length(dag))
(theta.p <- c(0, 0))

mcorr <- cov2cor(dag_covariance(dag, theta.p))
round(mcorr * 100)

(nc <- nrow(mcorr))
(theta.c <- (0.5:nc - nc/2)/2)

dd <- diag(exp(theta.c))
mcov <- dd %*% mcorr %*%dd

round(mcov, 1)

n <- 500

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

cGmodel <- cgeneric_dag_model(
    dag = dag,
    lambda = 5,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc),
    iprior = 3
)

ff <- y ~ 0 + factor(i) +
    f(i, model = cGmodel, replicate = r, vb.correct = FALSE)

fit <- inla(
    formula = ff,
    family = "poisson",
    data = dataf,
    control.mode = list(
        theta = rep(c(-2, 1), c(nc, np)), 
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE) ### if true prints looooooottttssss of details

fit$cpu

rbind(true = c(theta.c, theta.p), 
      cg = fit$mode$theta)

plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)

tail(fit$logfile, 30)

cc.fit <- cov2cor(dag_covariance(dag, fit$mode$theta[nc+1:np]))

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(fit$mode$theta[1:nc]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)


detach("package:corGraphs", unload = TRUE)
library(corGraphs)
