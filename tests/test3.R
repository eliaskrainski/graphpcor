
library(corGraphs)
library(INLA)

dagp <- list(
    p1 ~ p2 + p3 + c1 + c2,
    p2 ~ c3 + c4,
    p3 ~ c5 + c6)

d2plot <- GraphPlot(dagp, base=0)

par(mar = c(1, 1, 1, 1))
plot(d2plot$gr, nodeAttrs = d2plot$nAttrs)

dag <- list(
    p1 ~ p2 + p3 + c1 + c2,
    p2 ~ -c3 + c4,
    p3 ~ c5 - c6)

(np <- length(dag))
(theta.p <- c(0, 1, -1))

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
    y = as.vector(xx)##rpois(n*nc, exp(1 + xx))
)

cGmodel <- cgeneric_dag_model(
    dag = dag,
    lambda = 5,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc)
)

cff <- y ~ 0 + factor(i) +
    f(i, model = cGmodel, replicate = r, vb.correct = FALSE)

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

cc.fit <- cov2cor(dag_covariance(dag, cfit$mode$theta[nc+1:np]))

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(cfit$mode$theta[1:nc]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
