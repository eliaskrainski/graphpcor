
library(corGraphs)
library(INLA)

inla.setOption(safe = FALSE,
               num.threads = 6)

dag <- list(
    p1 ~ p2 + p3 + c1 + c2,
    p2 ~ c3 + c4,
    p3 ~ c5 - c6)
np <- length(dag)
nc <- 6

(theta.c <- (1:nc - nc/2))
(theta.p <- c(0, 1, 0))

qq <- dag_precision(dag, theta.p)
mcorr <- cov2cor(solve(qq)[1:nc, 1:nc])
dd <- diag(exp(-0.5 * theta.c))
mcov <- dd %*% mcorr %*%dd

round(mcorr * 100)
round(mcov, 1)


n <- 500

ll <- chol(mcov)
xx <- matrix(rnorm(n * nc), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(
    i = rep(1:nc, eac = n),
    r = rep(1:n, nc),
    y = c(rpois(n, exp(3 + xx[, 1])),
          rpois(n, exp(3 + xx[, 2])),
          rpois(n, exp(3 + xx[, 3])),
          rpois(n, exp(3 + xx[, 4])),
          rpois(n, exp(3 + xx[, 5])),
          rpois(n, exp(3 + xx[, 6])))
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
    family = "poisson",
    data = dataf,
    control.mode = list(
        theta = c(theta.c, theta.p), 
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE) ### if true prints looooooottttssss of details

cfit$cpu

rbind(true = c(theta.c, theta.p), 
      cg = cfit$mode$theta)

plot(cfit, F, F, F, F, F, F, plot.opt.trace = TRUE)

qq.fit <- dag_precision(dag, cfit$mode$theta[nc+1:np])
cc.fit <- cov2cor(solve(qq.fit)[1:nc, 1:nc])

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(-0.5*cfit$mode$theta[1:nc]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
