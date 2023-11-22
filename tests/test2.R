
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

SP <- GraphDens(dag)
names(SP)

(np <- length(dag))
(theta.p <- seq(0.5, -0.5, length = np))

mcorr <- cov2cor(dag_covariance(dag, theta.p))

(nc <- nrow(mcorr))
(theta.c <- (0.5:nc - nc/2)/2)

dd <- diag(exp(-1.0 * theta.c))
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
    y = rpois(n * nc, exp(3 + xx))
)

ff <- y ~ 0 + factor(i) +
    f(i, model = rGmodel, replicate = r, vb.correct = FALSE)

ArgsList <- list(
    S = dag,
    lambda = 7,
    SP = SP,
    Tdist = Tdist,
    GraphPrior = GraphPrior,
    init = 0,
    new = FALSE)

##str(ArgsList)

rGmodel <- inla.rgeneric.define(
    corGraphs_rgeneric,
    args = c(ArgsList, list(test = !FALSE)))

fit <- inla(
    formula = ff,
    family = "poisson",
    data = dataf,
    control.mode = list(
        theta = c(theta.c, theta.p), 
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE)

c0Gmodel <- cgeneric0_dag_model(
    dag = dag,
    lambda = 5,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc)
)

c0ff <- y ~ 0 + factor(i) +
    f(i, model = c0Gmodel, replicate = r, vb.correct = FALSE)

c0fit <- inla(
    formula = c0ff,
    family = "poisson",
    data = dataf,
    control.mode = list(
        theta = c(theta.c, theta.p), 
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE) ### if true prints looooooottttssss of details

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

rbind(fit$cpu, c0fit$cpu, cfit$cpu)

rbind(true = c(theta.c, theta.p), 
      rg = fit$mode$theta,
      c0g = c0fit$mode$theta,
      cg = cfit$mode$theta)

plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)
plot(c0fit, F, F, F, F, F, F, plot.opt.trace = TRUE)
plot(cfit, F, F, F, F, F, F, plot.opt.trace = TRUE)

tail(fit$logfile, 30)
tail(c0fit$logfile, 30)
tail(cfit$logfile, 30)

cc.fit <- cov2cor(dag_covariance(dag, cfit$mode$theta[nc+1:np]))

round(cor(xx)*100)
round(cc.fit*100)

ss.fit <- diag(exp(-1.0*cfit$mode$theta[1:nc]))
mcov.fit <- ss.fit %*% cc.fit %*% ss.fit

round(cov(xx), 2)
round(mcov.fit, 2)


detach("package:corGraphs", unload = TRUE)
library(corGraphs)
