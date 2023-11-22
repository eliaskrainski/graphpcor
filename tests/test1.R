
library(corGraphs)
library(INLA)

inla.setOption(safe = FALSE,
               num.threads = 6)

dag1 <- list(p1 ~ c1 + c2 + c3)
np <- 1
nc <- 3

theta.p <- 0
theta.c <- c(-1, 0, 1)

q1 <- dag_precision(dag1, theta.p)
q1

mcov0 <- solve(q1)[1:nc, 1:nc]
mcor <- cov2cor(mcov0)

round(100 * mcor)

mcov <- diag(exp(-0.5*theta.c)) %*% mcor %*% diag(exp(-0.5*theta.c))
mcov

n <- 300
m <- nrow(mcov)

ll <- chol(mcov)
xx <- matrix(rnorm(n * m), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(w1 = runif(n), w2 = runif(n))
dataf$y <- xx[,1] + xx[,2] * dataf$w1 + xx[,3] * dataf$w2

dataf$idx1 <- rep(1L, n)
dataf$idx2 <- rep(2L, n)
dataf$idx3 <- rep(3L, n)
dataf$repl <- 1:n

ff <- y ~ 0 +
    f(idx1, model = rGmodel, replicate = repl) +
    f(idx2, w1, copy = "idx1", replicate = repl) +
    f(idx3, w2, copy = "idx1", replicate = repl)

d1plot <- GraphPlot(dag1, base=0)

par(mar = c(1, 1, 1, 1))
plot(d1plot$gr, nodeAttrs = d1plot$nAttrs)

SP <- GraphDens(dag1)
names(SP)

ArgsList <- list(
    S = dag1,
    lambda = 7,
    SP = SP,
    Tdist = Tdist,
    GraphPrior = GraphPrior,
    init = 0,
    new = FALSE)

str(ArgsList)

rGmodel <- inla.rgeneric.define(
    corGraphs_rgeneric, 
    args = ArgsList)

hfix <- list(prec = list(initial = 10, fixed = TRUE))

fit <- inla(
    formula = ff,
    data = dataf,
    control.family = list(list(hyper = hfix)),
    control.mode = list(theta = rep(1, nc+np), restart = TRUE),
    verbose = !TRUE)

cGmodel <-cgeneric_dag_model(
    dag = dag1,
    lambda = 5, 
    sigma.prior.reference = rep(1, m),
    sigma.prior.probability = rep(0.1, m))

ffc <- y ~ 0 +
    f(idx1, model = cGmodel, replicate = repl) +
    f(idx2, w1, copy = "idx1", replicate = repl) +
    f(idx3, w2, copy = "idx1", replicate = repl)

cfit <- inla(
    formula = ffc, 
    data = dataf,
    control.family = list(list(hyper = hfix)),
    control.mode = list(theta = rep(1, nc+np), restart = TRUE, fixed = !TRUE),
    verbose = !TRUE)

rbind(fit$cpu, cfit$cpu)
c(max(fit$misc$nfunc), max(cfit$misc$nfunc))

rbind(true = c(theta.c, theta.p),
      rg = fit$mode$theta,
      cg = cfit$mode$theta)

plot(fit, F, F, F, F, F, F, plot.opt.trace = TRUE)
plot(cfit, F, F, F, F, F, F, plot.opt.trace = TRUE)

tail(fit$logfile, 30)
tail(cfit$logfile, 30)

ctest <- inla(
    formula = ffc, 
    data = dataf,
    control.family = list(list(hyper = hfix)),
    control.mode = list(theta = rep(0, nc+np), restart = !TRUE, fixed = TRUE),
    verbose = TRUE)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
