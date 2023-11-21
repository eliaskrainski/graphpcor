
library(corGraphs)
library(INLA)

inla.setOption(safe = FALSE,
               num.threads = 6)

S <- list(p1 ~ p2 + p3 + c1, p2 ~ c2 + c3 + p4, p3 ~ c4, p4 ~ c5)
SP_plot <- GraphPlot(S, base=0)

par(mar = c(1, 1, 1, 1))
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

SP <- GraphDens(S)
names(SP)

theta0 <- c(-0.5, 0, 0.5, 1)
theta1 <- -2:2
mcov <- ThetaCor(
    SP, c(theta1, theta0), 
    COV = TRUE)
round(mcov, 3)

n <- 500
m <- nrow(mcov)

ll <- chol(mcov)
xx <- matrix(rnorm(n * m), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(
    i = rep(1:m, eac = n),
    r = rep(1:n, m),
    y = c(rpois(n, exp(1 + xx[, 1])),
          rpois(n, exp(2 + xx[, 2])),
          rpois(n, exp(3 + xx[, 3])),
          rpois(n, exp(1 + xx[, 4])),
          rpois(n, exp(2 + xx[, 5])))
)

ff <- y ~ 0 + factor(i) +
    f(i, model = rGmodel, replicate = r, vb.correct = FALSE)

ArgsList <- list(
    S = S,
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
        theta = c(theta0, theta1), 
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    num.threads = 6,
    verbose = !TRUE)

rbind(c(theta1, theta0), 
      fit$mode$theta)

cGmodel <- cgeneric_dag_model(
    dag = S,
    lambda = 5,
    sigma.prior.reference = rep(1, m),
    sigma.prior.probability = rep(0.1, m)
)

cff <- y ~ 0 + factor(i) +
    f(i, model = cGmodel, replicate = r, vb.correct = FALSE)

cfit <- inla(
    formula = cff,
    family = "poisson",
    data = dataf,
    control.mode = list(
        theta = c(theta0, theta1), 
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    num.threads = 6,
    verbose = !TRUE)

rbind(fit$cpu, cfit$cpu)

rbind(c(theta1, theta0), 
      fit$mode$theta,
      cfit$mode$theta)

tail(fit$logfile, 30)
tail(cfit$logfile, 30)

round(fit$summary.fixed, 2)
round(fit$summary.hyperpar, 2)

ESTcov <- ThetaCor(
    SP, fit$summary.hyperpar$mean,
    COV = FALSE)
colnames(ESTcov) <- rownames(ESTcov) <-
    SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]

round(mcov, 2)
round(cov(xx), 2)
round(ESTcov, 2)

si <- diag(exp(-0.5 * cfit$mode$theta[1:m]))

dag_precision(S, rep(1, SP$NC))

qq <- dag_precision(S, cfit$mode$theta[-(1:m)])

round(qq, 3)

cc <- cov2cor(solve(qq)[1:m, 1:m])
vv <- si%*%cc%*%si

round(mcov, 2)
round(cc, 2)
round(vv, 2)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
