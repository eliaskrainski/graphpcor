
library(corGraphs)
library(INLA)

S <- list(p1 ~ p2 + p3 + c5, p2 ~ c1 + c2 + p4, p3 ~ c3, p4 ~ c4)
SP_plot <- GraphPlot(S, base=0)

par(mar = c(1, 1, 1, 1))
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

SP <- GraphDens(S)
names(SP)

mcov <- ThetaCor(
    SP, c(rep(0, 5), rep(1, 4)), 
    COV = TRUE)
round(mcov, 3)

mcor <- mcov

n <- 300
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
    init = 0)

str(ArgsList)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)

rGmodel <- inla.rgeneric.define(
    CorGraphs.rmodel,
    args = c(ArgsList, list(test = !FALSE)))

inla.setOption(safe = FALSE)

ofit <- fit

fit <- inla(
    formula = ff,
    family = "poisson",
    data = dataf,
    control.mode = list(
        theta = rep(c(0, 0),
                    c(5, 4)),
        restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    num.threads = 6,
    verbose = TRUE)

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

library(numDeriv)
round(hessian(SP$JD[[1]], rep(0, 9), SDev = rep(1, 5)), 4)

