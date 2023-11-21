
library(corGraphs)
library(INLA)

inla.setOption(safe = FALSE,
               num.threads = 6)

mcor <- matrix(c(1, 0.9, 0.9, 1), 2)
sds <- c(2, 0.7)
mcov <- t(mcor * sds) * sds
mcov

n <- 3000
m <- nrow(mcov)

ll <- chol(mcov)
xx <- matrix(rnorm(n * m), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(w = runif(n))
dataf$y <- xx[,1] + xx[,2] * dataf$w

dataf$idx1 <- rep(1L, n)
dataf$idx2 <- rep(2L, n)
dataf$repl <- 1:n

ff <- y ~ 0 +
    f(idx1, model = rGmodel, replicate = repl) +
    f(idx2, w, copy = "idx1", replicate = repl)

S <- list(p1 ~ c1 + c2)

SP_plot <- GraphPlot(S, base=0)

par(mar = c(1, 1, 1, 1))
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

SP <- GraphDens(S)
names(SP)

ArgsList <- list(
    S = S,
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
    control.mode = list(theta = c(1,1,1), restart = TRUE),
    control.family = list(list(hyper = hfix)),
    verbose = !TRUE)

sds
exp(-0.5 * fit$mode$theta[1:2])

cGmodel <-cgeneric_dag_model(
    dag = S,
    lambda = 5, 
    sigma.prior.reference = rep(1, m),
    sigma.prior.probability = rep(0.1, m))

ffc <- y ~ 0 +
    f(idx1, model = cGmodel, replicate = repl) +
    f(idx2, w, copy = "idx1", replicate = repl)

cfit <- inla(
    formula = ffc, 
    data = dataf,
    control.family = list(list(hyper = hfix)),
    verbose = !TRUE)

rbind(fit$cpu, cfit$cpu)

rbind(fit$mode$theta,
      cfit$mode$theta)

sds
exp(-0.5 * cfit$mode$theta[1:2])

##rbind(true = c(log(diag(cov(xx)))/2, theta0),
  ##    mode = res$mode$theta)


##summary(fit)

ESTcor <- ThetaCor(
    SP, fit$summary.hyperpar$mean,
    COV = FALSE)
colnames(ESTcor) <- rownames(ESTcor) <-
    SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]
print(ESTcor)

si <- diag(-0.5 * cfit$mode$theta[1:2])
qq <- dag_precision(S, cfit$mode$theta[3])
qq
cc <- cov2cor(solve(qq)[1:2, 1:2])
vv <- si%*%cc%*%si

mcov
vv[1:2, 1:2]

cc

mcov
cov2cor(mcov)[1,2]
