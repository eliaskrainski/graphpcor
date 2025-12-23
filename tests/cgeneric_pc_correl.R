library(INLA)

library(graphpcor)

n <- 4
m <- n*(n-1)/2

theta1 <- rnorm(m,-1)
base <- basecor(theta1, n)
base

lambda <- 5
Cmodel <- cgeneric(
    model = "pc_correl",
    n = n,
    base = theta1,
    lambda = lambda,
    useINLAprecomp = FALSE)

str(Cmodel)

graph(Cmodel, optimize = TRUE)

graph(Cmodel)

initial(Cmodel)

all.equal(as.matrix(solve(prec(Cmodel, theta = theta1))),
          base$base)

dat1 <- data.frame(
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fit0 <- inla(
    y ~ 0 + f(i, model = Cmodel),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode,
    verbose = !TRUE
)

all.equal(base$base, as.matrix(solve(prec(fit0))))

n
(sigmas <- c(2,1,0.5,0.1))
V <- diag(sigmas) %*% base$base %*% diag(sigmas)

ns <- 100
xx <- matrix(rnorm(ns * n), ns) %*% chol(V)

V
cov(xx)
cor(xx)

dataf <- data.frame(
    y = as.vector(xx),
    i = rep(1:n, each = ns),
    r = rep(1:ns, n)
)

Vmodel <- cgeneric(
    model = "pc_correl",
    n = n,
    base = theta1,
    lambda = lambda,
    sigma.prior.reference = rep(1, n), 
    sigma.prior.probability = rep(0.5, n),
    useINLAprecomp = FALSE)

fit <- inla(
    y ~ 0 + f(i, model = Vmodel, replicate = r),
    data = dataf,
    control.family = cfam,
    verbose = !TRUE
)

sigmas
exp(fit$mode$theta[1:n])
m
cor(xx)
(Rfit <- basecor(fit$mode$theta[n+1:m], p = n)$base)

cov(xx)
(Vfit <- diag(exp(fit$mode$theta[1:n])) %*% Rfit %*%
    diag(exp(fit$mode$theta[1:n])))
