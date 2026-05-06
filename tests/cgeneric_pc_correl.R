library(graphpcor)

library(INLA)

## n = 1, m = 1
C1 <- cgeneric(
    model = "pc_correl",
    n = 2L,
    lambda = 3.0
)

cgeneric_initial(C1)
cgeneric_graph(C1)
cgeneric_Q(C1, theta = 0)
cgeneric_Q(C1, theta = -1)
cgeneric_Q(C1, theta = 1)
mu(C1)

cgeneric_prior(C1, theta = 0.0)

integrate(function(x) exp(cgeneric_prior(C1, theta = matrix(x, 1))), -5, 5)

plot(function(x) exp(cgeneric_prior(C1, theta = matrix(x, 1))), -5, 5)

par(mfrow = c(3,3), mar=c(3,3,1,1), mgp = c(2,.5,0), bty = 'n')
for(l in 10^(-3:5)) {
    b <- .1 + 1/sqrt(l)
    cl <- cgeneric("pc_correl", n = 2L, base = 0, lambda = l)
    print(integrate(function(x) exp(cgeneric_prior(cl, theta = matrix(x, 1))), -b, b))
    plot(function(x) exp(cgeneric_prior(cl, theta = matrix(x, 1))), -b, b, n=1001)
}

## n =3, m = 3
c0 <- matrix(c( 1.0,  0.8, -0.5,
                0.8,  1.0, -0.4,
               -0.5, -0.4,  1.0), 3)
c0

(p <- ncol(c0))

b0 <- basecor(c0)
(th0b <- b0$theta)

h0 <- 0.02
th0 <- seq(-2+h0/2, 2-h0/2, h0)
(nth0 <- length(th0))
ncol(ths <- t(expand.grid(
         th1 = th0 + th0b[1],
         th2 = th0 + th0b[2],
         th3 = th0 + th0b[3])))

dim(p1ths <- array(exp(cgeneric_prior(C1, theta = ths)), rep(nth0, 3)))

sum(h0*apply(p1ths * (h0^2), 1:2, sum))
sum(h0*apply(p1ths * (h0^2), c(1,3), sum))
sum(h0*apply(p1ths * (h0^2), 2:3, sum))

## n = 4, m = 6
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
    lambda = lambda
)

Cmodel

cgeneric_graph(Cmodel, optimize = TRUE)

cgeneric_graph(Cmodel)

cgeneric_initial(Cmodel)

all.equal(as.matrix(solve(
  cgeneric_Q(Cmodel, theta = theta1))),
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

all.equal(base$base, as.matrix(solve(cgeneric_Q(fit0))))

n
(sigmas <- c(2,1,0.5,0.1))
V <- diag(sigmas) %*% base$base %*% diag(sigmas)

ns <- 200
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
    base = theta1, ##debug = 10,
    lambda = lambda,
    sigma.prior.reference = 1,#rep(1, n),
    sigma.prior.probability = 0.5
)

cgeneric_initial(Vmodel)
cgeneric_graph(Vmodel)
cgeneric_Q(Vmodel, theta = rep(0,10))
cgeneric_Q(Vmodel, theta = c(rep(0,4), rep(-1,6)))

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
