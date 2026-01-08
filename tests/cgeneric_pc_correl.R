library(graphpcor)
library(INLA)

## n = 1, m = 1
C1 <- cgeneric(
    model = "pc_correl",
    n = 2L,
    lambda = 3.0,
    useINLAprecomp = FALSE
)

integrate(function(x) exp(prior(C1, theta = matrix(x, 1))), -5, 5)

plot(function(x) exp(prior(C1, theta = matrix(x, 1))), -5, 5)

## n =3, m = 3
c0 <- matrix(c(1,.8,-.625, 0.8,1,-.5, -0.625,-.5,1), 3)
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

dim(p1ths <- array(exp(prior(M1, theta = ths)), rep(nth0, 3)))
p5ths <- array(exp(prior(M5, theta = ths)), rep(nth0, 3))

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
