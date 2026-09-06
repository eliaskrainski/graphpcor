library(INLA)
library(graphpcor)

## n = 2
stds2 <- cgeneric(
    model = "stds",
    n = 2L,
    sigma.prior.probability = rep(0.1,2),
    useINLAprecomp = FALSE
)

stds2
str(stds2)

cgeneric_initial(stds2)
cgeneric_graph(stds2)
cgeneric_Q(stds2, theta = 0:1)
cgeneric_Q(stds2, theta = -1:0)
cgeneric_Q(stds2, theta = 1:2)
cgeneric_mu(stds2)

cgeneric_prior(stds2, theta = cbind(0,0:1,-1:0,1:2))

cpc0s2 <- cgeneric(
    model = "pc_correl",
    n = 2L, lambda = 1,
    sigma.prior.probability = rep(0.1,2),
    cfixed = 1
)

cgeneric_initial(cpc0s2)
cgeneric_graph(cpc0s2)
cgeneric_Q(cpc0s2, theta = 0:1)
cgeneric_Q(cpc0s2, theta = -1:0)
cgeneric_Q(cpc0s2, theta = 1:2)
cgeneric_mu(cpc0s2)

cgeneric_prior(stds2, theta = cbind(0,0:1,-1:0,1:2))

n1 <- 1000
dataf <- data.frame(
    i = rep(1:2, each = n1),
    r = rep(1:n1, 2L),
    y = rnorm(2*n1, 0, rep(c(0.5,4), each = n1))
)
str(dataf)
tapply(dataf$y, dataf$i, sd)

fit <- inla(
###    formula = y ~ f(i, model = stds2, replicate = r),
    formula = y ~ f(i, model = cpc0s2, replicate = r),
    data = dataf,
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
)

exp(fit$mode$theta)

h0 <- 0.05
th0s <- seq(-5,5,h0); n0 <- length(th0s)
xx <- t(expand.grid(th0s, th0s))
lxx <- cgeneric_prior(stds2, xx)

image(th0s, th0s, matrix(lxx, n0))

sum(apply(matrix(exp(lxx)*h0, n0), 1, sum)*h0)

par(mfrow=c(1,2))
plot(function(x) exp(cgeneric_prior(stds2, theta = rbind(x, 0))), -5, 5)
plot(function(x) exp(cgeneric_prior(stds2, theta = rbind(0, x))), -5, 5)

## n = 4, graphpcor
stds4 <- cgeneric(
    model = 'stds', n = 4,
    prior.sigma.probability = rep(0.05, 4)
)

graphpcor(3L*4L, p = 4L)
Laplacian(graphpcor(3L*4L, p = 4L))

cgstds4 <- cgeneric(
    model = graphpcor(3L*4L, p = 4L),
    sigma.prior.probability = rep(0.1, 4),
    lambda = 1, cfixed = 1,
    d0 = rep(1, 4)
)

data4 <- data.frame(
    i = rep(1:4, each = n1),
    r = rep(1:n1, 4L),
    y = rnorm(4*n1, 0, rep(c(0.5,1,2,4), each = n1))
)
str(data4)
tapply(data4$y, data4$i, sd)

fit4 <- inla(
    formula = y ~ 0 + f(i, model = cgstds4, replicate = r),
    data = data4,
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
)

exp(fit4$mode$theta)

## n = 10, m = 4
np <- 4
ip <- c(1,1,2,2,2,3,3,4,4,4)
s10i4f2 <- cgeneric(
    model = 'stds',
    n = 10L,
    iparams = ip, 
    sigma.prior.probability = c(0.1,0,0.1,NA),
    useINLAprecomp = FALSE
)

s10i4f2
str(s10i4f2)

cgeneric_graph(s10i4f2)

cgeneric_initial(s10i4f2)

cgeneric_Q(s10i4f2, theta = 0:1)

cgeneric_Q(s10i4f2, theta = -1:0)

## fit some data
n <- s10i4f2$f$n
nd <- 1000
sds <- c(0.5,1,2,4)
dataf <- data.frame(
    i = sort(c(1:n, sample(1:n, nd-n, replace=TRUE)))
)
dataf$j <- ip[dataf$i]
dataf$y <- rnorm(nd, 0, sds[dataf$j])

(ni <- table(dataf$i))
dataf$repl <- unlist(lapply(ni, function(n) 1:n))

all(tapply(dataf$repl, dataf$i, max)==ni)

rbind(true = sds[ip],
      sample = tapply(dataf$y, dataf$i, sd))

fit <- inla(
    formula = y ~ 0 + f(i, model = s10i4f2, replicate = repl),
    data = dataf,
    control.family = list(hyper = list(prec = list(initial=10, fixed=TRUE)))
)

rbind(true = sds[c(1,3)],
      fitt = exp(fit$mode$theta))
