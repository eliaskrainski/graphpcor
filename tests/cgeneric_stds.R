library(INLA)
library(graphpcor)

## n = 2
stds2 <- cgeneric(
    model = "stds",
    n = 2L,
    sigma.prior.probability = rep(0.1,2)
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

h0 <- 0.05
th0s <- seq(-5,5,h0); n0 <- length(th0s)
xx <- t(expand.grid(th0s, th0s))
lxx <- cgeneric_prior(stds2, xx)

image(th0s, th0s, matrix(lxx, n0))

sum(apply(matrix(exp(lxx)*h0, n0), 1, sum)*h0)

par(mfrow=c(1,2))
plot(function(x) exp(cgeneric_prior(stds2, theta = rbind(x, 0))), -5, 5)
plot(function(x) exp(cgeneric_prior(stds2, theta = rbind(0, x))), -5, 5)

## n = 10, m = 4
s10i4f2 <- cgeneric(
    model = 'stds',
    n = 10L,
    iparams = c(1,1,2,2,2,3,3,4,4,4),
    sigma.prior.probability = c(0.1,0,0.1,NA)
)

stds10f3
str(stds10f3)

cgeneric_initial(s10i4f2)
cgeneric_graph(s10i4f2)

cgeneric_Q(s10i4f2, theta = 0:1)
cgeneric_Q(s10i4f2, theta = -1:0)

