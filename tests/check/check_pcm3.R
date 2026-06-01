library(graphpcor)

gph_model <- graphpcor(
    graph = matrix(1, 3, 3) ## dense
)

C0 <- matrix(c( 1.0,  0.8, -0.5,
                0.8,  1.0, -0.4,
               -0.5, -0.4,  1.0), 3)
C0

hh <- hessian(gph_model, x = C0)

str(hh)

param <- 5
c_model <- cgeneric(
    model = gph_model,
    lambda = param,
    base = C0
)

## the base model represented at 
th0 <- c_model$f$cgeneric$data$doubles$thetabase
th0

### check: evaluating theta0 gives the base model
vcov(gph_model, theta = th0)
solve(cgeneric_Q(c_model, theta = th0))

## H^0.5: "SQRT" of the Hessian of KLD around the base model

h.5 <- matrix(c_model$f$cgeneric$data$matrices[[1]][-(1:2)], 3)
h.5 
graphpcor:::dspd(hessian(gph_model, th0))$sqrt

## its determinant
(dh.5 <- det(h.5))
(lc <- c_model$f$cgeneric$data$doubles$lc)
all.equal(log(abs(dh.5)), lc)

th1 <- c(.1,.2,-.3)
cgeneric_prior(c_model, theta = th1)

## set up a grid of theta around theta0
h <- 0.01
th0s <- seq(-1.5+h/2, 1.5-h/2, h)
th0x <- t(expand.grid(th0s, th0s, th0s)) ## each column a m-dimentional theta
ncol(th0x)/1e6 ## size of the grid

## eval, in R
R <- sqrt(colSums( (h.5 %*% th0x)^2 ))
Sr <- 2 * (pi^(3/2)) * R^2 / gamma(3/2)
dthR <- param * exp(-param * R) * abs(dh.5) / Sr
head(log(dthR))

## eval, in C (through 'cgeneric' query method in INLAtools)
dth <- exp(
    cgeneric_prior(
        c_model, 
        theta = rbind(th0x[1,] + th0[1],
                      th0x[2,] + th0[2],
                      th0x[3,] + th0[3])
    )
)

all.equal(dthR, dth)

## organize into an array
d3th <- array(dth, dim = rep(length(th0s), 3))

d1th <- apply(d3th * h*h, 1, sum)
d2th <- apply(d3th * h*h, 2, sum)
d3th <- apply(d3th * h*h, 3, sum)

c(sum(d1th*h), sum(d2th*h), sum(d3th*h))

par(mfrow = c(1, 3), bty = 'n')
plot(th0s + th0[1], d1th); abline(v = th0[1])
plot(th0s + th0[2], d2th); abline(v = th0[2])
plot(th0s + th0[3], d3th); abline(v = th0[3])

## sampling from the p(theta|lambda)
fakedata <- data.frame(y = NA, i = 1:c_model$f$n)
library(INLA)
out.inla <- inla(
    y ~ 0 + f(i, model = c_model),
    control.family = list(hyper = list(prec = list(initial = 10, fixed = TRUE))),
    data = fakedata)

str(inla.hyperpar.sample(1000, out.inla))

th.samples <- inla.hyperpar.sample(n = 2000, result = out.inla)

iil <- which(lower.tri(diag(ncol(C0))))
th2corr <- function(th) {
    vcov(gph_model, theta = th)[iil]
}

r.corr <- t(sapply(1:nrow(th.samples), function(i)
    th2corr(th.samples[i, ])))
str(r.corr)

par(mfrow = c(1, 3), bty = 'n')
for(i in 1:3) {
    hist(r.corr[, i], 50, freq = FALSE)
    abline(v=C0[c(2,3,6)[i]], col = 2, lwd = 4, lty = 2)
}

## sampling direclty from the sphere with radius r
