library(INLA)

library(graphpcor)

n <- 4
lambda <- 5

model <- cgeneric(
    model = "pc_prec_correl", 
    n = n,
    lambda = lambda)

class(model)

graph(model, optimize = TRUE)

graph(model)

initial(model)
ith <- rep(3, n*(n-1)/2)

m <- n * (n-1)/2
theta1 <- rnorm(m)

(qq <- prec(model, theta = theta1))

(vv <- solve(qq))

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

## fit with no data
fit.fix <- inla(
    y ~ 0 + f(i, model = model),
    data = list(i = 1:n, y = rep(NA,n)),
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

all.equal(qq, prec(fit.fix))

### fit some data
nrep <- 200
xx <- matrix(rnorm(nrep * n), nrep) %*% as.matrix(chol(vv))
str(xx)

dat2 <- data.frame(
    i = rep(1:n, each = nrep),
    r = rep(1:nrep, n),
    y = as.vector(xx)
)
str(dat2)

fitr <- inla(
    y ~ 0 + f(i, model = model, replicate = r),
    data = dat2,
    control.family = cfam,
    control.inla = cinla
)

(Lfitted <- Lcorrel(fitr$mode$theta))
round(crossprod(Lfitted), 2)
round(vv, 2)

idxc <- which(lower.tri(diag(n)))
fncorr <- function(th) {
    crossprod(Lcorrel(th))[idxc]
}

scorrels <- t(inla.cgeneric.sample(
    n = 1000, result = fitr, name = 'i', 
    from.theta = fncorr, simplify = TRUE
))

par(mfrow = c(2, 3))
for(i in 1:6) {
    hist(scorrels[, i])
    abline(v = vv[idxc[i]])
}

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
