library(INLA)

library(graphpcor)

n <- 4
lambda <- 5

model <- cgeneric(
    model = "pc_correl",
    n = n,
    lambda = lambda)

graph(model, optimize = TRUE)

graph(model)

round(ith <- initial(model), 4)

m <- n * (n-1)/2
theta1 <- rnorm(m)

(qq <- prec(model, theta = theta1))

(vv <- solve(qq))

dat1 <- data.frame(
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fit <- inla(
    y ~ 0 + f(i, model = model),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

all.equal(qq, prec(fit))

pc1 <- prior(model, theta = theta1)
pc1

thetapi <- pi/(1+exp(-theta1))
pc0 <- INLA:::inla.pc.cormat.dtheta(thetapi, lambda = lambda, log = TRUE)
pc0 +sum(log(pi * exp(-theta1) / ( (1 + exp(-theta1))^2 ) )) ## 1st Jacobian

if(FALSE) {
    INLA:::inla.pc.cormat.dtheta
    INLA:::inla.pc.multvar.simplex.d
    INLA:::inla.pc.multvar.simplex.core
    INLA:::inla.pc.multvar.simplex.d.core
    INLA:::inla.pc.multvar.h.default
}

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
    control.inla = cinla,
    control.mode = cmode
)

(Lfitted <- Lprec(fitr$mode$theta))
round(tcrossprod(Lfitted), 2)
round(vv, 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
