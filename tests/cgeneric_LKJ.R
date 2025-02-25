library(INLA)

library(graphpcor)

n <- 4
(m <- n*(n-1)/2)
eta <- 10

model <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta)

graph(model, optimize = TRUE)

graph(model)

round(ith <- initial(model), 4)

theta1 <- rnorm(m)

(qq <- precision(model, theta = theta1))

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

all.equal(qq, precision(fit))

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

(Lfitted <- graphpcor:::theta2gamma2L(fitr$mode$theta))
round(tcrossprod(Lfitted), 2)
round(vv, 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
