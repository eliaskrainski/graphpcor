library(INLA)

library(graphpcor)

n <- 4
(m <- n*(n-1)/2)
eta <- 10

model <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta,
    debug = 9999L)

str(model)

graph(model, optimize = TRUE)

graph(model)

round(ith <- initial(model), 4)

theta1 <- rnorm(m)
x1 <- pi*plogis(theta1)

x1
cos(x1)
sin(x1)

(qq <- precision(model, theta = theta1))

(b1 <- graphpcor:::theta2gamma2L(theta1))
all.equal(as.matrix(qq), solve(tcrossprod(b1)))

(vv <- solve(qq))
all.equal(as.matrix(vv), tcrossprod(b1))

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

### 
nrep <- 2
xx <- matrix(rnorm(nrep * n), nrep) %*% chol(vv)
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
    ##, verbose = TRUE
)

round(precision(fitr), 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
