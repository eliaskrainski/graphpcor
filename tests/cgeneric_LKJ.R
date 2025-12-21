library(INLA)
library(graphpcor)

n <- 4
(m <- n*(n-1)/2)
eta <- 10

Cmodel <- cgeneric(
    model = "LKJ",
    n = n,
    eta = eta,
    useINLAprecomp = FALSE
)

str(Cmodel)

graph(Cmodel, optimize = TRUE)

graph(Cmodel)

initial(Cmodel)

theta1 <- rnorm(m)

(qc <- prec(Cmodel, theta = theta1))

(cc <- solve(qc))

all.equal(as.matrix(cc),
          basecor(theta1, p=n)$base)

## fake data
dat1 <- data.frame(
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fit <- inla(
    y ~ 0 + f(i, model = Cmodel),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

all.equal(qc, prec(fit))

### now consider variances as well (PC-prior for this)
Vmodel <- cgeneric(
    model = "LKJ",
    n = n,
    eta = eta,
    sigma.prior.reference = rep(1, n),
    sigma.prior.probability = rep(0.5, n),
    useINLAprecomp = FALSE
)

sigmas <- n:1/2
diag(sigmas) %*% cc %*% diag(sigmas)

QiV <- prec(Vmodel, theta = c(log(sigmas), theta1))
(V <- chol2inv(chol(QiV)))

### simulate some data
nrep <- 200
xx <- matrix(rnorm(nrep * n), nrep) %*% as.matrix(chol(V))
str(xx)

dat2 <- data.frame(
    i = rep(1:n, each = nrep),
    r = rep(1:nrep, n),
    y = as.vector(xx)
)
str(dat2)

fitr <- inla(
    y ~ 0 + f(i, model = Vmodel, replicate = r),
    data = dat2,
    control.family = cfam,
    control.inla = cinla
)

cc
round(solve(prec(Cmodel, theta = tail(fitr$mode$theta, m))), 4)
round(solve(prec(Vmodel, theta = c(rep(0, n), tail(fitr$mode$theta,m)))), 4)

V
round(solve(prec(Vmodel, theta = fitr$mode$theta)), 4)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
