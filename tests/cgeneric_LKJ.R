library(INLA)
library(graphpcor)

n <- 4
(m <- n*(n-1)/2)
eta <- 10

cmodel <- cgeneric(
    model = "LKJ",
    n = n,
    eta = eta
)

str(cmodel)

graph(cmodel, optimize = TRUE)

graph(cmodel)

initial(cmodel)

theta1 <- rnorm(m)

(qq <- prec(cmodel, theta = theta1))

(vv <- solve(qq))

all.equal(as.matrix(vv),
          basecor(theta1, n)$base)

## fake data
dat1 <- data.frame(
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fit <- inla(
    y ~ 0 + f(i, model = cmodel),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

all.equal(qq, prec(fit))

fit$summary.hy

par(mfrow = c(1,1), mar = c(4,4,1,1), mgp = c(2,0.5,0))
plot(fit$internal.marginals.hyperpar[[1]],
     type = 'n', bty = 'n',
     xlab = '', ylab = '', main = '')
for(i in 1:6) {
    lines(fit$internal.marginals.hyperpar[[i]], col = i)
    abline(v=theta1[i], col = i, lwd = 2)
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
    y ~ 0 + f(i, model = cmodel, replicate = r),
    data = dat2,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

(Cfitted <- basecor(fitr$mode$theta, n))
round(vv, 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
