library(INLA)
library(graphpcor)

## a cgeneric model for the LKJ prior
## for theta from the CPC parametrization
cglkj <- cgeneric(
    model = "LKJ", n = 3, eta = 2,
    useINLAprecomp = FALSE)

## correlation matrix, p = 3
cc <- matrix(c(1,.8,-.625, 0.8,1,-.5, -0.625,-.5,1), 3)

## CPC parametrization: C(theta)
(bb <- basecor(cc))
th <- bb$theta

## p(theta | eta)
prior(cglkj, theta = th)

## precision
(qq <- prec(cglkj, theta = th))
solve(qq)

## setup a grid around theta0
h <- 0.1
th0 <- seq(-3+h/2, 3-h/2, h)
n0 <- length(th0)
thg <- t(expand.grid(th0 + th[1], th0 + th[2], th0 + th[3]))
dim(thg)

ptheta <- array(exp(prior(cglkj, theta = thg)), c(n0, n0, n0))
sum(ptheta*(h^3))

par(mfrow = c(2, 2), mar = c(3, 3, 0.1, 0.1), mgp = c(2, 0.5, 0))
for(i in 1:3) {
    i2 <- setdiff(1:3, i)
    pp <- apply(ptheta * h, i2, sum)
    image(th0 + th[i2[1]], th0 + th[i2[2]], pp,
          xlab = as.expression(bquote(theta[.(i2[1])])),
          ylab = as.expression(bquote(theta[.(i2[2])])))
    contour(th0 + th[i2[1]], th0 + th[i2[2]], pp, add = TRUE)
}
ppp <- lapply(1:3, function(i) apply(ptheta * (h^2), i, sum))
plot(th0, xlim = c(min(th)-3, max(th)+3), ylim = range(unlist(ppp)))
for(i in 1:3) {
    lines(th0 + th[i], ppp[[i]], col = i, lwd = 3, lty = i)
}

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

all.equal(as.matrix(qc),
          as.matrix(prec(fit)))

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
