library(graphpcor)

## random drawn 
m <- 6
x <- rnorm(m)

to.x <- graphpcor:::rphi2x
to.rphi <- graphpcor:::x2rphi

table(replicate(1000, {x <- rnorm(10); all.equal(x, to.x(to.rphi(x)))}))

print(mean(abs(to.x(to.rphi(x)) - x)))

jacobian(to.rphi, x)
xx <- to.rphi(x)

## this is a test that forward-Jacobian is 1/backward-Jacobian
print(abs(det(jacobian(to.rphi, x))))
print(1/abs(det(jacobian(to.x, xx))))

## this is using formula for the 'a closed-form expression for the volume element in spherical
## coordinates' in https://en.wikipedia.org/wiki/N-sphere

r <- xx[1]
phi <- xx[-1]

ldJfn <- function(rphi) {
    m <- length(rphi)
    out <- (m-1)*log(rphi[1])
    if(m>2)
        out <- out + sum( (m-2):1 * log(sin(rphi[2:(m-1)])) )
    return(out)
}

ldJfn(xx)

print(ld.a <- (m-1) * log(r) +
          (m-2) * log(sin(phi[1])) +
          (m-3) * log(sin(phi[2])) +
          (m-4) * log(sin(phi[3])) +
          (m-5) * log(sin(phi[4])))
print(ld.num <- determinant(jacobian(to.x, xx))$modulus)

qcmodel <- cgeneric(
    model = "pc_prec_correl",
    n = 4,
    theta.base = rep(-1, 6),
    lambda = 1,
    debug = 1e9)

prior(qcmodel, theta = x)

Qx <- prec(qcmodel, theta = x)
n <- ncol(Qx)
Qx

solve(Qx)

hprec <- list(initial = 10, fixed = TRUE)

library(INLA)

fit0 <- inla(
    y ~ 0 + f(i, model = qcmodel),
    data = data.frame(i = 1:n, y = NA),
    control.family = list(hyper = list(prec = hprec)),
    control.mode = list(theta = x, fixed = TRUE),
    verbose = TRUE
)

(x - fit0$mode$theta)

q.fit0 <- prec(fit0)

all.equal(Qx, q.fit0)

## some data
cc <- graphpcor:::theta2C(x)
cc

nrep <- 100
xx <- matrix(rnorm(nrep * n), ncol = n) %*% chol(cc)
cov(xx)
cc

dataf <- data.frame(
    i = rep(1:n, each = nrep),
    r = rep(1:nrep, n),
    y = as.vector(xx)
)

str(dataf)

fit <- inla(
    formula = y ~ 0 + f(i, model = qcmodel, replicate = r),
    data = dataf,
    control.family = list(hyper = list(prec = hprec)), 
    verbose = TRUE
)

x-fit$mode$theta

Q1 <- prec(qcmodel, theta = fit$mode$theta)
Q1

cc
solve(Q1)
