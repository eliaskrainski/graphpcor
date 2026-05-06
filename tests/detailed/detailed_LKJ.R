library(INLA)
library(graphpcor)

n <- 4
(m <- n*(n-1)/2)

eta <- 3.0

cmodel <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta,
    debug = 1e9*0,
    useINLAprecomp = FALSE
)

cgeneric_graph(cmodel, optimize = TRUE)

cgeneric_graph(cmodel)

cgeneric_initial(cmodel)

theta1 <- rnorm(m)

(qq <- cgeneric_Q(cmodel, theta = theta1))

(vv <- chol2inv(chol(qq)))
(ll <- t(chol(vv)))

all.equal(graphpcor:::c4theta(theta1),
          t(as.matrix(ll)))

all.equal(as.matrix(vv),
          crossprod(graphpcor:::c4theta(theta1)))

lconst <- function(n, e) {
    e1 <- e - 1.0
    lres <- log(2)*n*(n-1) * (e1 + (n+n-1)/6)
    j <- 1:(n-1)
    a <- 2.0 * e1 + j+1
    lres <- lres + sum(2*j*lgamma(0.5*a)-j*lgamma(a))
    return(lres)
}

thU <- function(x) graphpcor:::c4theta(x)

thC <- function(th)
    crossprod(graphpcor:::c4theta(th))

thD <- function(th) {
    d1 <- -2*sum(log(cosh(th))) 
    x <- thC(th)
    p <- ncol(z)
    lrk <- matrix(0, p, p)
    for(i in 1:(p-2)) { ## p - (p-1) -1  == 0
    ##    cat(i, '')
        for(j in (i+1):p) {
            lrk[i, j] <- (p-i-1)*log(1 - x[i,j]^2)
  ##          cat(j, '')
        }
##        cat('\n')
    }
    lrk
}

thetaJacobian <- function(theta) {
    aux <- thD(theta) 
    return(-2*sum(log(cosh(theta))) +
           0.5 * sum(aux[upper.tri(aux)]))
}

th1 <- rnorm(n*(n-1)/2)
t(l1 <- graphpcor:::c4theta(th1))
c1 <- crossprod(l1)
c1

cgeneric_prior(cmodel, theta = th1)
c(sum(log(diag(l1))),
  (eta-1)*2*sum(log(diag(l1))),
  lconst(n, eta), 
  thetaJacobian(th1))
dLKJ(c1, eta, TRUE) + thetaJacobian(th1)

ths <- matrix(rnorm(1000 * m), m)

ds <- sapply(1:1000, function(j)
    dLKJ(crossprod(graphpcor:::c4theta(ths[, j])),
         eta = eta, log = TRUE))

lj <- cgeneric_prior(cmodel, theta = ths) - ds

ljn <- sapply(1:1000, function(j) 
    log(det(
        jacobian(function(x)
            crossprod(graphpcor:::c4theta(x))[lower.tri(diag(4))],
            x = ths[, j])))
    )

ljr <- sapply(1:1000, function(j) thetaJacobian(ths[, j]))

pairs(data.frame(cg=lj, R=ljr, num=ljn), upper.panel = NULL)

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

all.equal(qq, cgeneric_Q(fit))

fit$summary.hy

par(mfrow = c(1,1), mar = c(4,4,1,1), mgp = c(2,0.5,0))
plot(fit$internal.marginals.hyperpar[[1]],
     type = 'n', bty = 'n',
     xlab = '', ylab = '', main = '')
for(i in 1:m) {
    lines(fit$internal.marginals.hyperpar[[i]], col = i)
    abline(v=theta1[i], col = i, lwd = 2)
}

### fit some data
nrep <- 2
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

if(nrep<4)
    round(solve(cgeneric_Q(fitr)), 2)

(Lfitted <- graphpcor:::c4theta(fitr$mode$theta))
round(crossprod(Lfitted), 2)
round(vv, 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
