library(INLA)

library(graphpcor)

dlkj <- function(R, eta, verbose = FALSE) {
    lR <- chol(R)
    ldR <- 2*sum(log(diag(lR)))
    if(verbose) print(ldR)
    d <- ncol(R)
    k <- 1:(d-1)
    lbk <- lbeta(eta + (d-k-1)/2,
                 eta + (d-k-1)/2) * (d-k)
    if(verbose) print(lbk)
    p2 <- sum((2*eta -2 + d - k)*(d-k))
    if(verbose) print(p2)
    if(verbose) print(sum(lbk) + p2*log(2))
    o <- sum(lbk) + p2*log(2) + (eta-1)*ldR
    return(o)
}

################################################################
### n = 2
n <- 2
(m <- n*(n-1)/2)

(theta1 <- rnorm(m))

lR <- graphpcor:::theta2gamma2L(theta1)
log(attr(lR, 'determinant'))
R <- tcrossprod(lR)

dlkj(R,1)
dlkj(R,10)

eta <- 5.01

cmodel <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta,
    debug = 1e9)

str(cmodel)

p2 <- prior(cmodel, theta = theta1)
p2

log(attr(lR, 'determinant'))

dlkj(R, eta, TRUE)

j1 <- sum(log(pi*exp(-theta1)/((1+exp(-theta1))^2)))
j2fn <- function(x) {
    m <- length(x)
    r <- 0
    if(m>1) {
        r <- (m-1) * log(x[1])
        for(i in 1:(m-1)) {
            r <- r + (m-i-1) * log(sin(x[i+1]))
        }
    }
    return(r)
}
j2 <- j2fn(graphpcor:::x2rphi(pi/(1+exp(-theta1))))
c(j1 = j1, j2 = j2)

all.equal(p2, dlkj(R, eta) + j1 + j2)

################################################################
## n = 3
n <- 3
(m <- n*(n-1)/2)

(theta1 <- rnorm(m))

lR <- graphpcor:::theta2gamma2L(theta1)
log(attr(lR, 'determinant'))
R <- tcrossprod(lR)

gamma1 <- pi/(1+exp(-theta1))
gamma1

co <- cbind(cos(gamma1[1:2]),
            c(0, cos(gamma1[3])))
co

si <- cbind(sin(gamma1[1:2]),
            c(0, sin(gamma1[3])))
si
psi <- cbind(si[,1], si[,1]*si[,2])
psi

cbind(c(1, co[,1]),
      c(0, si[1,1], si[2,1]*co[2,2]),
      c(0,0,psi[2,2]))

lR
log(attr(lR, 'determinant'))


dlkj(R,1)
dlkj(R,10)

eta <- 5.01

cmodel <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta,
    debug = 1e9)

str(cmodel)

prec(cmodel)
solve(prec(cmodel))
solve(prec(cmodel, theta = theta1))

prec(cmodel, theta = theta1)
(Q1 <- chol2inv(t(graphpcor:::theta2gamma2L(theta1))))

itest <- inla(y ~ 0 + f(i, model = cmodel),
     data = list(y = rep(NA, n), i = 1:n),
     control.family = list(
         hyper = list(prec = list(initial = 10, fixed = TRUE))
     ),
     control.mode = list(theta = theta1, fixed = TRUE),
     control.compute = list(config = TRUE),
     verbose = !TRUE)

itest$mode$theta
solve(prec(itest))
solve(Q1)

p3 <- prior(cmodel, theta = theta1)
p3

gamma1
pi*exp(-theta1)/((1+exp(-theta1))^2)

j1 <- sum(log(pi*exp(-theta1)/((1+exp(-theta1))^2)))
rphi <- graphpcor:::x2rphi(pi/(1+exp(-theta1)))
c(rphi=rphi)
j2 <- j2fn(rphi)
c(j1 = j1, j2 = j2)

dlkj(R, eta) + c(no=0, ok=j1 + j2)
p3

all.equal(p3, dlkj(R, eta) + j1 + j2)

################################################################
## n = 4
n <- 4
(m <- n*(n-1)/2)

(theta1 <- rnorm(m))

lR <- graphpcor:::theta2gamma2L(theta1)
log(attr(lR, 'determinant'))
R <- tcrossprod(lR)

gamma1 <- pi/(1+exp(-theta1))

co <- cbind(cos(gamma1[1:3]),
            c(0, cos(gamma1[4:5])),
            c(0,0,cos(gamma1[6])))
co

si <- cbind(sin(gamma1[1:3]),
            c(0, sin(gamma1[4:5])),
            c(0,0,sin(gamma1[6])))
si
psi <- cbind(si[,1], si[,1]*si[,2], si[,1]*si[,2]*si[,3])
psi

cbind(c(1, co[,1]),
      c(0, si[1,1], si[2:3,1]*co[2:3,2]),
      c(0,0,psi[2,2], psi[3,2]*co[3,3]),
      c(0,0,0,psi[3,3]))

lR
log(attr(lR, 'determinant'))

dlkj(R,1)
dlkj(R,10)

eta <- 5.01

cmodel <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta,
    debug = 1e9)

str(cmodel)

p4 <- prior(cmodel, theta = theta1)
p4

gamma1
pi*exp(-theta1)/((1+exp(-theta1))^2)

j1 <- sum(log(pi*exp(-theta1)/((1+exp(-theta1))^2)))
rphi <- graphpcor:::x2rphi(pi/(1+exp(-theta1)))
c(rphi=rphi)
j2 <- j2fn(rphi)
c(j1 = j1, j2 = j2)

dlkj(R, eta) + c(no=0, ok=j1 + j2)
p4

all.equal(p4, dlkj(R, eta) + j1 + j2)


#######################################################################
### n = 5
n <- 5
(m <- n*(n-1)/2)

(theta1 <- rnorm(m))

lR <- graphpcor:::theta2gamma2L(theta1)
log(attr(lR, 'determinant'))
R <- tcrossprod(lR)

eta <- 5.01

cmodel <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta)

p5 <- prior(cmodel, theta = theta1)
p5

gamma1
pi*exp(-theta1)/((1+exp(-theta1))^2)

j1 <- sum(log(pi*exp(-theta1)/((1+exp(-theta1))^2)))
rphi <- graphpcor:::x2rphi(pi/(1+exp(-theta1)))
c(rphi=rphi)
j2 <- j2fn(rphi)
c(j1 = j1, j2 = j2)

dlkj(R, eta) + c(no=0, ok=j1 + j2)
p5

all.equal(p5, dlkj(R, eta) + j1 + j2)


## other tests
graph(cmodel, optimize = TRUE)

graph(cmodel)

round(ith <- initial(cmodel), 4)

theta1 <- rnorm(m)
x1 <- pi*plogis(theta1)

x1
cos(x1)
sin(x1)

(qq <- prec(cmodel, theta = theta1))

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
    y ~ 0 + f(i, model = cmodel),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

all.equal(qq, prec(fit))

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
    y ~ 0 + f(i, model = cmodel, replicate = r),
    data = dat2,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
    ##, verbose = TRUE
)

round(qq, 2)
round(prec(fitr), 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
