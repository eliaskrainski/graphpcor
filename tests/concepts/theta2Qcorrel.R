## this illustrates the theta2Qcorrel parametrization
## step 1: |     1                                    |
##         | \theta[1]      1         SYMMETRIC       |
##    Q0 = | \theta[2]   \theta[n]                    |
##         |    :           ...        ...            |
##         | \theta[n-1]    ...      \theta[m]    1   |
## step 2: V0 = Q0^{-1}
## step 3: C = diag(V0)^{-1/2} %*% V0 %*% diag(V0)^{-1/2}
## step 4: Q = C^{-1}

## 1. elements of Q for a random correlation matrix: summary properties

library(graphpcor)

lQCrandom <- function(p, lambda = NA, test = FALSE) {
    ## this do the opposite:
    ## specify C(\theta[1:m]) and then Q = C^{-1}
    ## return Cholesky (lower triangle) of Q
    m <- p * (p-1)/2
    if(is.na(lambda)) {
        theta <- pi/(1 + exp(-rnorm(m, 0, 1)))
    } else {
        stopifnot(lambda>0)
        theta <- INLA:::inla.pc.cormat.rtheta(n=1, p, lambda)
    }
    if(FALSE) {
        Q <- chol2inv(chol(INLA:::inla.pc.cormat.theta2R(theta)))
    } else {
        lC <- graphpcor:::theta2gamma2L(theta, fromR = FALSE)
        Q <- chol2inv(t(lC))
    }
    lQ <- t(chol(Q))
    if(test) {
        cc <- chol2inv(t(lQ))
        stopifnot(all.equal(diag(cc), rep(1,p)))
    }
    lQ    
}

lQCrandom(3)
lQCrandom(3)

lQCrandom(5) ## chol(rcorrel(5))
lQCrandom(5,10)

chol2inv(t(lQCrandom(3,0.5)))
chol2inv(t(lQCrandom(3,10)))

nrepl <- 1000

replicate(nrepl, lQCrandom(5,NA,TRUE)) -> p5l0
replicate(nrepl, lQCrandom(5,1,!TRUE)) -> p5l1
replicate(nrepl, lQCrandom(5,10,TRUE)) -> p5l10

replicate(nrepl, lQCrandom(25,NA,!TRUE)) -> p25l0
replicate(nrepl, lQCrandom(25,1,!TRUE)) -> p25l1
replicate(nrepl, lQCrandom(25,10,TRUE)) -> p25l10

par(mfrow = c(2, 3), mar = c(2,2,0,0), bty = 'n')
plot(1, xlim = c(1, 5), ylim = range(apply(p5l0, 3, diag)), log = 'y', type = 'n')
for(i in 1:nrepl)
    lines(diag(p5l0[,,i]))
plot(1, xlim = c(1, 5), ylim = range(apply(p5l1, 3, diag)), log = 'y', type = 'n')
for(i in 1:nrepl)
    lines(diag(p5l1[,,i]))
plot(1, xlim = c(1, 5), ylim = range(apply(p5l10, 3, diag)), log = 'y', type = 'n')
for(i in 1:nrepl)
    lines(diag(p5l10[,,i]))
plot(1, xlim = c(1, 25), ylim = range(apply(p25l0, 3, diag)), log = 'y', type = 'n')
for(i in 1:nrepl)
    lines(diag(p25l0[,,i]))
plot(1, xlim = c(1, 25), ylim = range(apply(p25l1, 3, diag)), log = 'y', type = 'n')
for(i in 1:nrepl)
    lines(diag(p25l1[,,i]))
plot(1, xlim = c(1, 25), ylim = range(apply(p25l10, 3, diag)), log = 'y', type = 'n')
for(i in 1:nrepl)
    lines(diag(p25l10[,,i]))


par(mfrow = c(2, 3), mar = c(2,2,0,0), bty = 'n')
plot(1, xlim = c(1, 10), ylim = range(apply(p5l0, 3, function(x) x[lower.tri(x)])), type = 'n')
for(i in 1:nrepl)
    lines(p5l0[,,i][lower.tri(diag(5))])
plot(1, xlim = c(1, 10), ylim = range(apply(p5l1, 3, function(x) x[lower.tri(x)])), type = 'n')
for(i in 1:nrepl)
    lines(p5l1[,,i][lower.tri(diag(5))])
plot(1, xlim = c(1, 10), ylim = range(apply(p5l10, 3, function(x) x[lower.tri(x)])), type = 'n')
for(i in 1:nrepl)
    lines(p5l10[,,i][lower.tri(diag(5))])
jj <- sample(which(lower.tri(diag(25))), 25)
plot(1, xlim = c(1, 20), ylim = range(apply(p25l0, 3, function(x) x[jj])), type = 'n')
for(i in 1:nrepl)
    lines(p25l0[,,i][jj])
plot(1, xlim = c(1, 20), ylim = range(apply(p25l1, 3, function(x) x[jj])), type = 'n')
for(i in 1:nrepl)
    lines(p25l1[,,i][jj])
plot(1, xlim = c(1, 20), ylim = range(apply(p25l10, 3, function(x) x[jj])), type = 'n')
for(i in 1:nrepl)
    lines(p25l10[,,i][jj])

par(mfrow = c(2, 3))
hist(apply(p5l0, 3, function(x) x[lower.tri(x)]), 100, main = '')
hist(apply(p5l1, 3, function(x) x[lower.tri(x)]), 100, main = '')
hist(apply(p5l10, 3, function(x) x[lower.tri(x)]), 100, main = '')
hist(apply(p25l0, 3, function(x) x[lower.tri(x)]), 100, main = '')
hist(apply(p25l1, 3, function(x) x[lower.tri(x)]), 100, main = '')
hist(apply(p25l10, 3, function(x) x[lower.tri(x)]), 100, main = '')
