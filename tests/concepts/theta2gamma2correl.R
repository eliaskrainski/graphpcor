## this illustrates the theta2gamma2R parametrization
## using the hypershere decomposition,
## Rapisarda, Brigo and Mercurio (2007)

## 1. random C (correlation matrix), summary properties
## 2. C samples considering PC-prior with different lambda
## 3. how to permute C?

library(graphpcor)

nrepl <- 1000
table(replicate(nrepl, {
    th1 <- pi/(1 + exp(rnorm(10)))
    all.equal(tcrossprod(graphpcor:::theta2gamma2L(th1, fromR=FALSE)),
              INLA:::inla.pc.cormat.theta2R(th1))   
}))

table(replicate(nrepl, {
    th1 <- pi/(1 + exp(rnorm(10,0,10)))
    all.equal(tcrossprod(graphpcor:::theta2gamma2L(th1, fromR=FALSE)),
              INLA:::inla.pc.cormat.theta2R(th1))   
}))

theta2correl(c(0))
theta2correl(c(-1))
theta2correl(c(1))
theta2correl(c(-3))
theta2correl(c(3))

theta2correl(c(0,0,0))

theta2correl(c(0,1,1))
theta2correl(c(1,0,1))
theta2correl(c(1,1,0))
theta2correl(c(1,1,-1))

rcorrel(2)
rcorrel(3)
rcorrel(5)

## 1.  random correlation matrices

## function to extract correlation matrix properties
fCsummary <- function(m) {
    stopifnot(all.equal(diag(m), rep(1, nrow(m))))
    out <- c(det=attr(m, "determinant"),
             kld=attr(m, "kld"))
    m <- m[lower.tri(m)]
    c(out, mean = mean(m), rmin = min(m), rmax = max(m), r = m)
}

s5 <- t(replicate(nrepl, fCsummary(rcorrel(5))))
s25 <- t(replicate(nrepl, fCsummary(rcorrel(25))))

apply(s5, 2, summary)
apply(s25, 2, summary)

stem(s5[,1])
stem(s5[,2])

stem(s25[,1])
stem(s25[,2])

s5[1,]

bkc <- -10:10/10
par(mfrow = c(2,5), mar = c(2, 2, 0,0), bty = "n")
for(k in 5+1:10)
    hist(s5[, k], bkc, main = '')

par(mfrow = c(15, 20), mar = c(0,0,0,0), bty = "n")
for(k in 6:ncol(s25))
    hist(s25[, k], bkc, main = '', axes = F)

### 2. from the PC-prior, draw random correlation matrix
set.seed(1)
rcorrel(5, 0.001)
set.seed(1)
rcorrel(5, 1)
set.seed(1)
rcorrel(5, nrepl)

s5l.01 <- t(replicate(nrepl, fCsummary(rcorrel(5, .01))))
s5l1 <- t(replicate(nrepl, fCsummary(rcorrel(5, 1))))
s5l100 <- t(replicate(nrepl, fCsummary(rcorrel(5, 100))))

apply(s5l.01, 2, summary)
apply(s5l1, 2, summary)
apply(s5l100, 2, summary)

par(mfcol = c(3,10), mar = c(2, 2, 0,0), bty = "n")
for(k in 6:ncol(s5l1)) {
    hist(s5l.01[, k], bkc, main = '')
    hist(s5l1[, k], bkc, main = '')
    hist(s5l100[, k], bkc, main = '')
}


s25l.01 <- t(replicate(nrepl, fCsummary(rcorrel(25, .01))))
s25l1 <- t(replicate(nrepl, fCsummary(rcorrel(25, 1))))
s25l100 <- t(replicate(nrepl, fCsummary(rcorrel(25, 100))))

apply(s25l.01[, 1:5], 2, summary)
apply(s25l1[, 1:5], 2, summary)
apply(s25l100[, 1:5], 2, summary)

dim(s25l1)

par(mfcol = c(3,10), mar = c(2, 2, 0,0), bty = "n")
for(k in sample(6:ncol(s25l1), 10)) {
    hist(s25l.01[, k], bkc, main = '')
    hist(s25l1[, k], bkc, main = '')
    hist(s25l100[, k], bkc, main = '')
}

## 3. how to permute C?

summary(replicate(nrepl, sum(log(diag(chol(
                            INLA:::inla.pc.cormat.permute(
                                       rcorrel(3, 1))))))))
summary(replicate(nrepl, sum(log(diag(chol(
                            INLA:::inla.pc.cormat.permute( ### ooops
                                       rcorrel(5, 1))))))))
summary(replicate(nrepl, sum(log(diag(chol(
                            INLA:::inla.pc.cormat.permute( ### ooops
                                       rcorrel(10, 1))))))))


