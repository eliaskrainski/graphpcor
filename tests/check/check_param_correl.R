## this illustrates the correlation parametrization

library(graphpcor)

args(theta2correl)

## compare the "SAP" in corgraph with the one in INLA: pi/(1+exp(-theta))
nrepl <- 5000
table(replicate(nrepl, {
    th1 <- rnorm(10)
    all.equal(theta2correl(th1, parametrization = 'sap'),
              INLA:::inla.pc.cormat.theta2R(pi/(1+exp(-th1))))   
}))

c0 <- matrix(c(1,0.7,0.1, 0.7,1,-0.3, 0.1,-0.3,1), 3)
c0

lc <- t(chol(c0))
lc

lq <- t(chol(chol2inv(t(lc))))
lq

findtheta <- function(cc, parametrization = 'cpc') {
    p <- nrow(cc)
    L <- t(chol(cc))
    il <- which(lower.tri(cc))
    m <- length(il)
    l <- L[il]
    if(parametrization=="itp") {
        lq <- t(chol(chol2inv(t(L))))
        r <- lq[il]
        attr(r, 'd0') <- diag(lq)
        return(r)
    }
    return(
        optim(rep(0, m), 
              function(x) mean((theta2L(x, p, parametrization)[il]-l)^2),
              method = 'BFGS')$par)
}

th0a <- findtheta(c0, "itp")
th0b <- findtheta(c0, "cpc")
th0c <- findtheta(c0, "sap")

c0
theta2correl(th0a, 3, 'itp')
theta2correl(th0b, 3, 'cpc')
theta2correl(th0c, 3, 'sap')

rbind(th0a, th0b, lq[lower.tri(lq)])

## 1.  random correlation matrices

## function to extract correlation matrix properties
fCsummary <- function(m) {
    stopifnot(all.equal(diag(m), rep(1, nrow(m))))
    out <- c(det=attr(m, "determinant"),
             kld=attr(m, "kld"))
    m <- m[lower.tri(m)]
    c(out, mean = mean(m), rmin = min(m), rmax = max(m), r = m)
}

system.time(s5 <- t(replicate(nrepl, fCsummary(rcorrel(5)))))
system.time(s25 <- t(replicate(nrepl, fCsummary(rcorrel(25)))))

apply(s5, 2, summary)
## apply(s25, 2, summary)

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


