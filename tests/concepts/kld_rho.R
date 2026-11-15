
## x ~ N(0, C)
##     |1  r  r|
## C = |r  1  r|
##     |r  r  1|
## r \in (-1, 1)

rfn <- function(r0, r1) {
    c0 <- toeplitz(c(1,r0,r0))
    c1 <- toeplitz(c(1,r1,r1))
    graphpcor:::KLD10(c1, c0)
}

rseq1 <- seq(-0.05,0.9,0.1)
kseq1 <- sapply(rseq1, function(r) rfn(r,r+.1))

range(rseq2 <- seq(0.99,0.997,0.001))
kseq2 <- sapply(rseq2, function(r) rfn(r,r+.001))

par(mfrow = c(1,2), mar = c(4,4,1,1), mgp = c(2,1,0), bty = "n")
plot(rseq1, kseq1, xlab = expression(rho), ylab = "KLD", pch = 19)
plot(rseq2, kseq2, xlab = expression(rho), ylab = "KLD", pch = 19)

## NOTICE the effect on samples, with AR(1) case
ar1fn <- function(n, r) {
    x <- numeric(n)
    x[1] <- rnorm(1, 0, 1/sqrt(1-r^2))
    for(i in 2:n)
        x[i] <- r * x[i-1] + rnorm(1, 0, 1)
    return(x)
}

rhos <- c(5:9/10, .945, 0.99, 0.9945, 0.999)
(drhos <- diff(rhos)) ## the last ones are ~2x closer

xx <- sapply(rhos, function(r) {
    set.seed(1) ## SAME SEED
    ##        arima.sim(model = list(ar = r), n = 1000, n.start = 100) ## has bad choices for starting
   ar1fn(300, r)
})

1/(1-rhos^2) ## marginal var
apply(xx, 2, var)

### now, look:
### sd(xx[, k] - xx[, k-1]): sd between realizations of close rho
sdnr <- apply(xx[, 2:ncol(xx)] - xx[, 1:(ncol(xx)-1)], 2, sd)
### corr(xx[, k], xx[, k-1]): cor between realizations of close rho
cxnr <- sapply(2:ncol(xx), function(k) cor(xx[, k], xx[, k-1]))
data.frame(diff_rho = drhos,
           sd_rho_near = sdnr,
           cx_rho_near = cxnr,
           row.names = paste(rhos[-length(rhos)], "to", rhos[-1])) 

par(mfrow = c(3,3), mar = c(4,4,1,1), mgp = c(2,1,0), bty = "n")
for(k in 1:ncol(xx))
    plot.ts(xx[,k], ylab = expression(x[t]),
            main = as.expression(bquote(rho==.(rhos[k]))))
