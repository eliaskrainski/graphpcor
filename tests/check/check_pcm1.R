
## Laplace
fm1 <- function(r, l) {
    l * exp(-l*abs(r))*0.5
}

lseq <- c(0.3, 1, 5)
nl <- length(lseq)
b <- 4; a <- -b

plot(function(x) fm1(x, 5), a, b, n = 1001,
     type = 'n', bty = 'n')
for(i in 1:nl) {
    plot(function(x) fm1(x, lseq[i]), a, b, n=1001,
         add = TRUE, col = i, lty = i)
    print(integrate(function(x) fm1(x, lseq[i]), a*5, b*5))
}

## Ex1: x ~ N(0, C(r))
##  C(r) = [1  r]
##         [r  1]
fcr <- function(r) {
    matrix(c(1, r, r, 1), 2)
}

fcr(0.9)

## internal parametrization
fcx <- function(x) {
  r <- tanh(x)
  matrix(c(1, r, r, 1), 2)
}

fcx(atanh(0.9))

## I(r)
Ir2f <- function(r) {
    aux <- 1-r^2
    return((1/aux) + 2*r^2/(aux^2))
}

## I(x)
Ix2f <- function(x) {
    return(tanh(x)^2 + 1)
}

xseq <- seq(-3, 3, 0.1)
rseq <- seq(-0.9, 0.9, 0.01)

library(graphpcor)

Ir.n <- sapply(rseq, function(r0) {
    C0 <- fcr(r0)
    hessian(function(x) graphpcor:::KLD10(fcr(x), C0),
            x = r0)
})

Ix.n <- sapply(xseq, function(x0) {
    C0 <- fcx(x0)
    hessian(function(x) graphpcor:::KLD10(fcx(x), C0),
            x = x0)
})

par(mfrow = c(1, 2), mar = c(4,4,1,1), mgp = c(3,1,0), las = 0, bty = 'n')
plot(Ir2f, -.99, .99, n=2001, log = 'y',
     xlab = expression(rho), ylab = 'I(r_0)')
points(rseq, Ir.n, pch = 19)
legend("top", c("Numeric", "Analytic"), bty = 'n',
       pch = 19, pt.cex = 1:0, lty = 0:1)
plot(Ix2f, -3, 3, n=2001, log = 'y',
     xlab = expression(xi), ylab = 'I(xi_0)')
points(xseq, Ix.n, pch = 19)

## density
drho <- function(rho, rho0, l) {
  x0 <- atanh(rho0)
  I0 <- 1-rho0^2
  x <- atanh(rho)
  sI0 <- sqrt(I0)
  xi <- sI0 * (x - x0)
  ld <- log(l*sI0/2) -l*abs(xi) -log(1-rho^2)
  exp(ld) 
}

rho0s <- c(-0.5, 0, 0.7)

lseq <- c(0.3, 1, 5)
nl <- length(lseq)
lcol <- c(1, 2, 4)
b <- 3; a <- -b
llab <- lapply(lseq, function(l) 
    bquote(lambda[rho] == .(l)))

for(k in 1:3) {
    for(i in 1:nl) {
        print(integrate(function(x) drho(x, rho0s[k], lseq[i]), -1, 1))
    }
}

par(mfrow = c(1, 3), mar = c(4,4,1,1), mgp = c(3,1,0), las = 1)
for(k in 1:3) {
  plot(function(x) drho(x,rho0s[k],lseq[nl]), 
       -1, 1, n=1001, type = "n", bty = "n", 
       xlab = expression(rho), 
       ylab = expression(pi(rho|lambda)))
  for(i in 1:nl) {
    plot(function(x) drho(x, rho0s[k], lseq[i]), 
         -1, 1, n=1001, add = TRUE, col = lcol[i], lty = i, lwd = 2)
  }
}
legend("top", bty = "n", 
       as.expression(llab), 
       ##title = expression(lambda), 
       col = lcol, lty = 1:nl, lwd = 2) 


