
## 1. sampling on the sphere
rsphere <- function(n, m) {
    z <- matrix(rnorm(n*m), m)
    r <- sqrt(colSums(z^2))
    t(z)/r
}

par(mfrow = c(2, 2), mar = c(2,2,0,0),
    mgp = c(1,0.5,0), bty = 'n')
plot(rsphere(1000, 2), xlab = '', ylab = '')
plot(rsphere(1000, 3)[,1:2], xlab = '', ylab = '')
plot(rsphere(1000, 3)[,2:3], xlab = '', ylab = '')
plot(rsphere(1000, 3)[,c(1,3)], xlab = '', ylab = '')

## 2. sampling within the ball, with r ~ Exp(\lambda)
rmpcp <- function(n, m, lambda = 1, Hneg.5, theta0) {
  out <- rsphere(n, m) * rexp(n, lambda)
  if(!missing(Hneg.5)) {
    out <- out %*% Hneg.5
  }
  if(!missing(theta0)) {
    out <- sweep(out, 2, -theta0)
  }
  return(out)
}

par(mfrow = c(2, 2), mar = c(2,2,0,0),
    mgp = c(1,0.5,0), bty = 'n')
plot(rmpcp(1000, 2, 1), xlab = '', ylab = '')
plot(rmpcp(1000, 3, 1)[,1:2], xlab = '', ylab = '')
plot(rmpcp(1000, 3, 1)[,2:3], xlab = '', ylab = '')
plot(rmpcp(1000, 3, 1)[,c(1,3)], xlab = '', ylab = '')

th2 <- c(-1,1)
H2 <- matrix(c(2,-1,-1,1), 2)/10
th3 <- (-1:1)
H3 <- matrix(c(3,-1,1, -1,2,0, 1,0,1), 3)/10
par(mfrow = c(2, 2), mar = c(2,2,0,0),
    mgp = c(1,0.5,0), bty = 'n')
plot(rmpcp(1000, 2, 1, H2, th2), xlab = '', ylab = '')
plot(rmpcp(1000, 3, 1, H3, th3)[,1:2], xlab = '', ylab = '')
plot(rmpcp(1000, 3, 1, H3, th3)[,2:3], xlab = '', ylab = '')
plot(rmpcp(1000, 3, 1, H3, th3)[,c(1,3)], xlab = '', ylab = '')

library(graphpcor)

c0 <- matrix(
    c( 1.0,  0.8, -0.5,
       0.8,  1.0, -0.4,
      -0.5, -0.4,  1.0), 3)
c0

bc0 <- basecor(c0)
I0 <- hessian(bc0)
I0d <- graphpcor:::dspd(I0)

str(I0d)

th1 <- rmpcp(10000, m = 3, lambda = 1,
             Hneg.5 = I0d$sqrtInv,
             theta0 = bc0$theta)
th5 <- rmpcp(10000, m = 3, lambda = 5,
             Hneg.5 = I0d$sqrtInv,
             theta0 = bc0$theta)

il3 <- c(2,3,6)
cs1 <- t(sapply(1:nrow(th1), function(i)
    tcrossprod(cholcor(th1[i,]))[il3]))
cs5 <- t(sapply(1:nrow(th5), function(i)
    tcrossprod(cholcor(th5[i,]))[il3]))

par(mfrow = c(2,3), mar = c(4,4,0.5,0.5), mgp = c(2,0.5,0),
    las = 1, bty = 'n')
hist(cs1[, 1], xlab = expression(rho[2~","~1]), ylab = '', main = '')
abline(v = c0[il3[1]], col = 2, lwd = 2, lty = 2)
legend("topleft", expression(lambda == 1))
hist(cs1[, 2], xlab = expression(rho[3~","~1]), ylab = '', main = '')
abline(v = c0[il3[2]], col = 2, lwd = 2, lty = 2)
hist(cs1[, 3], xlab = expression(rho[3~","~2]), ylab = '', main = '')
abline(v = c0[il3[3]], col = 2, lwd = 2, lty = 2)
hist(cs5[, 1], xlab = expression(rho[2~","~1]), ylab = '', main = '')
abline(v = c0[il3[1]], col = 2, lwd = 2, lty = 2)
legend("topleft", expression(lambda == 5))
hist(cs5[, 2], xlab = expression(rho[3~","~1]), ylab = '', main = '')
abline(v = c0[il3[2]], col = 2, lwd = 2, lty = 2)
hist(cs5[, 3], xlab = expression(rho[3~","~2]), ylab = '', main = '')
abline(v = c0[il3[3]], col = 2, lwd = 2, lty = 2)

## 3. set up two different ordering for the 'same' graph

p <- 5
gs <- graphpcor(paste0("X1~", paste0("X",2:p)))
gs

Ls <- Laplacian(gs)

oj <- c(2:p,1)
ogs <- graphpcor(Ls[oj, oj])
ogs

par(mfrow=c(1,2))
plot(gs)
plot(ogs)

## 4. define the SAME base model/correlation matrix from each one
thp <- c(-1, 1.0, 1.0, -0.5)

bp <- basepcor(thp, p, iLtheta = gs)
bp

round(100*bp$base[oj, oj])

obp <- basepcor(bp$base[oj, oj], iLtheta = ogs)
obp

all.equal(bp$base[oj, oj],
          obp$base)

## Hessian
Hs <- hessian(bp)
Hso <- hessian(obp)
sHs <- graphpcor:::dspd(Hs)
sHso <- graphpcor:::dspd(Hso)

## 5. Simulate from both and visualize
ns <- 5000
xxs <- rmpcp(ns, p-1, H = sHs$sqrt, theta0 = bp$theta)
xxso <- rmpcp(ns, p-1, H = sHso$sqrt, theta0 = obp$theta)

iil <- which(lower.tri(diag(p)))

str(xxs)

par(mfcol = c(2, p-1), mar = c(1.5,1.5,0,0), 
    mgp = c(1,0.5,0), bty = 'n')
for(k in 1:(p-1)) {
    hist(xxs[, k], 50, xlab = '', ylab = '', main = '', border = 'transparent')
    hist(xxso[, k], 50, xlab = '', ylab = '', main = '', border = 'transparent')
}

ccs <- t(sapply(
    1:nrow(xxs), function(i)
        basepcor(xxs[i, ], p, iLtheta = gs)$base[oj, oj][iil]
))

ccso <- t(sapply(
    1:nrow(xxso), function(i)
        basepcor(xxso[i, ], p, iLtheta = ogs)$base[iil]
))

str(ccs)

p*(p-1)/2
cfill <- rgb(1:0, 0.3, 0:1, 0.5)
par(mfrow = c(3, 4), mar = c(1.5,1.5,0,0), 
    mgp = c(1,0.5,0), bty = 'n')
for(k in 1:(p*(p-1)/2)) {
    d1 <- density(ccs[, k])
    d2 <- density(ccso[, k])
    h1 <- hist(ccs[, k], 50, plot = FALSE)
    h2 <- hist(ccso[, k], 50, plot = FALSE)
    plot(h1, xlab = '', ylab = '', main = '', freq = FALSE,
         xlim = range(h1$breaks, h2$breaks),
         ylim = range(h1$dens, h2$dens),
         col = cfill[1], border = 'transparent')
    plot(h2, add = TRUE, freq = FALSE,
         col = cfill[2], border = 'transparent')
    lines(d1, col = 2, lwd = 1)
    lines(d2, col = 4, lwd = 1)
    abline(v = obp$base[iil[k]], lty = 2)
}
legend("topright", c("initial", "ordered"), bty = "n",
       fill = cfill, border = 'transparent', cex = 1.5)
