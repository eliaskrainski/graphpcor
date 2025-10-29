library(graphpcor)

g1 <- treepcor(p1 ~ c1 + c2 + c3)
d1 <- dim(g1)

plot(g1)

(theta.p <- -0.5)

mcov0 <- vcov(g1, theta = theta.p)
mcov0

mcorr <- cov2cor(mcov0)
round(mcorr, 2)

(nc <- nrow(mcorr))
(theta.c <- seq(0.5, -0.5, length = nc))

ss <- diag(exp(theta.c))
mcov <- ss %*% mcorr %*% ss

round(100 * mcorr)
round(mcov, 1)

n <- 5000
m <- nrow(mcov)

ll <- chol(mcov)
xx <- matrix(rnorm(n * m), n) %*% ll

cov(xx)
cor(xx)

