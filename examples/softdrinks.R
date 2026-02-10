
library(INLA)
library(graphpcor)

inla.setOption(
    safe = FALSE,
    num.threads = 6L
)

data(softdrinks, package = "msos")

## 8 variables
set.seed(1)
jj <- sample(1:ncol(softdrinks))
dat <- softdrinks[, jj]
colnames(dat) <- gsub("7up", "SevenUp", colnames(dat))
(p <- ncol(V <- cov(dat)))
V

round(100 * cov2cor(V))

image(cov2cor(V))

lV <- chol(V)
lV

Q <- chol2inv(lV)
round(Q, 2)

## partial correlation matrix
pC <- cov2cor(cov(dat) + diag(p)*.0001)
dimnames(pC) <- dimnames(Q) <- dimnames(V) <-
    list(colnames(dat), colnames(dat))
round(pC*100)

## define a graphpcor from pC threshold
library(graphpcor)
rowSums(abs(pC)>0.5)
g <- graphpcor(abs(pC)>0.7)

g
c(p, p*(p-1)/2)

par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(g)

