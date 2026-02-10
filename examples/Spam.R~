
## A double matrix with 130 observations on the following 7 variables.
##  Baseball Baseball’s ranking out of seven sports.
## Football Football’s ranking out of seven sports.
## Basketball Basketball’s ranking out of seven sports.
## Tennis Tennis’ ranking out of seven sports.
## Cycling Cycling’s ranking out of seven sports.
## Swimming Swimming’s ranking out of seven sports.
## Jogging Jogging’s ranking out of seven sportssbp systolic blood pressure

data(sportsranks, package = "msos")

pairs(sportsranks)

sdat <- scale(sportsranks)

(n <- nrow(sdat))
(p <- ncol(sdat))

round((cc <- cov(sdat)) * 100)

image(cc)

lcc <- chol(cc + diag(p) * 0.01)
lcc

qc <- chol2inv(lcc)
round(qc, 2)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(colnames(sdat), colnames(sdat))
round(pC*100)

## define a graphpcor from pC threshold
rowSums(abs(pC)>0.94)

library(graphpcor)

g <- graphpcor(abs(pC)>0.93)
(dg <- dim(g))

c(n=n, p=p)
g
c(p, p*(p-1)/2)

par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(g)

c0 <- cgeneric(g, lambda = 1,
               base = rep(0, dg[2]), 
               useINLAprecomp = FALSE)

idat <- data.frame(
    i = rep(1:p, each = n),
    r = rep(1:n, p),
    y = as.vector(sdat)
)

head(sdat,3)
head(idat,3)

cfam <- list(hyper = list(
                 prec = list(intial = 10, fixed = TRUE)
             ))

library(INLA)

fit <- inla(
    y ~ 0 + f(i, model = c0, replicate = r),
    data = idat,
    control.family = cfam,
    verbose = !TRUE
)

fit$mode$theta
