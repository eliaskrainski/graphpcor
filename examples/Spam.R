
data(Spam, package = "msos")

sdat <- scale(Spam)

(n <- nrow(sdat))
(p <- ncol(sdat))

round((cc <- cov(sdat)) * 100)

image(cc)

length(jj <- 21:45)

image(cc[jj, jj])

cc <- cc[jj, jj]

lcc <- chol(cc + diag(ncol(cc)) * 0.0)
lcc

qc <- chol2inv(lcc)
round(qc, 2)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(colnames(sdat)[jj], colnames(sdat)[jj])
round(pC*100)

## define a graphpcor from pC threshold
rowSums(abs(pC)>0.1)

library(graphpcor)

g <- graphpcor(abs(pC)>0.05)

(dg <- dim(g))

c(p=p, length(jj), n=n)

g
c(p, p*(p-1)/2)

par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(g)

c0 <- cgeneric(g, lambda = 1,
               base = rep(0, dg[2]),
               debug = TRUE,
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
