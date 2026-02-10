
data(cars, package = "msos")

sdat <- scale(cars)

(n <- nrow(sdat))
(p <- ncol(sdat))

round((cc <- cov(sdat)) * 100)

image(cc)

lcc <- chol(cc + diag(p) * 0.0)
lcc

qc <- chol2inv(lcc)
round(qc, 2)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(colnames(sdat), colnames(sdat))
round(pC*100)

## define a graphpcor from pC threshold
rowSums(abs(pC)>0.5)

library(graphpcor)

g <- graphpcor(abs(pC)>0.3)
(dg <- dim(g))

c(n=n, p=p)
g
c(p, p*(p-1)/2)

plot(g)

ig <- graph_from_adjacency_matrix(attr(g, "graph"))
lypl <- layout.reingold.tilford(
    ig, circular = TRUE
)

set.seed(1)
par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(ig, arrow.mode = 0, layout = lypl)

plot(g, circular = TRUE)
plot(g, layout = lypl)

c0 <- cgeneric(g, lambda = 1,
               base = rep(0, dg[2]),
               useINLAprecomp = FALSE)

idat <- list(
    i = rep(1:p, each = n),
    r = rep(1:n, p),
    y = as.vector(sdat)
)

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
