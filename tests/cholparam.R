
library(graphpcor)

p <- 6
gs <- graphpcor(paste0("X1~", paste0("X",2:p)))
gs

Ls <- Laplacian(gs)
Ls

chol(Ls + diag(p))

oj <- c(2:p,1)
ogs <- graphpcor(Ls[oj, oj])
ogs
Laplacian(ogs)
chol(Laplacian(ogs) + diag(p))

par(mfrow=c(1,2))
plot(gs)
plot(ogs)

## A correlation matrix
thp <- rnorm(dim(gs)[2], -1)

bp <- basepcor(thp, p, iLtheta = gs)
bp

bp$base[oj, oj]

obp <- basepcor(bp$base[oj, oj], iLtheta = ogs)
obp

round(chol(bp$base), 4)
round(chol(obp$base), 4)

hessian(bp)
hessian(obp)

## simulate data
##thd <- rnorm(dim(gs)[2], -1)
##cc <- basepcor(thd, p, gs)
##cc$base
##Lc <- chol(cc$base)

Lc <- chol(bp$base)

n <- 50

xx <- scale(matrix(rnorm(n*p), n) %*% Lc)
cov(xx)

m1 <- cgeneric(gs, base = bp, lambda = 5, useINLAprecomp = FALSE)
m2 <- cgeneric(ogs, base = obp, lambda = 5, useINLAprecomp = FALSE)

idat <- data.frame(
    i = rep(1:p, each = n),
    r = rep(1:n, p),
    y1 = as.vector(xx),
    y2 = as.vector(xx[, oj])
)

cfam <- list(hyper = list(prec = list(initial = 20, fixed = TRUE)))

library(INLA)

fit1 <- inla(
    formula = y1 ~ 0 + f(i, model = m1),
    data = idat, control.family = cfam,
    control.mode = list(
        theta = rep(0, dim(gs)[2]),
        restart = TRUE)
)

fit2 <- inla(
    formula = y2 ~ 0 + f(i, model = m2),
    data = idat, control.family = cfam,
    control.mode = list(
        theta = rep(0, dim(gs)[2]),
        restart = TRUE)
)

c(fit1$misc$nfunc, fit2$misc$nfunc)

c(fit1$cpu.used[["Total"]], fit2$cpu.used[["Total"]])

c(fit1$misc$nfunc, fit2$misc$nfunc)/
    c(fit1$cpu.used[["Total"]], fit2$cpu.used[["Total"]]) 

list(
    obs = cor(xx)[oj, oj],
    fitted = vcov(ogs, theta = fit2$mode$theta)
)

all.equal(
    vcov(gs, theta = fit1$mode$theta)[oj, oj],
    vcov(ogs, theta = fit2$mode$theta))
