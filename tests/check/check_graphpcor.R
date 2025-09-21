library(graphpcor)

## graph in Example 2.6 of the GMRF book
g <- graphpcor(x1~x2+x3, x2~x4, x3~x4)
g
(ne <- dim(g))

## Laplacian
(G <- Laplacian(g))

all.equal(g, graphpcor(G))
all.equal(g, graphpcor(G!=0))

## base model (theta for L)
theta0l <- rep(-0.3, ne[2])

hessian(g, theta0l)

## build the cgeneric model
## Note: here 'model' is a 'graphpcor'
cmodel <- cgeneric(
    model = g, ## use the graphpcor
    lambda = 1,
    base = theta0l,
    sigma.prior.reference = rep(1, ne[1]),
    sigma.prior.probability = rep(0.5, ne[1]),
    debug = 1e9)

graph(cmodel)

## define initial L
theta.low <- rep(-1, ne[2])
idx.low <- (!is.zero(G) & lower.tri(G))
L1a <- diag(ne[1]); L1a[idx.low] <- theta.low
L1a

## the C code to fill.in
idx.fill <- which(is.zero(G) & (!is.zero(t(chol(G + diag(ne[1]))))))
nfi <- length(idx.fill)
nfi

if(nfi) {
    ij.fill <- cbind(
        row(G)[idx.fill],
        col(G)[idx.fill])
    print(ij.fill)
    Lf <- .C("fillL",
             as.integer(ne[1]),
             as.integer(nfi),
             as.integer(ij.fill[, 1]-1),
             as.integer(ij.fill[, 2]-1),
             l=L1a)$l
    print(Lf)
    Q0 <- crossprod(Lf)
} else {
    Q0 <- crossprod(L1a)
}
Q0

V0 <- chol2inv(t(Lf))
V0

C1 <- cov2cor(V0)
C1
##           [,1]      [,2]      [,3]      [,4]
## [1,] 1.0000000 0.9486833 0.9128709 0.7745967
## [2,] 0.9486833 1.0000000 0.8660254 0.8164966
## [3,] 0.9128709 0.8660254 1.0000000 0.7071068
## [4,] 0.7745967 0.8164966 0.7071068 1.0000000


## sigma parameters
logsigmas <- log(c(0.3, 0.7, 1.5, 0.9))

V1 <- diag(exp(logsigmas)) %*% C1 %*% diag(exp(logsigmas))
V1
##           [,1]      [,2]      [,3]      [,4]
## [1,] 0.0900000 0.1992235 0.4107919 0.2091411
## [2,] 0.1992235 0.4900000 0.9093267 0.5143928
## [3,] 0.4107919 0.9093267 2.2500000 0.9545942
## [4,] 0.2091411 0.5143928 0.9545942 0.8100000

stopifnot(all(sqrt(diag(V1)) == exp(logsigmas)))

sr <- exp(logsigmas)/sqrt(diag(V0))
sr
##[1] 0.07745967 0.28577380 1.06066017 0.90000000

diag(sr) %*% V0 %*% diag(sr)

round(chol(V1), 4)
##      [,1]   [,2]   [,3]   [,4]
## [1,]  0.3 0.6641 1.3693 0.6971
## [2,]  0.0 0.2214 0.0000 0.2324
## [3,]  0.0 0.0000 0.6124 0.0000
## [4,]  0.0 0.0000 0.0000 0.5196

theta1 <- c(d = logsigmas,
            l = theta.low)

cmodel$f$cgeneric$data$doubles$thetabasescaled
## [1] -0.2190589 -0.2256071 -0.2690858 -0.2829668
drop(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]) %*% theta0l)

round(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]), 3)
###        [,1]   [,2]   [,3]   [,4]
### [1,]  2.855 -0.126 -0.052 -0.002
### [2,] -0.126  1.916 -0.018 -0.029
### [3,] -0.052 -0.018  1.005 -0.001
### [4,] -0.002 -0.029 -0.001  1.003

drop(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]) %*% theta.low)
##[1] -2.6759597 -1.7432539 -0.9338943 -0.9716843

drop(matrix(cmodel$f$cgeneric$data$matrices$h[-(1:2)], ne[2]) %*% theta.low) -
    cmodel$f$cgeneric$data$doubles$thetabase
## [1] -1.873172 -1.220278 -0.653726 -0.680179

Q1c <- prec(cmodel, theta = theta1)

round(Q1c, 4)
round(Q1 <- chol2inv(chol(V1)), 4)
##          [,1]     [,2]     [,3]    [,4]
## [1,] 166.6667 -45.1754 -12.1716  0.0000
## [2,] -45.1754  24.4898   0.0000 -3.8881
## [3,] -12.1716   0.0000   2.6667  0.0000
## [4,]   0.0000  -3.8881   0.0000  3.7037

all.equal(Q1, as.matrix(Q1c))

prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = theta1)

dataf <- list(
    i = 1:ne[1],
    y = rep(NA, ne[1]))

library(INLA)

fit1 <- inla(
    formula = y ~ 0 + f(i, model = cmodel),
    family = 'poisson',
    data = dataf,
    control.mode = list(theta = theta1, fixed = TRUE)
)

all.equal(Q1c, prec(fit1))


## PRIOR
library(INLA)
library(graphpcor)

## graph in Example 2.6 of the GMRF book
g <- graphpcor(x1~x2+x3, x2~x4, x3~x4)
(ne <- dim(g))
Laplacian(g)

th.base <- c(-1,0.5,1,0)
c.base <- vcov(g, theta = th.base)
c.base

lambs <- c(0.1, 0.5, 2, 10
           ); names(lambs) <- paste0("l", lambs)
lmodels <- lapply(lambs, function(l)
    cgeneric(
        model = g, ## use the graphpcor
        lambda = l, base = th.base,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.0, ne[1]),
        useINLAprecomp = FALSE)
    )

sapply(lmodels, function(x) x$f$cgeneric$data$doubles$lconst)

sapply(lmodels, initial)

ifits <- lapply(
    lmodels, function(mc) {
        inla(y ~ 0 + f(i, model = mc),
             data = data.frame(i=1:ne[1], y = NA),
             control.family = list(
                 hyper = list(prec = list(initial = 10, fixed = TRUE))
             )
             )
    }
)
sapply(ifits, function(r) r$mode$theta)

lapply(ifits, function(x) grep("aluations =", x$logfile, value = TRUE))

vcov(g, theta = ifits[[1]]$mode$theta)
g2c <- function(th) vcov(g, theta=th)[lower.tri(diag(ne[1]))]
g2c(ifits[[1]]$mode$theta)

inla.cgeneric.sample(n = 2, result = ifits[[1]], name = 'i', from.theta = g2c)

sr <- lapply(ifits, function(r)
    inla.cgeneric.sample(
        n = 2000, result = r, name = 'i',
        from.theta = g2c, simplify = TRUE))
str(sr,1)

cnams <- c("c[2,1]", "c[3,1]", "c[4,1]",
           "c[3,2]", "c[4,2]", "c[4,3]")

png("corr4priors.png", width = 1600, height=1000, res = 100)
par(mfrow = c(length(sr), 6), mar = c(4,4,2,1), mgp = c(3,1,0))
for(i in 1:length(sr)) {
    for(k in 1:6) {
        ck0 <- c.base[lower.tri(diag(ne[1]))][k]
        r <- range(ck0, sr[[i]][k,])+c(-0.1, 0.1)
        r[r<(-1)] <- -1
        r[r>1] <- 1
##        r <- c(-1,1)
        hist(sr[[i]][k,], seq(r[1], r[2], length = 30),
             freq=FALSE,
        main = paste0("lambda = ", lambs[i], " : ", cnams[k]), xlab = '')
        abline(v=ck0, lwd = 2, col = 2)
    }
}
dev.off()

system("eog corr4priors.png &")

###

## marginalizing the cgneric prior, not through inla
m <- dim(g)[2]
hth <- 0.2
th0 <- seq(-4, 4, hth)
ths0 <- t(as.matrix(do.call("expand.grid", lapply(1:(m-1), function(x)
    th0))))
str(ths0)

## number of evals /M
length(th0) * ncol(ths0) * m /1e6

str(prior(lmodels[[1]], theta = matrix(rnorm(m*3), m)))

library(parallel)
system.time(pth <- mclapply(lmodels, function(cm) {
    sapply(1:m, function(i) {
        thr <- rbind(0, ths0)
        thr[1, ] <- thr[i, ]
        thr[i, ] <- 0.0
        sapply(th0, function(th) {
            thr[i, ] <- th
            sum(exp(prior(cm, theta = thr) + log(hth)*(m-1)))
        })
    })
}, mc.cores = length(lmodels)))

str(pth)

sapply(pth, apply, 2, function(x)
       sum(hth * x, na.rm = TRUE))

par(mfrow = c(length(lambs),m), mar = c(3,3,1,1),
    mgp = c(2,1,0), las = 1, bty = "n")
for(i in 1:length(lambs)) {
    for(k in 1:m) {
        plot(th0, pth[[i]][, k], type = 'o', xlab = '', ylab = '',
             main = paste0("lambda = ", lambs[i]))
        abline(v=th.base[k])
    }
}



#############################################################################
####

rho0 <- c(-0.3, 0.7)
C0 <- matrix(c(1, rho0, 
               rho0[1], 1, prod(rho0), 
               rho0[2], prod(rho0), 1), 3)

library(graphpcor)
gpc <- graphpcor(x1 ~ x2 + x3)

H <- hessian(gpc, x = C0)
(theta0 <- attr(H, "base"))

ptheta <- function(theta, H) {
    xi <- attr(H, "h.5") %*% (theta - attr(H, "base"))
    rphi <- graphpcor:::x2rphi(xi)
    
}

all.equal(C0,
          vcov(gpc, theta = theta0))
prec(gpc, theta = theta0)

cpc <- cgeneric(gpc, base = theta0, lambda = 5,
          ##      debug = 100000,
                useINLAprecomp = FALSE)

graph(cpc)

all.equal(prec(cpc, theta = theta0),
          prec(gpc, theta = theta0))

all.equal(prec(cpc, theta = theta0+c(0.33, -0.1)),
          prec(gpc, theta = theta0+c(0.33, -0.1)))

prior(cpc, theta = cbind(c(0,0), c(1,1)))

h0 <- 0.05
x0 <- seq(-10+h0/2, 10-h0/2, h0)
nx0 <- length(x0)
xx <- t(expand.grid(x0 + theta0[1],
                    x0 + theta0[2]))

summary(t(xx))

dxx <- prior(cpc, theta = xx)

dxx2 <- matrix(exp(dxx), nx0)

sum(h0 * rowSums(dxx2*h0))
sum(h0 * colSums(dxx2*h0))

par(mfrow = c(1, 2))
plot(x0+theta0[1], rowSums(dxx2*h0))
plot(x0+theta0[2], colSums(dxx2*h0))


library(plot3D)

surf3D(matrix(xx[1, ], nx0),
       matrix(xx[2, ], nx0),
       matrix(exp(dxx), nx0), 
       colkey = !TRUE,
       xlab = '', ylab  = '', zlab = '',
       box = TRUE, bty = "b")

  text(-0.15, -0.43, expression(xi[1]), adj = 1)
  text(0.28, -0.37, expression(xi[2]), adj = 1)
  legend("topleft", bty = "n", as.expression(llab[[i]]), cex = 2)


ribbon3D(x0+theta0[1], x0+theta0[2],
         matrix(exp(dxx), nx0),
         ticktype = 'detailed',
         curtain = TRUE)


