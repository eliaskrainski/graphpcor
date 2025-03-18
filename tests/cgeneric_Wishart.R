library(INLA)

library(graphpcor)

(V <- matrix(c(2,-1,1, -1,7,3, 1,3,3), 3))
(cc <- cov2cor(V))
(lV <- chol(V))
(Q <- chol2inv(lV))
(lQ <- chol(Q))

n <- nrow(V)
dof <- 10
RR <- toeplitz(n:1)
R <- c(diag(RR), RR[lower.tri(RR)])

W <- cgeneric(
    model = "Wishart",
    n = n,
    dof = dof,
    R = R
)

graph(W, optimize = TRUE)

graph(W)

round(ini <- initial(W), 4)

theta1 <- c(log(diag(lQ)), lQ[upper.tri(lQ)])
theta1

(theta0 <- c(-log(diag(V)), log((1 + cc[upper.tri(cc)]) / (1 -cc[upper.tri(cc)]))))

dat1 <- data.frame(    
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fit0 <- inla(
    y ~ 0 + f(i, model = 'iid3d', order = n, n = n*1),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta0, fixed = TRUE)
)

fit1 <- inla(
    y ~ 0 + f(i, model = 'iidkd', order = n, n = n*1,
              hyper = list(theta1 = list(param = c(dof, rep(1, n), rep(0, n*(n-1)/2))))),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

fit2 <- inla(
    y ~ 0 + f(i, model = W),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

pp3 <- list(prec(fit0),
            prec(fit1),
            prec(fit2))
pp3[[2]]

all.equal(pp3[[1]], pp3[[2]])
all.equal(pp3[[1]], pp3[[3]])

cbind(fit0$mlik, fit1$mlik, fit2$mlik)

### 
nrep <- 300
xx <- matrix(rnorm(nrep * n), nrep) %*% lV
str(xx)

dat2 <- data.frame(
    i = rep(1:n, each = nrep),
    r = rep(1:nrep, n),
    y = as.vector(xx)
)
str(dat2)

fit0r <- inla(
    y ~ 0 + f(i, model = 'iid3d', order = n, n = n, replicate = r), 
    data = dat2,
    control.family = cfam,
    control.inla = cinla
)

fit1r <- inla(
    y ~ 0 + f(i, model = 'iidkd', order = n, n = n, replicate = r, 
              hyper = list(theta1 = list(
                               param = c(dof, rep(1, n), rep(0, n*(n-1)/2))))),
    data = dat2,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta1)
)

fit2r <- inla(
    y ~ 0 + f(i, model = W, replicate = r),
    data = dat2,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta1)
)

cbind(fit0r$mlik, fit1r$mlik, fit2r$mlik)

rbind(fit2r$mode$theta, fit1r$mode$theta)

(fit2r$mode$theta - fit1r$mode$theta)/abs(fit1r$mode$theta)

-1 + fit1r$summary.hyperpar[, 2] / fit2r$summary.hyperpar[, 2]

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
