
library(INLA)

library(graphpcor)

sigmas <- c(2, 1, 0.5)
cc <- matrix(c( 1.0,  0.8, -0.5,
               0.8,  1.0, -0.4,
               -0.5, -0.4,  1.0), 3)
(V <- diag(sigmas) %*% cc %*% diag(sigmas))

(lV <- chol(V))
(Q <- chol2inv(lV))
(lQ <- chol(Q))

n <- nrow(V)
dof <- n+10
##RR <- toeplitz(n:1)
RR <- diag(n)
R <- c(diag(RR), RR[lower.tri(RR)])

W <- cgeneric(
    model = "Wishart",
    n = n,
    dof = dof,
    R = R,
##    debug = 1e9,
    useINLAprecomp = FALSE
)

W

str(W)

graph(W, optimize = TRUE)

graph(W)

round(ini <- initial(W), 4)

theta1 <- c(log(diag(lQ)), lQ[upper.tri(lQ)])
theta1

prec(W, theta = theta1)


(theta0 <- c(-log(diag(V)), log((1 + cc[upper.tri(cc)]) / (1 -cc[upper.tri(cc)]))))

dat0 <- data.frame(    
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fit0 <- inla(
    y ~ 0 + f(i, model = paste0('iid', n, 'd'), order = n, n = n*1),
    data = dat0,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta0, fixed = TRUE)
)

fit1 <- inla(
    y ~ 0 + f(i, model = 'iidkd', order = n, n = n*1,
              hyper = list(theta1 = list(param = c(dof, rep(1, n), rep(0, n*(n-1)/2))))),
    data = dat0,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)


fit2 <- inla(
    y ~ 0 + f(i, model = W),
    data = dat0,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

pp3 <- list(prec(fit0),
            prec(fit1),
            prec(fit2))
pp3[[1]]

c(q12=all.equal(pp3[[1]], pp3[[2]]),
  q13=all.equal(pp3[[1]], pp3[[3]]))

cbind(fit0$mlik, fit1$mlik, fit2$mlik)

### 
nrep <- 1000
xx <- matrix(rnorm(nrep * n), nrep) %*% lV
str(xx)

dat1 <- data.frame(
    i = rep(1:n, each = nrep),
    r = rep(1:nrep, n),
    y = as.vector(xx)
)
str(dat1)

fit0r <- inla(
    y ~ 0 + f(i, model = 'iid3d', order = n, n = n, replicate = r), 
    data = dat1,
    control.family = cfam,
    control.inla = cinla
)

fit1r <- inla(
    y ~ 0 + f(i, model = 'iidkd', order = n, n = n, replicate = r, 
              hyper = list(theta1 = list(
                               param = c(dof, rep(1, n), rep(0, n*(n-1)/2))))),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta1)
)

fit2r <- inla(
    y ~ 0 + f(i, model = W, replicate = r),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta1)
)

fne <- c(grep("evaluations = ", fit0r$logfile, value = TRUE),
         grep("evaluations = ", fit1r$logfile, value = TRUE),
         grep("evaluations = ", fit2r$logfile, value = TRUE))
fne

fn <- as.integer(sapply(strsplit(fne, "= "), tail, 1))
fn

c(fit0r$cpu.used["Total"],
  fit1r$cpu.used["Total"],
  fit2r$cpu.used["Total"]) / fn

cbind(fit0r$mlik, fit1r$mlik, fit2r$mlik)

rbind(fit2r$mode$theta, fit1r$mode$theta)

(fit2r$mode$theta - fit1r$mode$theta)/abs(fit1r$mode$theta)

-1 + fit1r$summary.hyperpar[, 2] / fit2r$summary.hyperpar[, 2]

### p = 10
p <- 10
m.p <- p * (p-1)/2
dof.p <- 10 + p
Wp <- cgeneric("Wishart", n = p, dof = dof.p,
               R = rep(1:0, c(p,m.p)), useINLAprecomp = FALSE)

datp <- data.frame(
    i = rep(1:p, each = nrep),
    r = rep(1:nrep, p),
    y = rnorm(nrep * p)
)

cmdp <- list(theta = rep(0, p+m.p), restart = TRUE)

fitp <- inla(
    y ~ 0 + f(i, model = 'iidkd', order = p, n = p, replicate = r, 
              hyper = list(theta1 = list(
                               param = c(dof.p, rep(1, p), rep(0,m.p))))), 
    data = datp,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmdp
)

fitWp <- inla(
    y ~ 0 + f(i, model = Wp, replicate = r),
    data = datp,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmdp
)

fnp0 <- c(grep("evaluations = ", fitp$logfile, value = TRUE),
          grep("evaluations = ", fitWp$logfile, value = TRUE))
fnp0
fnp <- as.integer(sapply(strsplit(fnp0, "= "), tail, 1))
fnp

c(fitp$cpu.used["Total"],
  fitWp$cpu.used["Total"]) / fnp

cbind(fitp$mlik, fitWp$mlik)

summary((fitp$mode$theta - fitWp$mode$theta)/abs(fitWp$mode$theta))
summary(-1 + fitp$summary.hyperpar[, 2] / fitWp$summary.hyperpar[, 2])

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
