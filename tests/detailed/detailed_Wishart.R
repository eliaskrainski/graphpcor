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
    R = R,
    debug = 1e8 * 1## print lots of details
)

cgeneric_graph(W, optimize = TRUE)

cgeneric_graph(W)

round(ini <- cgeneric_initial(W), 4)

theta1 <- c(log(diag(lQ)), lQ[upper.tri(lQ)])
theta1

l2q <- function(l) {
    m <- length(l)
    n <- (sqrt(1+8*m)-1)/2
    q <- diag(exp(l[1:n]))
    q[upper.tri(q, diag = FALSE)] <- l[-(1:n)]
    crossprod(q)
}

myJ <- function(l, h = 0.005, verbose = FALSE) {
    h.2 <- h * 2
    m <- length(l)
    J <- matrix(0, m, m)
    for(i in 1:m) {
        l2 <- l1 <- l
        l1[i] <- l[i] - h
        l2[i] <- l[i] + h
        q1 <- l2q(l1)
        dd <- (l2q(l2) - q1)/h.2
        J[i, ]  <- c(diag(dd), dd[lower.tri(dd, diag = FALSE)])
    }
    if(verbose) {
        print(J)
        print(diag(qr(J)$qr))     
    }
    det(J)
}

myJ(rnorm(length(R)), verbose = TRUE)

cgeneric_prior(W, theta = theta1)
myJ(theta1, 1e-3, verbose = TRUE)
myJ(theta1, 1e-5, verbose = TRUE)

cgeneric_Q(W, theta = theta1)

theta2iidkd <- function(th, old = FALSE, covar = TRUE, corr = FALSE) {
    m <- length(th)
    n <- (sqrt(1+8*m)-1)/2
    if(old) {
        V <- diag(rep(0.5, n), nrow = n, ncol = n)
        V[upper.tri(V, diag = FALSE)] <- 2/(1+exp(-th[-(1:n)])) -1
        V <- V + t(V)
        ss <- exp(-th[1:n]/2)
        V <- t(V * ss) * ss
        if(!covar) return(chol2inv(chol(V)))
        if(corr) {
            V <- cov2cor(V)
        }
    } else {
        V <- diag(exp(th[1:n]))
        V[lower.tri(V)] <- th[-(1:n)]
        if(covar|corr) {
            V <- chol2inv(t(V))
            if(corr) {
                V <- cov2cor(V)
            }            
        } else {
            return(tcrossprod(V))
        }
    }
    return(V)
}

(theta0 <- c(-log(diag(V)), log((1 + cc[upper.tri(cc)]) / (1 -cc[upper.tri(cc)]))))

V
theta2iidkd(theta0, old = TRUE)
theta2iidkd(theta1, old = FALSE)
solve(cgeneric_Q(W, theta = theta1))

theta2iidkd(theta0, old = TRUE, covar = FALSE)
theta2iidkd(theta1, old = FALSE, covar = FALSE)

dW <- function(Q, R, d) {
    Q <- as.matrix(Q)
    R <- as.matrix(R)
    lq <- chol(Q)
    lr <- chol(R)
    p <- nrow(Q)
    tr <- sum(diag(crossprod(R,Q)))
    hldr <- sum(log(diag(lr)))
    hldq <- sum(log(diag(lq)))
    n1 <- hldq * (d - p -1) - tr/2
    d1 <- 0.5 * d * p * log(2) - d * hldr
    d2 <- 0.25 * (p * (p-1)) * log(pi) + sum(lgamma((d + 1 - (1:p))/2)) 
    cat('hldr =', hldr, 'hldq = ', hldq, 'tr = ', tr,
        'n = ', n1, 'd1 = ', d1, 'd2 = ', d2,
        'cprior = ', d1+d2, '\n')
    return(n1 - d1 - d2)
}

theta2prior <- function(th, R, d) {
    Q <- theta2iidkd(th, old = FALSE, covar = FALSE)
    dW(Q, R, d) + log(abs(myJ(th)))
}

cgeneric_prior(W, theta = theta1) -log(abs(myJ(theta1)))
dW(Q, RR, dof) 
theta2cgeneric_prior(theta1, RR, dof) -log(abs(myJ(theta1)))

c(qq=CholWishart::dWishart(
                      as.matrix(Q),
                      dof,
                      as.matrix(RR)),
  qv=CholWishart::dWishart(
                      as.matrix(Q),
                      dof,
                      solve(as.matrix(RR))),
  vq=CholWishart::dWishart(
                      solve(as.matrix(Q)),
                      dof,
                      as.matrix(RR)),
  vv=CholWishart::dWishart(
                      solve(as.matrix(Q)),
                      dof,
                      solve(as.matrix(RR)))
  )

dat1 <- data.frame(    
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)
ccomp <- list(config = TRUE)

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
    control.mode = cmode,
    control.compute = ccomp
)

fit2 <- inla(
    y ~ 0 + f(i, model = W),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode,
    control.compute = ccomp
)

c(all.equal(Q, as.matrix(cgeneric_Q(fit0))),
  all.equal(cgeneric_Q(fit0),
            cgeneric_Q(fit1)),
  all.equal(cgeneric_Q(fit0),
            cgeneric_Q(fit2)))

cbind(fit0$mlik, fit1$mlik, fit2$mlik)

### 
nrep <- 2

dat2 <- data.frame(
    i = rep(1:n, each = nrep),
    r = rep(1:nrep, n),
    y = NA
)
str(dat2)

fit0r <- inla(
    y ~ 0 + f(i, model = 'iid3d', order = n, n = n, replicate = r), 
    data = dat2,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta0, fixed = TRUE)
)

fit1r <- inla(
    y ~ 0 + f(i, model = 'iidkd', order = n, n = n, replicate = r, 
              hyper = list(theta1 = list(
                               param = c(dof, rep(1, n), rep(0, n*(n-1)/2))))),
    data = dat2,
    control.family = cfam,
    control.inla = cinla,
    control.mode = list(theta = theta1, fixed = TRUE)
)

fit2r <- inla(
    y ~ 0 + f(i, model = W, replicate = r),
    data = dat2,
    control.family = cfam,
    control.inla = cinla, ##verbose =T,
    control.mode = list(theta = theta1, fixed = TRUE)
)

Q
round(cgeneric_Q(fit0r), 3)
round(cgeneric_Q(fit1r), 3)
round(cgeneric_Q(fit2r), 3)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
