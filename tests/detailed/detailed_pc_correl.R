library(INLA)

library(graphpcor)

n <- 4
lambda <- 5

model <- cgeneric(
    model = "pc_correl", 
    n = n,
    lambda = lambda,
    debug = 1e9)

graph(model, optimize = TRUE)

graph(model)

round(ith <- initial(model), 4)

m <- n * (n-1)/2
theta1 <- rnorm(m)

theta1
(qq <- prec(model, theta = theta1))

sum(diag(chol(qq)))

(vv <- solve(qq))

if(FALSE) {
    ## compare with inla.pc.cormat.dtheta + 1st Jacobian

    pc1 <- prior(model, theta = theta1)
    pc1
    
    ldJ1 <- sum(log(pi * exp(-theta1) / ( (1 + exp(-theta1))^2 ))) ## 1st Jacobian
    
    if(FALSE) {
        INLA:::inla.pc.cormat.dtheta
        INLA:::inla.pc.multvar.simplex.d
        INLA:::inla.pc.multvar.simplex.core
        INLA:::inla.pc.multvar.simplex.d.core
        INLA:::inla.pc.multvar.h.default
    }
    
    thetapi <- pi/(1+exp(-theta1))
    pc0 <- INLA:::inla.pc.cormat.dtheta(thetapi, lambda = lambda, log = TRUE)
    
    c(pc0 + ldJ1, pc1)

}

## call inla with NA to compare 
dat1 <- data.frame(    
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fit <- inla(
    y ~ 0 + f(i, model = model),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode
)

all.equal(qq, prec(fit))


m2 <- model; m2$f$cgeneric$debug <-
                 m2$f$cgeneric$data$ints$debug <- FALSE
th2corr <- function(th)
    chol2inv(chol(prec(m2, theta = th)))[lower.tri(diag(n))]

th2corr(rep(0,6))

s1 <- inla.cgeneric.sample(
    n=1000, result = fit, name = 'i',
    from.theta = th2corr,
    simplify = TRUE
)

str(s1)

par(mfrow = c(2, 3), mar = c(3,3,1,1), mgp = c(2,1,0), las = 1, bty = "n")
for(k in 1:6)
    hist(s1[k, ], 30, main = '', xlab = '', ylab = '', freq = FALSE)

### now fit some data

nrep <- 200
xx <- matrix(rnorm(nrep * n), nrep) %*% as.matrix(chol(vv))
str(xx)
cor(xx)

dat2 <- data.frame(
    i = rep(1:n, each = nrep),
    r = rep(1:nrep, n),
    y = as.vector(xx)
)
str(dat2)

fitr <- inla(
    y ~ 0 + f(i, model = model, replicate = r),
    data = dat2,
    control.family = cfam,
    control.inla = cinla
)

(Lfitted <- graphpcor:::theta2gamma2L(fitr$mode$theta))
round(tcrossprod(Lfitted), 2)
round(vv, 2)

ths0 <- t(as.matrix(do.call("expand.grid", lapply(1:5, function(x) seq(-2,2,.5)))))
str(ths0)

length(th0 <- seq(-3, 3, 0.1)+0.0001)
str(ths)

prior(m2, theta=rnorm(6))

system.time(pth <- sapply(1:6, function(i) {
    thr <- rbind(0, ths0)[c(i, setdiff(1:6, i)), ]
    sapply(th0, function(th) {
        thr[i, ] <- th
        mean(exp(prior(m2, theta = thr)))
    })
})
)

str(pth)

par(mfrow = c(2, 3), mar = c(3,3,1,1), mgp = c(2,1,0), las = 1, bty = "n")
for(k in 1:6)
    plot(th0, pth[, 1], type = 'o')

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
