library(INLA)

library(graphpcor)

n <- 4
lambda <- 5

model <- cgeneric(
    model = "pc_prec_correl", 
    n = n,
    lambda = lambda,
##    debug = 1e9,
    useINLAprecomp = FALSE)

graph(model, optimize = TRUE)

graph(model)

round(ith <- initial(model), 4)

m <- n * (n-1)/2
theta1 <- rnorm(m)

theta1
(qq <- prec(model, theta = theta1))

sum(diag(chol(qq)))

(vv <- solve(qq))

prior(model, theta = theta1)

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

pc1 <- prior(model, theta = theta1)
pc1

thb <- rep(0, length(theta1))
Hb <- graphpcor:::theta2H(thb)
graphpcor:::dtheta(theta1, lambda=lambda,
                   theta.base = thb, 
                   H.elements = Hb)

lseq <- c(0.2, 0.5, 1, 2, 5, 20); names(lseq) <- paste0("lambda", lseq)
lsm <- lapply(lseq, function(l) {
    cgeneric(
        model = "pc_prec_correl", 
        n = n,
        lambda = l)    
})

th2corr <- function(th, m)
    chol2inv(chol(as.matrix(prec(m, theta = th))))[lower.tri(diag(n))]

th2corr(rep(0,6), lsm[[1]])

lfit <- lapply(lsm, function(cm) 
    fit <- inla(
        y ~ 0 + f(i, model = cm),
        data = dat1,
        control.family = cfam,
        control.inla = cinla,
        control.mode = cmode
    )
    )

ssl <- lapply(1:length(lsm), function(i)
    inla.cgeneric.sample(
        n=1000, result = lfit[[i]], name = 'i',
        from.theta = function(th) th2corr(th, lsm[[i]]),
        simplify = TRUE
    )
)

str(ssl)

par(mfrow = c(6,6), mar = c(3,3,1,1), mgp = c(2,1,0), las = 1, bty = "n")
for(m in 1:6) {
    for(k in 1:6) {
        hist(ssl[[m]][k, ], 30, main = '', xlab = '', ylab = '', freq = FALSE)
    }
}

## marginals
sapply(lfit, function(fit) 
    sapply(fit$internal.marginals.hyperpar, function(m)
           inla.pmarginal(10, m)))

par(mfrow = c(6,6), mar = c(3,3,1,1), mgp = c(2,1,0), las = 1, bty = "n")
for(m in 1:length(lsm)) {
    for(k in 1:6) {
        plot(lfit[[m]]$internal.marginals.hyperpar[[k]])
    }
}

## marginalizing the cgneric prior, not through inla
ths0 <- t(as.matrix(do.call("expand.grid", lapply(1:5, function(x) seq(-2,2,.5)))))
str(ths0)

length(th0 <- seq(-3, 3, 0.1)+0.0001)

system.time(pth <- lapply(lsm, function(cm) {
    sapply(1:6, function(i) {
        thr <- rbind(0, ths0)
        thr[1, ] <- thr[i, ]
        thr[i, ] <- 0.0
        sapply(th0, function(th) {
            thr[i, ] <- th
            mean(exp(prior(cm, theta = thr)))
        })
    })
}))

str(pth)

sapply(pth, apply, 2, function(x)
       sum(0.1 * x))

par(mfrow = c(6,6), mar = c(3,3,1,1), mgp = c(2,1,0), las = 1, bty = "n")
for(m in 1:6) {
    for(k in 1:6)
        plot(th0, pth[[m]][, 1], type = 'o', xlab = '', ylab = '',
             main = paste0("lambda = ", lseq[m]))
}

### fit some data
nrep <- 200
xx <- matrix(rnorm(nrep * n), nrep) %*% as.matrix(chol(vv))
str(xx)

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
    control.inla = cinla,
    control.mode = cmode
)

(Lfitted <- graphpcor:::theta2gamma2L(fitr$mode$theta))
round(tcrossprod(Lfitted), 2)
round(vv, 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
