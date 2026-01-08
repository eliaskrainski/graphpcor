library(INLA)

library(graphpcor)

n <- 4
lambda <- 5

model <- cgeneric(
    model = "pc_correl",
    n = n,
    lambda = lambda,
    debug = 1e9,
    useINLAprecomp = FALSE)
model

graph(model, optimize = TRUE)

graph(model)

round(ith <- initial(model), 4)

m <- n * (n-1)/2
theta1 <- rnorm(m)

theta1
(qq <- prec(model, theta = theta1))

sum(diag(chol(qq)))

(vv <- solve(qq))

basecor(theta1, n)

## call inla with NA to compare
dat1 <- data.frame(
    i = 1:n,
    y = rep(NA, n)
)

cinla <- list(int.strategy = 'eb')
cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
cmode <- list(theta = theta1, fixed = TRUE)

fitfix <- inla(
    y ~ 0 + f(i, model = model),
    data = dat1,
    control.family = cfam,
    control.inla = cinla,
    control.mode = cmode,
    verbose = !TRUE
)

all.equal(qq, prec(fitfix))

iil <- which(lower.tri(diag(n)))
iil

m2 <- cgeneric(
    model = "pc_correl",
    n = n,
    lambda = lambda,
    debug = FALSE,
    useINLAprecomp = FALSE)

th2corr <- function(th) {
    chol2inv(chol(prec(m2, theta = th)))[iil]
}

th2corr(rep(0,6))
th2corr(rep(-1,6))
th2corr(rep(1,6))

summary(t(sapply(1:1000, function(i)
                 th2corr(rnorm(m)))))


s1 <- inla.cgeneric.sample(
    n=1000, result = fitfix, name = 'i',
    from.theta = th2corr,
    simplify = TRUE
)

par(mfrow = c(2, 3), mar = c(3,3,1,1), mgp = c(2,1,0), las = 1, bty = "n")
for(k in 1:6)
    hist(s1[k, ], 30, main = '', xlab = '', ylab = '', freq = FALSE)


## fit the prior for different lambda

lambdas <- c(0.5, 1, 10, 100, 10000)

(th0 <- (1:m)-m/2)
b0 <- basepcor(th0, n)
b0
b0$base[iil]

par(mfrow = c(5, 6), mar = c(3,3,.1,.1), mgp = c(2,1,0), las = 1, bty = "n")
for(i in 1:length(lambdas)) {
    Cmodel <- cgeneric(
        model = "pc_correl",
        n = n, base = th0,
        lambda = lambdas[i],
        useINLAprecomp = FALSE)
    fit <- inla(
        y ~ 0 + f(i, model = Cmodel),
        data = dat1,
        control.family = cfam,
        control.inla = cinla
    )
    hs <- inla.hyperpar.sample(1000, fit)
    ccs <- sapply(1:nrow(hs), function(i) {
        L <- graphpcor:::Lprec0(
                             theta = hs[i, ],
                             p = n,
                             itheta = iil,
                             d0 = n:1)
        cov2cor(chol2inv(t(L)))[iil]
    })
    for(k in 1:6) {
        h <- hist(ccs[k, ], (-100:100)/100, plot = FALSE)
        ck <- b0$base[iil[k]]
        plot(h, xlim = c(max(-1, ck-0.2), min(ck+0.2,1)),
             main = '', xlab = '', ylab = '', freq = FALSE)
        abline(v = ck, col = 2, lty = 2, lwd = 2)
    }
}


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

(Cfitted <- basecor(fitr$mode$theta, p = n)$base)
round(vv, 2)

ths0 <- t(as.matrix(do.call(
    "expand.grid", lapply(1:5, function(x) seq(-2,2,.5)))))
str(ths0)

h0 <- 0.01
length(th0 <- c(seq(-1+h0/2, -h0/2, h0),
                seq(h0/2, 1-h0/2, h0)))
str(th0)

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
