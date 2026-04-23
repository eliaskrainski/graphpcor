## 1. setup a random correlation matrix: C
## 2. drawn samples from N(0, C)
## 3. set PC and LKJ priors for C
## 4. fit the model considering each prior
## 5. visualize the results

near.zero.corr <- FALSE

library(graphpcor)
library(INLA)

n <- 50
p <- 15 ## with p>10, 'ccd' is not available, but 'eb' is ok
(m <- p*(p-1)/2)

## 1. setup a random correlation matrix: C
if(near.zero.corr) {
    ## lower than base model (favor LKJ)
    theta.true <- rep(0, m) 
} else {
    ## above base model, favor PC but still far from the base model!
    theta.true <- log(rep(p:1, p:1-1))/rep((0.25+1:p)^2, p:1-1) 
}
L.true <- cholcor(theta.true)
C.true <- tcrossprod(L.true)
round(C.true * 100)

## 2. drawn samples from N(0, C), (but n*10 samples)
x <- matrix(rnorm(n*p*10), ncol = p) %*% t(L.true)
round(cor(x) * 100) ## for all samples

## take the first n samples to be used
cc.obs.n <- cor(x[1:n, ])
th.obs.n <- graphpcor:::corr2theta(cc.obs.n)
dataf <- data.frame(
    r = rep(1:n, p),
    i = rep(1:p, each = n),
    y = as.vector(x[1:n, ]) ## first n samples
)

## 3. set PC and LKJ priors for C
## base model: exchangeable
r <- 0.7
round(Cbase <- toeplitz(c(1,rep(r,p-1))), 2)
b <- basecor(Cbase)
b

## PC-prior at the base model
m1 <- cgeneric(b, lambda = 5)
## LKJ prior model
m0 <- cgeneric("LKJ", n = p, eta = 5)

## 4. fit the model considering each prior

theta.ini <- rep(ifelse(near.zero.corr,0.1,1), length(theta.true))
##theta.ini <- theta.true + rnorm(length(theta.true), 0, 0.2)

fit <- inla(
    formula = y ~ 0 + f(i, model = m1, replicate = r),
    data = dataf,
    control.family = list(
        hyper = list(
            prec = list(initial = 10, fixed = TRUE)
        )
    ),
    control.inla = list(
        int.strategy = ifelse(p>10,'eb','ccd') 
    ),
    control.mode = list(
        theta = theta.ini,
        restart = TRUE
    )
)

fit0 <- inla(
    formula = y ~ 0 + f(i, model = m0, replicate = r),
    data = dataf,
    control.family = list(
        hyper = list(
            prec = list(initial = 10, fixed = TRUE)
        )
    ),
    control.inla = list(
        int.strategy = ifelse(p>10,'eb','ccd') 
    ),
    control.mode = list(
        theta = theta.ini, 
        restart = TRUE
    )
)

fit$misc$nfunc ## lower it providing start around theta.true
fit0$misc$nfunc ## lower it providing start around theta.true

## sample from the joint posterior for the hyperparamters
hsamples <- inla.hyperpar.sample(
    n = 10000,
    result = fit,
    intern = TRUE
)
h0samples <- inla.hyperpar.sample(
    n = 10000,
    result = fit0,
    intern = TRUE
)

iiu <- which(upper.tri(C.true))
## compute the correlation at each hyperpar sample
ccsamples <- t(sapply(1:nrow(hsamples), function(i) {
    ll <- cholcor(hsamples[i, ], p = p, parametrization = 'cpc')
    tcrossprod(ll)[iiu]
}))
cc0samples <- t(sapply(1:nrow(h0samples), function(i) {
    ll <- cholcor(h0samples[i, ], p = p, parametrization = 'cpc')
    tcrossprod(ll)[iiu]
}))

library(coda)
ci1 <- HPDinterval(as.mcmc(ccsamples))
ci0 <- HPDinterval(as.mcmc(cc0samples))

c(mean(ci1[,2]-ci1[,1]), mean(ci0[,2]-ci0[,1])) ## average interval length 

cc.true <- C.true[iiu]
mean((cc.true >= ci1[,1]) & (cc.true <= ci1[,2])) ## truth coverage
mean((cc.true >= ci0[,1]) & (cc.true <= ci0[,2])) ## truth coverage

mean((cc.true >= r) & (cc.true <= r)) ## prior coverage
mean((cc.true >= 0) & (cc.true <= 0)) ## prior coverate

## 5. visualize the results
k2 <- k1 <- 1; 
par(mfrow = c(p, p), mar = c(1.5,1.5,0.1,0.1), mgp = c(0.5,0.5,0), bty = 'n')
for(i in 1:p) {
    for(j in 1:p) {
        if (i==j) {
            if(i==1) {
                plot(0, 0, type = 'n', axes = FALSE, xlab = '', ylab = '')
                legend("topleft", bty = 'n',
                       c("TRUE", "sample", "base", "PC-prec", "LKJ"), 
                       col = c(2, 3, 1, 4, 6), lty = c(3,3,3,1,1), lwd = 2)
            } else {
                plot(0, 0, pch = paste(i), axes = FALSE, xlab = '', ylab = '')
            }
        }
        if(j>i) {
#            h <- hist(ccsamples[, k1], 50, plot = FALSE)
 #           h0 <- hist(cc0samples[, k1], 50, plot = FALSE)
            d1 <- density(ccsamples[, k1])
            d0 <- density(cc0samples[, k1])
            xlm <- range(r, cc.obs.n[iiu[k1]],
                         #h$breaks, h0$breaks,
                         d1$x, d0$x)
            ylm <- range(#h$dens, h0$dens,
                         d1$y, d0$y)
  #          plot(h, xlim = xlm, ylim = ylm, border = 'transparent',
   #              freq = FALSE, xlab = '', ylab = '', main = '')
    #        plot(h0, add = TRUE, freq = FALSE,
     #            col = rgb(0.3,0.7,1,0.5), border = 'transparent')
            plot(d1, lwd = 2, col = 4,
                 xlim = xlm, ylim = ylm,
                 xlab = '', ylab = '', main = '')
            lines(d0, lwd = 2, col = 6)
            abline(v = c(C.true[iiu[k1]], cc.obs.n[iiu[k1]], r),
                   lty = 3, lwd = 2, col = c(2,3,1))
            k1 <- k1 + 1
        }
        if(i>j) {
            mg1 <- inla.smarginal(fit$internal.marginals.hyperpar[[k2]])
            mg0 <- inla.smarginal(fit0$internal.marginals.hyperpar[[k2]])
            xlm <- range(b$theta[k2], mg1$x, mg0$x)
            ylm <- range(mg1$y, mg0$y)
            plot(mg1, type = 'l', xlim = xlm, ylim = ylm,
                 xlab = '', ylab = '', main = '', lwd = 2, col = 4)
            lines(mg0, lwd = 2, col = 6)
            abline(v = c(theta.true[k2], th.obs.n[k2], b$theta[k2]),
                   lty = 3, lwd = 2, col = c(2,3,1))
            k2 <- k2 + 1
        }
    }
}
