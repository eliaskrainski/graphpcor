library(graphpcor)

## A correlation matrix
c1 <- matrix(c( 1.0,  0.8, -0.5,
                0.8,  1.0, -0.4,
               -0.5, -0.4,  1.0), 3)
(p <- ncol(c1))

## The cgeneric model for correlation matrix
## considering the CPC parametrization
c0model <- cgeneric(
    model = "pc_correl", n = p,
    lambda = 3,
    useINLAprecomp = FALSE)
c0model

c1model <- cgeneric(
    model = "pc_correl", n = p,
    base = c1, lambda = 3,
    useINLAprecomp = FALSE)
c1model

dataf0 <- data.frame(
    i = 1:p,
    y = NA
)

library(INLA)

fh <- list(prec = list(initial = 20, fixed = TRUE))

fit0 <- inla(
    y ~ 0 + f(i, model = c0model),
    data = dataf0,
    control.family = list(hyper = fh)
)
fit1 <- inla(
    y ~ 0 + f(i, model = c1model),
    data = dataf0,
    control.family = list(hyper = fh)
)

th0 <- fit0$mode$theta
th1 <- fit1$mode$theta

cholcor(th1)
tcrossprod(cholcor(th1))

basecor(base = th0, p = p)
basecor(base = th1, p = p)

## sample theta from its joint posterior
th0post <- inla.hyperpar.sample(
    n = 3000, result = fit0)
summary(th0post)
th1post <- inla.hyperpar.sample(
    n = 3000, result = fit1)
summary(th1post)

## lower.tri index
iil <- which(lower.tri(c1))
iil

## transform from theta to correlations
cor0post <- t(sapply(1:nrow(th0post), function(i)
    tcrossprod(cholcor(th0post[i, ]))[iil]))
summary(cor0post)
cor1post <- t(sapply(1:nrow(th1post), function(i)
    tcrossprod(cholcor(th1post[i, ]))[iil]))
summary(cor1post)

## labels for plots
thlabs <- lapply(1:3, function(i)
    as.expression(bquote(theta[.(i)])))
clabs <- list(expression(rho[1~","~2]),
              expression(rho[1~","~3]),
              expression(rho[2~","~3]))

## visualize marginals for theta (lower) and correlations (upper)
k1 <- k2 <- 0
par(mfrow = c(3, 3), mar = c(3,3,0.5,0.5),
    mgp = c(2,0.5,0), bty = "n")
for(i in 1:p) {
    for(j in 1:p) {
        if (i==j) {
            plot(0, type = "n", xlab = "",
                 ylab = "", axes = FALSE)
        }
        if(i<j) {
            k1 <- k1 + 1
            h0 <- hist(cor0post[,k1], 100, plot = FALSE)
            h1 <- hist(cor1post[,k1], 100, plot = FALSE)
            plot(h0, xlim = range(h0$breaks, h1$breaks), 
                 col = gray(0.5), border = "transparent",
                 main = "", xlab = clabs[[k1]], freq = FALSE)
            plot(h1, add = TRUE, freq = FALSE,
                 col = rgb(1,0.5,0,0.5), border = "transparent")
            abline(v = c(0, c1[iil[k1]]), lty = 2, col = c(1, 2))
        }
        if(i>j) {
            k2 <- k2 + 1
            h0 <- hist(th0post[,k2], 100, plot = FALSE)
            h1 <- hist(th1post[,k2], 100, plot = FALSE)
            plot(h0, xlim = range(h0$breaks, h1$breaks), 
                 col = gray(0.5), border = "transparent",
                 main = "", xlab = thlabs[[k2]], freq = FALSE)
            plot(h1, add = TRUE, freq = FALSE,
                 col = rgb(1,0.5,0,0.5), border = "transparent")
            lines(fit0$internal.marginals.hyperpar[[k2]])
            lines(fit1$internal.marginals.hyperpar[[k2]], col = 2)
        } 
    }
}
