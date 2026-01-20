
library(graphpcor)

## A base correlation matrix
c0 <- matrix(c( 1.0,  0.8, -0.5,
                0.8,  1.0, -0.4,
               -0.5, -0.4,  1.0), 3)
p <- ncol(c0)

## The cgeneric model for correlation matrix
## considering the CPC parametrization
c0model <- cgeneric(
    model = "pc_correl", n = p,
    base = c0, lambda = 10,
    useINLAprecomp = FALSE)
c0model

dataf0 <- data.frame(
    i = 1:p,
    y = NA
)

library(INLA)

fh <- list(prec = list(initial = 20, fixed = TRUE))

fit <- inla(
    y ~ 0 + f(i, model = c0model),
    data = dataf0,
    control.family = list(hyper = fh)
)

th0 <- fit$mode$theta

cholcor(th0)
tcrossprod(cholcor(th0))

basecor(base = th0, p = p)

## sample theta from its joint posterior
thpost <- inla.hyperpar.sample(
    n = 3000, result = fit)
summary(thpost)

## lower.tri index
iil <- which(lower.tri(c0))
iil

## transform from theta to correlations
corpost <- t(sapply(1:nrow(thpost), function(i)
    tcrossprod(cholcor(thpost[i, ]))[iil]))
summary(corpost)

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
            hist(corpost[,k1], 100, main = '',
                 col = gray(0.5), border = "transparent",
                 xlab = clabs[[k1]], freq = FALSE)
        }
        if(i>j) {
            k2 <- k2 + 1
            hist(thpost[,k2], 100, main = '',
                 col = gray(0.5), border = "transparent",
                 xlab = thlabs[[k2]], freq = FALSE)
            lines(fit$internal.marginals.hyperpar[[k2]])
        } 
    }
}
