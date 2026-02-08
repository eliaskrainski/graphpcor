library(INLA)
library(graphpcor)

### for concepts: test/concepts/concepts_graphpcor.R
### for details: tests/check/check_graphpcor.R
### for details: tests/detailed/detailed_graphpcor.R

## the graph in Example 2.6 of the GMRF book
g <- graphpcor(x1 ~ x2+x3, x2~x4, x3~x4)
class(g)

g

summary(g)

(ne <- dim(g))

plot(g)

(Lg <- Laplacian(g))

## define the cgeneric model (see test/concepts/graphpcor.R for other options)
theta.base <- rnorm(ne[2])
s0 <- 3:0+0.5

cmodel <- cgeneric(
    model = g,
    lambda = 10,
    base = theta.base,
    useINLAprecomp = FALSE, 
    sigma.prior.reference = s0, 
    sigma.prior.probability = rep(0.01, ne[1]))

cmodel

c(log(s0), theta.base)

mu(cmodel)

initial(cmodel)

sigmas <- c(5, 1, 0.5, 0.1)
thetaL <- c(-5, 1, 2, -0.1)
theta1 <- c(log(sigmas), thetaL)

Vg <- vcov(g, theta = theta1)
Vg

cov2cor(vcov(g, theta = theta.base)) ## base model correlation
cov2cor(Vg) ## correlation to be used to sample data

Q1 <- prec(cmodel, theta = theta1)
Q1

all.equal(Vg, solve(as.matrix(Q1)), check.attributes = FALSE)

## some data
nrep <- 100
nd <- nrep * ne[1]

xx <- matrix(rnorm(nd), nrep) %*% chol(Vg)
(Vxx <- cov(xx))

theta.y <- 10
datar <- data.frame(
    r = rep(1:nrep, ne[1]),
    i = rep(1:ne[1], each = nrep),
    y = rnorm(nd, 1 + xx, exp(-2*theta.y))
)

m1 <- y ~ f(i, model = cmodel, replicate = r)
fit <- inla(
    formula = m1,
    data = datar,
    control.family = list(hyper = list(prec = list(initial = theta.y, fixed = TRUE)))
)
fit$cpu.used

grep("function evaluations =", fit$logfile, value = TRUE)
grep("fn-calls=", fit$logfile, value = TRUE)

rbind(true = c(theta1),
      fit = fit$mode$theta)
fit$mode$theta - c(theta1)

prec(cmodel, theta = fit$mode$theta)

round(Vg, 2)

round(vcov(g, theta = thetaL), 2)
round(cor(xx), 2)
round(Vfit <- vcov(g, theta = fit$mode$theta[5:8]), 2)

ptheta.samples <- t(inla.hyperpar.sample(
    n = 10000, result = fit, intern = TRUE))

il4 <- which(lower.tri(diag(4)))

system.time(
    csamples <- apply(ptheta.samples, 2, function(th) {
        v <- vcov(g, theta = th)
        c(sqrt(diag(v)), cov2cor(v)[il4])
    }))

il4
ii4 <- c(1, 5:7, 2, 8:9, 3, 10, 4)

true.p <- c(sqrt(diag(Vg)), cov2cor(Vg)[il4])
xx.p <- c(sqrt(diag(Vxx)), cov2cor(Vxx)[il4])

par(mfrow = c(4, 4), mar = c(3,3,0.3,0.3), mgp = c(2, 0.5, 0), las = 1)
k <- 0; k.th <- 0
for(i in 1:4) {
    for(j in 1:4) {
        if(i>j) {
            if(is.zero(Lg[i,j])) {
                plot(0, type = "n", xlab = "", ylab = "", bty = "n", axes = FALSE)
            } else {
                k.th <- k.th + 1
                plot(fit$internal.marginals.hyperpar[[4+k.th]],
                     bty = "n", type = "l", main = '', ylab = "Density",
                     xlab = as.expression(bquote(theta[.(k.th)])))
                abline(v = thetaL[k.th], lty = 1, lwd = 3)
            }
            if((i==4) & (j==1)) {
                legend("bottomleft",
                       c("true", "sample", "posterior", "post. sample"),
                       lwd = c(3,3,1,0), lty = c(1,3,1,0), bty = "n",
                       fill = rep(c("transparent", gray(0.4)), c(3,1)), 
                       border = rep(c("transparent", "black"), c(3,1)))
            }
        } else {
            k <- k + 1
            if(i==j) {
                hist(csamples[i, ], 30, freq = FALSE, main = '',
                     xlab = as.expression(bquote(sigma[.(i)])),
                     col = gray(0.3), border = 'transparent')
                abline(v = c(true.p[i], xx.p[i]), lty = c(1,3), lwd = 3)
            }
            if(i<j) {
                hist(csamples[ii4[k], ], 30, freq = FALSE, main = '',
                     xlab = as.expression(bquote(rho[.(i)~.(j)])),
                     col = gray(0.5), border = 'transparent')
                abline(v = c(true.p[ii4[k]], xx.p[ii4[k]]), lty = c(1,3), lwd = 3)
            }
        } 
    }
}
