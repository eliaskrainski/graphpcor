library(INLA)
library(graphpcor)

################################################################
### n = 2, m = 1
n <- 2
(m <- n*(n-1)/2)

(theta1 <- rnorm(m))

lR <- cholcor(theta1, 2)
log(attr(lR, 'determinant'))
R <- tcrossprod(lR)
R
basecor(theta1, p = 2)

dLKJ(R, 1, log=TRUE)
dLKJ(R, 10, log=TRUE)

seq(-.9, .9, length = 19)
head(seq(-.99, .99, length = 199))
2e2-1

for(e in c(0.1, 0.2, 0.5, 1, 2, 10, 50))
    print(integrate(function(x) sapply(x, function(xx)
        dLKJ(matrix(c(1,xx,xx,1), 2), e)),
        -.999, .999, subdivisions = 2e3-1))

etas <- c(0.1, 0.5, 1, 3, 30)

par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0),
    las = 1, bty = "n")
plot(function(x) sapply(x, function(xx)
    dLKJ(matrix(c(1,xx,xx,1), 2), max(etas))),
    -.99, .99, n = 199, type = "n",
    xlab = expression(rho), ylab = 'p(R|eta)')
for(i in 1:length(etas)) {
    plot(function(x) sapply(x, function(xx)
        dLKJ(matrix(c(1,xx,xx,1), 2), etas[i])),
        -.99, .99, n = 199,
        add = TRUE, col = i+1, lwd = 2)
    plot(function(x)
        dbeta((x+1)/2, etas[i] + (n-2)/2, etas[i] + (n-2)/2)/2,
        -1, 1, n = 2001, 
        add = TRUE, col = 1, lty = 2)
}

jacobian(function(x) tcrossprod(cholcor(x))[2], theta1)

ptheta <- function(th, eta, p, iL = which(lower.tri(diag(p)))) {
    L <- cholcor(th, p, itheta = iL)
    R <- tcrossprod(L)
    J <- jacobian(function(x) tcrossprod(cholcor(x))[iL], th)
    exp(dLKJ(R, eta, log = TRUE) +
        new("numeric", determinant(J)$modulus))
}

for(e in etas)
    print(integrate(function(x) sapply(x, function(xx)
        ptheta(xx, e, 2)), -5, 5, subdivisions = 2e3+1))

hth <- 0.05
length(sth <- seq(-4, 4, hth))
(nsth <- length(sth))

par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0),
    las = 1, bty = "n")
plot(function(x) sapply(x, function(xx)
    ptheta(xx, max(etas), 2)),
    -3, 3, n = 601, type = "n",
    xlab = expression(theta), ylab = 'p(theta|eta)')
for(i in 1:length(etas)) {
    cm <- cgeneric("LKJ", n = n, eta = etas[i], useINLAprecomp = FALSE)
    plot(function(x) sapply(x, function(xx)
        ptheta(xx, etas[i], 2)), -3, 3, n = 601, 
        add = TRUE, lwd = 3, col = i+1)
    lines(sth, exp(prior(cm, theta = matrix(sth,1))), lty = 2)
}

################################################################
## n = 3, m = 3
n <- 3
(m <- n*(n-1)/2)

(theta1 <- rnorm(m))

lR <- cholcor(theta1, n)
lR
log(attr(lR, 'determinant'))
R <- tcrossprod(lR)
R
basecor(theta1, p = n)

### evaluate  fore one parameter
eta <- 1 + rgamma(1, 10, 2)
eta

cmodel <- cgeneric(
    model = "LKJ", 
    n = n,
    eta = eta,
    useINLAprecomp = FALSE)

##str(cmodel)

graph(cmodel)

initial(cmodel)

prec(cmodel, theta = theta1)

all.equal(as.matrix(solve(prec(cmodel, theta = theta1))),
          R)

prior(cmodel, theta = theta1)

th3 <- t(expand.grid(th1 = sth, th2 = sth, th3 = sth))
str(th3)

p3 <- array(prior(cmodel, theta = th3), rep(nsth, 3))

par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5),
    mgp = c(2, 0.5, 0), bty = "n")
image(sth, sth, p3[,,nsth/2])
contour(sth, sth, p3[,,nsth/2], add = TRUE, nlevels = 5)
image(sth, sth, p3[,nsth/2,])
contour(sth, sth, p3[,nsth/2,], add = TRUE, nlevels = 5)
image(sth, sth, p3[nsth/2,,])
contour(sth, sth, p3[nsth/2,,], add = TRUE, nlevels = 5)

sum(exp(p3)*(hth^3))

hth2 <- hth * hth

array(2:30, c(3,5,2))
apply(array(2:30, c(3,5,2)), 3, sum)

sum(apply(exp(p3)*hth2, 2, sum) * hth)
sum(apply(exp(p3)*hth2, 2, sum) * hth)
sum(apply(exp(p3)*hth2, 3, sum) * hth)

plot(sth, apply(exp(p3)*hth2, 1, sum),
     xlab = expression(theta), ylab = "density",
     type = "l", lwd = 2, lty = 2)
lines(sth, apply(exp(p3)*hth2, 2, sum), col = 2)
lines(sth, apply(exp(p3)*hth2, 3, sum), col = 3)

## eval 3 different parameters
## fit with no data: get back the prior

pthlabs <- c(
    expression(pi[theta[1]](theta[1]~"|"~eta)),
    expression(pi[theta[2]](theta[2]~"|"~eta)),
    expression(pi[theta[3]](theta[3]~"|"~eta))
)
rholabs <- c(
    expression(rho[2~","~1]),
    expression(rho[3~""~1]),
    expression(rho[3~""~2])
    )
prholabs <- c(
    expression(pi[rho[2~","~1]](rho[2~","~1]~"|"~eta)),
    expression(pi[rho[3~","~1]](rho[3~","~1]~"|"~eta)),
    expression(pi[rho[3~","~2]](rho[3~","~2]~"|"~eta))
)

par(mfrow = c(3, 6), mar = c(4,4,0,0), mgp = c(2,0.5,0))
for(e in c(1/2, 3, 30)) {
    cm <- cgeneric("LKJ", n = n, eta = e)
    p3a <- array(exp(prior(cm, theta = th3)), rep(nsth, 3))
    itest <- inla(
        y ~ 0 + f(i, model = cm),
        data = list(y = rep(NA, n), i = 1:n),
        control.family = list(
            hyper = list(prec = list(initial = 10, fixed = TRUE))
        )
    )
    for(i in 1:3) {
        plot(itest$marginals.hyperpar[[i]], pch = 19, type = 'b',
             xlab = as.expression(bquote(theta[.(i)])),
             ylab = pthlabs[[i]], main = "")
        lines(sth, apply(p3a*hth2, i, sum), col = i, lwd = 3)
        if(i==1)
            legend("topleft", bty = "n",
                   as.expression(bquote(eta == .(e))))
    }
    hprior3 <- inla.hyperpar.sample(10000, itest)
    il3 <- which(lower.tri(diag(n)))
    cprior3 <- apply(hprior3, 1, function(x)
        tcrossprod(cholcor(x, n))[il3])
    for(i in 1:3) {
        hist(cprior3[i, ], seq(-1, 1, 0.05), freq = FALSE,
             main = "",
             xlab = rholabs[[i]],
             ylab = prholabs[[i]])
        plot(function(x)
            dbeta((x+1)/2, e + (n-2)/2, e + (n-2)/2)/2,
             -1, 1, n = 2001, 
             add = TRUE, col = i, lwd = 3)
    }
}

## draw samples, for two different sample size
nsims <- c(50, 500)
xxs <- lapply(nsims, function(ns) 
    tcrossprod(matrix(rnorm(ns * n), ns), lR))

R
Rs <- lapply(xxs, cor)
Rs

## base model at the sample correlation (to get observed theta)
Bs <- lapply(Rs, basecor, p = n)
(oths <- sapply(Bs, function(b) b$theta))

datf3 <- lapply(1:2, function(i)
    data.frame(
        y = as.vector(xxs[[i]]),
        i = rep(1:n, each = nsims[i]),
        r = rep(1:nsims[i], n)))

fr <- y ~ 0 + f(i, model = cmodel, replicate = r)
cfam <- list(
    hyper = list(prec = list(initial = 10, fixed = TRUE))
)

## consider different priors for these 2 data
par(mfrow = c(3, 6), mar = c(4,4,0,0), mgp = c(2,0.5,0))
for(e in c(1/2, 3, 30)) {
    cmodel <- cgeneric("LKJ", n = n, eta = e)
    p3a <- array(exp(prior(cmodel, theta = th3)), rep(nsth, 3))
    itests <- lapply(datf3, function(ddf)
        inla(formula = fr, 
             data = ddf, 
             control.family = cfam
             )
        )
    for(i in 1:3) {
        pm1 <- itests[[1]]$marginals.hyperpar[[i]]
        pm2 <- itests[[2]]$marginals.hyperpar[[i]]
        xlm <- range(theta1[i],
                     Bs[[1]]$theta[i],
                     Bs[[2]]$theta[i],
                     pm1[, 1], pm2[, 1])
        ylm <- c(0, max(pm1[, 2], pm2[, 2]))
        plot(pm1, type = "l", col = 2, lwd = 2,
             xlim = xlm, ylim = ylm, 
             xlab = as.expression(bquote(theta[.(i)])),
             ylab = pthlabs[[i]], main = "")
        lines(pm2, col = 3, lwd = 2)
        lines(sth, apply(p3a*hth2, i, sum), lty = 3, lwd = 3)
        abline(v = theta1[i], col = 1, lty = 1)
        abline(v = Bs[[1]]$theta[i], lty = 2, col = 2, lwd = 2)
        abline(v = Bs[[2]]$theta[i], lty = 2, col = 3, lwd = 2)
        if(i==1)
            legend("topleft", bty = "n",
                   c("true", "prior",
                     as.expression(
                         lapply(nsims, function(ns)
                             bquote("Obs."~theta~", n="~.(ns)))),
                     as.expression(
                         lapply(nsims, function(ns)
                             bquote("Post."~theta~", n="~.(ns))))), 
                   lty = c(1, 3, 2, 2, 1, 1),
                   col = c(1, 1, 2:3, 2:3),
                   lwd = c(1, 3, 2, 2, 2, 2),
                   title = as.expression(bquote(eta == .(e))))
    }
    hpriors3 <- lapply(itests, function(r)
                       inla.hyperpar.sample(10000, r))
    cpriors3 <- lapply(hpriors3, apply, 1, function(x)
        tcrossprod(cholcor(x, n))[il3])
    hh <- 0.01
    for(i in 1:3) {
        h1 <- hist(cpriors3[[1]][i, ], seq(-1, 1, hh), plot = FALSE)
        h2 <- hist(cpriors3[[2]][i, ], seq(-1, 1, hh), plot = FALSE)
        ih <- which((h1$counts>0) | (h2$counts>0))
        plot(h1, freq = FALSE, xlab = rholabs[[i]],
             ylab = prholabs[[i]],
             xlim = range(h1$mids[ih], h2$mids[ih]) + c(-1,1)*hh,
             ylim = range(h1$dens[ih], h2$dens[ih]),
             main = "", col = rgb(1,0,0,0.5))
        plot(h2, add = TRUE, freq = FALSE, col = rgb(0,1,0,0.5))
        plot(function(x)
            dbeta((x+1)/2, e + (n-2)/2, e + (n-2)/2)/2,
             -1, 1, n = 2001, 
            add = TRUE, col = 1, lty = 3, lwd = 3)
        abline(v = R[il3[i]], lty = 1, col = 1)
        rug(R[il3[i]], col = 1, lwd = 1)
        rug(Rs[[1]][il3[i]], col = 2, lwd = 5, lty = 2)
        rug(Rs[[2]][il3[i]], col = 3, lwd = 5, lty = 2)
        if(i==1)
            legend("topleft", bty = "n",
                   paste("n = ", nsims), 
                   fill = c(2,3), border = 2:3,
                   title = as.expression(bquote(eta == .(e))))
    }
}


################################################################
## n = 4
n <- 4
(m <- n*(n-1)/2)

(theta1 <- rnorm(m))

lR <- cholcor(theta1, p = n)
log(attr(lR, 'determinant'))
R <- tcrossprod(lR)
il4 <- which(lower.tri(R))

## simulate two datasets
xxs <- lapply(nsims, function(ns) 
    tcrossprod(matrix(rnorm(ns * n), ns), lR))

R
Rs <- lapply(xxs, cor)
Rs

## base model at the sample correlation (to get observed theta)
Bs <- lapply(Rs, basecor, p = n)
(oths <- sapply(Bs, function(b) b$theta))

datfs <- lapply(1:2, function(i)
    data.frame(
        y = as.vector(xxs[[i]]),
        i = rep(1:n, each = nsims[i]),
        r = rep(1:nsims[i], n)))

rholabs4 <- c(
    expression(rho[2~","~1]),
    expression(rho[3~""~1]),
    expression(rho[4~""~1]),
    expression(rho[3~""~2]),
    expression(rho[4~""~2]),
    expression(rho[4~""~3])
)
prholabs4 <- c(
    expression(pi[rho[2~","~1]](rho[2~","~1]~"|"~eta)),
    expression(pi[rho[3~","~1]](rho[3~","~1]~"|"~eta)),
    expression(pi[rho[4~","~1]](rho[4~","~1]~"|"~eta)),
    expression(pi[rho[3~","~2]](rho[3~","~2]~"|"~eta)),
    expression(pi[rho[4~","~2]](rho[4~","~2]~"|"~eta)),
    expression(pi[rho[4~","~3]](rho[4~","~3]~"|"~eta))
)

## consider different priors for these 2 data
## visualize the posterior for each correlation
par(mfrow = c(3, m), mar = c(3,3,0,0), mgp = c(1.5,0.5,0))
for(e in c(1/2, 3, 30)) {
    cmodel <- cgeneric("LKJ", n = n, eta = e)
    itests <- lapply(datfs, function(ddf)
        inla(formula = fr, 
             data = ddf, 
             control.family = cfam
             )
        )
    hpriors <- lapply(itests, function(r)
        inla.hyperpar.sample(10000, r))
    cpriors <- lapply(hpriors, apply, 1, function(x)
        tcrossprod(cholcor(x, n))[il4])
    hh <- 0.01
    for(i in 1:m) {
        h1 <- hist(cpriors[[1]][i, ], seq(-1, 1, hh), plot = FALSE)
        h2 <- hist(cpriors[[2]][i, ], seq(-1, 1, hh), plot = FALSE)
        ih <- which((h1$counts>0) | (h2$counts>0))
        plot(h1, freq = FALSE, xlab = rholabs4[[i]],
             ylab = prholabs4[[i]],
             xlim = range(h1$mids[ih], h2$mids[ih]) + c(-1,1)*hh,
             ylim = range(h1$dens[ih], h2$dens[ih]),
             main = "", col = rgb(1,0,0,0.5))
        plot(h2, add = TRUE, freq = FALSE, col = rgb(0,1,0,0.5))
        plot(function(x)
            dbeta((x+1)/2, e + (n-2)/2, e + (n-2)/2)/2,
             -1, 1, n = 2001, 
            add = TRUE, col = 1, lty = 3, lwd = 3)
        abline(v = R[il4[i]], lty = 1, col = 1)
        rug(R[il4[i]], col = 1, lwd = 1)
        rug(Rs[[1]][il4[i]], col = 2, lwd = 5, lty = 2)
        rug(Rs[[2]][il4[i]], col = 3, lwd = 5, lty = 2)
        if(i==1)
            legend("topleft", bty = "n",
                   paste("n = ", nsims), 
                   fill = c(2,3), border = 2:3,
                   title = as.expression(bquote(eta == .(e))))
    }
}


detach("package:graphpcor", unload = TRUE)
library(graphpcor)
