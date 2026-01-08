library(graphpcor)
library(INLA)

################################################################
### n = 2, m = 1
n <- 2
(m <- n*(n-1)/2)

## evaluate some lambda
lambdas <- c(0.5, 1, 2, 10, 50, 1000)
nlambs <- length(lambdas)

for(i in 1:nlambs) {
    Cmodel <- cgeneric(
        model = "pc_correl",
        n = n,
        lambda = lambdas[i],
        base = 0,
        useINLAprecomp = FALSE)
    print(integrate(function(th) exp(prior(Cmodel, theta = th)), -3, 3))
}


cfam <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))
fi <- y ~ 0 + f(i, model = Cmodel)
d1n <- data.frame(y = NA, i = 1:n)

hth <- 0.005
ths <- c(seq(-1.5+hth/2, 1.5-hth/2, hth))
(nths <- length(ths))

par(mfrow = c(2, 3), mar = c(3,3,0.1,0.1), mgp = c(1.5, 0.5, 0), bty = "n")
for(i in 1:nlambs) {
    Cmodel <- cgeneric(
        model = "pc_correl",
        n = n,
        lambda = lambdas[i],
        useINLAprecomp = FALSE)
    ifit <- inla(
        formula = fi, data = d1n, control.family = cfam
    )
    pm <- ifit$marginals.hyperpar[[1]]
    plot(ths, exp(prior(Cmodel, theta = matrix(ths, 1))),
         xlim = range(pm[, 1]) * 3, type = "l")
    lines(pm, col = 2, lty = 3)
}

basecor(-1, 2)
solve(prec(Cmodel, theta = -1))

################################################################
## n = 3, m = 3
c0 <- matrix(c(1,     .8,  -.5,
               0.8,    1,  -.4,
               -.5,  -.4,   1), 3)
c0

(n <- ncol(c0))

b0 <- basecor(c0)
(theta3 <- b0$theta)
(m <- n*(n-1)/2)

lR <- cholcor(theta3, n)
lR

R <- tcrossprod(lR)
R

cmodel <- cgeneric(
    model = "pc_correl", 
    base = c0, 
    lambda = 5,
    useINLAprecomp = FALSE)

cmodel

graph(cmodel)

initial(cmodel)

prec(cmodel, theta = theta3)

all.equal(as.matrix(solve(prec(cmodel, theta = theta3))),
          R)

prior(cmodel, theta = theta3)

hth <- 0.02; hth2 <- hth^2
str(sths <- lapply(1:3, function(i)
    seq(-2+hth/2, 2-hth/2, hth) + theta3[i]))
(nsth <- length(sths[[1]]))
ncol(th3 <- t(expand.grid(sths)))/1e6

p3 <- array(prior(cmodel, theta = th3), rep(nsth, 3))
str(p3)

par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5),
    mgp = c(2, 0.5, 0), bty = "n")
image(sths[[1]], sths[[2]], p3[,,nsth/2])
contour(sths[[1]], sths[[2]], p3[,,nsth/2], add = TRUE, nlevels = 5)
image(sths[[1]], sths[[3]], p3[,nsth/2,])
contour(sths[[1]], sths[[3]], p3[,nsth/2,], add = TRUE, nlevels = 5)
image(sths[[2]], sths[[3]], p3[nsth/2,,])
contour(sths[[2]], sths[[3]], p3[nsth/2,,], add = TRUE, nlevels = 5)

exp(-8:-6)

sum(exp(p3)*(hth^3))

hth2 <- hth * hth

array(2:30, c(3,5,2))
apply(array(2:30, c(3,5,2)), 3, sum)

sum(apply(exp(p3)*hth2, 2, sum) * hth)
sum(apply(exp(p3)*hth2, 2, sum) * hth)
sum(apply(exp(p3)*hth2, 3, sum) * hth)

plot(sths[[1]], apply(exp(p3)*hth2, 1, sum),
     type = "l", xlab = "", ylab = "density", 
     xlim = range(unlist(sths)))
lines(sths[[2]], apply(exp(p3)*hth2, 2, sum), col = 2)
lines(sths[[3]], apply(exp(p3)*hth2, 3, sum), col = 3)
legend("topleft", bty = "n", lty = 1, col = 1:3,
       c(expression(theta[1]), expression(theta[2]), expression(theta[3])))

## consider different prior parameters
## fit with no data: get back the prior

pthlabs <- c(
    expression(pi[theta[1]](theta[1]~"|"~lambda)),
    expression(pi[theta[2]](theta[2]~"|"~lambda)),
    expression(pi[theta[3]](theta[3]~"|"~lambda))
)
rholabs <- c(
    expression(rho[2~","~1]),
    expression(rho[3~""~1]),
    expression(rho[3~""~2])
    )
prholabs <- c(
    expression(pi[rho[2~","~1]](rho[2~","~1]~"|"~lambda)),
    expression(pi[rho[3~","~1]](rho[3~","~1]~"|"~lambda)),
    expression(pi[rho[3~","~2]](rho[3~","~2]~"|"~lambda))
)

il3 <- which(lower.tri(diag(n)))
il3

par(mfrow = c(nlambs, 6), mar = c(4,4,0,0), mgp = c(2,0.5,0))
for(il in 1:nlambs) {
    lamb <- lambdas[il]
    Cmodel <- cgeneric(
        model = "pc_correl", n = n, lambda = lamb,
        base = theta3, useINLAprecomp = FALSE)
    itest <- inla(
        formula = fi, control.family = cfam,
        data = list(y = rep(NA, n), i = 1:n)
    )
    for(i in 1:3) {
        plot(itest$marginals.hyperpar[[i]], pch = 19, type = 'b',
             xlab = as.expression(bquote(theta[.(i)])),
             ylab = pthlabs[[i]], main = "")
        if(i==1)
            legend("topleft", bty = "n",
                   as.expression(bquote(lamda == .(lamb))))
    }
    hprior3 <- inla.hyperpar.sample(10000, itest)
    cprior3 <- apply(hprior3, 1, function(x)
        tcrossprod(cholcor(x, n))[il3])
    for(i in 1:3) {
        hist(cprior3[i, ], seq(-.1, .1, 0.01)+R[il3[i]], freq = FALSE,
             main = "",
             xlab = rholabs[[i]],
             ylab = prholabs[[i]])
        abline(v = R[il3[i]], col = 2, lty = 2, lwd = 2)
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
par(mfrow = c(nlambs, 6), mar = c(4,4,0,0), mgp = c(2,0.5,0))
for(lamb in lambdas) {
    cmodel <- cgeneric("pc_correl", n = n, lambda = lamb,
                       useINLAprecomp = FALSE)
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
        xlm <- range(theta3[i],
                     Bs[[1]]$theta[i],
                     Bs[[2]]$theta[i],
                     pm1[, 1], pm2[, 1])
        ylm <- c(0, max(pm1[, 2], pm2[, 2]))
        plot(pm1, type = "l", col = 2, lwd = 2,
             xlim = xlm, ylim = ylm, 
             xlab = as.expression(bquote(theta[.(i)])),
             ylab = pthlabs[[i]], main = "")
        lines(pm2, col = 3, lwd = 2)
        lines(sths[[i]], apply(p3a*hth2, i, sum), lty = 3, lwd = 3)
        abline(v = theta3[i], col = 1, lty = 1)
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
                   title = as.expression(bquote(lambda == .(lamb))))
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
        abline(v = R[il3[i]], lty = 1, col = 1)
        rug(R[il3[i]], col = 1, lwd = 1)
        rug(Rs[[1]][il3[i]], col = 2, lwd = 5, lty = 2)
        rug(Rs[[2]][il3[i]], col = 3, lwd = 5, lty = 2)
        if(i==1)
            legend("topleft", bty = "n",
                   paste("n = ", nsims), 
                   fill = c(2,3), border = 2:3,
                   title = as.expression(bquote(lambda == .(lamb))))
    }
}


################################################################
## n = 4
n <- 4
(m <- n*(n-1)/2)

(theta0 <- rnorm(m))

lR <- cholcor(theta0, p = n)
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
    expression(pi[rho[2~","~1]](rho[2~","~1]~"|"~lambda)),
    expression(pi[rho[3~","~1]](rho[3~","~1]~"|"~lambda)),
    expression(pi[rho[4~","~1]](rho[4~","~1]~"|"~lambda)),
    expression(pi[rho[3~","~2]](rho[3~","~2]~"|"~lambda)),
    expression(pi[rho[4~","~2]](rho[4~","~2]~"|"~lambda)),
    expression(pi[rho[4~","~3]](rho[4~","~3]~"|"~lambda))
)

## consider different priors for these 2 data
## visualize the posterior for each correlation
par(mfrow = c(nlambs, m), mar = c(3,3,0,0), mgp = c(1.5,0.5,0))
for(il in 1:nlambs) {
    lamb <- lambdas[il]
    cmodel <- cgeneric("pc_correl", n = n, lambda = lamb, useINLAprecomp = FALSE)
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
        abline(v = R[il4[i]], lty = 1, col = 1)
        rug(R[il4[i]], col = 1, lwd = 1)
        rug(Rs[[1]][il4[i]], col = 2, lwd = 5, lty = 2)
        rug(Rs[[2]][il4[i]], col = 3, lwd = 5, lty = 2)
        if(i==1)
            legend("topleft", bty = "n",
                   paste("n = ", nsims), 
                   fill = c(2,3), border = 2:3,
                   title = as.expression(bquote(lambda == .(lamb))))
    }
}


detach("package:graphpcor", unload = TRUE)
library(graphpcor)
