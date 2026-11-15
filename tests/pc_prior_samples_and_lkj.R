### show the prior for different lambda, compare with LKJ

library(graphpcor)

c4 <- matrix(c( 1.0, 0.7, 0.5,-0.2,
                0.7, 1.0,-0.0, 0.4,
                0.5,-0.0, 1.0,-0.9,
               -0.2, 0.4,-0.9, 1.0), 4)

c4
all.equal(c4,t(c4))
chol(c4)

cbase <- basecor(c4)
cbase

diag(hessian(cbase))

sample.basecor(cbase, size = 2, lambda = 1)
sample.basecor(cbase, size = 2, lambda = 100)

## index for the the correlations
p <- ncol(c4)
iil <- which(lower.tri(diag(p)))
ijl <- lapply(iil, function(il)
    c(row(diag(p))[il],
      col(diag(p))[il]))
ijl

lambdas <- c(.01, 0.1, 0.3, 1)
etas <- lambdas * 10
names(lambdas) <- paste0("lambda=", lambdas)
names(etas) <- paste0("eta=", etas)
nsim <- 1000

### define a base model with base = I
c0base <- basecor(diag(p))

prior_c0_samples <- lapply(lambdas, function(l)
    sample.basecor(c0base, nsim, l))
prior_corr_samples <- lapply(lambdas, function(l)
    sample.basecor(cbase, nsim, l))

str(prior_corr_samples)

cols <- c(1:6)
ltys <- rep(1, 6)
lwds <- ltys

par(mfrow = c(length(lambdas), length(iil)),
    mar = c(2,2,1.5,0.5),
    mgp = c(2, 1, 0.5), bty = "n")
for(i1 in 1:length(lambdas)) {
    for (i2 in 1:length(iil)) {
        i <- ijl[[i2]][1]
        j <- ijl[[i2]][2]
        cij <- cbase$base[i, j]
        c0s <- prior_c0_samples[[i1]][i, j, ]
        ccs <- prior_corr_samples[[i1]][i, j, ]
        d0 <- density(c0s, bw = 0.1, from = -1, to = 1)
        d1 <- density(ccs, bw = 0.1, from = -1, to = 1)
        ymax <- max(d0$y, d1$y, dbeta(0.5, etas[i1], etas[i1])/2)
        plot(function(x) dbeta((1+x)/2, etas[i1], etas[i1])/2, -1, 1,
             ylim = c(0, ymax), xlab = '', ylab = '',
             main = bquote(rho[.(i)~","~.(j)]))
        lines(d0, col = 2)
        lines(d1, col = 4)        
        rug(x = cij, 0.1, lwd = 3)
    }
    lg <- list(bquote("LKJ : "~eta ~ "=" ~.(etas[i1])),
               bquote("PC_I : "~lambda~"="~.(lambdas[i1])),
               bquote("PC_C : "~lambda~"="~.(lambdas[i1])))
    legend("topright", bty = "n",
           as.expression(lg), lty = 1, col = c(1,2,4))
}
