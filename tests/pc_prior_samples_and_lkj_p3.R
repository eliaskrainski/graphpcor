
### case specific comparison of pc-correlation prior and LKJ

library(graphpcor)

## lambda values 
lambdas <- c(.2, 1, 5)
names(lambdas) <- paste0("lambda=", sapply(lambdas, format, digits = 2))
etas <- lambdas * 5
names(etas) <- paste0("eta=",etas)
nsim <- 5000

### define a base model with base = I
p <- 3
iil <- which(lower.tri(diag(p)))
ijl <- lapply(iil, function(il)
    c(row(diag(p))[il], col(diag(p))[il]))
npars <- length(iil)

c0base <- basecor(diag(p))

round(hessian(c0base), 4)

sample.basecor(c0base, size = 2, lambda = 1)

prior_c0_samples <- lapply(lambdas, function(l)
    sample.basecor(c0base, nsim, l))

## case specific correlation matrix
c3 <- matrix(c(1.0, 0.9, 0,
               0.9, 1.0,-0.4,
               0.0, -0.4, 1.0), 3)

c3
all.equal(c3,t(c3))
chol(c3)

## the base
cbase <- basecor(c3)
cbase

diag(hessian(cbase))

sample.basecor(cbase, size = 2, lambda = 1)
sample.basecor(cbase, size = 2, lambda = 100)

prior_corr_samples <- lapply(lambdas, function(l)
    sample.basecor(cbase, nsim, l))

cols <- c(1:6)
ltys <- rep(1, 6)
lwds <- ltys
ym <- 2.7

if(FALSE)
    png("lkj_pc_correl_priors_p3.png", 5000, 3000, pointsize = 6, res = 600)
par(mfrow = c(length(lambdas), npars), 
    mar = c(2.5, 2.5, 1.5, 0.5), mgp = c(2, 1, 0),
    bty = "n", las = 1)
for(i1 in 1:length(lambdas)) {
    for (i2 in 1:npars) {
        i <- ijl[[i2]][1]
        j <- ijl[[i2]][2]
        cij <- cbase$base[i, j]
        c0s <- prior_c0_samples[[i1]][i, j, ]
        ccs <- prior_corr_samples[[i1]][i, j, ]
        d0 <- density(c0s, bw = 0.1,
                      from = -1, to = 1)
        d1 <- density(ccs, bw = 0.1,
                      from = -1, to = 1)
        ymax <- max(d0$y, d1$y, ym, dbeta(0.5, etas[i1], etas[i1])/2)
        plot(function(x) dbeta((1+x)/2, etas[i1], etas[i1])/2, -1, 1, n=1001,
             ylim = c(0, ymax), xlab = '', ylab = '',
             main = bquote(rho[.(i)~","~.(j)]))
        lines(d0, col = 2)
        lines(d1, col = 4)
##        h0 <- hist(c0s, 100, plot = FALSE)
  ##      h1 <- hist(ccs, 100, plot = FALSE)
    ##    plot(h0, add = TRUE, freq = FALSE, border = 'transparent', col = rgb(1,0,0,.5))
      ##  plot(h1, add = TRUE, freq = FALSE, border = 'transparent', col = rgb(0,0,1,.5))
        rug(x = cij, 0.1, lwd = 3)
        if(i2==1) {
            lg <- list(bquote("LKJ ("~eta ~ "=" ~.(format(etas[i1], digits = 2)) ~ ")"),
                       bquote("PC_I ("~lambda~"="~.(lambdas[i1]) ~ ")"),
                       bquote("PC_C ("~lambda~"="~.(lambdas[i1]) ~ ")"))
            legend(-0.95-(i1>1)/10, ymax*.9, ##"topright",
                   bty = "n", seg.len = 1, 
                   as.expression(lg), lty = 1, col = c(1,2,4))            
        }
    }
}
if(FALSE) {
    dev.off()
    system("eog lkj_pc_correl_priors_p3.png &")
}
