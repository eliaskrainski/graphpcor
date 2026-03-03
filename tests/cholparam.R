

rcgeneric <- function(n, cg) {
    p <- length(initial(cg))
    m <- 0; s = 1
    ns <- 0
    out <- matrix(0, n, p)
    repeat{
        x0 <- matrix(rnorm(1000*p, m, s), nrow = p)
        px0 <- exp(prior(cg, theta = x0))
        jj0 <- which(runif(1000) < (px0 / dnorm(0,0,1)))
        if(length(jj0)>0) {
            jj <- head(jj0, n-ns)
            out[ns + 1:length(jj), ] <- t(x0[, jj, drop = FALSE])
            ns <- ns + length(jj)
            if(ns>=n) break            
        } else {
            if(ns>0) {
                m <- colMeans(out[1:ns, , drop = FALSE])
                s <- sd(out)
            }
        }
    }
    return(out)
}


library(graphpcor)

p <- 5
gs <- graphpcor(paste0("X1~", paste0("X",2:p)))
gs

Ls <- Laplacian(gs)
Ls

chol(Ls + diag(p))

oj <- c(2:p,1)
ogs <- graphpcor(Ls[oj, oj])
ogs
Laplacian(ogs)
chol(Laplacian(ogs) + diag(p))

par(mfrow=c(1,2))
plot(gs)
plot(ogs)

## A correlation matrix
thp <- c(-1,1,-3,0.3)

bp <- basepcor(thp, p, iLtheta = gs)
bp

bp$base[oj, oj]

obp <- basepcor(bp$base[oj, oj], iLtheta = ogs)
obp

round(chol(bp$base), 4)
round(chol(obp$base), 4)

hessian(bp)
hessian(obp)

## simulate data
##thd <- rnorm(dim(gs)[2], -1)
##cc <- basepcor(thd, p, gs)
##cc$base
##Lc <- chol(cc$base)

##Lc <- chol(bp$base)

lseq <- c(0.1, 1, 10); names(lseq) <- paste0("l", lseq)
nl <- length(lseq)

mm <- list(
    m1 = lapply(lseq, function(l)
        cgeneric(gs, base = bp, lambda = l, useINLAprecomp = FALSE)),
    m2 = lapply(lseq, function(l)
        cgeneric(ogs, base = obp, lambda = l, useINLAprecomp = FALSE)
        )
)

rr <- rcgeneric(10000, mm$m1[[1]])

par(mfrow = c(2, 2))
for(j in 1:4) {
    hist(rr[, j], freq = FALSE)
    abline(v = mm[[1]][[1]]$f$cgeneric$data$doubles$thetabase[j],
           col = 2, lwd = 2, lty = 2)
}

mm[[1]][[1]]$f$cgeneric$data$ints
mm[[1]][[1]]$f$cgeneric$data$doubles
mm[[1]][[1]]$f$cgeneric$data$doubles$thetabase

par(mfrow = c(2,4))
plot(function(x) prior(mm[[1]][[1]], theta = rbind(x,0,0,0)),
     -10, 10, n = 1001, ylim = c(-100, 0))
for(j in 2:nl)
    plot(function(x) prior(mm[[1]][[j]], theta = rbind(x,0,0,0)),
         -10, 10, n = 1001, add = TRUE, col = j)
plot(function(x) prior(mm[[1]][[1]], theta = rbind(0,x,0,0)),
     -10, 10, n = 1001)
plot(function(x) prior(mm[[1]][[1]], theta = rbind(0,0,x,0)),
     -10, 10, n = 1001)
plot(function(x) prior(mm[[1]][[1]], theta = rbind(0,0,0,x)),
     -10, 10, n = 1001)
plot(function(x) prior(mm[[2]][[1]], theta = rbind(x,0,0,0)),
     -10, 10, n = 1001)
plot(function(x) prior(mm[[2]][[1]], theta = rbind(0, x,0,0)),
     -10, 10, n = 1001)
plot(function(x) prior(mm[[2]][[1]], theta = rbind(0,0,x,0)),
     -10, 10, n = 1001)
plot(function(x) prior(mm[[2]][[1]], theta = rbind(0,0,0,x)),
     -10, 10, n = 1001)

fits[[1]][[1]]$mode$theta

mm[[1]][[1]]$f$cgeneric$data$doubles$lambda
sapply(mm, sapply, function(x) x$f$cgeneric$data$doubles$lambda)

idat <- data.frame(
    i = 1:p, y = NA)

#    i = rep(1:p, each = n),
 #   r = rep(1:n, p),
  #  y1 = as.vector(xx),
   # y2 = as.vector(xx[, oj])

cfam <- list(hyper = list(prec = list(initial = 20, fixed = TRUE)))

library(INLA)

fits <- lapply(mm, lapply, function(m)
    inla(formula = y ~ 0 + f(i, model = m),
         data = idat, control.family = cfam,
         control.mode = list(
             theta = rep(-1, dim(gs)[2]),
             restart = TRUE)
         )
    )

lapply(fits, sapply, function(x) x$mode$theta)

sapply(fits, sapply, function(x) mean(x$summary.hy$sd))

sapply(fits, sapply, function(r) r$misc$nfunc)
sapply(fits, sapply, function(r) r$cpu.used[["Total"]])

sapply(fits, sapply, function(r) r$misc$nfunc)/
sapply(fits, sapply, function(r) r$cpu.used[["Total"]])

## posterior hyperpar samples
ill <- which(lower.tri(Ls))
cc1 <- lapply(fits[[1]], function(r) {
    hh <- inla.hyperpar.sample(5000, r, TRUE)
    t(sapply(1:nrow(hh), function(i)
        vcov(gs, theta = hh[i, ])[ill]))
})
cc2 <- lapply(fits[[2]], function(r) {
    hh <- inla.hyperpar.sample(5000, r, TRUE)
    t(sapply(1:nrow(hh), function(i)
        vcov(ogs, theta = hh[i, ])[ill]))
})


i2 <- c(4,7,9,10, 1,2,3, 5,6, 8)

cols <- c(gray(0.1,0.5), rgb(.3,.5,1,.7))

k1 <- 0
par(mfrow = c(p,p), mar = c(3,3,0.5,0.5),
    mgp = c(1.5,0.5,0), bty = "n")
for(i in 1:p) {
    for(j in 1:p) {
        if(i<j) {
            k1 <- k1 + 1
            ld1 <- lapply(cc1, function(cc)
                density(cc[, k1])) ##, 100, plot = FALSE))
            ld2 <- lapply(cc2, function(cc)
                density(cc[, i2[k1]]))##, 100, plot = FALSE))
            xlm <- range(unlist(lapply(ld1, function(x) x$x)), 
                         unlist(lapply(ld2, function(x) x$x))
                         ) + c(-.1,.1)
            ylm <- range(unlist(lapply(ld1, function(x) x$y)),
                         unlist(lapply(ld2, function(x) x$y)))
            xlm[xlm<(-1)] <- -1
            xlm[xlm>1] <- 1
            plot(ld1[[1]], xlim = xlm, ylim = ylm, freq = FALSE,
                 col = cols[1], ##border = 'transparent',
                 main = '', ylab = "Density",
                 xlab = bquote(rho[.(i)~","~.(j)]))
            for(i in 1:length(lseq)) {
                lines(ld1[[i]], col = 1, lty =i)
                lines(ld2[[i]], col = 2, lty =i)
            }            
#            plot(h2, add = TRUE, freq = FALSE,
 #                col = cols[2], border = 'transparent')
            abline(v = bp$base[ill[k1]], lty = 2, lwd = 2)
        } else {
            #if(i==j) {
                plot(0, type = 'n', xlab = '', ylab = '', axes = FALSE)
            #} else {                
            #}
        }
    }
}

if(FALSE) {
    
        if(i<j) {
            
        }
        if(j==1) {
            if(i>1) {
                plot(fit1$marginals.hyperpar[[i-1]], type = "l")
            }
        }
        if(i==p) {
            if(j==1) {
                lines(fit2$marginal.hyperpar[[j]], lty = 2, lwd = 2)
            } else {
                if(j<p)
                    plot(fit2$marginals.hyperpar[[j]], type = "l")
            }
        }

}

