
library(graphpcor)
library(sf)

setwd(here::here("examples", "cancer12br4"))
getwd()

## making the base models to sample from the prior
source("cancer12prior.R")
base0 <- basepcor(diag(nc), iLtheta = cGraph)
base9 <- basepcor(rep(-7, dim(cGraph)[2]), p = nc, iLtheta = cGraph)
mean(base9$base[lower.tri(diag(nc))])

## read the dataset and map
dataf <- readRDS("rds/southBR12cancer.rds")
smap <- readRDS("rds/sbrmap.rds")

stopifnot(nrow(dataf) == nrow(smap))
stopifnot(all.equal(rownames(dataf),
                    paste(smap$code_muni)))

## number of areas selected
(n.s <- nrow(dataf))

causes_r <- causes[ireord]
jcause <- match(causes_r, colnames(dataf))
jcause

stopifnot(all(causes_r==names(dataf)[jcause]))

## build the long data for INLA
nc
ldat <- data.frame(
    cause = rep(causes_r, each = n.s),
    icause = rep(1:nc, each = n.s),
    iarea = rep(1:n.s, nc),
    ijoint = 1:(nc*n.s),
    E = as.vector(as.matrix(dataf[, jcause-nc])),
    y = as.vector(as.matrix(dataf[, jcause]))
)
ldat$i2cause <- ldat$icause
ldat$i2area <- ldat$iarea

nDeaths <- tapply(ldat$y, factor(ldat$cause, causes_r), sum)
stopifnot(all(nDeaths == colSums(dataf[,jcause])))
odeaths <- order(nDeaths)

## components of the fitted models
ffl <- list(
    IndC = y ~ f(icause, model = cstds, replicate = iarea) + cause-1,
    CICC = y ~ f(icause, model = cstds, replicate = iarea) +
        f(i2cause, model = cgrpc0, replicate = i2area) + cause -1,
    CICCi = y ~ f(icause, model = cstds, replicate = iarea) +
        f(i2cause, model = cgrpci, replicate = i2area) + cause -1,
    CICC9 = y ~ f(icause, model = cstds, replicate = iarea) +
        f(i2cause, model = cgrpc9, replicate = i2area) + cause -1,
    SStr = y ~ f(ijoint, model = replBesag) + cause - 1,
    RBym = y ~ f(icause, model = cstds, replicate = iarea) +
        f(ijoint, model = replBesag) + cause-1,
    ICSC = y ~ f(icause, model = cstds, replicate = iarea) +
        f(ijoint, model = grpc0Besag) + cause -1,
    ICSCi = y ~ f(icause, model = cstds, replicate = iarea) +
        f(ijoint, model = grpciBesag) + cause -1,
    ICSC9 = y ~ f(icause, model = cstds, replicate = iarea) +
        f(ijoint, model = grpc9Besag) + cause -1
)
length(ffl)

### file names for the fitted models
fitflk <- paste0("rds/fit_c", nc, names(ffl), "_n", n.s, ".rds")

fits <- vector("list", length(ffl))
names(fits) <- names(ffl)
for(k in 1:length(ffl)) {
    flk <- fitflk[k]
    if(file.exists(flk)) {
        cat("Loading model", names(ffl)[k], "\n")
        fits[[k]] <- readRDS(flk)
    }
}

sapply(fits, function(fit) length(fit$mode$theta))

sapply(fits, function(fit) fit$cpu.used[["Total"]])/60

nfn <- data.frame(
    fntot = sapply(fits, function(r) r$misc$nfunc))
nfn

sapply(fits, function(fit) fit$cpu.used["Total"])/nfn$fntot

tabf <- sapply(fits, function(r)
    c(DIC=sum(r$dic$local.dic),
      WAIC=sum(r$waic$local.waic),
      CPO=-sum(log(r$cpo$cpo)),
      GCPO=-sum(log(r$gcpo$cv))))

tabf

apply(tabf,1,min)

tabf-apply(tabf,1,min)

xtable::xtable(tabf-apply(tabf,1,min), digits = 1)

par(mfrow = c(2,2), mar = c(4,4,1,1), bty = "n")
plot(fits[[2]]$mode$theta, fits[[3]]$mode$theta); abline(0:1)
plot(fits[[2]]$mode$theta, fits[[4]]$mode$theta); abline(0:1)
plot(fits[[7]]$mode$theta, fits[[8]]$mode$theta); abline(0:1)
plot(fits[[7]]$mode$theta, fits[[9]]$mode$theta); abline(0:1)

cor(cbind(fits[[7]]$summary.random$ijoint$mean,
          fits[[8]]$summary.random$ijoint$mean,
          fits[[9]]$summary.random$ijoint$mean))

nc
cpars <- sapply(fits[c(2:4,7:9)], function(m)
    m$mode$theta[nc*2+1:dg[2]])
cor(cpars, baseC$theta)

lapply(fits[7:9], function(fit)
    data.frame(iid=exp(fit$mode$theta[1:nc]),
               besag=exp(fit$mode$theta[nc+1:nc]),
               frac=exp(fit$mode$theta[nc+1:nc]*2)/
                   (exp(fit$mode$theta[1:nc] * 2) +
                    exp(fit$mode$theta[nc+1:nc]*2)),
               row.names = causes))



par(mfrow = c(5,5), mar = c(3,3,0.5,0.5),
    mgp = c(2,0.5,0), bty = "n")
for(k in 1:dg[2]) {
    marg0 <- fits[[7]]$internal.marginals.hyperpar[[2*nc+k]]
    marg1 <- fits[[8]]$internal.marginals.hyperpar[[2*nc+k]]
    marg2 <- fits[[9]]$internal.marginals.hyperpar[[2*nc+k]]
    xlm <- range(marg0[,1], marg1[,1], marg2[,1], 0, baseC$theta[k])
    ylm <- range(marg0[,2], marg1[,2], marg2[,2])
    plot(marg0, type = 'l', xlim = xlm, ylim = ylm, las = 1,
         xlab = bquote(theta[.(k)]), ylab = 'Density',
         col = gray(0.5), lwd = 4)
    lines(marg1, col = 2)
    lines(marg2, col = 4)
    rug(baseC$theta[k], 0.1, lwd = 2, col = 2)
##    rug(mctheta[k], 0.1, lwd = 2, col = c(6,1)[agree[k]+1])
}

## prior samples
nsampls <- 5000
sprior.5 <- list(
    sample.basepcor(base0, nsampls, lambda = 1),
    sample.basepcor(baseC, nsampls, lambda = 1),
    sample.basepcor(base9, nsampls, lambda = 1/2)
)

## posterior samples
library(INLA)
hsampls <- lapply(fits[7:9], function(r)
    inla.hyperpar.sample(
        n = nsampls, result = r,
        intern = TRUE, improve.marginals = TRUE))

tau.sampls <- exp(hsampls[[3]][, 1:nc])
sig.sampls <- exp(hsampls[[3]][, nc+1:nc])
phi.sampls <- sig.sampls^2 / (sig.sampls^2 + tau.sampls^2)
colnames(phi.sampls) <- causes_r

micfn <- function(x) {
    c(mean(x), quantile(x, c(0.025,0.975)))
}
fn1 <- function(x) {
    paste0(x[1], " (", x[2], "; ", x[3], ")")
}
fn2 <- function(x,d=2) {
    fn1(format(micfn(x), digits = d))
}

mtaus <- t(sapply(fits[[9]]$internal.marginals.hyperpar[1:12], function(m) {
    m <- inla.tmarginal(exp, m)
    unlist(inla.zmarginal(m, TRUE))[c(1,3,7)]
}))
msigmas <- t(sapply(fits[[9]]$internal.marginals.hyperpar[12+1:12], function(m) {
    m <- inla.tmarginal(exp, m)
    unlist(inla.zmarginal(m, TRUE))[c(1,3,7)]
}))

stab <- data.frame(
    Deaths = nDeaths, 
    tau = apply(format(mtaus, digits = 2), 1, fn1),
    sigma = apply(format(msigmas, digits = 3), 1, fn1), 
    phi = apply(phi.sampls, 2, fn2, d = 3))
##    row.names = causes_r)

stab[odeaths,]

xtable::xtable(stab[odeaths, ])

th2c <- function(theta) {
    cov2cor(chol2inv(t(graphpcor:::Lprec0(theta, p=12, iLtheta=baseC$iLtheta, d0=baseC$d0))))
}

csampls <- lapply(hsampls, function(th) {
    ccs <- sapply(1:nrow(th), function(i)
        th2c(th[i, 12*2+1:24]))
    dim(ccs) <- c(12,12,nrow(th))
    ccs
})

ill <- which(lower.tri(diag(12)))

csmeans <- sapply(csampls, function(cs) {
    sapply(ill, function(ii) {
        i = row(diag(12L))[ii]
        j = col(diag(12L))[ii]
        mean(cs[i,j,])
    })
})

summary(csmeans)
table(csmeans[,2]>=csmeans[,1])

dmplot <- function(dd, xlab = '', ylab = '', main = '', xlim, ylim,
                   lty=1:length(dd), lwd=1:length(dd), col=1:length(dd)) {
    if(missing(xlim))
        xlim <- range(unlist(lapply(dd, function(d) d$x)))
    if(missing(ylim))
        ylim = range(unlist(lapply(dd, function(d) d$y)))   
    plot(dd[[1]]$x, dd[[1]]$y, xlim = xlim, ylim = ylim,
         type = "l", xlab = xlab, ylab = ylab, main = main,
         lty = lty[1], lwd = lwd[1], col = col[1], axes = FALSE)
    for(k in 2:length(dd))
        lines(dd[[k]]$x, dd[[k]]$y, lty = lty[k], lwd = lwd[k], col = col[k])
}

cols6 <- c(gray(.5), 2, 4, gray(.5), 2, 4)

dir.create("figures")
flc1 <- "figures/besag_Cposteriors.png"
png(flc1, 7000, 5000, res = 400)
par(mfrow = c(12-1, 12-1), mar = c(.1,.1,.1,.1),
    mgp = c(1.5, 0.5, 0), bty = 'n')
for(i0 in 1:11) {
    for(j0 in 1:12) {
        iedge <- ((i0-1)*nc + j0) %in% baseC$iLtheta
        if(j0>i0) {
            i <- odeaths[i0]
            j <- odeaths[j0]
            a <- max(-1, min(C0r[i,j], baseC$base[i,j], 0,
                             quantile(csampls[[2]][i,j,],0.001)))
            b <- min(1, max(C0r[i,j], baseC$base[i,j],0,1,
                            quantile(csampls[[2]][i,j,],0.999)))
            a = -0.3; b = 1
            dd <- list(pr1=density(sprior.5[[1]][i,j,], 0.05, from = a, to = b))
            dd$pr2 <- density(sprior.5[[2]][i,j,], 0.05, from = a, to = b)
            dd$pr3 <- density(sprior.5[[3]][i,j,], 0.05, from = a, to = b)
            dd$po1 <- density(csampls[[1]][i,j,], 0.05, from = a, to = b)
            dd$po2 <- density(csampls[[2]][i,j,], 0.05, from = a, to = b)
            dd$po3 <- density(csampls[[3]][i,j,], 0.05, from = a, to = b)
            dmplot(dd[c(1,2,3:6)], xlab = '', ylab = '', main = '',
                   lty = c(3,3,3,1,1,1),
                   lwd = c(2,2,2,2,2,2),
                   col = cols6)
            rug(seq(1,a,-.1), 0.05)
            rug(seq(1,a,-.2), 0.1)
            legend("top", "", title = paste0(causes_r[i], "\n", causes_r[j]),
                   title.col = 1+iedge, bty = "n")
            yh <- sapply(dd, function(d) max(d$y)*0.99)
            segments(0, 0, 0, yh/2, lty = 3, col = gray(0.5))
            segments(1, 0, 1, yh, lty = 1, col = 1)
            rug(x = C0r[i,j], ticksize = 0.15, col = 'black', lwd = 2)
            rug(x = baseC$base[i,j], ticksize=0.1, col = 2, lwd = 2, lty = 3)
        } else {
            if(j0<i0) {
                plot(0,type='n',axes=F,xlab='',ylab='')
            }
        }
        if((i0==4) & (j0==1)) {
            par(xpd=NA)
            legend("top", ncol = 1, bty = "n", lwd = 2, lty = c(3,3,3,1,1,1), 
                   c("Prior P0", "Prior P1", "Prior P2",
                     "Posterior|P0", "Posterior|P1", "Posterior|P2"),
                   col = cols6, cex = 1.2)
            par(xpd=TRUE)
        }
    }
}
dev.off()

if(FALSE)
    system(paste("eog", flc1, "&"))

### visualize selected pairs
ijsel <- list(c(1,8), c(1,9), c(2,3), c(3,6), c(3,7), c(3,12),
              c(5,11), c(6,12), c(10,12))
flc1s <- "figures/besag_Cposteriors_sel.png"
png(flc1s, 4000, 2500, res = 600)
par(mfrow = c(3,3), mar = c(1.5,1.5,1,0),
    mgp = c(1.5, 0.5, 0), bty = 'n')
jjp <- 1:length(ijsel)
for(k in jjp) {
    i = ijsel[[k]][1]
    j = ijsel[[k]][2]    
        iedge <- ((i-1)*nc + j) %in% baseC$iLtheta
            a <- max(-1, min(C0r[i,j], baseC$base[i,j], 0,
                     quantile(csampls[[2]][i,j,],0.001)))
            b <- min(1, max(C0r[i,j], baseC$base[i,j],0,1,
                            quantile(csampls[[2]][i,j,],0.999)))
            dd <- list(pr1=density(sprior.5[[1]][i,j,], 0.04, from = a, to = b))
            dd$pr2 <- density(sprior.5[[2]][i,j,], 0.04, from = a, to = b)
            dd$pr3 <- density(sprior.5[[3]][i,j,], 0.04, from = a, to = b)
            dd$po1 <- density(csampls[[1]][i,j,], 0.04, from = a, to = b)
            dd$po2 <- density(csampls[[2]][i,j,], 0.04, from = a, to = b)
            dd$po3 <- density(csampls[[3]][i,j,], 0.04, from = a, to = b)
            dmplot(dd[c(1,2,3:6)], xlab = '', ylab = '', main = '',
                   lty = c(3,3,3,1,1,1),
                   lwd = c(2,2,2,2,2,2),
                   col = cols6)
            title(main = paste0(causes_r[i], ", ", causes_r[j]),
                  cex.main = 0.95, col.main = 1+iedge)
            abline(v = 0, lty = 2, col = gray(0.05))
            rug(x = C0r[i,j], ticksize = 0.15, col = 'blue', lwd = 2)
            rug(x = baseC$base[i,j], ticksize=0.1, col = 2, lwd = 2, lty = 3)        
        if(k==1) {
            legend("top", ncol = 2, bty = "n", lwd = 2, lty = c(3,3,3,1,1,1), 
                   c("Prior P0", "Prior P1", "Prior P2",
                     "Posterior|P0", "Posterior|P1", "Posterior|P2"),
                   col = cols6, cex = 0.7)
        }
}
dev.off()

system(paste("eog", flc1s, "&"))

head(smap)
table(substr(smap$CD_MUN,1,2))

## aggregate municipalities into states (for cosmetic map, not really need)
library(tidyverse)
ufmap <- smap %>% group_by(abbrev_state) %>%
    summarize() %>% ungroup()

range(log(smro <- (0.1+dataf[,nc+ireord])/(0.1+dataf[,ireord])))
brko <- seq(-1,1,length=51)*max(abs(log(smro)))*1.001
range(brko)

cor(sapply(fits[7:9], function(r) r$summary.random$icause$mean))

fit <- fits[[length(fits)]]

unstr.m <- t(matrix(fit$summary.random$icause$mean, nc))
str.m <- matrix(fit$summary.random$ijoint$mean, n.s)
range(str.m)

brku <- seq(-1,1,length=21)*max(abs(unstr.m))*1.001
brks <- seq(-1,1,length=21)*max(abs(str.m))*1.001

library(ggplot2)

gg00 <- ggplot(smap) + theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.8))

gg0o <- gg00 +
    scale_fill_distiller(palette = "RdBu", minor_breaks = brko)
gg0u <- gg00 + 
    scale_fill_distiller(palette = "RdBu", minor_breaks = brku)
gg0s <- gg00 + 
    scale_fill_distiller(palette = "RdBu", minor_breaks = brks)

ggobs <- lapply(1:nc, function(i)
    gg0o +
    geom_sf(aes(fill = log(smro[, i])), 
            color = 'transparent') +
    labs(fill = causes_r[i]) +
    geom_sf(data=ufmap, fill = 'transparent'))
gguns <- lapply(1:nc, function(i)
    gg0u +
    geom_sf(aes(fill = unstr.m[, i]), color = 'transparent') +
    labs(fill = causes_r[i]) +
    geom_sf(data=ufmap, fill = 'transparent'))
ggstr <- lapply(1:nc, function(i)
    gg0s +
    geom_sf(aes(fill = str.m[, i]), color = 'transparent') +
    labs(fill = causes_r[i]) +
    geom_sf(data=ufmap, fill = 'transparent'))

library(ggpubr)

##ggarrange(plotlist = ggobs)#c(ggobs, gguns, ggstr),
##ggarrange(plotlist = ggobs)#c(ggobs, gguns, ggstr),
##          nrow = 3, ncol = nc)

apply(matrix(st_bbox(map), 2), 1, diff)


flmaps1 <- "figures/smoothedSMR12c.png"
png(flmaps1, 3600, 3600, res = 200)
ggarrange(plotlist = ggstr, nrow = 3, ncol = 4)
dev.off()

flmaps2 <- "figures/maps_SMR_all.png"
png(flmaps2, 4800, 3600, res = 200)
ggarrange(plotlist = c(ggobs, gguns, ggstr),
          nrow = 3, ncol = nc)
dev.off()

system(paste("eog", flmaps1, "&"))
