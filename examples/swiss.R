
## A data frame with 47 observations on 6 variables,
## _each_ of which is in percent, i.e., in [0, 100].
## ‘Fertility’         Ig, common standardized fertility measure
## ‘Agriculture’       % of males involved in agriculture as occupation
## ‘Examination’       % draftees receiving highest mark on army examination
## ‘Education’         % education beyond primary school for draftees.
## ‘Catholic’          % ‘catholic’ (as opposed to ‘protestant’).
## ‘Infant.Mortality’  live births who live less than 1 year.
## NOTE: All variables but ‘Fertility’ give proportions of the population.

data(swiss)

str(swiss)

par(mfrow = c(1, 1), mar = c(0,0,0,0), bty = "n")
pairs(swiss)

sdat <- scale(swiss)

(n <- nrow(sdat))
(p <- ncol(sdat))

round((cc <- cov(sdat)) * 100)

vnams <- colnames(sdat)
ii <- 1:p
ilp <- which(lower.tri(cc))

cobs.v <- cc; cobs.v[ilp] <- NA

library(fields)

par(mfrow = c(1, 1), mar = c(0,0,0,0), bty = "n")
image.plot(ii, ii, cobs.v, xlab = "", ylab = "")
text(ii, ii, vnams)

lcc <- chol(cc + diag(p) * 0.0)
qc <- chol2inv(lcc)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
  list(vnams, vnams)
round(pC*100)

## define a graphpcor from a minimum spanning tree
nb <- lapply(1:p, function(i)
  setdiff(1:p,i)); class(nb) <- 'nb'
nbc <- lapply(1:p, function(i) 1-abs(pC[i, -i]))

library(spdep)
nbw <- nb2listw(nb, nbc, style="B")
mst <- mstree(nbw)

G0 <- matrix(0, p, p, dimnames = dimnames(pC))
for(i in 1:nrow(mst)) {
  G0[mst[i,1], mst[i,2]] <- 1
  G0[mst[i,2], mst[i,1]] <- 1
}
G0

library(graphpcor)
g0 <- graphpcor(G0)
(dg0 <- dim(g0))

p*(p-1)/2

graphpcor(abs(pC)>0.2)
graphpcor(abs(pC)>0.25)
graphpcor(abs(pC)>0.3)
graphpcor(abs(pC)>0.35)
graphpcor(abs(pC)>0.4)

g1 <- graphpcor(abs(pC)>0.3)
g1

G1 <- attr(g1, "graph")
(dg1 <- dim(g1))

c(p=p, n=n)
g0
g1

attr(g1, 'graph')

png("g01swiss.png", width = 1800, height = 900, res = 300)
par(mfrow = c(1, 2), mar = c(0,2,0,2))
plot(g0, Rgraphviz = TRUE)
plot(g1, Rgraphviz = TRUE)
dev.off()

if(FALSE)
    system("eog g01swiss.png &")

c0 <- cgeneric(g0, lambda = 5,
     base = rep(0, dg0[2]))
c1 <- cgeneric(g1, lambda = 5,
     base = rep(0, dg1[2]))
cpc <- cgeneric("pc_correl", lambda = 5, n = p,
                base = rep(0, p*(p-1)/2))
lkj <- cgeneric("LKJ", eta = 10, n = p)

idat <- data.frame(
    i = rep(1:p, each = n),
    r = rep(1:n, p),
    y = as.vector(sdat)
)

head(sdat,3)
head(idat,3)

cfam <- list(
    hyper = list(
        prec = list(intial = 10, fixed = TRUE)
    ))

library(INLA)

fit0 <- inla(
    y ~ 0 + f(i, model = c0, replicate = r),
    data = idat,
    control.family = cfam
)
fit1 <- inla(
    y ~ 0 + f(i, model = c1, replicate = r),
    data = idat,
    control.family = cfam
)
fitcpc <- inla(
    y ~ 0 + f(i, model = cpc, replicate = r),
    data = idat,
    control.family = cfam
)
fitlkj <- inla(
    y ~ 0 + f(i, model = lkj, replicate = r),
    data = idat,
    control.family = cfam
)

rbind(fit0$cpu.used, fit1$cpu.used,
      fitcpc$cpu.used, fitlkj$cpu.used)

c(fit0$misc$nfunc, fit1$misc$nfunc,
  fitcpc$misc$nfunc, fitlkj$misc$nfunc)

c0fit <- vcov(g0, theta = fit0$mode$theta)
c1fit <- vcov(g1, theta = fit1$mode$theta)
cpcfit <- tcrossprod(cholcor(fitcpc$mode$theta, p = p))

c0fit.v <- c0fit; c0fit.v[upper.tri(cc, TRUE)] <- NA
c1fit.v <- c1fit; c1fit.v[upper.tri(cc, TRUE)] <- NA
cpcfit.v <- cpcfit; cpcfit.v[upper.tri(cc, TRUE)] <- NA

rgbcfn <- function(x, rank = FALSE, ab = range(x, na.rm = TRUE)) {
    if(rank) {
        u <- rank(x)/length(x)
    } else {
        u <- (x-ab[1])/diff(ab)
    }
    rgb(dnorm(u,1.0,0.5)/dnorm(1.0,1.0,0.5),
        dnorm(u,0.5,1.0)/dnorm(0.5,0.5,1.0/2),
        dnorm(u,0.0,0.5)/dnorm(0.0,0.0,0.5))
}

range(cc[ilp])
(cm <- max(abs(cc[ilp]), abs(c0fit[ilp]), abs(c1fit[ilp]), abs(cpcfit[ilp])))
bkc <- seq(-1, 1, 0.1)*cm*1.001
ncols <- length(bkc)-1
cols <- rgbcfn((bkc[-1] + bkc[1:ncols])/2)

c(length(bkc), length(cols))

ig0v <- graph_from_adjacency_matrix(
    upperPadding(attr(g0, "graph")))
ig0v$layout <- layout.fruchterman.reingold ##layout.kamada.kawai
ig1v <- graph_from_adjacency_matrix(
    upperPadding(attr(g1, "graph")))
ig1v$layout <- layout.fruchterman.reingold ##layout.kamada.kawai

par(mfrow = c(2,2), mar = c(0,0,0,0), bty = "n")
image.plot(ii, ii, cobs.v, breaks = bkc, col = cols, 
           axes = FALSE, xlab = "", ylab = "")
text(ii, ii, vnams)
text(row(cc), col(cc), 
     ifelse(is.na(t(c0fit.v)), "", format(cc*100, digits = 2)))
image.plot(ii, ii, cpcfit.v, breaks = bkc, col = cols, add = TRUE,
           axes = FALSE, xlab = "", ylab = "")
text(row(cc), col(cc),
     ifelse(is.na(cpcfit.v), "", format(cpcfit*100, digits = 1)))
plot(ig0v, edge.arrow.mode = 0)
plot(ig1v, edge.arrow.mode = 0)
image.plot(ii, ii, t(c0fit.v), breaks = bkc, col = cols, 
           axes = FALSE, xlab = "", ylab = "")
image.plot(ii, ii, c1fit.v, breaks = bkc, col = cols, 
           axes = FALSE, xlab = "", ylab = "", add = TRUE)
text(ii, ii, vnams)
text(row(cc), col(cc),
     ifelse(is.na(c0fit.v), "", format(c0fit*100, digits = 1)))
text(row(cc), col(cc),
     ifelse(is.na(t(c1fit.v)), "", format(t(c1fit)*100, digits = 1)))


nps <- 10000
h0sampls <- inla.hyperpar.sample(
    n = nps, result = fit0, intern = TRUE
)
h1sampls <- inla.hyperpar.sample(
    n = nps, result = fit1, intern = TRUE
)
hcpcsampls <- inla.hyperpar.sample(
    n = nps, result = fitcpc, intern = TRUE
)
hlkjsampls <- inla.hyperpar.sample(
    n = nps, result = fitlkj, intern = TRUE
)

iup <- which(upper.tri(cc))
c0sampls <- t(sapply(1:nps, function(i) {
    vcov(g0, theta = h0sampls[i, ])[iup]
}))
c1sampls <- t(sapply(1:nps, function(i) {
    vcov(g1, theta = h1sampls[i, ])[iup]
}))
cpcsampls <- t(sapply(1:nps, function(i) {
    tcrossprod(cholcor(hcpcsampls[i, ]))[iup]
}))
clkjsampls <- t(sapply(1:nps, function(i) {
    tcrossprod(cholcor(hlkjsampls[i, ]))[iup]
}))


lG0 <- upperPadding(G0)
lG1 <- upperPadding(G1)

fcols <- c(gray(0.35), rgb(0,0,1,0.7))

##plot(1:2, col = fcols, pch = 19, cex = 5)

png("swissResult.png", width = 4000, height = 3000, res = 300)
par(mfrow = c(p, p), mar = c(2.1, 2.1, 0.1, 0.1),
    mgp = c(1.5, 0.5, 0), bty = 'n')
kc <- k2 <- k1 <- k0 <- 0
for(i in 1:p) {
    for(j in 1:p) {
        if(i==j) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, vnams[j], cex = 2)
        }
        if(j>i) {
            iij0 <- which((i == (lG0@i+1) & (j == (lG0@j+1))))
            iij1 <- which((i == (lG1@i+1) & (j == (lG1@j+1))))
            if(length(iij0)>0) {
                k0 <- k0 + 1
                m0 <- inla.smarginal(fit0$internal.marginals.hyperpar[[k0]])
                h0 <- hist(h0sampls[, k0], 100, plot = FALSE)
                h0$ok <- TRUE
            } else {
                h0 <- list(ok=FALSE)
            }
            if(length(iij1)>0) {
                k1 <- k1 + 1
                m1 <- inla.smarginal(fit1$internal.marginals.hyperpar[[k1]])
                h1 <- hist(h1sampls[, k1], 100, plot = FALSE)
                h1$ok <- TRUE
            } else {
                h1 <- list(ok = FALSE)
            }
            if(h0$ok) {
                plot(m0, type = "l", ##h0, freq = FALSE,
                     main = '', xlab = '', ylab = '',
                     xlim = range(0, h0$breaks, h1$breaks),
                     ylim = range(h0$dens, h1$dens),
                     col = fcols[1], border = 'transparent')
                if(h1$ok) {
                    lines(m1, ##plot(h1, freq = FALSE, add = TRUE,
                         col = fcols[2], border = 'transparent')
                    legend("topleft", bty = "n",
                           as.expression(lapply(c(k0,k1), function(i) bquote(theta[.(i)]))),
                           fill = fcols, border = 'transparent')
                } else {
                    legend("topleft", bty = "n",
                           as.expression(lapply(c(k0), function(i) bquote(theta[.(i)]))),
                           fill = c(fcols[1]), border = 'transparent')
                }
                rug(0, 0.1, lty = 3, lwd = 2)
            } else {
                if(h1$ok) {
                    plot(m1, type = "l", ##h1, freq = FALSE,
                         xlim = range(0, h1$breaks),
                         ylim = range(h0$dens, h1$dens),
                         main = '', xlab = '', ylab = '',
                         col = fcols[2], border = 'transparent')
                    legend("topleft", bty = "n",
                           as.expression(lapply(c(k1), function(i) bquote(theta[.(i)]))),
                           fill = fcols[2], border = 'transparent')                    
                    rug(0, 0.1, lty = 3, lwd = 2)
                } else {
                    plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')                    
                    rug(0, 0.1, lty = 3, lwd = 2)
                }
            }
        }
        if(j<i) {
            k2 <- k2 + 1
            c_obs <- cc[iup[k2]]
            ic_obs <- tanh(c_obs + qnorm(c(0.025, 0.975)) / sqrt(n-3))
            pits <- c(mean(c0sampls[, k2] <= c_obs),
                      mean(c1sampls[, k2] <= c_obs),
                      mean(cpcsampls[, k2] <= c_obs),
                      mean(clkjsampls[, k2] <= c_obs))
            dcpc <- density(cpcsampls[, k2])
            dlkj <- density(clkjsampls[, k2])
            h0 <- hist(c0sampls[, k2], 100, plot = FALSE)
            h1 <- hist(c1sampls[, k2], 100, plot = FALSE)
            plot(h0, freq = FALSE, 
                 xlim = range(ic_obs, h0$breaks, h1$breaks, dcpc$x, dlkj$x),
                 ylim = range(h0$dens, h1$dens), 
                 main = '', xlab = '', ylab = '',
                 col = fcols[1], border = 'transparent')
            plot(h1, freq = FALSE, add = TRUE,
                 col = fcols[2], border = 'transparent')
            lines(dcpc, lty = 2)
            lines(dlkj, lwd = 2, lty = 3, col = 6)
            rug(c(c_obs, ic_obs), 0.2, lty = 1, lwd = 1, col = 'red')
            if(FALSE)
                legend('topleft', title = "PIT", bty = 'n',
                       paste(c("mst", "G1", "pc-cor", "LKJ"), ":", 
                             format(pits, digits = 2)), ncol = 1)
            if(k2==1) {
                legend("left", title = "", bty = "n",
                       c("Obs. IC"), lty = c(1), lwd = c(1), col = c(2))
                legend("topleft", title = "prior", bty = "n",
                       c("PC with MST", "PC with G1", "PC-correl", "LKJ"),
                       lty = c(0, 0, 2, 3), lwd = c(0,0, 1, 2), col = c(0,0,1,6),
                       fill = c(fcols,rep("transparent",2)), border = 'transparent')
            }
        }
    }
}
dev.off()

if(FALSE)
    system("eog swissResult.png &")
