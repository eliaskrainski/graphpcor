
##   A data frame with 462 observations on the following 10 variables.
##  sbp systolic blood pressure
##  tobacco cumulative tobacco (kg)
##  ldl low density lipoprotein cholesterol
##  adiposity a numeric vector
##  famhist family history of heart disease,
##       a factor with levels ‘"Absent"’ and ‘"Present"’
##  typea type-A behavior
##  obesity a numeric vector
##  alcohol current alcohol consumption
##  age age at onset
##  chd response, coronary heart disease


data(SAheart, package = "msos")

str(SAheart)

sdat <- scale(within(SAheart, {
    famhist <- (famhist=="Present")
}))

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
lcc

qc <- chol2inv(lcc)
round(qc, 2)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(vnams, vnams)
round(pC*100)

## define a graphpcor from a minimum spanning tree
library(spdep)
nb <- lapply(1:p, function(i)
    setdiff(1:p,i)); class(nb) <- 'nb'
nbc <- nbcosts(nb, scale(mtcars))
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
g1 <- graphpcor(abs(pC)>0.125)
G1 <- attr(g1, "graph")
(dg1 <- dim(g1))

c(p=p, n=n)
g0
g1

attr(g1, 'graph')

par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(g0)
plot(g1)

c0 <- cgeneric(g0, lambda = 1,
               base = rep(0, dg0[2]), 
               useINLAprecomp = FALSE)
c1 <- cgeneric(g1, lambda = 1,
               base = rep(0, dg1[2]), 
               useINLAprecomp = FALSE)

idat <- data.frame(
    i = rep(1:p, each = n),
    r = rep(1:n, p),
    y = as.vector(sdat)
)

head(sdat,3)
head(idat,3)

cfam <- list(hyper = list(
                 prec = list(intial = 10, fixed = TRUE)
             ))

library(INLA)

fit0 <- inla(
    y ~ 0 + f(i, model = c0, replicate = r),
    data = idat,
    control.family = cfam,
    control.mode = list(
        theta = rnorm(dg0[2]),
        restart = TRUE
    )
)
fit1 <- inla(
    y ~ 0 + f(i, model = c1, replicate = r),
    data = idat,
    control.family = cfam,
    control.mode = list(
        theta = rnorm(dg1[2]),
        restart = TRUE
    )
)

fit0$mode$theta
fit1$mode$theta

c0fit <- vcov(g0, theta = fit0$mode$theta)
c1fit <- vcov(g1, theta = fit1$mode$theta)

round(cc*100)
round(c0fit*100)
round(c1fit*100)

c1fit.v <- c1fit; c1fit.v[upper.tri(cc, TRUE)] <- NA

range(cc[ilp])
(cm <- max(abs(cc[ilp]), abs(c1fit[ilp])))
bkc <- seq(-1, 1, 0.1)*cm*1.001
ncols <- length(bkc)-1
cols <- rgb(ncols:1/ncols, 0.3, 1:ncols/ncols)

c(length(bkc), length(cols))

par(mfrow = c(1, 1), mar = c(0,0,0,0), bty = "n")
image.plot(ii, ii, cobs.v, breaks = bkc, col = cols, 
           axes = FALSE, xlab = "", ylab = "")
image.plot(ii, ii, t(c1fit.v), breaks = bkc, col = cols, 
           axes = FALSE, xlab = "", ylab = "", add = TRUE)
text(ii, ii, vnams)
text(row(cc), col(cc),
     ifelse(is.na(cobs.v), "", format(cc*100, digits = 2)))
text(row(cc), col(cc),
     ifelse(is.na(t(cobs.v)), "", format(t(c1fit)*100, digits = 1)))


hsampls <- inla.hyperpar.sample(
    n = 10000, result = fit1, intern = TRUE
)
iup <- which(upper.tri(cc))
csampls <- t(sapply(1:nrow(hsampls), function(i) {
    vcov(g1, theta = hsampls[i, ])[iup]
}))

dim(hsampls)
dim(csampls)

lG1 <- t(upperPadding(G1))
str(lG1)
lG1

par(mfrow = c(10, 10), mar = c(2.1, 2.1, 0.1, 0.1),
    mgp = c(1.5, 0.5, 0), bty = 'n')
k1 <- 0; k2 <- 0
for(i in 1:p) {
    for(j in 1:p) {
        if(i==j) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, vnams[j], cex = 2)
        }
        if(j<i) {
            iij <- which((i == (lG1@i+1) & (j == (lG1@j+1))))
            if(length(iij)>0) {
                k1 <- k1 + 1
                h <- hist(hsampls[, k1], 100, plot = FALSE)
                plot(h, xlim = range(0, h$breaks),
                     main = '', xlab = '', ylab = '',
                     col = gray(0.35), border = 'transparent')
                abline(v = 0, lty = 3, lwd = 2)
                ##legend('topleft', paste(i,j), bty = 'n')                
            } else {
                plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
            }
        }
        if(j>i) {
            k2 <- k2 + 1
            h <- hist(csampls[, k2], 100, plot = FALSE)
            plot(h, xlim = range(range(h$breaks), cc[iup[k2]]),
                 main = '', xlab = '', ylab = '',
                 col = gray(0.35), border = 'transparent')
            abline(v = cc[iup[k2]], lty = 2, lwd = 2, col = 'red')
            abline(v = 0, lty = 3, lwd = 2)
            ##legend('topleft', paste(i,j), bty = 'n')
        }
    }
}

png('graph1SAheart.png', 2000, 2000, res = 300)
par(mar = c(0,0,0,0))
plot(g1)
dev.off()

system("eog graph1SAheart.png &")
