library(graphpcor)

## STEP 1: data model definition

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

## model idea:
##   model chd and famhist = y1[1:n, 1:2] as
##     y1[i,j] ~ bernoulli(p_{ij}), p_{ij} = inv_logit(mu[j] + x_sigmas[j] * x[i,j])
##   model the other p2=p-p2=p-2 variables as
##     y2[i,j] = Normal(mu[k] + x_sigmas[k] * x[i,k], y2_sigmas[j]), k = 2+j
##   model the latent x as
##     x ~ multi_normal(0, R), R a correlation matrix

## STAN model base code
Scode0 <- '
data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> p1;
  int<lower=1> p2;
  array[n,p1] int<lower=0, upper=1> y1;
  matrix[n,p2] y2;
  real<lower=0> x_sigmas_lambda;
  real<lower=0> y2_sigmas_lambda;
}
parameters {
  vector[p] mu;
  vector<lower=0>[p] x_sigmas;
  vector<lower=0>[p2] y2_sigmas;
}
transformed parameters {
  corr_matrix[p] rho;
}
model {
  vector[p] x[n];
  x ~ multi_normal(rep_vector(0, p), rho);
  for(j in 1:p1) {
    for(i in 1:n) {
      y1[i,j] ~ bernoulli(inv_logit(mu[j] + x_sigmas[j] * x[i,j]));
    }
  } 
  for(j in 1:p2) {
    for(i in 1:n) {
      y2[i,j] ~ normal(mu[p1+j] + x_sigmas[p1+j], y2_sigmas[j]);
    }
  }
  for(j in 1:p) {
    x_sigmas[j] ~ exponential(x_sigmas_lambda);
  }
  for(j in 1:p2) {
    y2_sigmas[j] ~ exponential(y2_sigmas_lambda);
  }
}
'

## STEP 2: update the STAN code with code for the
## graphpcor prior for 'rho'
Sgrpc <- stan_add(
    Scode0, 'graphpcor',
    lambda = 1, name = "rho")

cat(Sgrpc)

## STEP 3:  compile STAN code
library(rstan)
options(mc.cores = 4L)

system.time(
    Sgpc_cmpld <- stan_model(
        model_code = gsub("\n  \n", "\n", Sgrpc, fixed = TRUE),
        model_name = "graphpcor"
    )
)

## STEP 4: prepare the data
data(SAheart, package = "msos")

str(SAheart)

## re-order the columns having 'chd' at first and 'famhist'
## at the second column, converted to binary with
## 0 for 'Absent', 1 for 'Present'
dataf <- data.frame(
    chd = SAheart$chd,
    famhist = (SAheart$famhist=="Present")+0L,
    SAheart[, c(1:4, 6:9)]
)
str(dataf)

(n <- nrow(dataf))
(p <- ncol(dataf))

## STAN (initial) data
Sdata0 <- list(
    n = as.integer(n),
    p = as.integer(p),
    p1 = 2L,
    p2 = as.integer(p-2),
    y1 = as.matrix(dataf[, 1:2]),
    y2 = as.matrix(dataf[, 3:p]),
    x_sigmas_lambda = 1,
    y2_sigmas_lambda = 1
)

str(Sdata0)

round((cc <- cor(dataf)) * 100)

vnams <- colnames(dataf)

lcc <- chol(cc)
qc <- chol2inv(lcc)

## partial correlation matrix
pC <- cov2cor(qc)
dimnames(pC) <- dimnames(qc) <- dimnames(cc) <-
    list(vnams, vnams)
round(pC*100, 2)

## define a minimum spanning tree considering the functions
## that I have implemented (>20y ago!) in the spdep package
## this MST consider the weights w_{ij} = 1-|\rho_{ij}|
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

## define a graph consideing partial correlations (abs) > 0.09
tanh(qnorm(c(0.025, 0.975)) / sqrt(n-3))

g1 <- graphpcor(abs(pC)>0.09)
G1 <- attr(g1, "graph")
(dg1 <- dim(g1))

c(p=p, n=n)
g0
g1

if(FALSE)
    png("g01SAheart.png", width = 1800, height = 900, res = 300)
par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(g0, Rgraphviz = TRUE)
plot(g1, Rgraphviz = TRUE)
if(FALSE)
    dev.off()

if(FALSE)
    system("eog g01SAheart.png &")

## base model
baseM0 <- basepcor(diag(p), iLtheta = g0)
baseM1 <- basepcor(diag(p), iLtheta = g1)

##   update STAN data
SdataM0 <- stan_add(Sdata0, baseM0, lambda = 1, name = 'rho')
SdataM1 <- stan_add(Sdata0, baseM1, lambda = 1, name = 'rho')

## STAN sampling
Samples1 <- sampling(
    Sgpc_cmpld, 
    data = SdataM0,
    iter = 30000,
    warmup = 5000,
    chains = 4
)

## PLOTS
library(coda)

thnams <- paste0("grpc_theta[", 1:m, "]")
rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

th.base
th.true
mcmc_intervals(Samples, thnams)

baseC$base
corr
ycorr
mcmc_intervals(Samples, rhonams)

library(ggplot2)

iil <- which(lower.tri(corr))

ggcc <- mcmc_areas(
  Samples,
  pars = rhonams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) + geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = ycorr[iil],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'black', lty = 2) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = corr[iil],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'red', lty = 2) +
    geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = baseC$base[iil],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'blue', lty = 2) 

ggcc


## my own plots
Sth <- Reduce("cbind", extract(Samples, thnams))
colnames(Sth) <- thnams
Scorr <- Reduce("cbind", extract(Samples, rhonams))
colnames(Scorr) <- rhonams

summary(Scorr)
(ycorr <- cor(y))
postCmean <- colMeans(Scorr)
postCmean


par(mfrow = c(p, p), mar = c(3.5,1.5,0.5,0.5),
    mgp = c(1.5,0.5,0), bty = "n")
k2 <- k1 <- 0
for(i in 1:p) {
    for(j in 1:p) {
        if(i==j) {
            plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
            if(j==1)
                legend("topleft", c("Observed", "TRUE", "base", "Posterior"),
                       bty = "n", col = c(1, 2, 4, 0),
                       lty = 2, lwd = 2, fill = c(0,0,0,gray(0.5)),
                       border = 'transparent')
        }
        if(j>i) {
            if(length(Sdata$grpc_ii)>k1 &&
               (i==Sdata$grpc_jj[k1+1]) &
               (j==Sdata$grpc_ii[k1+1])) {
                k1 <- k1 + 1
                thk1 <- c(NA, th.true[k1], th.base[k1])
                h <- hist(Sth[, k1], 100, plot = FALSE)
                plot(h, main = '', freq = FALSE, border = 'transparent',
                     xlim = range(h$breaks, thk1, na.rm = TRUE),
                     xlab = as.expression(bquote(theta[.(k1)])))
                abline(v = thk1, col = c(1,2,4),
                       lty = 2, lwd = 2)                
            } else {
                plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
            }
        }
        if(j<i) {
            k2 <- k2 + 1
            cc <- c(ycorr[iil[k2]], corr[iil[k2]],
                    baseC$base[iil[k2]])
            hS <- hist(Scorr[,k2], 100, plot = FALSE)
            plot(hS, freq = FALSE, main = '', border = 'transparent',
                 xlab = as.expression(bquote(rho[.(i)~","~.(j)])),
                 xlim = range(hS$breaks, cc, na.rm = TRUE))
            abline(v = cc, col = c(1,2,4), lty = 2, lwd = 2)
        }
    }
}



c0 <- cgeneric(g0, lambda = 10,
               base = rep(0, dg0[2]), 
               useINLAprecomp = FALSE)
c1 <- cgeneric(g1, lambda = 10,
               base = rep(0, dg1[2]), 
               useINLAprecomp = FALSE)

## two dense priors
##  pc_correl
cpc <- cgeneric("pc_correl", lambda = 10, n = p,
                base = rep(0, p*(p-1)/2), 
                useINLAprecomp = FALSE)
## LKJ
lkj <- cgeneric("LKJ", eta = 30, n = p,
                useINLAprecomp = FALSE)

idat <- data.frame(
    i = rep(1:p, each = n),
    r = rep(1:n, p),
    y = as.vector(dataf)
)

head(dataf,3)
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

if(FALSE)
    all.equal(chol2inv(chol(as.matrix(prec(lkj, theta = fitcpc$mode$theta)))),
              chol2inv(chol(as.matrix(prec(cpc, theta = fitcpc$mode$theta)))))
if(FALSE)
    all.equal(chol2inv(chol(as.matrix(prec(lkj, theta = fitlkj$mode$theta)))),
              chol2inv(chol(as.matrix(prec(cpc, theta = fitlkj$mode$theta)))))

all.equal(chol2inv(t(cholcor(fitcpc$mode$theta, p = p))),
          as.matrix(prec(cpc, theta = fitcpc$mode$theta)))

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

png("SAheartResult.png", width = 4000, height = 3000, res = 300)
par(mfrow = c(p,p), mar = c(1.6, 1.6, 0.1, 0.1),
    mgp = c(1, 0.5, 0), bty = 'n')
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
                 xlim = range(ic_obs, h0$breaks, h1$breaks),#, dcpc$x, dlkj$x),
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
        }
        if((i==1) & (j==2)) {
            legend("bottom", title = "", bty = "n",
                   c("Obs. IC"), lty = c(1), lwd = c(1), col = c(2))
        }
        if((i==1) & (j==2)) {
            legend("top", title = "prior", bty = "n",
                   c("PC with MST", "PC with G1", "PC-correl"), ##"LKJ"),
                   lty = c(0, 0, 2, 3), lwd = c(0,0, 1, 2), col = c(0,0,1,6),
                   fill = c(fcols,rep("transparent",2)),border = 'transparent')
        }
    }
}
dev.off()


if(FALSE)
    system("eog SAheartResult.png &")

}
