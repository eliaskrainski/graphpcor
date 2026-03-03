## EXAMPLE definition
## STAN model definition
##   initial code
##   update code
##   compile code
## DATA sample
##   CORRELATION matrix definition
##   SAMPLE drawn
##   initial STAN data
##   base model definition
##   update STAN data with PC-prior parameters
## MCMC sampling
## PLOTS

## EXAMPLE definition
## y ~ N(\mu, \Sigma)
##  where \mu = [0,...,0] and \Sigma = C(\theta) 
## theta[1:m]: vector of parameters to be placed in L[iil]
##       where iil is an index vector for the lower triangle

## MODEL definition
##   initial STAN code
Scode0 <- '
data {
  int n;
  int p;
  vector[p] y[n];
  vector[p] ymu;
}
model {
  y ~ multi_normal_cholesky(ymu, LC);
}
generated quantities {
  corr_matrix[p] rho = tcrossprod(LC);
}
'

##   update STAN code 
library(graphpcor)
## add the code for the PC prior for a correlation 
## matrix through its Cholesky decomposition,
## which was named "LC" in the initial code
Scode <- stan_add(
    Scode0, "pc_correl",
    lambda = 1, name = "LC"
)

##   compile STAN code
library(rstan)
options(mc.cores = 4L)

system.time(
    Spccor <- stan_model(
        model_code = Scode, 
        model_name = "pc_dense"
    )
)

## DATA sample
##   CORRELATION matrix definition
corr <- matrix(c( 1.0, 0.7, 0.5,-0.3,
                  0.7, 1.0,-0.2, 0.4,
                  0.5,-0.2, 1.0,-0.8,
                 -0.3, 0.4,-0.8, 1.0), 4)
all.equal(corr,t(corr))
Uc <- chol(corr)

p <- ncol(corr)

## DATA sample
n <- 100
set.seed(1)
y <- matrix(rnorm(n*p), n) %*% Uc

## initial STAN data
Sdata0 <- list(
    n = as.integer(n),
    p = as.integer(p),
    y = y,
    ymu = rep(0,p))


## compute the correlation parameters at the observed data
## to add into the plots later
corr
(ycorr <- cor(y))
th.obs <- basecor(ycorr)$theta
(th.true <- basecor(corr)$theta)

## base model definition
m <- p * (p-1)/2
set.seed(2)
th.base <- rnorm(m)
baseC <- basecor(th.base, p = p)

## update STAN data
Sdata <- stan_add(
    Sdata0, baseC, lambda = 1, name = "LC")

str(Sdata)

## STAN sampling
Samples <- sampling(
    Spccor,
    data = Sdata,
    iter = 30000,
    warmup = 5000,
    chains = 4
)

## PLOTS
library(coda)

thnams <- paste0("LC_theta[", 1:m, "]")

rhonams <- paste0("rho[",
                  unlist(lapply(2:p, function(i) i:p)), ",",
                  unlist(lapply(2:p, function(i) rep(i-1, p-i+1))),
                  "]")
rhonams

library("bayesplot")

th.true
th.base

mcmc_intervals(Samples, thnams)

mcmc_intervals(Samples, rhonams)

library(ggplot2)

ggcc <- mcmc_areas(
  Samples,
  pars = rhonams, 
  prob = 0.9, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) + geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = corr[lower.tri(corr)],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'red', lty = 2) +
 geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = ycorr[lower.tri(corr)],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'black', lty = 2) +
 geom_segment(aes(x = x, y = ya, xend = x, yend = yb),
                  data.frame(
                     x = baseC$base[lower.tri(corr)],
                     ya = 6:1-0.1,
                     yb = 6:1+0.5), color = 'blue', lty = 2)

ggcc

## my own plots
Sth <- Reduce("cbind", extract(Samples, thnams))
colnames(Sth) <- thnams
Scorr <- Reduce("cbind", extract(Samples, rhonams))
colnames(Scorr) <- rhonams

summary(Scorr)
postCmean <- colMeans(Scorr)
postCmean

iil <- which(lower.tri(corr))

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
            k1 <- k1 + 1
            thk1 <- c(th.obs[k1], th.true[k1], th.base[k1])
            h <- hist(Sth[, k1], 100, plot = FALSE)
            plot(h, main = '', freq = FALSE, border = 'transparent',
                 xlim = range(h$breaks, thk1),
                 xlab = as.expression(bquote(theta[.(k1)])))
            abline(v = thk1, col = c(1,2,4),
                   lty = 2, lwd = 2)
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

