### Focus only on C, for scaled data

## load the packages for the data analysis
options(mc.cores = 5L)
library(coda)
library(rstan)
library(graphpcor)

## set the working directory
setwd(here::here("uci_heart"))
getwd()

## get the data
if((!file.exists("dat1.rds")) &&
   (!file.exists("dat2.rds"))) {
    source("data_prepare.R")
    rm(list = setdiff(ls(all = TRUE), c("dat1", "dat2")))
} else {
    dat1 <- scale(readRDS("dat1.rds")) ## scaled data1
    dat2 <- scale(readRDS("dat2.rds")) ## scaled data1
}
stopifnot(all(colnames(dat1)==colnames(dat2)))

(n1 <- nrow(dat1))
(n2 <- nrow(dat2))
(p <- ncol(dat1))

ynames <- colnames(dat1)
C1obs <- cor(dat1)
round(C1obs*100)

###############################################################
### 1) Define the model to be fitted with STAN              ###
### The dense correlation models consider the correlation   ###
### matrix from its Cholesky, which is directly             ###
### parametrized in the LKJ model.                          ###
### The sparse model consider the Cholesky of the precision ###
### matrix instead and then map to the correlation matrix.  ###  
###############################################################

### First define a model to be completed later,
###  without the correlation matrix or its Cholesky.
Smodel0 <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  vector[p] y[n];
  vector[p] mu;
}
model {
  y ~ multi_normal_cholesky(mu, LCorr);
}
"

Smodel_LKJ <- stan_add_code(
    x = Smodel0,
    to_add = list(
        data = "  real<lower=0> eta;\n",
        parameters = "  cholesky_factor_corr[p] LCorr;\n",
        "generated quantities" = 
            "  corr_matrix[p] Corr;
  Corr = multiply_lower_tri_self_transpose(LCorr);\n",
        model = "  LCorr ~ lkj_corr_cholesky(eta);\n")
)

cat(Smodel_LKJ)

## Compile the LKJ model
system.time(
    cmpl_LKJ <- 
        stan_model(
            model_code = Smodel_LKJ, 
            model_name = "LKJ"
        )
)

## use STAN to drawn samples form the LKJ prior
## by using a zero correlated fake dataset
dat0 <- matrix(rnorm(n1 * p), n1)
(sIdat0 <- graphpcor:::dspd(cov(dat0))$sqrtInv)
cov(Dat0 <- dat0 %*% sIdat0)
Ssampl0 <- sampling(
     object = cmpl_LKJ, 
     data = list(y=Dat0,n=nrow(Dat0),p=p,mu=rep(0,p),eta=1),##(1+n1)/2),
     iter = 20500,
     warmup = 500,
     thin = 10,
     chains = 5
)

## prepare dat1 for STAN (without the correlation prior parameter)
Sdat1 <- list(
    y = dat1,      ## the 1st data
    n = n1,         ## number of observations in 1st data
    p = p,           ## number of variables
    mu = rep(0,p)     ## prior mean for the intercepts
)

## draw samples from the posterior
Ssampl1 <- sampling(
    object = cmpl_LKJ, 
    data = c(Sdat1,
             list(eta = 1)), ## LKJ prior parameter
    iter = 20500,
    warmup = 500,
    thin = 10,
    chains = 5
)

## names to extract the correlations (lower diag)
rhonams <- unlist(lapply(1:(p-1), function(j)
    paste0("Corr[", (j+1):p, ",", j, "]")))
rhonams

if(FALSE) {
    stan_trace(Ssampl1, rhonams)
}

## collect and organize the 4 parallel sample chains
C0Samples <- mcmc(Reduce("cbind", extract(Ssampl0, rhonams)))
C1Samples <- mcmc(Reduce("cbind", extract(Ssampl1, rhonams)))

##########################################################
### Define dense and sparse base models 
### based on the correlation with data 1
cc1mode <- sapply(hc1dens, function(d)
    d$x[which.max(d$y)])

C1mode <- diag(p)/2
C1mode[lower.tri(C1mode)] <- cc1mode
C1mode <- C1mode + t(C1mode)

round(C1obs*100)
round(C1mode*100)

BaseDens <- basecor(C1mode)
BaseGrph <- basepcor(C1mode, n = n1)
BaseGrph

round(abs((BaseGrph$base - C1mode)*100),1)

## build the graphpcor to visualization
Grph <- graphpcor(BaseGrph$iLtheta, p = p, nodes = ynames)

gfl <- "~/github/corr-priors/examples/figures/uci_graph.png"
png(gfl, 1000, 3000, res = 600)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot(Grph, Rgraphviz = TRUE)
dev.off()

system(paste("eog", gfl, "&"))

## sample from the prior and compare with the posterior from dat1
nsamples <- 1e4
cLambda <- 3
baseDsamples <- sample.basecor(BaseDens, nsamples, lambda = cLambda)
baseGsamples <- sample.basepcor(BaseGrph, nsamples, lambda = cLambda)

dmplot <- function(dd, xlab = '', ylab = '', main = '', xlim, ylim,
                   lty=1:length(dd), lwd=1:length(dd), col=1:length(dd)) {
    if(missing(xlim))
        xlim <- range(unlist(lapply(dd, function(d) d$x)))
    if(missing(ylim))
        ylim = range(unlist(lapply(dd, function(d) d$y)))   
    plot(dd[[1]]$x, dd[[1]]$y, xlim = xlim, ylim = ylim,
         type = "l", xlab = xlab, ylab = ylab, main = main,
         lty = lty[1], lwd = lwd[1], col = col[1])
    for(k in 2:length(dd))
        lines(dd[[k]]$x, dd[[k]]$y, lty = lty[k], lwd = lwd[k], col = col[k])
}

db2 <- function(x,a,b) {
   exp( (a-1)*log(1+x) + (b-1)*log(1-x) - lbeta(a,b) - (a+b-1)*log(2))
}

## visualize
flC1prior <- "~/github/corr-priors/examples/figures/uci_C1prior.png"
png(flC1prior, 4500, 3000, res = 600)
k = 0
par(mfcol = c(p-1, p-1), mar = c(3,3,0.5,0.5),
    mgp = c(1.5,0.7,0), las = 1, bty = 'n', xpd = NA)
for(j in 1:(p-1)) {
    for(i in 2:p) {
        il <- (j-1)*p + i
        if(j>=i) {
            plot(0, type = 'n', axes = FALSE,
                 xlab = '', ylab = '', main = '')
        } else {
            k <- k + 1
            a <- max(-1, min(0, C1mode[i,j]-.2))
            b <- min(1, max(0, C1mode[i,j]+.2))
            d0 <- density(C0Samples[,k], 0.05, from = a, to = b)
            d1 <- density(C1Samples[,k], 0.05, from = a, to = b)
            dd <- density(baseDsamples[i,j,], 0.02, from = a, to = b)
            ds <- density(baseGsamples[i,j,], 0.02, from = a, to = b)
            if(FALSE) {
                dmplot(list(d0,d1,dd,ds), main = '', ylim = c(0, 12),
                   lty=c(2,1,1,1), lwd=c(2,2,2,2), col = c(gray(.5),1,2,4),
                   xlab = paste0("C(", ynames[i], ", ", ynames[j], ")"),
                   ylab = "Density")
            } else {
                dmplot(list(d1,dd,ds), main = '', ylim = c(0, 12),
                   lty=c(1,1,1), lwd=c(2,2,2), col = c(1,2,4),
                   xlab = paste0("C(", ynames[i], ", ", ynames[j], ")"),
                   ylab = "Density") 
            }
            rug(C1obs[i,j], 0.25, lwd = 2, lty = 2,
                col = c(2,4)[1+(il %in% BaseGrph$iLtheta)])
            par(xpd = FALSE)
            abline(v = 0, lty = 2, col = gray(0.5))
            abline(h = 7, lty = 3, col = gray(0.3))
            par(xpd = NA)
        }
        if(j==(p-1)) {
            if(i==2) {
                legend("topleft", c(as.expression(bquote(C(theta)~'|'~D[1])),
                                    "Prior, dense", "Prior, graph"), 
                       lty = 1, lwd = 2, col = c(1,2,4), bty = "n")
            }
        }
    }
}
dev.off()

system(paste("eog", flC1prior, "&"))

## use the other hospital data
## update data for LKJ model
Sdat12 <- Sdat2 <- Sdat1  
Sdat2$y <- as.matrix(dat2)  ## the 2nd dataset
Sdat2$n <- n2      ## number of observations on 2nd dataset
Sdat12$y <- rbind(dat1, dat2)
Sdat12$n <- n1 + n2

### sample from the LKJ (non-informative prior)
Ssampl2 <- sampling(
    object = cmpl_LKJ, 
    data = c(Sdat2, list(eta = 1)), 
    iter = 20500,
    warmup = 500,
    thin = 10,
    chains = 5
)

Ssampl12 <- sampling(
    object = cmpl_LKJ, 
    data = c(Sdat12, list(eta = 1)), 
    iter = 20500,
    warmup = 500,
    thin = 10,
    chains = 5
)

## STAN model with informative prior to
## set the prior for 'Lcorr'
Smodel0_PCdense <- stan_add(
    x = Smodel0,
    model = "pc_correl",
    name = "LCorr" ## the Cholesky of Corr
)
## add how to get 'Corr' from 'Lcorr'
Smodel_PCdense <- stan_add_code(
    x = Smodel0_PCdense,
    to_add = list(
        "generated quantities" = 
            "  corr_matrix[p] Corr;
  Corr = multiply_lower_tri_self_transpose(LCorr);\n")
)

cat(Smodel_PCdense)

## Compile the 'PC_dense' model
system.time(
    cmpl_dense <- stan_model(
        model_code = Smodel_PCdense, 
        model_name = "PC_dense"
    )
)


## Sparse model: add the prior definition for 'Corr'
Smodel0_sparse <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  vector[p] y[n];
  vector[p] mu;
}
transformed parameters {
  corr_matrix[p] Corr;
}
model {
  y ~ multi_normal(mu, Corr);
}
"
Smodel_sparse <- stan_add(
    x = Smodel0_sparse, 
    model = "graphpcor", 
    name = "Corr"
)

cat(Smodel_sparse)

## Compile the 'PC_sparse' model
system.time(
    cmpl_sparse <- 
        stan_model(
            model_code = Smodel_sparse, 
            model_name = "PC_sparse"
        )
)

## add the data for the dense case
Sdat2dense <- stan_add(
    x = Sdat2, model = BaseDens,
    lambda = cLambda, name = 'LCorr')

str(Sdat2dense)

Ssampl2dense <- sampling(
    object = cmpl_dense,
    data = Sdat2dense,
    iter = 20500,
    warmup = 500,
    thin = 10,
    chains = 5
)

## add the data for the sparse case
Sdat2sparse <- stan_add(
    x = Sdat2, model = BaseGrph,
    lambda = cLambda, name = "Corr")

Ssampl2sparse <- sampling(
    object = cmpl_sparse,
    data = Sdat2sparse,
    iter = 20500,
    warmup = 500,
    thin = 10,
    chains = 5
)

if(FALSE) {
    stan_trace(Ssampl2sparse, paste0("grpc_theta[", 1:5, "]"))
    stan_dens(Ssampl2sparse, paste0("grpc_theta[", 1:5, "]"))
    
    stan_trace(Ssampl2dense, rhonams)
    stan_trace(Ssampl2sparse, rhonams)
}

CLKJsamples <- mcmc(Reduce("cbind", extract(Ssampl2, rhonams)))
CLKJsamples12 <- mcmc(Reduce("cbind", extract(Ssampl12, rhonams)))
CDensSamples <- mcmc(Reduce("cbind", extract(Ssampl2dense, rhonams)))
CGrphSamples <- mcmc(Reduce("cbind", extract(Ssampl2sparse, rhonams)))

C2obs <- cor(dat2)
round(C2obs*100)

## visualize
flC1post <- "~/github/corr-priors/examples/figures/uci_C1posterior.png"
png(flC1post, 4500, 3000, res = 600)
k = 0; 
par(mfcol = c(p-1, p-1), mar = c(3,3,0.5,0.5),
    mgp = c(1.5,0.7,0), las = 1, bty = 'n')
for(j in 1:(p-1)) {
    for(i in 2:p) {
        il <- (j-1)*p + i
        if(j>=i) {
            plot(0, type = 'n', axes = FALSE,
                 xlab = '', ylab = '', main = '')
        } else {
            k <- k + 1
            mk <- c(C1obs[i,j], C2obs[i,j])
            a <- max(-1, min(0, mk-0.15))
            b <- min(1, max(0, mk+0.15))
            d01 <- density(C1Samples[, k], 0.02, from = a, to = b)
            d02 <- density(CLKJsamples[,k], 0.02, from = a, to = b)
            d012 <- density(CLKJsamples12[,k], 0.02, from = a, to = b)
            dd <- density(baseDsamples[i,j,], 0.02, from = a, to = b)
            d2d <- density(CDensSamples[,k], 0.02, from = a, to = b)
            ds <- density(baseGsamples[i,j,], 0.02, from = a, to = b)
            d2s <- density(CGrphSamples[,k], 0.02, from = a, to = b)
            if(FALSE) {
                dmplot(list(d01,d02,d012,dd,d2d,ds,d2s), main = '',
                       lty = c(2,1,3,2,1,2,1), lwd=c(2,2,2,2,2,2,2),
                       col = c(1,1,6,2,2,4,4),
                       xlab = paste0("C(", ynames[i], ", ", ynames[j], ")"))
            } else {
                dmplot(list(dd,d2d, ds,d2s), main = '',
                       lty = c(3,1,3,1), lwd=c(2,2,2,2),
                       col = c(2,2,4,4),
                       xlab = paste0("C(", ynames[i], ", ", ynames[j], ")"))
            }
            abline(v = 0, lty = 2, col = gray(0.5))
            rug(BaseDens$base[i,j], -0.05, lwd = 4, col = 2)
            rug(BaseGrph$base[i,j], -0.05, lwd = 2, col = 4)
            rug(C2obs[i,j], -0.1, lwd = 3, col = 6)
        }
        if(j==(p-1)) {
            if(i==2) {
                legend("topleft", c("Prior, dense", "Prior, graph",
                                    "Post., dense", "Post., graph"),
                       lty = c(3,3,1,1), lwd = 2, col = c(2,4,2,4),
                       bty = "n")
            }
        }
    }
}
dev.off()

system(paste("eog", flC1post, "&"))
