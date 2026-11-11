########################################################
#### Example for the 12 types of cancer considering
####  - prior from 540 Valencian municipalities
####  - data from 1191 municipalities from South Brazil

library(sf)
library(INLA)
library(graphpcor)

setwd(here::here("examples", "cancer12br4"))
getwd()

## read the dataset and map
dataf <- readRDS("rds/southBR12cancer.rds")

### the municipalities graph
graphfl <- paste0("graph", n.s)
if(!file.exists(graphfl)) {
    source("mapnb.R")
}

### the south Brazil municipalities map
smap <- readRDS("rds/sbrmap.rds")

head(dataf,2)
head(smap,2)

## prior specification
### build a base correlation matrix considering the
### - graph built from the Valencian data
### - baseC: base at the Valencian data correlation
source("cancer12prior.R")
ls()

### - base0:  base at I (no-correlation a priori)
base0 <- basepcor(diag(nc), iLtheta = cGraph)

### - base9: base with correlation near 11'
base9 <- basepcor(rep(-7, dim(cGraph)[2]), p = nc, iLtheta = cGraph)
mean(base9$base[lower.tri(diag(nc))])

### set the graphpcor cgeneric model for each of the 3 base models
### correlations equal zero, under graph
cgrpc0 <- cgeneric(
    model = base0,
    lambda = 1, 
    sigma.prior.reference = rep(0.5, nc),
    sigma.prior.probability = rep(0.05, nc))

### the full information from the Valencian data
cgrpci <- cgeneric(
    model = baseC,
    lambda = 1, 
    sigma.prior.reference = rep(0.5, nc),
    sigma.prior.probability = rep(0.05, nc))

### correlations near 1, under the graph
cgrpc9 <- cgeneric(
    model = base9,
    lambda = 1/2, 
    sigma.prior.reference = rep(0.5, nc),
    sigma.prior.probability = rep(0.05, nc))

### set a independency model but different standard deviations
cstds <- cgeneric(
    model = 'stds', n = nc,
    sigma.prior.reference = rep(0.5,nc),
    sigma.prior.probability = rep(0.05,nc),
    developing = TRUE)

## number of areas selected
(n.s <- nrow(dataf))

## build the spatial structure matrix
ireadgraph <- inla.read.graph(graphfl)
sGraph <- Sparse(sparseMatrix(
    i = rep(1:ireadgraph$n, ireadgraph$nnbs),
    j = unlist(ireadgraph$nbs), x = 1L
))
R0 <- Diagonal(n.s, colSums(sGraph)) - sGraph

## scale the spatial structure matrix
ctrn <- list(A = matrix(1, 1, n.s), e = 0)
Rscaled <- inla.scale.model(
    Q = R0, constr = ctrn)

## build the cgeneric Besag model
cbesag1 <- cgeneric(
    model = "generic0",
    R = Rscaled, scale = FALSE,
    constr = FALSE, extraconstr = ctrn,
    param = c(1, NA))

## build the Kronecker product model: independent (x) Besag
## for a replicated Besag model for each disease
replBesag <- kronecker(cstds, cbesag1, developing = TRUE)
stopifnot(length(cgeneric_initial(replBesag))==
          length(cgeneric_initial(cstds)))

## build the Kronecker product model: graphpcor (x) Besag
## for the graph based cross-correlation model and spatial Besag
grpc0Besag <- kronecker(cgrpc0, cbesag1, developing = TRUE)
grpciBesag <- kronecker(cgrpci, cbesag1, developing = TRUE)
grpc9Besag <- kronecker(cgrpc9, cbesag1, developing = TRUE)
stopifnot(length(cgeneric_initial(grpc0Besag))==
          length(cgeneric_initial(cgrpc0)))

## index for the disease in the data
head(dataf, 2)
causes_r <- causes[ireord]
jcause <- match(causes_r, colnames(dataf))
jcause

stopifnot(all(causes_r==names(dataf)[jcause]))

## build the long data for INLA
ldat <- data.frame(
    cause = rep(causes_r, each = n.s),
    icause = rep(1:nc, each = n.s),
    iarea = rep(1:n.s, nc),
    ijoint = 1:(nc*n.s), 
    E = as.vector(as.matrix(dataf[, -nc + jcause])),
    y = as.vector(as.matrix(dataf[, jcause]))
)
ldat$i2cause <- ldat$icause
ldat$i2area <- ldat$iarea

head(ldat,2)

(nthetas <- length(cgeneric_initial(cstds)) +
     length(cgeneric_initial(grpc0Besag)))
(nth <- c(available = max(parallel::detectCores()-4,4),
          optimal = nthetas*2, ## optimal considering biggest model
          ok = nthetas+1)) ## ok, if optimal is not possible

### the formulae for each of the models and priors
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

ccomp <- list(dic = TRUE, waic = TRUE, cpo = TRUE)

### file names for the fitted models
fitflk <- paste0("rds/fit_c", nc, names(ffl), "_n", n.s, ".rds")

fits <- vector("list", length(ffl))
names(fits) <- names(ffl)
for(k in 1:length(ffl)) {
    flk <- fitflk[k]
    if(file.exists(flk)) {
        cat("Loading model", names(ffl)[k], "\n")
        fits[[k]] <- readRDS(flk)
    } else {
        cat("Fitting model", names(ffl)[k], "\n")
        print(ffl[[k]])
        t0 <- Sys.time()
        fits[[k]] <- inla(
            formula = ffl[[k]],
            family = 'poisson',
            data = ldat, E = E,
            control.compute = ccomp,
            verbose = !TRUE,
            num.threads = nth[1]
##            num.threads = nth[3], inla.call = "remote"
        )
        print(Sys.time()-t0)
        cat("Evaluating model", names(ffl)[k], "\n")
        print(ffl[[k]])
        tmp <- fits[[k]]
        tmp$.args$verbose <- FALSE
        tmp$.args$control.compute$smtp <- "taucs"
        t0 <- Sys.time()
        fits[[k]]$gcpo <- inla.group.cv(
            tmp, num.level.sets = 10
        )
        print(Sys.time()-t0)
        saveRDS(fits[[k]], flk)
    }
}

sapply(fits, function(fit) length(fit$mode$theta))
## IndC CICC SStr RBym ICSC 
##   12   48   12   24   48 

sapply(fits, function(fit) fit$cpu.used)
##n.s=1191: Total  23.0754509 783.9336100 25.3827488 119.3854756 2894.2066107

sapply(fits, function(r) grep("valuations =", r$logfile, val = TRUE))
sapply(fits, function(fit) fit$misc$nfunc) ## Total
##n.s=1191:   751 18858   579  2047 23541

sapply(fits, function(fit) fit$cpu.used["Total"])/
    sapply(fits, function(fit) fit$misc$nfunc)
## 0.03072630 0.04157035 0.04383894 0.05832217 0.12294323 

names(fits[[1]]$gcpo)

tabf <- sapply(fits, function(r)
    c(DIC=sum(r$dic$local.dic),
      WAIC=sum(r$waic$local.waic),
      lCPO=-sum(log(r$cpo$cpo)),
      lGCPO=-sum(log(r$gcpo$cv))))
tabf
tabf-apply(tabf,1,min)

nc
cpars <- sapply(fits[c(2,3,6,7)], function(m)
    m$mode$theta[nc*2+1:dg[2]])
cor(cpars, baseC$theta)

fit <- fits[[length(fits)]]

data.frame(iid=exp(fit$mode$theta[1:nc]),
           besag=exp(fit$mode$theta[nc+1:nc]),
           frac=exp(fit$mode$theta[nc+1:nc]*2)/
               (exp(fit$mode$theta[1:nc] * 2) +exp(fit$mode$theta[nc+1:nc]*2)),
           row.names = causes)

vmode <- vcov(cGraph, theta = fit$mode$theta[-(1:nc)])
mctheta <- fit$mode$theta[nc*2 + 1:dg[2]]
cmode <- vcov(cGraph, theta = mctheta)
round(cmode*100)

sum(agree <- ((baseC$theta<0) & (mctheta<0)) |
        ((baseC$theta>=0) & (mctheta>=0)) )
which(!agree)

pbase <- sapply(1:nc, function(k)
    inla.pmarginal(baseC$theta[k],
                   fit$internal.marginals.hyperpar[[dg[1]+k]]))
round(pbase,2)
round(1-2*abs(pbase-0.5), 2)

dg

