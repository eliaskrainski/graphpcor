### fit some models to model 4 diseaes in Germany

library(graphpcor)
library(INLA)

### function to scale covariance matrix as
##    sqrt(variance) in the diagonal
###   correlation at the off diagonal
scc.fn <- function(x, mat = FALSE) {
    if(mat) {
        r <- x
    } else {
        r <- cov(x)
    }
    r[upper.tri(r)] <- cov2cor(r)[upper.tri(r)]
    diag(r) <- sqrt(diag(r))
    return(r)
}

### inla setup controls
ctrc <- list(
    dic = TRUE, waic = TRUE, cpo = TRUE,
    control.gcpo = list(
        enable = TRUE,
        num.level.sets = 5,
        strategy = "posterior")
)
ctri <- list(
    int.strategy = "ccd"
)

### Load the data
data('Germany4')
names(Germany4)

## extract the data from geometry
dataf <- sf::st_drop_geometry(Germany4)

## the observed smr
smr.obs <- dataf[, c(2,4,6,8)]/
    dataf[, c(3,5,7,9)]
summary(smr.obs)

## covariance, correlation and std for the observed smr
svcor.obs <- scc.fn(smr.obs)
round(svcor.obs, 2) ## upper.tri as in Table 2 of Held et. al. (2005)

## size of the correlation model
(K <- ncol(smr.obs))
(n.areas <- nrow(smr.obs))

### Models 1: scaled Besag model for each disease
m1f <- Obs ~
    f(id.area, model = "besag", graph = graphGermany,
      scale.model = TRUE,
      hyper = list(theta = list(prior = "pc.prec", param = c(0.5, 0.05))))

vnames <- c("oral", "osph", "lary", "lung"
            ); names(vnames) <- vnames

lres.m1 <- lapply(vnames, function(v) {
    inla(formula = m1f, 
         family = "poisson",
         data = data.frame(
             id.area = Germany4[['id']],
             Obs = Germany4[[paste0(v, "_obs")]],
             expected = Germany4[[paste0(v, "_exp")]]), 
         E = expected, 
         control.inla = ctri,
         control.compute = ctrc
         )
})

## covariance, correlation and std for the fitted log(smr)
smr_fit1 <- sapply(lres.m1, function(r)
    r$summary.random$id$mean)
svcor.m1 <- scc.fn(smr_fit1)
round(svcor.m1, 2)

### long data format
ldata <- data.frame(
    iddata = 1:(K * n.areas),
    disease = factor(rep(vnames, each = n.areas), labels = vnames),
    id.spatial = rep(Germany4$id, K),
    id.disease = rep(1:K, each = n.areas),
    Obs = unlist(dataf[paste0(vnames, "_obs")]),
    E = unlist(dataf[paste0(vnames, "_exp")])
)

##########################################################
### Model 2: correlated disease models 
##########################################################
cgWishart <- cgeneric(
    model = "Wishart", n = K, dof = K + 2,
    R = c(rep(1, K), diag(K)[lower.tri(diag(K))]))
cgLKJ <- cgeneric(
    model = "LKJ", n = K, eta = 2,
    sigma.prior.reference = rep(0.5, K),
    sigma.prior.probability = rep(0.05, K))
cgPC <- cgeneric(
    model = "pc_correl", n = K, lambda = 1, 
    sigma.prior.reference = rep(0.5, K),
    sigma.prior.probability = rep(0.05, K))

### The cgeneric Besag model definition 
R.spatial <- inla.as.sparse(
    Diagonal(n.areas, x = colSums(graphGermany)) - graphGermany)

## define the cgeneric 'besag' model with R.spatial as the precision structure
cBesag <- cgeneric(
    model = "generic0",
    R = R.spatial,  
    param = c(1, NA) ## Fix it to one (let variances in the graphical model)
)

str(cBesag$f$extraconstr)

cstds <- cgeneric("stds", n=K, lambda = 1,
                  sigma.prior.reference = rep(0.5, K),
                  sigma.prior.probability = rep(0.05, K))
brepl <- kronecker(cstds, cBesag)
bWishart <- kronecker(cgWishart, cBesag)
bLKJ <- kronecker(cgLKJ, cBesag)
bPC <- kronecker(cgPC, cBesag)

ffl <- list(
    rpl = Obs ~ f(iddata, model = brepl) + disease-1,
   ## Wish = Obs ~ f(iddata, model = bWishart) + disease-1,
    LKJ = Obs ~ f(iddata, model = bLKJ) + disease-1,
    PC = Obs ~ f(iddata, model = bPC) + disease-1)

inla.setOption(inla.call = "/home/eliask/.cache/R/INLA/stiles-binary/v26.08.20/bin/inla.run")
inla.setOption(smtp = "stiles")

fits <- lapply(ffl, function(ff) {
    print(ff)
    t0 <- Sys.time()
    o <- inla(
        formula = ff, 
        family = "poisson",
        data = ldata, 
        E = E,
        verbose = TRUE,
        control.compute = ctrc,
       control.inla = ctri
    )
    print(Sys.time()-t0)
    o
})

sapply(fits, function(r) r$cpu.used)

sapply(fits, function(r) r$misc$nfunc)

c(sum(sapply(lres.m1, function(x) x$dic$dic)),
  sapply(fits, function(x) x$dic$dic)
  )

sapply(fits, function(x) sum(x$gcv$gcv))

## covariance, correlation and std for the fitted log(smr)
svcor.m2 <- scc.fn(
    vcov(tree3, theta = res.m2$mode$theta),
    mat = TRUE)
colnames(svcor.m2) <- rownames(svcor.m2) <-
    gsub("obs", "m2", rownames(svcor.obs))


####################################################
### Model 3: besag with graphpcor
####################################################

## graphical model for correlation
graph4 <- graphpcor(
    Oral ~ Osph + Lary,
    Lung ~ Osph + Lary
)
graph5 <- graphpcor(
    Oral ~ Osph + Lary + Lung,
    Lung ~ Osph + Lary
)
graph5b <- graphpcor(
    Oral ~ Osph + Lary,
    Lung ~ Osph + Lary,
    Osph ~ Lary
)

par(mfrow = c(2, 2), mar = c(0,0,0,0))
plot(tree3)
plot(graph4)
plot(graph5)
plot(graph5b)

nsim1 <- 500
il4 <- which(lower.tri(diag(4)))
names(il4) <- c(paste0("Or_", c("Os", "La", "Lu")),
                paste0("Os_", c("La", "Lu")), "La_Lu")
vc3 <- sapply(1:nsim1, function(x)
    cov2cor(vcov(tree3, theta = runif(3, -3, 3)))[il4])
vc4 <- sapply(1:nsim1, function(x)
    cov2cor(vcov(graph4, theta = runif(4, -3, 3)))[il4])
vc5 <- sapply(1:nsim1, function(x)
    cov2cor(vcov(graph5, theta = runif(5, -3, 3)))[il4])
vc5b <- sapply(1:nsim1, function(x)
    cov2cor(vcov(graph5b, theta = runif(5, -3, 3)))[il4])

par(mfcol = c(4, 6), mar = c(3.5,3.5,1.5,0.5), mgp = c(2,0.5,0), bty = "n")
for(k in 1:6) {
    hist(vc3[k, ], main = paste("Tree:", names(il4)[k]),
         xlab = 'Correlation', freq = FALSE)
    hist(vc4[k, ], main = paste("Graph 4:", names(il4)[k]),
         xlab = 'Correlation', freq = FALSE)
    hist(vc5[k, ], main = paste("Graph 5:", names(il4)[k]),
         xlab = 'Correlation', freq = FALSE)
    hist(vc5b[k, ], main = paste("Graph 5b:", names(il4)[k]),
         xlab = 'Correlation', freq = FALSE)
}

### Two cgeneric graphpcor models
g4model <- cgeneric(
    model = graph4, 
    lambda = 1,
    sigma.prior.reference = rep(1, K),
    sigma.prior.probability = rep(0.1, K)
)
g5model <- cgeneric(
    model = graph5, 
    lambda = 1,
    sigma.prior.reference = rep(1, K),
    sigma.prior.probability = rep(0.1, K)
)
g5bmodel <- cgeneric(
    model = graph5b, 
    lambda = 1,
    sigma.prior.reference = rep(1, K),
    sigma.prior.probability = rep(0.1, K)
)

### The kronecker product model
graph4Besag <- kronecker(g4model, cBesag)
graph5Besag <- kronecker(g5model, cBesag)
graph5bBesag <- kronecker(g5bmodel, cBesag)

str(graph5Besag$f$extraconstr)

m3f <- Obs ~ 0 + disease + 
    f(iddata, model = graph4Besag)
m4f <- Obs ~ 0 + disease + 
    f(iddata, model = graph5Besag)
m5f <- Obs ~ 0 + disease + 
    f(iddata, model = graph5bBesag)

### fit the model with graph4
res.m3 <- inla(
    formula = m3f, 
    family = "poisson",
    data = ldata,
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

## fit the model with graph5
res.m4 <- inla(
    formula = m4f, 
    family = "poisson",
    data = ldata, 
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

## fit the model with graph5b
res.m5 <- inla(
    formula = m5f, 
    family = "poisson",
    data = ldata, 
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

## similar marginal std (see next plot)
res.m3$summary.hyperpar[1:K, c(1,2,3,5)]
res.m4$summary.hyperpar[1:K, c(1,2,3,5)]
res.m5$summary.hyperpar[1:K, c(1,2,3,5)]

edg.all <- matrix(
    paste(rep(vnames, each = K), "~", vnames),
    K)
edg.l <- edg.all[lower.tri(diag(K))]
edg.all
edg.l

lg4 <- Laplacian(graph4)
lg4
lg5 <- Laplacian(graph5)
lg5
lg5b <- Laplacian(graph5b)
lg5b

edg.g4 <- edg.all[lower.tri(lg4) & (!is.zero(lg4))]
edg.g5 <- edg.all[lower.tri(lg5) & (!is.zero(lg5))]
edg.g5b <- edg.all[lower.tri(lg5b) & (!is.zero(lg5b))]

edg.g4
res.m3$summary.hyperpar[-(1:K), c(1,2,3,5)]
edg.g5
res.m4$summary.hyperpar[-(1:K), c(1,2,3,5)] ## edge (Lary ~ Osph) "non-significant"
edg.g5b
res.m5$summary.hyperpar[-(1:K), c(1,2,3,5)] ## !!!

ii <- list(1:4, c(1:2,4:5,3), c(1:2, 4:5,3))
ii
ije <- list(pmatch(edg.g4, edg.all[lower.tri(lg4)]),
           pmatch(edg.g5, edg.all[lower.tri(lg5)]),
           pmatch(edg.g5b, edg.all[lower.tri(lg5b)]))
ije

par(mfrow = c(2,5), mar = c(3,3,1,0.5), mgp = c(2,0.5,0))
for(i in 1:K) {
    nami <- vnames[i]
    m3.i <- inla.tmarginal(exp, res.m3$internal.marginals.hyperpar[[i]])
    m4.i <- inla.tmarginal(exp, res.m4$internal.marginals.hyperpar[[i]])
    m5.i <- inla.tmarginal(exp, res.m5$internal.marginals.hyperpar[[i]]) 
    plot(m3.i, type = "l", bty = 'n', lwd = 2,
         xlim = range(m3.i[, 1], m4.i[, 1], m5.i[,1]),
         ylim = range(m3.i[, 2], m4.i[, 2], m5.i[,2]),
         xlab = bquote(sigma[.(nami)]), ylab = 'Density')
    lines(m4.i, lwd = 3, lty = 2)
    lines(m5.i, lwd = 4, lty = 3)
}
for(i in 1:4) {
    edjnam <- edg.l[ije[[1]][i]]
    m3.i <- inla.smarginal(res.m3$internal.marginals.hyperpar[[K+ii[[1]][i]]])
    m4.i <- inla.smarginal(res.m4$internal.marginals.hyperpar[[K+ii[[2]][i]]])
    m5.i <- inla.smarginal(res.m5$internal.marginals.hyperpar[[K+ii[[3]][i]]])
    xlm <- range(m3.i$x, m4.i$x, m5.i$x)
    ylm <- range(m3.i$y, m4.i$y, m5.i$y)
    plot(m3.i, type = "l", bty = 'n', lwd = 2,
         xlim = xlm, ylim = ylm, 
         xlab = edjnam, ylab = 'Density')
    lines(m4.i, lwd = 3, lty = 2)
    lines(m5.i, lwd = 4, lty = 3)
    abline(v=0, col=gray(0.5), lwd = 2, lty = 2)
}
plot(inla.smarginal(res.m4$internal.marginals.hyperpar[[K+ii[[2]][5]]]),
     type = "l", bty = "n", lwd = 3, lty = 2,
     xlab = edg.l[ije[[2]][3]], ylab = "Density")
abline(v=0, col=gray(0.5), lwd = 2, lty = 2)
plot(inla.smarginal(res.m5$internal.marginals.hyperpar[[K+ii[[3]][5]]]),
     type = "l", bty = "n", lwd = 4, lty = 3,
     xlab = edg.l[ije[[3]][3]], ylab = "Density")
abline(v=0, col=gray(0.5), lwd = 2, lty = 2)


## covariance, correlation and std for the fitted log(smr)
svcor.m3 <- scc.fn(
    vcov(graph4, theta = res.m3$mode$theta),
    mat = TRUE)
colnames(svcor.m3) <- rownames(svcor.m3) <-
    gsub("obs", "m3", rownames(svcor.obs))
svcor.m4 <- scc.fn(
    vcov(graph5, theta = res.m4$mode$theta),
    mat = TRUE)
colnames(svcor.m4) <- rownames(svcor.m4) <-
    gsub("obs", "m4", rownames(svcor.obs))
svcor.m5 <- scc.fn(
    vcov(graph5b, theta = res.m5$mode$theta),
    mat = TRUE)
colnames(svcor.m5) <- rownames(svcor.m5) <-
    gsub("obs", "m5", rownames(svcor.obs))

round(svcor.obs, 2) ## upper.tri as in Table 2 of Held et. al. (2005)
round(svcor.m1, 2)
round(svcor.m2, 2) 
round(svcor.m3, 2)
round(svcor.m4, 2)
round(svcor.m5, 2)

if(FALSE) {
    
    tbs <- list(
        m0 = svcor.obs,
        m1 = svcor.m1,
        m2 = svcor.m2,
        m3 = svcor.m3,
        m4 = svcor.m4,
        m5 = svcor.m5)
    for(i in 1:length(tbs))
        dimnames(tbs[[i]])  <- list(
            c("Oral", "Oesophagus", "Larynx", "Lung"),
            c("Oral", "Oesophagus", "Larynx", "Lung")
        )
    for(i in 1:length(tbs))
        tbs[[i]][lower.tri(tbs[[i]])] <- NA

    knitr::kable(tbs[[1]], digits = 2)

    knitr::kable(tbs[[2]], digits = 2)

    knitr::kable(tbs[[3]], digits = 2)

    knitr::kable(tbs[[4]], digits = 2)
    knitr::kable(tbs[[5]], digits = 2)
    knitr::kable(tbs[[6]], digits = 2)
    
}

gcpo.tab <- data.frame(
    M1 = sapply(lres.m1, function(r)
        -sum(log(r$gcpo$gcpo))), 
    apply(cbind(M2 = res.m2$gcpo$gcpo,
                M3 = res.m3$gcpo$gcpo,
                M4 = res.m4$gcpo$gcpo,
                M5 = res.m5$gcpo$gcpo), 2,
    tapply, rep(vnames, each = n.areas),
          function(x) -sum(log(x)))
)

gcpo.tab

## gcpo difference to the lowest
gcpo.tab - apply(gcpo.tab, 1, min)

### correlations uncertainty
nsampl <- 3e3
th.m2.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m2, intern = TRUE, improve = TRUE))
th.m3.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m3, intern = TRUE, improve = TRUE))
th.m4.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m4, intern = TRUE, improve = TRUE))
th.m5.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m5, intern = TRUE, improve = TRUE))

scc.m2 <- t(apply(th.m2.samples, 2, function(x)
    scc.fn(vcov(tree3, theta=x), mat = TRUE)))
scc.m3 <- t(apply(th.m3.samples, 2, function(x)
    scc.fn(vcov(graph4, theta=x), mat = TRUE)))
scc.m4 <- t(apply(th.m4.samples, 2, function(x)
    scc.fn(vcov(graph5, theta=x), mat = TRUE)))
scc.m5 <- t(apply(th.m5.samples, 2, function(x)
    scc.fn(vcov(graph5b, theta=x), mat = TRUE)))


par(mfrow = c(K, K), mar = c(3,2,2,0.5), mgp = c(2,0.5,0), las = 1)
for(i in 1:K) {
    nami <- vnames[i]
    for (j in 1:K) {
        if(j<i) {
            if(i==4) {
                if(j==1) {
                    plot(graph4)
                    legend("topleft", "graph4", bty = "n")
                }
                if(j==2) {
                    plot(graph5)
                    legend("topleft", "graph 5", bty = "n")
                }
                if(j==3) {
                    plot(graph5b)
                    legend("topleft", "graph 5b", bty = "n")
                }
            } else {
                plot(0, type = "h", bty = "n", axes = FALSE, xlab = "", ylab = "")
            }                
        } else {
            namj <- vnames[j]
            k <- (j-1)*K + i
            bk <- pretty(c(scc.m3[, k], scc.m4[, k], scc.m5[, k]), 20)
            hist(scc.m3[, k], bk, freq = FALSE, ylab = '',
                 xlab = ifelse(i==j,
                               as.expression(bquote(sigma[.(nami)])),
                        ifelse(i<j,
                               paste0("Cor(", nami, ", ", namj, ")"), 
                               paste0("Cov(", nami, ", ", namj, ")"))), 
                 dens = 20, #col = gray(0.4, 0.5), border = 'transparent',
                 main = ifelse(i==j, vnames[i], ""))
            hist(scc.m4[, k], bk, freq = FALSE, add = TRUE, 
                 col = rgb(1, 0.5, 0.3, 0.5), border = 'transparent')
            hist(scc.m5[, k], bk, freq = FALSE, add = TRUE, 
                 col = rgb(0.3, 0.5, 1, 0.5), border = 'transparent')
        }
        if((i==3) & (j==2)) {
            legend("topleft", c("graph 4", "graph 5", "graph 5b"),
                   dens = c(20,NA,NA),
                   fill = c(1, rgb(1,0.5,0.3), rgb(0.3,0.5,1)),
                   bty = 'n', border = 'transparent', cex = 1.5)
        }
    }
}


nrepl <- 3e3
prec4 <- replicate(nrepl, cgeneric_Q(graph4, theta = rnorm(4))[il4])
prec5 <- replicate(nrepl, cgeneric_Q(graph5, theta = rnorm(5))[il4])
prec5b <- replicate(nrepl, cgeneric_Q(graph5b, theta = rnorm(5))[il4])

par(mfrow = c(2, 3), mar = c(4,4,1,1), mgp = c(3, 0.5, 0))
for(i in 1:6) {
    bk <- pretty(c(prec4[i,], prec5[i,], prec5b[i,]), n = 50)
    hist(prec4[i,], bk, xlab = '', ylab = '', main = edg.l[i], dens = 20)
    hist(prec5[i,], bk, add = TRUE, col = rgb(1,0.5,0.3,0.5), border = 'transparent')
    hist(prec5b[i,], bk, add = TRUE, col = rgb(0.3,0.5,1,0.5), border = 'transparent')
}
legend("topleft", c("graph 4", "graph 5", "graph 5b"),
       dens = c(20,NA,NA),
       fill = c(1, rgb(1,0.5,0.3), rgb(0.3,0.5,1)),
       bty = 'n', border = 'transparent', cex = 1.5)

sessionInfo()

