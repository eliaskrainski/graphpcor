### useful packages

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

### INLA setup
inla.setOption(
    num.threads = "6:1",
    safe = FALSE
)

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
pcprec <- list(prior = "pc.prec", param = c(1, 0.01))
m1f <- Obs ~
    f(id.area, model = "besag", graph = graphGermany)

vnames <- c("oral", "osph", "lary", "lung")

lres.m1 <-lapply(vnames, function(v) {
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
svcor.m1 <- scc.fn(sapply(
    lres.m1, function(r) r$summary.random$id$mean))
colnames(svcor.m1) <- rownames(svcor.m1) <-
    gsub("obs", "m1", rownames(svcor.obs))
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
### Model 2: tree model (tree in Fig. 8 of paper1) 
##########################################################

tree3 <- treepcor(
    p1 ~ p2 + c4,
    p2 ~ p3 + c2 + c3,
    p3 ~ c1)

### cgeneric tree model
tree3model <- cgeneric(
    model = tree3, 
    lambda = 5,
    sigma.prior.reference = c(1, 1, 1, 1),
    sigma.prior.probability = c(0.1, 0.1, 0.1, 0.1)
)

### The kronecker product model: besag with treepcor
R.spatial <- inla.as.sparse(
    Diagonal(n.areas, x = colSums(graphGermany)) - graphGermany)
## define the cgeneric 'besag' model
cBesag <- cgeneric(
    model = "generic0",
    R = R.spatial,  
    param = c(1, NA) ## Fix it to one (let variances in the graphical model)
)

str(cBesag$f$extraconstr)

tree3Besag <- kronecker(tree3model, cBesag)
str(tree3Besag$f$extraconstr)

m2f <- Obs ~ 0 + disease + 
    f(iddata, model = tree3Besag)

### fit the model
res.m2 <- inla(
    formula = m2f, 
    family = "poisson",
    data = ldata, 
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

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
graph5 <- graphpcor(
    Oral_1 ~ Osph_2 + Lary_3,
    Lung_4 ~ Osph_2 + Lary_3 + Oral_1
)
graph4 <- graphpcor(
    Oral_1 ~ Osph_2 + Lary_3,
    Lung_4 ~ Osph_2 + Lary_3
)

par(mfrow = c(1, 3), mar = c(0,0,0,0))
plot(tree3)
plot(graph5)
plot(graph4)

### Two cgeneric graphpcor models
g5model <- cgeneric(
    model = graph5, 
    lambda = 5,
    sigma.prior.reference = rep(1, K),
    sigma.prior.probability = rep(0.1, K)
)
g4model <- cgeneric(
    model = graph4, 
    lambda = 5,
    sigma.prior.reference = rep(1, K),
    sigma.prior.probability = rep(0.1, K)
)

### The kronecker product model
graph5Besag <- kronecker(g5model, cBesag)
graph4Besag <- kronecker(g4model, cBesag)

str(graph5Besag$f$extraconstr)

m3f <- Obs ~ 0 + disease + 
    f(iddata, model = graph5Besag)
m4f <- Obs ~ 0 + disease + 
    f(iddata, model = graph4Besag)

### fit the model with graph5
res.m3 <- inla(
    formula = m3f, 
    family = "poisson",
    data = ldata,
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

## fit the model with graph4
res.m4 <- inla(
    formula = m4f, 
    family = "poisson",
    data = ldata, 
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

## similar marginal std (see next plot)
res.m3$summary.hyperpar[1:K, c(1,2,3,5)]
res.m4$summary.hyperpar[1:K, c(1,2,3,5)]

edgl.all <- matrix(paste(rep(vnames, each = K), "~", vnames), K)
edgl.all

lg5 <- Laplacian(graph5)
lg5
lg4 <- Laplacian(graph4)
lg4

edg.g5 <- edgl.all[lower.tri(lg5) & (!is.zero(lg5))]
edg.g4 <- edgl.all[lower.tri(lg4) & (!is.zero(lg4))]

edg.g5
res.m3$summary.hyperpar[-(1:K), c(1,2,3,5)] ## edge (Lary ~ Osph) "non-significant"
edg.g4
res.m4$summary.hyperpar[-(1:K), c(1,2,3,5)]

par(mfrow = c(3,3), mar = c(3,3,1,0.5), mgp = c(2,0.5,0))
for(i in 1:K) {
    nami <- vnames[i]
    m3.i <- inla.tmarginal(exp, res.m3$internal.marginals.hyperpar[[i]])
    m4.i <- inla.tmarginal(exp, res.m4$internal.marginals.hyperpar[[i]]) 
    plot(m3.i, type = "l", bty = 'n', lwd = 2,
         xlim = range(m3.i[, 1], m4.i[, 1]),
         ylim = range(m3.i[, 2], m4.i[, 2]),
         xlab = bquote(sigma[.(nami)]), ylab = 'Density')
    lines(m4.i, lwd = 3, lty = 2)
}
for(i in (K+1):(K+5)) {
    edjnam <- edg.g5[i-K]
    m3.i <- inla.smarginal(res.m3$internal.marginals.hyperpar[[i]])
    k <- i - (i>(K+2))
    cat(c(i=i, k=k), "\n")
    if(i!=(K+3)) {
        m4.i <- inla.smarginal(res.m4$internal.marginals.hyperpar[[k]])
        xlm <- range(m3.i$x, m4.i$x)
        ylm <- range(m3.i$y, m4.i$y)
    } else {
        xlm <- range(m3.i$x)
        ylm <- range(m3.i$y)
    }
    plot(m3.i, type = "l", bty = 'n', lwd = 2,
         xlim = xlm, ylim = ylm, 
         xlab = edjnam, ylab = 'Density')
    if(i!=(K+3)) {
        lines(m4.i, lwd = 3, lty = 2)
    }
    abline(v=0, col=gray(0.5), lwd = 2, lty = 2)
}

## covariance, correlation and std for the fitted log(smr)
svcor.m3 <- scc.fn(
    vcov(graph5, theta = res.m3$mode$theta),
    mat = TRUE)
colnames(svcor.m3) <- rownames(svcor.m3) <-
    gsub("obs", "m3", rownames(svcor.obs))
svcor.m4 <- scc.fn(
    vcov(graph4, theta = res.m4$mode$theta),
    mat = TRUE)
colnames(svcor.m4) <- rownames(svcor.m4) <-
    gsub("obs", "m4", rownames(svcor.obs))

round(svcor.obs, 2) ## upper.tri as in Table 2 of Held et. al. (2005)
round(svcor.m1, 2)
round(svcor.m2, 2) 
round(svcor.m3, 2)
round(svcor.m4, 2)

if(FALSE) {
    
    tbs <- list(
        m0 = svcor.obs,
        m1 = svcor.m1,
        m2 = svcor.m2,
        m3 = svcor.m3,
        m4 = svcor.m4)
    for(i in 1:length(tbs))
        dimnames(tbs[[i]])  <- list(
            c("Oral", "Oesophagus", "Larynx", "Lung"),
            c("Oral", "Oesophagus", "Larynx", "Lung")
        )

    for(i in 1:length(tbs))
        tbs[[i]][lower.tri(tbs[[i]])] <- NA

    knitr::kable(tbs[[1]], digits = 2)

    knitr::kable(tb1, digits = 2)

    knitr::kable(tb2, digits = 2)

    knitr::kable(tb3, digits = 2)
    knitr::kable(tb4, digits = 2)
    
}

gcpo.tab <- data.frame(
    M1 = sapply(lres.m1, function(r)
        -sum(log(r$gcpo$gcpo))), 
    apply(cbind(M2 = res.m2$gcpo$gcpo,
                M3 = res.m3$gcpo$gcpo,
                M4 = res.m4$gcpo$gcpo), 2,
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

scc.m2 <- t(apply(th.m2.samples, 2, function(x)
    scc.fn(vcov(tree3,theta=x), mat = TRUE)))
scc.m3 <- t(apply(th.m3.samples, 2, function(x)
    scc.fn(vcov(graph5,theta=x), mat = TRUE)))
scc.m4 <- t(apply(th.m4.samples, 2, function(x)
    scc.fn(vcov(graph4,theta=x), mat = TRUE)))

par(mfrow = c(K, K), mar = c(3,2,2,0.5), mgp = c(2,0.5,0), las = 1)
for(i in 1:K) {
    nami <- vnames[i]
    for (j in 1:K) {
        namj <- vnames[j]
        k <- (i-1)*K + j
        bk <- pretty(c(scc.m2[, k], scc.m3[, k], scc.m4[, k]), 20)
        hist(scc.m2[, k], bk, freq = FALSE, ylab = '',
             xlab = ifelse(i==j,
                           as.expression(bquote(sigma[.(nami)])),
                           ifelse(i<j,
                                  paste0("Cor(", nami, ", ", namj, ")"), 
                                  paste0("Cov(", nami, ", ", namj, ")"))), 
             col = gray(0.4, 0.5), border = 'transparent',
             main = ifelse(i==j, vnames[i], ""))
        hist(scc.m3[, k], bk, freq = FALSE, add = TRUE, 
             col = rgb(1, 0.5, 0.3, 0.5), border = 'transparent')
        hist(scc.m4[, k], bk, freq = FALSE, add = TRUE, 
             col = rgb(0.3, 0.5, 1, 0.5), border = 'transparent')
    }
}
legend("topright", c("tree", "graph 5", "graph 4"),
       fill = c(gray(0.4), rgb(1,0.5,0.3), rgb(0.3,0.5,1)),
       bty = 'n', border = 'transparent', cex = 2)
