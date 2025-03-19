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
    disease = rep(vnames, each = n.areas),
    id.spatial = rep(Germany4$id, K),
    id.disease = rep(1:K, each = n.areas),
    Obs = unlist(dataf[paste0(vnames, "_obs")]),
    E = unlist(dataf[paste0(vnames, "_exp")])
)


##########################################################
### Model 2: tree model (tree in Fig. 8 of paper1) 
##########################################################

tree <- treepcor(
    p1 ~ p2 + c4,
    p2 ~ p3 + c2 + c3,
    p3 ~ c1)

### cgeneric tree model
treemodel <- cgeneric(
    model = tree, 
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

treeBesag <- kronecker(treemodel, cBesag)
str(treeBesag$f$extraconstr)

m2f <- Obs ~ 0 + disease + 
    f(iddata, model = treeBesag)

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
    vcov(tree, theta = res.m2$mode$theta),
    mat = TRUE)
colnames(svcor.m2) <- rownames(svcor.m2) <-
    gsub("obs", "m2", rownames(svcor.obs))


####################################################
### Model 3: besag with graphpcor
####################################################

## graphical model for correlation
graph <- graphpcor(
    Oral_1 ~ Osph_2 + Lary_3,
    Osph_2 ~ Lary_3,
    Lary_3 ~ Lung_4
)

par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(tree)
plot(graph)

### cgeneric graphpcor model
gmodel <- cgeneric(
    model = graph, 
    lambda = 5,
    sigma.prior.reference = c(1, 1, 1, 1),
    sigma.prior.probability = c(0.1, 0.1, 0.1, 0.1)
)

### The kronecker product model
graphBesag <- kronecker(gmodel, cBesag)
str(graphBesag$f$extraconstr)

m3f <- Obs ~ 0 + disease + 
    f(iddata, model = graphBesag)

### fit the model
res.m3 <- inla(
    formula = m3f, 
    family = "poisson",
    data = ldata, 
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

## covariance, correlation and std for the fitted log(smr)
svcor.m3 <- scc.fn(
    vcov(graph, theta = res.m3$mode$theta),
    mat = TRUE)
colnames(svcor.m3) <- rownames(svcor.m3) <-
    gsub("obs", "m3", rownames(svcor.obs))

round(svcor.obs, 2) ## upper.tri as in Table 2 of Held et. al. (2005)
round(svcor.m1, 2)
round(svcor.m2, 2) 
round(svcor.m3, 2)

if(FALSE) {
    
    tb0 <- svcor.obs
    tb1 <- svcor.m1
    tb2 <- svcor.m2
    tb3 <- svcor.m3
    dimnames(tb0) <- dimnames(tb1) <- dimnames(tb2) <-  
        dimnames(tb3) <- list(
            c("Oral", "Oesophagus", "Larynx", "Lung"),
            c("Oral", "Oesophagus", "Larynx", "Lung")
        )

    knitr::kable(tb0, digits = 2)

    knitr::kable(tb1, digits = 2)

    knitr::kable(tb2, digits = 2)

    knitr::kable(tb3, digits = 2)
    
}

gcpo.tab <- data.frame(
    M1 = sapply(lres.m1, function(r)
        -sum(log(r$gcpo$gcpo))), 
    apply(cbind(M2 = res.m2$gcpo$gcpo,
                M3 = res.m3$gcpo$gcpo), 2,
          tapply, rep(vnames, each = n.areas),
          function(x) -sum(log(x)))
)

gcpo.tab
gcpo.tab - apply(gcpo.tab, 1, min)

### correlations uncertainty
nsampl <- 3e3
th.m2.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m2, intern = TRUE, improve = TRUE))
th.m3.samples <- t(inla.hyperpar.sample(
    n = nsampl, res.m3, intern = TRUE, improve = TRUE))

scc.fn(vcov(tree, theta = th.m2.samples[, 1]), mat = TRUE)
scc.fn(vcov(graph, theta = th.m3.samples[, 1]), mat = TRUE)

scc.m3 <- scc.m2 <- matrix(NA, nsampl, 4*4)

for(i in 1:nsampl) {
    scc.m2[i, ] <- scc.fn(
        vcov(
            object = tree,
            theta = th.m2.samples[,i]),
        mat = TRUE)
}

for(i in 1:nsampl) {
    scc.m3[i, ] <- scc.fn(
        vcov(
            object = graph,
            theta = th.m3.samples[,i]),
        mat = TRUE)
}

scc2.li <- matrix(apply(scc.m2, 2, quantile, 0.025), 4)
scc2.ls <- matrix(apply(scc.m2, 2, quantile, 0.975), 4)

scc3.li <- matrix(apply(scc.m3, 2, quantile, 0.025), 4)
scc3.ls <- matrix(apply(scc.m3, 2, quantile, 0.975), 4)

round(scc2.li, 2)
round(svcor.m2, 2)
round(scc2.ls, 2)

round(scc3.li, 2)
round(svcor.m3, 2)
round(scc3.ls, 2)
