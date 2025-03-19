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

## size of the correlation model
(K <- ncol(smr.obs))
(n.areas <- nrow(smr.obs))

## long data format
ldata <- data.frame(
    iddata = 1:(K * n.areas),
    disease = rep(vnames, each = n.areas),
    id.spatial = rep(Germany4$id, K),
    id.disease = rep(1:K, each = n.areas),
    Obs = unlist(dataf[paste0(vnames, "_obs")]),
    E = unlist(dataf[paste0(vnames, "_exp")])
)


svcor.obs <- scc.fn(smr.obs)
round(svcor.obs, 2) ## upper.tri as in Table 2 of Held et. al. (2005)

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
dimnames(svcor.m1) <- dimnames(svcor.obs)
round(svcor.m1, 2)

## graphical model for correlation
graph <- graphpcor(
    Oral1 ~ Osph2 + Lary3,
    Lung4 ~ Osph2 + Lary3
)

plot(graph)

gmodel <- cgeneric(
    model = graph, 
    lambda = 5,
    sigma.prior.reference = c(1, 1, 1, 1),
    sigma.prior.probability = c(0.1, 0.1, 0.1, 0.1)
)

## Model 2: graph model without spatial correlation
m2f <- Obs ~ 0 + disease +
    f(id.disease, model = gmodel, replicate = id.spatial) 

res.m2 <- inla(
    formula = m2f,
    family = "poisson",
    data = ldata, 
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

res.m2$summary.hy

## covariance, correlation and std for the fitted log(smr)
svcor.m2 <- scc.fn(vcov(graph, theta = res.m2$mode$theta), mat = TRUE)
dimnames(svcor.m2) <- dimnames(svcor.obs)
round(svcor.m2, 2)

### Model 3: multivariate Besag model
## Laplacian
R.spatial <- inla.as.sparse(
    Diagonal(n.areas, x = colSums(graphGermany)) - graphGermany)

## define the cgeneric 'besag' model
## besag cgeneric model
cBesag <- cgeneric(
    model = "generic0",
    R = R.spatial,  
    param = c(1, NA) ## Fix it to one (let variances in the graphical model)
)

str(cBesag$f$extraconstr)

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
    control.mode = list(theta = res.m2$mode$theta, restart = TRUE),
    control.compute = ctrc,
    control.inla = ctri
)

## covariance, correlation and std for the fitted log(smr)
svcor.m3 <- scc.fn(vcov(graph, theta = res.m3$mode$theta), mat = TRUE)
dimnames(svcor.m3) <- dimnames(svcor.obs)
round(svcor.m3, 2)

#############################
### Model 4: tree model (tree in Fig. 8 of paper1) 
#############################

tree <- treepcor(p1 ~ p2 + c4, p2 ~ p3 + c2 + c3, p3 ~ c1)

par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(tree)
plot(graph)

### cgeneric tree model
treemodel <- cgeneric(
    model = tree, 
    lambda = 5,
    sigma.prior.reference = c(1, 1, 1, 1),
    sigma.prior.probability = c(0.1, 0.1, 0.1, 0.1)
)

### The kronecker product model
treeBesag <- kronecker(treemodel, cBesag)

m4f <- Obs ~ 0 + disease + 
    f(iddata, model = treeBesag)

### fit the model
res.m4 <- inla(
    formula = m4f, 
    family = "poisson",
    data = ldata, 
    E = E,
    control.compute = ctrc,
    control.inla = ctri
)

## covariance, correlation and std for the fitted log(smr)
svcor.m4 <- scc.fn(vcov(tree, theta = res.m4$mode$theta), mat = TRUE)
dimnames(svcor.m4) <- dimnames(svcor.obs)

round(svcor.obs, 2) ## upper.tri as in Table 2 of Held et. al. (2005)
round(svcor.m1, 2)
round(svcor.m2, 2) 
round(svcor.m3, 2)
round(svcor.m4, 2)

