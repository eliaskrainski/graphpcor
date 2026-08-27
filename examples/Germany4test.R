### packages 
library(graphpcor)
library(INLA)

### Load the data
data('Germany4')
names(Germany4)

## extract the data from geometry
dataf <- sf::st_drop_geometry(Germany4)

## size of the correlation model
(K <- 4)
(n.areas <- nrow(dataf))

### Models 1: scaled Besag model for each disease
pcprec <- list(prior = "pc.prec", param = c(1, 0.01))
m1f <- Obs ~
    f(id.area, model = "besag", graph = graphGermany,
      scale.model = TRUE, hyper = list(theta = pcprec))

vnames <- c("oral", "osph", "lary", "lung")

lres.m1 <- lapply(vnames, function(v) {
    inla(formula = m1f, 
         family = "poisson",
         data = data.frame(
             id.area = Germany4[['id']],
             Obs = Germany4[[paste0(v, "_obs")]],
             expected = Germany4[[paste0(v, "_exp")]]), 
         E = expected
         )
})

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
    lambda = 5, ##useINLAprecomp = FALSE,
    sigma.prior.reference = c(1, 1, 1, 1),
    sigma.prior.probability = c(0.1, 0.1, 0.1, 0.1)
)

### The kronecker product model: besag with treepcor
R.spatial <- inla.as.sparse(
    Diagonal(n.areas, x = colSums(graphGermany)) - graphGermany)

## define the cgeneric 'besag' model with R.spatial as the precision structure
cBesag <- cgeneric(
    model = "generic0",
    R = R.spatial, ## useINLAprecomp = FALSE,
    param = c(1, NA) ## Fix it to one (let variances in the graphical model)
)

str(cBesag$f$extraconstr)

## the Kroneker between the tree model and the Besag
tree3Besag <- kronecker(tree3model, cBesag)##, useINLAprecomp = FALSE)

str(tree3Besag$f$extraconstr)

m2f <- Obs ~ 0 + disease + 
    f(iddata, model = tree3Besag)

res.m2 <- inla(
    formula = m2f, 
    family = "poisson",
    data = ldata, 
    E = E,
    verbose = TRUE)

##inla.setOption(inla.call = "/home/eliask/.cache/R/INLA/stiles-binary/v26.08.20/bin/inla.run")

##inla.setOption(inla.call = "/home/eliask/.cache/R/INLA/stiles-binary/latest/bin/inla.run")

inla.setOption(smtp = "stiles")

### fit the model
res.m2 <- inla(
    formula = m2f, 
    family = "poisson",
    data = ldata, 
    E = E,
    verbose = TRUE,
    ##inla.call = "remote"##
    inla.call ="/home/eliask/.cache/R/INLA/stiles-binary/latest/bin/inla.run"
)

