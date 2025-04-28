
library(INLA)
library(INLAjoint)
library(graphpcor)
library(JM) # This package contains the dataset

data(pbc2) # dataset

library(ggplot2)
library(ggpubr)

gg0 <- ggplot(data = pbc2)

ggarrange(
    gg0 + geom_line(aes(x=year, y = serBilir, group = id)),
    gg0 + geom_line(aes(x=year, y = serChol, group = id)),
    gg0 + geom_line(aes(x=year, y = albumin, group = id)),
    gg0 + geom_line(aes(x=year, y = alkaline, group = id)),
    gg0 + geom_line(aes(x=year, y = SGOT, group = id)),
    gg0 + geom_line(aes(x=year, y = platelets, group = id)),
    gg0 + geom_line(aes(x=year, y = prothrombin, group = id)),
    gg0 + geom_line(aes(x=year, y = histologic, group = id))
)

## extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age",
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])

f1 <- function(x) x^2

### first prepare data with correct format using dataOnly and full sample
run0 <- joint(
    formLong = list(serBilir ~ (1 + year + f1(year)) * drug + (1 + year + f1(year)| id),
                    platelets ~ (1 + year + f1(year)) * drug + (1 + year + f1(year)| id)),
    dataLong = Longi, id = "id", timeVar = "year",
    family = c("lognormal", "poisson"),
    control=list(int.strategy = "eb", cfg = TRUE),
    corLong = TRUE, run = FALSE)
run.kd <- joint.run(run0)
vars.kd <- summary(run.kd, sdcor=T)$ReffList[[1]]

## graph-based structure model for partial correlations
G <- graphpcor(a1 ~ a2 + a3 + b1, ## 1st intercept with 1st slope and 2nd intercept
               b1 ~ b2 + b3) ## 2nd intercept with 2nd slope
plot(G)

## the correlation model
gmodel1 <- cgeneric(
  model = G,
  sigma.prior.reference = rep(5, dim(G)[1]),
  sigma.prior.probability = rep(0.2, dim(G)[1]),
  lambda = 3
)

## because of the indexing we need to do a "Kronecker":
## within individual correlation model (X) between individuals "iid" model
iidmodel <- cgeneric(
    model = "iid",
    n = 312, ## number of individuals
    param = c(1, 0.0),
    useINLAprecomp=FALSE
)

## the model to actually use
gmodel <- kronecker(gmodel1, iidmodel, useINLAprecomp=TRUE)

## update the formula with the model
run.G <- run0
run0$.args$formula
run.G$.args$formula <-
    Yjoint ~ -1 + Intercept_L1 + year_L1 + f1year_L1 + drugDpenicil_L1 + 
    year.X.drugDpenicil_L1 + f1year.X.drugDpenicil_L1 + Intercept_L2 + 
    year_L2 + f1year_L2 + drugDpenicil_L2 + year.X.drugDpenicil_L2 + 
    f1year.X.drugDpenicil_L2 +
    f(IDIntercept_L1, WIntercept_L1, model = gmodel) +
    f(IDyear_L1, Wyear_L1, copy = "IDIntercept_L1") + 
    f(IDf1year_L1, Wf1year_L1, copy = "IDIntercept_L1") +
    f(IDIntercept_L2, WIntercept_L2, copy = "IDIntercept_L1") +
    f(IDyear_L2, Wyear_L2, copy = "IDIntercept_L1") +
    f(IDf1year_L2, Wf1year_L2, copy = "IDIntercept_L1")

##    Yjoint ~ -1 + Intercept_L1 + year_L1 + drugDpenicil_L1 + year.X.drugDpenicil_L1 +
  ##      Intercept_L2 + year_L2 + drugDpenicil_L2 + year.X.drugDpenicil_L2 +
    ##    f(IDIntercept_L1, WIntercept_L1, model = gmodel) +
      ##  f(IDyear_L1, Wyear_L1, copy = "IDIntercept_L1") +
        ##f(IDIntercept_L2, WIntercept_L2, copy = "IDIntercept_L1") +
        ##f(IDyear_L2, Wyear_L2, copy = "IDIntercept_L1")

run.G <- joint.run(run.G)

# save(JM_INLA_GRAPH, file="JM_INLA_GRAPH.RData")

vfit <- vcov(G, theta = run.G$mode$theta[-1])
round(cov2cor(vfit), 2)

hsamples <- inla.hyperpar.sample(n=1e4, result = run.G, intern = TRUE)
scsamples <- t(sapply(1:nrow(hsamples), function(i) {
    v <- vcov(G, theta = hsamples[i, -1])
    return(c(sqrt(diag(v)), cov2cor(v)[lower.tri(v)]))
}))

gsummary <- data.frame(
    mean = apply(scsamples, 2, mean),
    sd = apply(scsamples, 2, sd),
    q0.025 = apply(scsamples, 2, quantile, 0.025),
    q0.975 = apply(scsamples, 2, quantile, 0.975)
)

round(cbind(vars.kd[, 1:2], gsummary[, 1:2]), 4)

## "12", "13", "14", "15", "16",
##   "23", "24", "25", "26", 
##     "34", "35", "36"
##      "45", "46", "56"
ijc <- list(c(1,2), c(1,3), c(1,4), c(1,5), c(1,6),
            c(2,3), c(2,4), c(2,5), c(2,6),
            c(3,4), c(3,5), c(3,6), c(4,5), c(4,6), c(5,6))
nvars <- 6
ncors <- length(ijc)
npars <- nvars + ncors

par(mfrow = c(1, 2), mar = c(4,4,1,1), mgp = c(3,1.5,0), las = 1, bty = 'n')
plot(1:nvars-0.1, vars.kd[1:nvars, 1], pch = 19, axes = FALSE, log = "y",
     xlim = c(0.5, 4.5), ylim = range(vars.kd[1:nvars, 3:5], gsummary[1:nvars, 3:4]),
     xlab = '', ylab = expression(sigma))
axis(2)
axis(1, 1:nvars, as.expression(lapply(1:nvars, function(i) bquote(sigma[.(i)]))))
segments(1:nvars-0.1, vars.kd[1:nvars, 3], 1:nvars-0.1, vars.kd[1:nvars, 5])
points(1:nvars+0.1, gsummary[1:nvars, 1], col = 2, pch = 8)
segments(1:nvars+0.1, gsummary[1:nvars, 3], 1:nvars+0.1, gsummary[1:nvars, 4], col = 2)
plot((nvars+1):npars-0.1, vars.kd[(nvars+1):npars, 1], pch = 19, axes = FALSE,
     ylim = range(vars.kd[(nvars+1):npars, 3:4]), xlim = c(nvars-1, ncors+1), 
     xlab = '', ylab = 'Correlation')
axis(2)
axis(1, (nvars+1):npars, as.expression(lapply(ijc, function(ij)
    bquote(rho[.(ij[1])~.(ij[2])]))))
segments((nvars+1):npars-0.1, vars.kd[(nvars+1):npars, 3],
   (nvars+1):npars-0.1, vars.kd[(nvars+1):npars, 5])
points((nvars+1):npars+0.1,gsummary[(nvars+1):npars, 1], col = 2, pch = 8)
segments((nvars+1):npars+0.1, gsummary[(nvars+1):npars, 3],
         (nvars+1):npars+0.1, gsummary[(nvars+1):npars, 4], col = 2)
abline(h=0)
legend("bottomleft", c("IIDKD", "GraphPCor"), bty = "n",
       lty = 1, col = 1:2, pch = c(19,8))
