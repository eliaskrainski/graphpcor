
if(FALSE) { ### just to avoid having to add dependencies we don't really need
    
    library(dplyr)
    library(correlation)
    library(see)
    library(ggraph)
    library(ggpubr)
    
    
    mc <- mtcars %>%
        correlation(partial = FALSE)
    pc <- mtcars %>%
        correlation(partial = TRUE)
    
    ggarrange(
        mc %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("blue", "green")),
        pc %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("blue", "green"))
    )

}

library(corGraphs)
library(INLA)

inla.setOption(
    num.threads = 6L)

n <- nrow(mtcars)
nc <- length(jjy <- 1:7)
jjx <- 8:11

xx <- as.matrix(mtcars[, 8:11])

mypc <- sapply(jjy, function(i) {
    sapply(jjy, function(j) {
        yi <- mtcars[, i]
        yj <- mtcars[, j]
        cor(resid(lm(yi ~ xx)),
            resid(lm(yj ~ xx)))
    })
}); rownames(mypc) <- colnames(mypc) <-
        colnames(mtcars)[jjy]

round(mypc * 100)

str(mtcars)

datax <- data.frame(
    kronecker(diag(nc), as.matrix(mtcars[, jjx])))
str(datax)

data0 <- data.frame(
    datax,
    iv = factor(rep(1:nc, each = n)),
    y = unlist(mtcars[, jjy]))

ff0 <- update(y ~ iv, paste(".~.+", paste(colnames(datax), collapse = "+")))
ff0

pprc <- list(prec = list(initial = 10, fixed = TRUE))
fit0 <- inla(
    formula = ff0,
    ##    family = "poisson",
    control.family=list(hyper = pprc),
    data = data0,
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE) ### if true prints looooooottttssss of details

round(fit0$summary.fix[, c(1, 2, 3, 5)], 2)

data1 <- data.frame(
    y = data0$y - fit0$summary.fitted.values$mean
)

summary(data1$y)
summary(matrix(data1$y, n))
round(100 * cov2cor(cov(matrix(data1$y, n))))
round(100 * mypc)

data1$iv <- data0$iv
data1$i <- rep(1:nc, each = n)
data1$r <- rep(1:n, nc)

summary(data1)

ff1 <- y ~ iv + f(i, model = cGmodel, replicate = r, vb.correct = FALSE)
ff1

## dag to plot 
mtcd1 <- list(
    p1 ~ p2 + p3 + p4 + c1,
    p2 ~ c7,
    p3 ~ p5 + p6 + c6,
    p4 ~ c5,
    p5 ~ c2,
    p6 ~ c3 + c4
)

d2plot <- GraphPlot(mtcd1, base=0)

par(mar = c(1, 1, 1, 1))
plot(d2plot$gr, nodeAttrs = d2plot$nAttrs)

mtcd <- list(
    p1 ~ p2 + p3 + p4 - c1,
    p2 ~ - c7,
    p3 ~ p5 + p6 + c6,
    p4 ~ -c5,
    p5 ~ c2,
    p6 ~ c3 + c4
)
(np <- length(mtcd))

cGmodel <- cgeneric_dag_model(
    dag = mtcd,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc),
    lambda = 1,
    iprior = 3,
    debug = 0
)

str(cGmodel, 5)

fit1 <- inla(
    formula = ff1,
    control.family=list(hyper = pprc), 
    data = data1,
    control.inla = list(int.strategy = "eb"),
    control.mode = list(theta = rep(c(-2, 0), c(nc, np)), restart = TRUE),
##    control.mode = list(theta = rep(0, nc+np), restart = FALSE, fixed = TRUE),
    verbose = !TRUE) ### if true prints looooooottttssss of details

fit1$cpu

fit1$mode$theta
fit1$mode$theta[nc+1:np]

plot(fit1, F, F, F, F, F, F, plot.opt.trace = TRUE)

cc.fit1 <- cov2cor(dag_covariance(
    mtcd, fit1$mode$theta[nc+1:np]))

round(mypc*100)
round(cc.fit1*100)

### joint
data2 <- data.frame(
    datax,
    data1
)

ff2 <- update(ff0, .~.+f(i, model = cGmodel, replicate = r, vb.correct = FALSE))

fit2 <- inla(
    formula = ff2,
    control.family=list(hyper = pprc), 
    data = data2,
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE) ### if true prints looooooottttssss of details

fit2$cpu

fit2$mode$theta
fit2$mode$theta[nc+1:np]

plot(fit2, F, F, F, F, F, F, plot.opt.trace = TRUE)

cc.fit2 <- cov2cor(dag_covariance(
    mtcd, fit2$mode$theta[nc+1:np]))

round(mypc*100)
round(cc.fit1*100)
round(cc.fit2*100)

ff3 <- update(ff0, .~.+f(i, model = "iidkd", order = nc, n = n*nc, vb.correct = FALSE)) 

lcc <- t(chol(solve(mypc)))
ini3 <- c(log(diag(lcc)), lcc[lower.tri(lcc)])
length(ini3)

fit3 <- inla(
    formula = ff3,
    control.family=list(hyper = pprc), 
    data = data2,
    control.mode = list(theta = ini3, restart = TRUE),
    control.inla = list(int.strategy = "eb"),
    control.compute = list(smtp = "pardiso"),
    inla.call = "remote",
    num.threads = 2*length(ini3),
    verbose = !TRUE
)

xxsamples <- inla.iidkd.sample(
    n = 3000,
    result = fit3,
    name = "i",
    return.cov = TRUE)

cvfit3 <- Reduce("+", xxsamples)/length(xxsamples)
cc.fit3 <- cov2cor(cvfit3)

round(mypc*100)
round(cc.fit1*100)
round(cc.fit2*100)
round(cc.fit3*100)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
