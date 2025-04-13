
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

    mc.f <- mc[abs(mc$r)>0.7, ]
    pc.f <- pc[abs(pc$r)>0.3, ]

    mc.c <- mc
    mc.c$r <- mc$r*(abs(mc$r)>0.7)
    mc.c$CI_low <- mc$CI_low*(abs(mc$r)>0.7)
    mc.c$CI_high <- mc$CI_high*(abs(mc$r)>0.7)
    mc.c$t <- mc$t*(abs(mc$r)>0.7)

    pc.c <- pc
    pc.c$r <- pc$r*(abs(pc$r)>0.3)
    pc.c$CI_low <- pc$CI_low*(abs(pc$r)>0.3)
    pc.c$CI_high <- pc$CI_high*(abs(pc$r)>0.3)
    pc.c$t <- pc$t*(abs(pc$r)>0.3)

    ggarrange(
        mc %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("gray", "red")),
        mc.c %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("gray", "red")),
        mc.f %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("gray", "red")),
        pc %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("gray", "red")),
        pc.c %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("gray", "red")),
        pc.f %>%
        plot() +
        scale_edge_colour_gradientn(
            limits = c(-1, 1),
            colors = c("gray", "red"))
    )

}

library(INLA)
library(graphpcor)

inla.setOption(
    safe = FALSE,
    num.threads = 6L
)

mtcars <- as.data.frame(
    kronecker(matrix(1, 10, 1),
              as.matrix(mtcars))
)

n <- nrow(mtcars)
(nc <- length(jjy <- c(1, 3, 4, 6, 7)))
length(jjx <- c(2, 5, 8:11))

names(mtcars)[jjx]
names(mtcars)[jjy]

xx <- as.matrix(mtcars[, jjx])

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
    control.family=list(hyper = pprc),
    data = data0,
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE)

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

ff1 <- y ~ iv + f(i, model = gmodel, replicate = r, vb.correct = FALSE)
ff1

## graph to plot
mtcd1 <- list(
    p1 ~ p2 + p3 + c1,
    p2 ~ p4 + c4,
    p3 ~ c3,
    p4 ~ p5 + c2,
    p5 ~ c5
)
d2plot <- GraphPlot(mtcd1, base=0)

par(mar = c(1, 1, 1, 1))
plot(d2plot$gr, nodeAttrs = d2plot$nAttrs)

mtcd <- list(
    p1 ~ p2 + p3 - c1,
    p2 ~ p4 + c4,
    p3 ~ c3,
    p4 ~ p5 + c2,
    p5 ~ c5
)
(np <- length(mtcd))

gmodel <- dcg_model(
    dcg = mtcd,
    sigma.prior.reference = rep(1, nc),
    sigma.prior.probability = rep(0.1, nc),
    lambda = 1,
    debug = 0
)

str(gmodel, 5)

fit1 <- inla(
    formula = ff1,
    control.family=list(hyper = pprc),
    data = data1,
    control.inla = list(int.strategy = "eb"),
    control.mode = list(theta = rep(c(-2, 0), c(nc, np)), restart = TRUE),
##    control.mode = list(theta = rep(0, nc+np), restart = FALSE, fixed = TRUE),
    verbose = !TRUE) ### if true prints looooooottttssss of details

fit1$cpu.used

fit1$mode$theta
fit1$mode$theta[nc+1:np]

plot(fit1, F, F, F, F, F, F, plot.opt.trace = TRUE)

cc.fit1 <- cov2cor(dcg_covcov(
    mtcd, fit1$mode$theta[nc+1:np]))

round(mypc*100)
round(cc.fit1*100)

### joint
data2 <- data.frame(
    datax,
    data1
)

ff2 <- update(ff0, .~.+f(i, model = gmodel, replicate = r, vb.correct = FALSE))

fit2 <- inla(
    formula = ff2,
    control.family=list(hyper = pprc),
    data = data2,
    control.inla = list(int.strategy = "eb"),
    verbose = !TRUE) ### if true prints looooooottttssss of details

fit2$cpu.used

fit2$mode$theta
fit2$mode$theta[nc+1:np]

plot(fit2, F, F, F, F, F, F, plot.opt.trace = TRUE)

cc.fit2 <- cov2cor(dcg_covcov(
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
##    control.mode = list(theta = ini3, restart = TRUE),
    control.inla = list(int.strategy = "eb"),
#    control.compute = list(smtp = "pardiso"),
 #    inla.call = "remote",
  #   num.threads = 2*length(ini3),
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

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
