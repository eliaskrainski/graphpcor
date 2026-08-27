
## Install bnlearn package, if not yet installed
if (!requireNamespace("bnlearn", quietly = TRUE)) {
  install.packages("bnlearn")
}

library(bnlearn)

## 1. Load a real available continuous dataset 
## The 'marks' dataset contains exam scores (0-100)
##  for 88 students across 5 subjects
data(marks)
head(marks)

n <- nrow(marks)
p <- ncol(marks)

library(ggplot2)
library(GGally)

cor_ci <- function(data, mapping, method = "pearson",
                   use = "complete.obs", ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  test <- cor.test(x, y, method = method)
  
  label <- sprintf("r = %.2f\nCI = [%.2f, %.2f]",
                   test$estimate,
                   test$conf.int[1],
                   test$conf.int[2])
  scol <- ifelse(test$conf.int[2]<0, "blue",
                 ifelse(test$conf.int[1]>0, "red", "gray"))
  ggplot(data = data, mapping = mapping) +
    annotate("text", 
             x = mean(range(x, na.rm = TRUE)),
             y = mean(range(y, na.rm = TRUE)),
             label = label,
             color = scol,
             size = 4) +
    theme_void()
}

## Create customized pairs plot
ggpairs(
  marks, 
  upper = list(continuous = wrap("points", alpha = 0.6, size = 0.7)),
  diag  = list(continuous = wrap("barDiag", bins = 15, fill = "gray70")),
  lower = list(continuous = cor_ci)
)


## 2. Structure Learning
## Learn the network structure using the Hill Climbing (hc) algorithm
## For continuous data, this optimizes the
##  Bayesian Information Criterion (BIC) score

continuous_structure <- hc(marks)

## View the learned structure
print(continuous_structure)

plot(continuous_structure,
     main = "Continuous Bayesian Network Structure (Exam Marks)")

library(graphpcor)

arcs <- continuous_structure$arcs
arcs

gpc_arcs <- graphpcor(
    paste(arcs[, 1], arcs[, 2], sep = "~")
)
gpc_arcs
g0 <- attr(gpc_arcs, "graph")
g0

## have it in the same order as the columns of marks
xnames <- colnames(marks)
ordn <- match(xnames, colnames(g0))
g0[ordn, ordn]

gpc <- graphpcor(g0[ordn, ordn])
summary(gpc)

par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(continuous_structure,
     main = "Continuous Bayesian Network Structure (Exam Marks)")
plot(gpc, Rgraphviz = TRUE)

dgpc <- dim(gpc)
dgpc

c0a <- basepcor(rep(-1, dgpc[2]), iLtheta = gpc)
c0b <- basepcor(rep(-2, dgpc[2]), iLtheta = gpc)

round((cc <- cor(marks))*100)
round(c0a$base*100)
round(c0b$base*100)

clkjmodel0 <- cgeneric(
    model = "LKJ", n = p, eta = 1
)
clkjmodel1 <- cgeneric(
    model = "LKJ", n = p, eta = 5
)
cpcmodel0 <- cgeneric(
    model = "pc_correl", n = p,
    base = rep(0, p*(p-1)/2),
    lambda = 1
)
cpcmodel1 <- cgeneric(
    model = basecor(c0b$base),
    lambda = 1
)
cgmodel0 <- cgeneric(
    model = gpc,
    base = rep(0, dgpc[2]),
    lambda = 1
)
cgmodel1 <- cgeneric(
    model = c0b,
    lambda = 1
)

library(INLA)

ldat <- data.frame(
    y = as.vector(scale(marks)),
    i = rep(1:p, each = n),
    s = rep(1:n, p)
)

hlik <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))

fit_lkj0 <- inla(
    formula = y ~ 0 + f(i, model = clkjmodel0, replicate = s),
    data = ldat,
    control.family = hlik,
    control.mode = list(theta = rep(0, p*(p-1)/2))
)
fit_lkj1 <- inla(
    formula = y ~ 0 + f(i, model = clkjmodel1, replicate = s),
    data = ldat,
    control.family = hlik,
    control.mode = list(theta = rep(0, p*(p-1)/2))
)

fit_c0 <- inla(
    formula = y ~ 0 + f(i, model = cpcmodel0, replicate = s),
    data = ldat,
    control.family = hlik
)
fit_c1 <- inla(
    formula = y ~ 0 + f(i, model = cpcmodel1, replicate = s),
    data = ldat,
    control.family = hlik
)

fit_g0 <- inla(
    formula = y ~ 0 + f(i, model = cgmodel0, replicate = s),
    data = ldat,
    control.family = hlik
)
fit_g1 <- inla(
    formula = y ~ 0 + f(i, model = cgmodel1, replicate = s),
    data = ldat,
    control.family = hlik
)

nmc <- 10000
hlkj0sampls <- inla.hyperpar.sample(
    n = nmc, result = fit_lkj0, intern = TRUE
)
hlkj1sampls <- inla.hyperpar.sample(
    n = nmc, result = fit_lkj1, intern = TRUE
)
hc0sampls <- inla.hyperpar.sample(
    n = nmc, result = fit_c0, intern = TRUE
)
hc1sampls <- inla.hyperpar.sample(
    n = nmc, result = fit_c1, intern = TRUE
)
hg0sampls <- inla.hyperpar.sample(
    n = nmc, result = fit_g0, intern = TRUE
)
hg1sampls <- inla.hyperpar.sample(
    n = nmc, result = fit_g1, intern = TRUE
)

iup <- which(upper.tri(cc))

clkj0sampls <- t(sapply(1:nmc, function(i) {
    tcrossprod(cholcor(hlkj0sampls[i, ]))[iup]
}))
clkj1sampls <- t(sapply(1:nmc, function(i) {
    tcrossprod(cholcor(hlkj1sampls[i, ]))[iup]
}))
cpc0sampls <- t(sapply(1:nmc, function(i) {
    tcrossprod(cholcor(hc0sampls[i, ]))[iup]
}))
cpc1sampls <- t(sapply(1:nmc, function(i) {
    tcrossprod(cholcor(hc1sampls[i, ]))[iup]
}))
cg0sampls <- t(sapply(1:nmc, function(i) {
    vcov(gpc, theta = hg0sampls[i,])[iup]
}))
cg1sampls <- t(sapply(1:nmc, function(i) {
    vcov(gpc, theta = hg1sampls[i,])[iup]
}))


lG0 <- upperPadding(matrix(1,p,p))
lG1 <- upperPadding(attr(gpc, "graph"))

fcols <- c(
    gray(0:1/2),
    rgb(1, c(.5,0.1), c(.1,.5)),
    rgb(c(.5,.1), c(0.1,0.5), 1)
)

#png("exam_marks.png", width = 4000, height = 2500, res = 300)
par(mfcol = c(p,p), mar = c(1.6, 1.6, 0.1, 0.1),
    mgp = c(1, 0.5, 0), bty = 'n')
kc <- k2 <- k1 <- k0 <- 0
for(i in 1:p) {
    for(j in 1:p) {
        if(i==j) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, xnames[j], cex = 2)
        }
        if(j>i) {
            iij0 <- which((i == (lG0@i+1) & (j == (lG0@j+1))))
            iij1 <- which((i == (lG1@i+1) & (j == (lG1@j+1))))
            if(length(iij0)>0) {
                k0 <- k0 + 1
                ml0 <- inla.smarginal(fit_lkj0$internal.marginals.hy[[k0]])
                ml1 <- inla.smarginal(fit_lkj1$internal.marginals.hy[[k0]])
                mc0 <- inla.smarginal(fit_c0$internal.marginals.hy[[k0]])
                mc1 <- inla.smarginal(fit_c1$internal.marginals.hy[[k0]])
                h0 <- TRUE
            } else {
                h0 <- FALSE
            }
            if(length(iij1)>0) {
                k1 <- k1 + 1
                mg0 <- inla.smarginal(fit_g0$internal.marginals.hy[[k1]])
                mg1 <- inla.smarginal(fit_g1$internal.marginals.hy[[k1]])
                h1 <- TRUE
            } else {
                h1 <- FALSE
            }
            if(h0) {
                plot(ml0, type = "l", 
                     main = '', xlab = '', ylab = '',
                     xlim = range(0, ml0$x, ml1$x, mc0$x, mg0$x, mg1$x), 
                     ylim = range(ml0$y, ml1$y, mc0$y, mc1$y, mg0$y, mg1$y),
                     col = fcols[1], lwd = 2)
                lines(ml1, col = fcols[2], lwd = 2)
                lines(mc0, col = fcols[3], lwd = 2)
                lines(mc1, col = fcols[4], lwd = 2)
                if(h1) {
                    lines(mg0, col = fcols[5], lwd = 2)
                    lines(mg1, col = fcols[6], lwd = 2)
                    legend("topleft", bty = "n",
                           title = c("Internal parameter"),
                           as.expression(lapply(
                               rep(c(k0,k1), c(4,2)), function(i)
                                   bquote(theta[.(i)]))),
                           lwd = 2, col = fcols)
                } else {
                    legend("topleft", bty = "n",
                           title = "Internal parameter",
                           as.expression(lapply(rep(k0,4), function(i)
                               bquote(theta[.(i)]))),
                           lwd = 2, col = c(fcols[1:4]))
                }
                rug(0, 0.1, lty = 3, lwd = 2)
            } 
        }
        if(j<i) {
            k2 <- k2 + 1
            c_obs <- cc[iup[k2]]
            ic_obs <- tanh(c_obs + qnorm(c(0.025, 0.975)) / sqrt(n-3))
            dl0 <- density(clkj0sampls[, k2])
            dl1 <- density(clkj1sampls[, k2])
            dc0 <- density(cpc0sampls[, k2])
            dc1 <- density(cpc1sampls[, k2])
            dg0 <- density(cg0sampls[, k2])
            dg1 <- density(cg1sampls[, k2])
            plot(dl0, freq = FALSE, 
                 xlim = range(ic_obs, dl0$x, dl1$x, dc0$x, dg0$x, dg1$x),
                 ylim = range(dl0$y, dl1$y, dc0$y, dc1$y, dg0$y, dg1$y), 
                 main = '', xlab = '', ylab = '',
                 col = fcols[1], lwd = 2)
            lines(dl1, col = fcols[2], lwd = 2)
            lines(dc0, col = fcols[3], lwd = 2)
            lines(dc1, col = fcols[4], lwd = 2)
            lines(dg0, col = fcols[5], lwd = 2)
            lines(dg1, col = fcols[6], lwd = 2)
            rug(c(c_obs, ic_obs), 0.1, lty = 1, lwd = 2, col = 6)
        }
        if((i==2) & (j==1)) {
            legend("topleft", "|", title = "Obs. IC", bty = "n",
                   lty = c(0), lwd = c(2), col = c(1), text.col = 6)
        }
        if((i==2) & (j==1)) {
            legend("topright", bty = "n", col = fcols, lwd = 2,
                   c("LKJ0", "LKJ1", "PC0", "PC1", "G0", "G1"))
        }
    }
}
#dev.off()

#system("eog exam_marks.png &")
