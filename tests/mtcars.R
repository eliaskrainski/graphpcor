
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

    mcl <- structure(mc$r, class = 'dist', Size = ncol(mtcars),
                     Diag = FALSE, Upper = FALSE, method='marginal')
    pcl <- structure(pc$r, class = 'dist', Size = ncol(mtcars),
                       Diag = FALSE, Upper = FALSE, method='partial')    

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

library(fields)
library(spdep)
library(graphpcor)
library(INLA)

## 11 variables
sdat <- scale(mtcars)
(p <- ncol(cc <- cov(sdat)))
lc <- chol(cc)
Qc <- chol2inv(lc)

ii <- 1:p
image.plot(ii, ii, Qc)

## partial correlation matrix
pC <- cov2cor(Q)
dimnames(pC) <- dimnames(Q) <- dimnames(V) <-
    list(colnames(mtcars), colnames(mtcars))
round(pC*100)

## define a graphpcor from pC threshold
library(spdep)
nb <- lapply(1:p, function(i) setdiff(1:p,i)); class(nb) <- 'nb'
nbc <- nbcosts(nb, scale(mtcars))
nbw <- nb2listw(nb, nbc, style="B")
mst <- mstree(nbw)

G <- matrix(0, p, p, dimnames = dimnames(pC))
G
for(i in 1:nrow(mst)) {
    G[mst[i,1], mst[i,2]] <- 1
    G[mst[i,2], mst[i,1]] <- 1
}
Sparse(G)
g <- graphpcor(G)

g

c(p, p*(p-1)/2)

par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(g)

lp <- Laplacian(g)

Lplot <- function(L, ij) {
    xy <- eigen(L)$vectors[,ij]    
    plot(xy[, 1:2], cex = 10, xlab = "", ylab = "", axes = FALSE)
    text(xy[, 1], xy[, 2], colnames(L))
    for(i in 1:nrow(lp)) {
        jj <- which(abs(lp[i,])>0)
        if(length(jj)>0) {
            segments(xy[i, 1], xy[i, 2],
                     xy[jj, 1], xy[jj, 2])
        }
    }
    return(NULL)
}

gg <- graph_from_adjacency_matrix(attr(g, "graph"))

par(mfrow = c(3,3), mar = c(0,0,0,0))
plot(g)
plot(gg)
Lplot(lp, 1:2)
Lplot(lp, c(1,3))
Lplot(lp, c(1,4))
Lplot(lp, -2:-1 + ncol(lp))
Lplot(lp, c(1, ncol(lp)-1))
Lplot(lp, c(2, ncol(lp)-1))

library(igraph)


plot(gg)

attr(g, "graph")
g2graph <- function(g) {
    g <- attr(g, "graph")
    
}

round(sigmas <- data.frame(obs = sqrt(diag(V))), 1)

## define the cgeneric model (including variances)
cmodel0 <- cgeneric(
    "LKJ", n = p, eta = 5, 
    sigma.prior.reference = sigmas$obs, 
    sigma.prior.probability = rep(0.5, p),
    useINLAprecomp = FALSE)
cmodel1 <- cgeneric(
    "pc_correl", n = p, lambda = 2, 
    sigma.prior.reference = sigmas$obs, 
    sigma.prior.probability = rep(0.5, p),
    useINLAprecomp = FALSE)
cmodel2 <- cgeneric(
    g, lambda = 1, 
    sigma.prior.reference = sigmas$obs, 
    sigma.prior.probability = rep(0.5, p),
    useINLAprecomp = FALSE)

## prepare the (long format) data for INLA
n <- nrow(mtcars)
nc <- c(nc=1, nc=5)
datac <- lapply(nc, function(a) {
    d <- data.frame(
        i = rep(rep(1:p, each = n), a), ## id for variables
        r = rep(rep(1:n, p), a),       ## id for replicates
        y = rep(as.vector(as.matrix(mtcars)), a))
    d$a <- factor(d$i) ## intercept for each variable
    d
})

head(mtcars)
str(datac)

## model formula
f0 <- y ~ a + f(i, model = cmodel0, replicate = r)
f1 <- y ~ a + f(i, model = cmodel1, replicate = r)
f2 <- y ~ a + f(i, model = cmodel2, replicate = r)

## fix the likelihood precision to a high value
pprc <- list(prec = list(initial = 10, fixed = TRUE))

## model fitting
fits <- lapply(datac, function(datax) {
    list(fit0 = inla(
             formula = f0,
             control.family=list(hyper = pprc),
             data = datax,
             control.inla = list(int.strategy = "eb")),
         fit1 = inla(
             formula = f1,
             control.family=list(hyper = pprc),
             data = datax,
             control.inla = list(int.strategy = "eb")),
         fit2 = inla(
             formula = f2,
             control.family=list(hyper = pprc),
             data = datax,
             control.inla = list(int.strategy = "eb"))
         )
})
    
cor(colMeans(mtcars),
    sapply(fits[[1]], function(x)
        x$summary.fix[, 1]))

iil <- which(lower.tri(V))

sc1f <- function(theta) {
    t(sapply(1:nrow(theta), function(i) {
        cc <- tcrossprod(cholcor(theta[i, -(1:p)]))
        c(s = exp(theta[i, 1:p]), c = cc[iil])       
    }))
}
sc2f <- function(theta) {
    t(sapply(1:nrow(theta), function(i) {
        v <- vcov(g, theta = theta[i, ])
        c(s=sqrt(diag(v)), c=cov2cor(v)[iil])
    }))
}

vcsamples <- lapply(fits, function(lr)
    list(
        m0 = sc1f(inla.hyperpar.sample(
            n = 10000, lr[[1]], intern = TRUE)),
        m1 = sc1f(inla.hyperpar.sample(
            n = 10000, lr[[2]], intern = TRUE)),
        m1 = sc2f(inla.hyperpar.sample(
            n = 10000, lr[[3]], intern = TRUE))
    ))

sigmas$m0 <- exp(fit0$mode$theta[1:p])
sigmas$m1 <- exp(fit1$mode$theta[1:p])
sigmas$m2 <- exp(fit2$mode$theta[1:p])

cols <- c(gray(.5), rgb(c(1,.5), .7, c(.5,1), .5))

par(mfrow = c(1,2), mar = c(2.5,2.5,0.5,0.5), mgp = c(1.5,0.5,0))
plot(sigmas[, 1:2], ylim = range(sigmas), asp = 1,
     pch = 19, col = cols[1], cex = 3, log = 'xy')
points(sigmas[, c(1, 3)], pch = 19, col = cols[2], cex = 3)
points(sigmas[, c(1, 4)], pch = 19, col = cols[3], cex = 3)
abline(0:1, lty = 2)
plot(sigmas[, 1], sigmas[, 2]/sigmas[, 1],
     ylim = range(sigmas[, -1]/sigmas[, 1]),
     xlab = 'Obs sigma', ylab = 'Relative fit',
     pch = 19, col = cols[1], cex = 3, log = 'x')
points(sigmas[, 1], sigmas[, 3]/sigmas[, 1],
     pch = 19, col = cols[2], cex = 3)
points(sigmas[, 1], sigmas[, 4]/sigmas[, 1],
       pch = 19, col = cols[3], cex = 3)
abline(h=1, lty = 2)

for(i in 1:(p-1))
    cat(cumsum(p:1)[i]-i + (i+1):p, "\n")

ij <- lapply(1:(p-1), function(i)
    cumsum(p:1)[i]-i + (i+1):p)
ij

for(i in 2:p) {
    cat(sapply(1:(i-1), function(j) ij[[j]][i-j]), "\n")
}
ji <- lapply(2:p, function(i) sapply(1:(i-1), function(j) ij[[j]][i-j]))
ji

par(mfrow = c(p, p), mar = c(3,3,0.5,0.5), mgp = c(1,.5,0), bty = "n")
for(i in 1:p) {
    for(j in 1:p) {
        if(j<i) {
            if(TRUE) {
                cij <- cor(mtcars[, i], mtcars[, j])
                jj <- ji[[i-1]][j]
                h0 <- hist(vcsamples[[2]][[1]][, jj],
                           -10:10/10, plot = FALSE)
                h1 <- hist(vcsamples[[2]][[2]][, jj],
                           -10:10/10, plot = FALSE)
                h2 <- hist(vcsamples[[2]][[3]][, jj],
                           -10:10/10, plot = FALSE)
                plot(h0, main = '', xlab = "", col = cols[1],
                     border = 'transparent')
                plot(h1, add = TRUE, col = cols[2], border = 'transparent')
                plot(h2, add = TRUE, col = cols[3], border = 'transparent')
                legend("topright", bty = "n", 
                       title = paste(i,j, collapse = ","),
                       format(cij*100, digits = 1))
                abline(v = cij, lty = 2, lwd = 2)
                ##plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
            } else {
                cij <- cor(mtcars[, i], mtcars[, j])
                jj <- ji[[i-1]][j]
                hist(vc0samples[, jj], -10:10/10,
                     col = cols[1], border = 'transparent',
                     main = '', xlab = "")
                legend("topright", bty = "n", 
                       title = paste(i,j, collapse = ","),
                       format(cij*100, digits = 1))
                abline(v = cij, lty = 2, lwd = 2, col = 2)
            }
        }
        if(i==j) {
            h0 <- hist(vcsamples[[1]][[1]][, i], plot = FALSE)
            h1 <- hist(vcsamples[[1]][[2]][, i], plot = FALSE)
            h2 <- hist(vcsamples[[1]][[3]][, i], plot = FALSE)
            plot(h0, xlim = range(h0$breaks, h1$breaks, h2$breaks),
                 main = '', xlab = paste("s", i),
                 col = cols[1], border = 'transparent')
            plot(h1, add = TRUE, col = cols[2], border = 'transparent')
            plot(h2, add = TRUE, col = cols[3], border = 'transparent')
            abline(v = sqrt(diag(V))[i], lty = 2, lwd = 2, col = 3)
            if(j==1)
                legend("topright", c("LKJ", "PC1", "PC2"), bty = "n",
                       col = cols, lty = 2, lwd = 2, fill = cols,
                       border = cols)
        }
        if(j>i) {
            cij <- cor(mtcars[, i], mtcars[, j])
            jj <- ij[[i]][j-i]
            h0 <- hist(vcsamples[[1]][[1]][, jj],
                       -10:10/10, plot = FALSE)
            h1 <- hist(vcsamples[[1]][[2]][, jj],
                       -10:10/10, plot = FALSE)
            h2 <- hist(vcsamples[[1]][[3]][, jj],
                       -10:10/10, plot = FALSE)
            plot(h0, main = '', xlab = "", col = cols[1],
                 border = 'transparent')
            plot(h1, add = TRUE, col = cols[2], border = 'transparent')
            plot(h2, add = TRUE, col = cols[3], border = 'transparent')
            legend("topright", bty = "n", 
                   title = paste(i,j, collapse = ","),
                   format(cij*100, digits = 1))
            abline(v = cij, lty = 2, lwd = 2)
        }
    }
}
legend("topright", paste("p", 0:2), bty = "n", fill = cols,
       border = 'transparent')


detach("package:graphpcor", unload = TRUE)
library(graphpcor)


tree2 <- treepcor(
    p1 ~ p2 + c1 + c2,
    p2 ~ -c3 + c4)

str(tree2)

attr(tree2, "relationship")


tp2a(tree2)

gt <- graph_from_adjacency_matrix(abs(tp2a(tree2)))

plot(gt)

tree2
dim(tree2)
summary(tree2)

plot(tree2)
