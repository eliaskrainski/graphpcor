
library(graphpcor)

getwd()

## correlations from Figure 3 of
## Botella-Rocamora, P. et. al. (2015)
## 540 municipalities, 20 years (1987-2006)
nc <- 12
C0 <- diag(0.5,nc,nc);
C0[upper.tri(C0)] <- c(
    41,
    36,38,
    27,37,45,
    46,33,35,35,
    46,36,48,46,37,
    69,51,22,20,32,25,
    70,52,48,39,48,49,57,
    20,37,22,20,06,37,28,30,
    61,43,45,59,45,45,51,76,31,
    51,34,45,48,28,48,40,48,18,46,
    39,30,51,42,54,48,34,40,37,49,38
)/100; C0 <- C0 + t(C0)

causes <- c(
    "Oral", "Stomach", "Colon", "Rectum",
    "Liver", "Pancreatic", "Larynx", "Lung",
    "Prostate", "Bladder", "Lymphatic", "Leukemia")
caus <- c(
    "Oral", "Stom", "Coln", "Rect", "Livr", "Panc",
    "Larx", "Lung", "Prst", "Blad", "Lymp", "Leuk")
dimnames(C0) <- list(caus, caus)
round(100*C0)

baseC0 <- basepcor(
    base = C0, d0 = (nc:1)^2,
    n = 540, gamma = 1.0, start = "empty")
length(baseC0$iLtheta) ## number of edges

round((baseC0$base-C0)*100)
mean(abs((baseC0$base - C0)*100))

cGraph0 <- graphpcor(baseC0)

ireord <- Cholesky(
    Laplacian(cGraph0) + Diagonal(nc),
    perm = TRUE)@perm + 1L
ireord

C0r <- C0[ireord, ireord]

G0 <- attr(cGraph0, "graph")
cGraph <- graphpcor(G0[ireord, ireord])

sum(abs(chol(Laplacian(cGraph0)+diag(nc)))>0)
sum(abs(chol(Laplacian(cGraph)+diag(nc)))>0)

(dg <- dim(cGraph))

Greord <- G0[ireord, ireord]
Greord
iLr <- which((abs(Greord)>0) & lower.tri(C0))
length(iLr)

## base at the reordered 
baseC <- basepcor(
    base = baseC0$base[ireord, ireord], iLtheta = iLr)

round((baseC$base - C0[ireord, ireord])*100)
mean(abs((baseC$base - C0[ireord, ireord])*100))
baseC$theta

if(FALSE) {
    
    H0 <- hessian(cGraph0, baseC0$theta)
    H <- hessian(cGraph, baseC$theta)
    length(bkH <- pretty(range(H0,H),65))
    imgcols <- topo.colors(length(bkH)-1)

    ## visualize
    par(mfrow = c(2,2), mar = c(0,0,0,0))
    plot(cGraph0, Rgraphviz = TRUE)
    plot(cGraph, Rgraphviz = TRUE)
    par(mar = c(0.5,0.5,0.5,0.5))
    image(H0, breaks = bkH, col = imgcols, axes = FALSE)
    image(H, breaks = bkH, col = imgcols, axes = FALSE)
    
    png("figures/graphpcor12cancer.png", 3000, 2000, res = 300)
    par(mar = c(0,0,0,0))
    plot(cGraph, Rgraphviz = TRUE)
    dev.off()

}
