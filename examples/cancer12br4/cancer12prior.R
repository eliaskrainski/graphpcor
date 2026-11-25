
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

    ## Drawn samples and visualize
    plot_corr_samples <- function(A, R, labels = NULL, n_overlay = 0,
                                  pal = colorRampPalette(c("firebrick", "white", "steelblue"))(201)) {
        p <- dim(A)[1]
        if (is.null(labels)) labels <- paste0("V", seq_len(p))
        if(missing(R) || is.null(R)) R <- apply(A, c(1, 2), mean)
        ## ellipse glyph for a correlation r, centred at 0, fitting in [-1, 1]
        ell <- function(r, scale = 0.9, k = 100) {
            t <- seq(0, 2 * pi, length.out = k)
            d <- acos(max(min(r, 1), -1))
            cbind(scale * cos(t + d / 2) / sqrt(2),
                  scale * cos(t - d / 2) / sqrt(2)) * sqrt(2)
        }
        col_of <- function(r) pal[round((r + 1) / 2 * 200) + 1]
        op <- par(mfrow = c(p, p), mar = c(0.3, 0.3, 0.3, 0.3)*0, oma = c(1, 1, 1, 1))
        on.exit(par(op))
        for (i in seq_len(p)) {
            for (j in seq_len(p)) {
                if (i == j) {       ## diagonal: labels
                    plot.new(); box()
                    text(0.5, 0.5, labels[i], cex = 1.5, font = 2)
                } else if (i > j) { ## lower: ellipses
                    plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1,
                         axes = FALSE, xlab = "", ylab = "")
                    box()
                    polygon(ell(R[i, j]), col = col_of(R[i, j]), border = "black")
                    ## optional: faint ellipses for a few individual samples
                    if (n_overlay > 0) {
                        idx <- sample(dim(A)[3], min(n_overlay, dim(A)[3]))
                        for (s in idx) lines(ell(A[i, j, s]), col = adjustcolor("grey40", 0.3))
                    }
                    text(-0.95*sign(R[i,j]), 0.9, sprintf("%.2f", R[i, j]),
                         adj = c(0.5,0.5), cex = 0.8)   
                } else {      ## upper: boxplots
                    boxplot(A[i, j, ], ylim = c(-1, 1), axes = FALSE,## outline = FALSE,
                            col = col_of(R[i,j]), border = "grey50",
                            outcol = "grey40", pch = 20, cex = 0.5)
                    abline(h = 0, lty = 2, col = "grey60")
                    box()
                }
            }
        }
    }

    set.seed(1)
    scorrels <- sample.basepcor(baseC0, 2000, lambda = 0.7)

    ## compare this plot with Figure 3 of
    ## Botella-Rocamora, P. et. al. (2015)
    plot_corr_samples(scorrels, baseC0$base, labels = causes)
    ##plot_corr_samples(scorrels, baseC0$base, n_overlay = 50)   ## overlay 20 sample ellipses

}
