
library(INLA)
library(graphpcor)

setwd(here::here("examples"))

###########################################################
## 1) get the UCI Heart Diseasex data from four hospitals
###########################################################

url <- paste0(
    "https://archive.ics.uci.edu/ml/",
    "machine-learning-databases/heart-disease/")

files <- c(
    cleveland = "processed.cleveland.data",
    hungary = "processed.hungarian.data",
    switzerland = "processed.switzerland.data",
    long_beach = "processed.va.data"
)

for(fl in files)
    if(!file.exists(fl))
        download.file(paste0(url, fl), fl)

## read each one as element in a list
ldat0 <- lapply(files, function(fl)
    data.frame(
        local = gsub("processed.", "",
                     gsub(".data", "", fl, fixed = TRUE),
                     fixed = TRUE),
        read.csv(
            file = fl, 
            header = FALSE,
            na.strings = "?",
            col.names = c(
                "age","sex","cp","trestbps","chol","fbs","restecg",
                "thalach","exang","oldpeak","slope","ca","thal","target"
            )
        )
    )
    )

sapply(ldat0, dim)

sapply(do.call("rbind", ldat0), summary)

par(mfrow = c(3,5), mar = c(4,4,1,1), mgp = c(3,1.5,0.5), bty = "n")
for(k in 2:15)
    hist(c(ldat0[[1]][, k], ldat0[[2]][, k],
           ldat0[[3]][, k], ldat0[[4]][, k]),
         main = names(ldat0[[1]])[k], xlab = '')

xsel <- c("age", "trestbps", "chol", "thalach", "oldpeak")
names(xsel) <- xsel
(p <- length(xsel))

alldx <- do.call("rbind", lapply(ldat0, function(d) d[xsel]))

par(mfrow = c(2,3), mar = c(3,3,0,0), mgp = c(2,0.5,0))
for(v in xsel) hist(alldx[,v], main = '', xlab = v)

##########################################################
### 2. Use first hospital data to define a graph by
## fitting a glm to each variable considering others
## as covariates, and select significant ones

## first we consider each variable's family
G0 <- matrix(0L, p, p)
for(i in 1:p) {
    jj <- which(!(xsel %in% xsel[i]))
    xxf <- paste(xsel[jj], collapse = "+")
    ff <- as.formula(paste(xsel[i], "~", xxf))
##    r <- inla(ff, data = ldat0[[1]][,xsel])
 ##   j <- (r$summary.fixed[-1, 3]>0) |
    ##      (r$summary.fixed[-1, 5]<0)
    r <- glm(ff, data=ldat0[[1]][, xsel])
    s <- coef(summary(r))
    j <- abs(s[-1,3])>1.5
    G0[jj[j], i] <- 1
}
dimnames(G0) <- list(xsel, xsel)

gph0 <- graphpcor(G0|t(G0))
gph0
G0 <- attr(gph0, "graph")

par(mfrow = c(1,1), mar = c(0,0,0,0))
plot(gph0, Rgraphviz = TRUE)

## find the best ordering
ch0 <- Cholesky(Laplacian(gph0) + Diagonal(ncol(G0)), perm = TRUE)
(ord0 <- ch0@perm + 1L)

## build the graph to minimizes the fill-in
gph <- graphpcor(G0[ord0, ord0])
dim(gph)

sum(abs(chol(Laplacian(gph0) + diag(p)))>0)
sum(abs(chol(Laplacian(gph) + diag(p)))>0)

par(mfrow = c(1,2), mar = c(0,0,0,0))
plot(gph0, Rgraphviz = TRUE)
plot(gph, Rgraphviz = TRUE)

## organize the data accounting for the order
summary(ldat0[[1]][xsel[ord0]])

## join the other 3 hospital data and order columns
vnams <- xsel[ord0]
dat3 <- do.call("rbind", lapply(ldat0[-1], function(x) x[vnams]))
(n3 <- nrow(dat3))
str(dat3)

## prepare terms
allnams <- c(vnams, "target")
b0v <- factor(rep(allnams, each = n3), allnams, allnams)
head(ic <- matrix(rep(1:p, each = n3), nrow = n3))
head(rc <- matrix(rep(1:n3, p), nrow = n3))
colnames(ic) <- paste0("i", 1:p)
colnames(rc) <- paste0("r", 1:p)
target3 <- unlist(lapply(ldat0[-1], function(d) d$target))

## data stack
dstack <- do.call(
    "inla.stack",
    c(lapply(vnams, function(y) {
        iy <- which(vnams %in% y)
        eff <- data.frame(b0=b0v[(iy-1)*n3 + 1:n3],
                          i = iy, r=1:n3)
        d <- list(lnk = rep(iy, n3), y = dat3[, iy])
        names(d)[2] <- y
        inla.stack(
            tag = y,
            data = d,
            effects = list(eff),
            A = list(1)
        )
    }),
    list(inla.stack(
        tag = "target",
        data = list(lnk = rep(p+1, n3),
                    target = target3),
        effects = list(data.frame(b0=b0v[p*n3+1:n3],
                                  ic, rc)),
        A = list(1)
    )))
)

str(dstack)

vnams
ff0 <- update(
    list(chol, oldpeak, trestbps, thalach, age, target) ~ 0 + b0,
    paste(".~.+",
          paste(paste0("f(i", 1:p, ", copy = 'i', replicate = r", 1:p,
                       ", fixed = FALSE)"), collapse = "+")))
ff0

ffd <- update(ff0, .~.+f(i, model = dmodel, replicate = r))
ffg <- update(ff0, .~.+f(i, model = gmodel, replicate = r))

## define the dense model
dmodel <- cgeneric(
    model = "pc_correl", n = p, lambda = 1,
    sigma.prior.probability = rep(0.01, p)
)

## define the sparse model 
gmodel <- cgeneric(
    model = gph, lambda = 1,
    sigma.prior.probability = rep(0.1, p)
)

## compare the number of parameters
c(M0 <- length(cgeneric_initial(dmodel)),
  M1 <- length(cgeneric_initial(gmodel)))
c(m0 <- p*(p-1)/2,
  m1 <- dim(gph)[2])

pfix <- list(prec = list(initial = 10, fixed = TRUE))
clk1 <- list(hyper = pfix)

fitd <- inla(
    formula = ffd, family = c(rep("gaussian", p), "poisson"),
    control.family = list(clk1, clk1, clk1, clk1, clk1, list()),
    data = inla.stack.data(dstack),
    control.predictor = list(link = lnk)
)

fitg <- inla(
    formula = ffg, family = c(rep("gaussian", p), "poisson"),
    control.family = list(clk1, clk1, clk1, clk1, clk1, list()),
    data = inla.stack.data(dstack),
    control.predictor = list(link = lnk)
)

## extract/transform marginals (std parameters)
sd.y <- lapply(fitd$internal.marginals.hyperpar[1:p], function(m)
    inla.tmarginal(function(x) exp(-x/2), m))
sd.x <- lapply(fitd$internal.marginals.hyperpar[p+1:p], function(m)
    inla.tmarginal(exp, m))

ss.y <- lapply(fitg$internal.marginals.hyperpar[1:p], function(m)
    inla.tmarginal(function(x) exp(-x/2), m))
ss.x <- lapply(fitg$internal.marginals.hyperpar[p+1:p], function(m)
    inla.tmarginal(exp, m))

par(mfrow = c(3,p), mar = c(3,3,0.5,0.5), mgp = c(1.5,0.5,0), bty = "n")
for(i in 1:p) {
    s1 <- inla.smarginal(fitd$marginals.fixed[[i]])
    s2 <- inla.smarginal(fitg$marginals.fixed[[i]])
    plot(s1$x, s1$y, xlim = range(s1$x, s2$x),
         ylim = range(s1$y, s2$y), type = "l", lwd = 2)
    lines(s2, col = 2, lty = 2, lwd = 2)
}
for(i in 1:p) {
    plot(sd.x[[i]][,1], sd.x[[i]][,2],
         xlim = range(sd.x[[i]][,1], ss.x[[i]][,1]),
         ylim = range(sd.x[[i]][,2], ss.x[[i]][,2]),
         type = "l", lwd = 2)
    lines(ss.x[[i]][,1], ss.x[[i]][,2],
          col = 2, lty = 2, lwd = 2)
}
plot(0,0,bty='n',axes=F,xlab='',ylab='')
for(i in 1:(p-1)) {
    plot(sd.y[[i]][,1], sd.y[[i]][,2],
         xlim = range(sd.y[[i]][,1], ss.y[[i]][,1]),
         ylim = range(sd.y[[i]][,2], ss.y[[i]][,2]),
         type = "l", lwd = 2)
    lines(ss.y[[i]][,1], ss.y[[i]][,2],
          col = 2, lty = 2, lwd = 2)
}

rbind(##fiti$cpu.used,
      fitd$cpu.used, fitg$cpu.used)

stds <- data.frame(
    obs = apply(dat3, 2, sd, na.rm = TRUE),
##    mi = exp(fiti$mode$theta[1:p]),
    m0x = exp(fitd$mode$theta[p+1:p]),
    m1x = exp(fitg$mode$theta[p+1:p]),
    m0y = exp(-0.5*fitd$mode$theta[1:p]),
    m1y = exp(-0.5*fitg$mode$theta[1:p]))
round(stds, 2)

nhsamples <- 5000
hsamples0 <- inla.hyperpar.sample(nhsamples, fit0)
hsamples1 <- inla.hyperpar.sample(nhsamples, fit1)

dim(hsamples0)
dim(hsamples1)

ill <- which(lower.tri(G1))
il1 <- which(lower.tri(G1) & (abs(G1)>0))

cholcor(hsamples0[1, p+1:m0], p)

ccsamples0 <- t(sapply(1:nhsamples, function(i) {
    tcrossprod(cholcor(hsamples0[i, p + 1:m0], p))[ill]
}))

th2c <- function(th)
    cov2cor(chol2inv(
        t(graphpcor:::Lprec0(th, p = p, iLtheta = il1, d0 = p:1))))

ccsamples1 <- t(sapply(1:nhsamples, function(i) {
    th2c(hsamples1[i, p + 1:m1])[ill]
}))

par(mfrow=c(4,4))
for(k in 1:15) hist(ccsamples0[, k], 100)

par(mfrow=c(4,4))
for(k in 1:15) hist(ccsamples1[, k], 100)


par(mfrow = c(p, p), mar = c(3,3,0.5,0.5), mgp = c(2,0.5,0), bty = "n")
k2 <- k1 <- 0
for(i in 1:p) {
    for(j in 1:p) {
        if(i==j) {
            mg0 <- inla.tmarginal(exp, fit0$internal.marginals.hyperpar[[i]])
            mg1 <- inla.tmarginal(exp, fit1$internal.marginals.hyperpar[[i]])
            plot(mg0, type = "l", ylab = "Density",
                 xlab = as.expression(bquote(sigma[.(i)])),
                 xlim = range(stds$obs[i], mg0[,1], mg1[,1]),
                 ylim = range(mg0[,2], mg1[,2]))
            lines(mg1, lwd = 2, lty = 2, col = 2)
            abline(v = stds$obs[i], lty = 2, col = 2)
        }
        if(i>j) {
            if(G1[i,j]!=0) {
                k1 <- k1 + 1
                mg1 <- inla.smarginal(
                    fit1$internal.marginals.hyperpar[[p+k1]])
                plot(mg1, type = "l", lwd = 2, ylab = "Density",
                     xlab = as.expression(bquote(theta[.(k1)])))
                abline(v = 0, lty = 2, col = 2)
            } else {
                plot(0, type = 'n', axes = FALSE, xlab = '', ylab = '')
            }
        }
        if(i<j) {
            k2 <- k2 + 1
            h0 <- hist(ccsamples0[, k2], 100, plot = FALSE)
            h1 <- hist(ccsamples1[, k2], 100, plot = FALSE)
            d0 <- density(ccsamples0[, k2])
            d1 <- density(ccsamples1[, k2])
            plot(h0, freq = FALSE,
                 col = gray(0.5, 0.5), border = 'transparent',
                 main = '', ylim = range(d0$y, d1$y),
                 xlim = range(0, d0$x, d1$x),
                 xlab = as.expression(bquote(rho[.(i)~","~.(j)])))
            plot(h1, add = TRUE, freq = FALSE)
            lines(d0$x, d0$y, lwd = 2)
            lines(d1$x, d1$y, lty = 2, col = 2, lwd = 2)
        }
    }
}
