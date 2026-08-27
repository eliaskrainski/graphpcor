
### automatic setup working directory
setwd(here::here("examples"))
getwd()

library(graphpcor)
library(data.table)
library(INLA)

## build the graph from Fig. 3 in 
## https://www.nature.com/articles/s41598-024-65845-0
## library(data.table)

## from arriving arrows in 3rd graph
g0 <- graphpcor(
    cough ~ diarrhea + sob + muscle_sore + fever,
    cough ~ runny_nose + loss_of_taste + fatigue,
    cough ~ sore_throat + headache,
    muscle_sore ~ diarrhea + fever + fatigue,
    fever ~ sob,
    loss_of_smell ~ loss_of_taste,
    runny_nose ~ headache,
    loss_of_taste ~ runny_nose + fatigue,
    sore_throat ~ diarrhea + loss_of_smell + runny_nose + fatigue,
    headache ~ diarrhea + muscle_sore + fatigue
)
g0

(dd <- dim(g0))

(G0 <- attr(g0, "graph"))

reord <- inla.qreordering(G0)
str(reord)

g1 <- graphpcor(G0[reord$ireordering, reord$ireordering])
g1

dd
sum((abs(chol(diag(dd[1]) + Laplacian(g0)))>0) & upper.tri(G0))
sum((abs(chol(diag(dd[1]) + Laplacian(g1)))>0) & upper.tri(G0))

plot(g1, Rgraphviz = TRUE)

(G1 <- attr(g1, "graph"))
onames <- rownames(G1)
onames

## data from the paper
symptoms <- read.csv("symptoms.csv")
str(symptoms)

ods <- match(onames, names(symptoms))
ods

nd <- nrow(symptoms)
lsympt <- data.frame(
    i = rep(1:dd[1], each =  nd),
    r = rep(1:nd, dd[1]),
    symptom = factor(rep(onames, each = nd), onames, onames),
    y = unlist(symptoms[, ods])
)

cgprior0 <- cgeneric(g1, lambda = 1) ## \sigma_i = 1 (fixed)
ff0 <- y ~ 0 + symptom +
    f(i, model = cgprior0, replicate = r)

fit0 <- inla(
    formula = ff0, family = 'binomial',
    control.family = list(link = 'probit'),
    data = lsympt,
    control.inla = list(int.strategy = 'eb'),
    control.compute = list(config = TRUE)
)

fit0$misc$nfunc

image(fit0$misc$config$config[[1]]$Q)
image(fit0$misc$config$config[[1]]$Q[1:110, 1:110])

cmode <- vcov(g1, theta = fit0$mode$theta)
round(cmode * 100)

cc0 <- cor(symptoms[, ods])
cor(cc0[lower.tri(cc0)],
    cmode[lower.tri(cc0)])

bb <- cc0
bb[lower.tri(bb)] <- cmode[lower.tri(bb)]

fields::image.plot(bb, breaks = -32:32/32)

## use this as a base model to set the prior for the Brazilian data
cgprior1 <- cgeneric(g1, base = cmode, lambda = 1)

## The hospitalized 2021 SRAG data and its dictionary are available at:
url0 <- "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/"

## Download the dictionary file
dicfl <- "dicionario-de-dados-2019-a-2025.pdf"
if(!file.exists(dicfl))
   download.file(paste0(url0, dicfl), dicfl)

## file source and file name
csvfl <- "INFLUD21-23-03-2026.csv" ## last available version
## download (if not yet done)
if(!file.exists(csvfl)) 
    download.file(paste0(url0, "2021/", csvfl), csvfl)

## read the data (fread is fast)
##   (Here I set as.data.frame instead of data.table format)
##   (However operating with data.table can be muuuuch faster)
srag0 <- as.data.frame(fread(csvfl))

## dimension and names
dim(srag0)

## The disease classification 
## 1 : Infl., 2 : Other Resp, 3 : Other Etiol.,
## 4 : N. Specified, 5 : COVID19
##   (set NA as 0) :
na0 <- function(x) ifelse(is.na(x), 0, x)
table(na0(srag0$CLASSI_FIN))

## index for COVID
icov19 <- (na0(srag0$CLASSI_FIN)==5) ## logical (useful for crosstable)
iicov19 <- which(icov19)             ## index   (use this to select!!!)

fn1 <- function(x) {
### need this due to missing: which(x==1) as [x==1] does not work
    r <- integer(length(x))
    r[which(x==1)] <- 1L
    return(r)
}

fn2 <- function(x, patt) {
### need this as some symptoms are in "OTHERS"
    r <- integer(length(x))
    for(pat in patt)
        r[grep(pat, x)] <- 1L
    return(r)
}

onames

set.seed(1)
iicov19 <- sort(sample(iicov19, size = 1e5))
(n <- length(iicov19))

## select/organize the Brazilian data
data1 <- data.frame(
    sob      = fn1(srag0$DISPNEIA[iicov19]),
    diarrhea = fn1(srag0$DIARREIA[iicov19]),
    fever    = fn1(srag0$FEBRE[iicov19]),
    muscle_sore = fn2(srag0$OUTRO_DES[iicov19], "MUS"),
    loss_of_smell = fn1(srag0$PERD_OLFT[iicov19]),
    runny_nose = fn2(srag0$OUTRO_DES[iicov19], "CORIZA"),
    loss_of_taste = fn1(srag0$PERD_PALA[iicov19]),
    sore_throat = fn1(srag0$GARGANTA[iicov19]),
    cough = fn1(srag0$TOSSE[iicov19]),
    fatigue = fn1(srag0$FADIGA[iicov19]),
    headache = fn2(srag0$OUTRO_DES[iicov19],
                   c("CABECA", "CABEÇA", "CEFAL"))
)

(dsumm <- data.frame(
     data0 = colMeans(symptoms[, ods]),
     data1 = colMeans(data1)))

(nd1 <- nrow(data1))
ldat1 <- data.frame(
    i = rep(1:dd[1], each = nd1),
    r = rep(1:nd1, dd[1]),
    symptom = factor(rep(onames, each = nd1), onames, onames),
    y = unlist(data1)
)

ff1 <- y ~ 0 + symptom +
    f(i, model = cgprior0, replicate = r)

th.ini <- c(7.3352, -8.5898, -3.7433, -4.2434, 4.0964, 0.6043,
            -6.1895, -7.3269, 1.3144, 1.6973, -0.5722, -18.3429,
            -0.8781, -0.3527, -0.5694, 1.6964, -4.0465, 4.0544,
            -1.5749, -4.8009, -1.1551, -1.5412, -0.4514, 0.1599)

fit1 <- inla(
    formula = ff1, family = 'binomial',
    control.family = list(link = 'probit'),
    data = ldat1,
    control.mode = list(
        theta = rep(0, length(fit0$mode$theta)),
        restart = TRUE),
    num.threads = 16,
##    inla.call = 'remote',
    verbose = TRUE
)

par(mfrow = c(4,6), mar = c(3,3,0.5,0.5), mgp = c(1.5,0.5,0), bty = 'n')
for(i in 1:dd[2]) {
    d0 <- inla.smarginal(fit0$internal.marginals.hyperpar[[i]])
    d1 <- inla.smarginal(fit1$internal.marginals.hyperpar[[i]])
    plot(d0, xlim = range(0, d0$x, d1$x),
         ylim = range(d0$y, d1$y),
         lwd = 2, lty = 2)
    lines(d1, lwd = 2)
    abline(h = 0, lwd = 2, lty = 2, col = "red")
}
