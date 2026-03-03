
library(INLA)
library(graphpcor)
library(coda)
library(rstan)
options(mc.cores = 4L)

## STAN model code without the prior for L, as LCorr
Scode0 <- "
data {
  int<lower=1> p;
  int<lower=1> n;
  vector[p] y[n];
  vector[p] mu;
}
transformed parameters {
  matrix[p,p] rho;
}
model {
  y ~ multi_normal(mu, rho);
}
"

## add the pc_correl code
Scode <- stan_add(Scode0, 'graphpcor', lambda = 1, name = "rho")
Scode

## it takes time to compile the model
system.time(
    stan_cmpld <- 
        stan_model(
            model_code = Scode,
            model_name = "pc_graphpcor"
        )
)

## Correlation model definition
g <- graphpcor(x1~x2+x3, x4~x2+x3)
p <- dim(g)[1]
m <- dim(g)[2]
c(p, m)

##   PARAMETER definition
set.seed(2)
th.true <- rnorm(m)
corr <- basepcor(th.true, p, g)$base
Uc <- chol(corr)

##   SAMPLE drawn
n <- 100
set.seed(3)
y <- matrix(rnorm(n*p), n) %*% Uc

(ycorr <- cor(y))

##   initial STAN data
Sdata0 <- list(
    n = as.integer(n),
    p = as.integer(p),
    y = y,
    mu = rep(0,p)
)
str(Sdata0)

##   base model definition
## prior parameters definition
set.seed(1)
th.base <- rnorm(m)*0
baseC <- basepcor(th.base, p = p, iLtheta = g)

## choose lambda>1 (as inla will only fit one mode)
lambdas <- c(0.1, 1, 10, 100)
names(lambdas) <- sprintf("%1.1f", lambdas)

stan.samples <- vector("list", length(lambdas))

for(i in 1:length(lambdas)) {
    cat("SAMPLING with lambda = ", lambdas[i], "\n")
    stan.samples[[i]] <- sampling(
        stan_cmpld,
        data = stan_add(
            x = Sdata0, base = baseC,
            lambda = lambdas[i], name = "rho"),
        iter = 30000,
        warmup = 5000,
        chains = 4
    )
    cat("Finish SAMPLING with lambda = ", lambdas[i], "\n")
}

## collect and organize the 4 parallel sample chains
thnams <- paste0("grpc_theta[", 1:m, "]")
stan.th.samples <- lapply(stan.samples, function(s)
    mcmc(Reduce("cbind", extract(s, thnams))))

sapply(stan.th.samples, dim)

## for the correlations
rhonams <- c("rho[2,1]", "rho[3,1]", "rho[3,2]")
stan.cor.samples <- lapply(stan.samples, function(s)
    mcmc(Reduce("cbind", extract(s, rhonams))))

## define the cgeneric models for each lambda
cgLambdas <- lapply(lambdas, function(l)
    cgeneric(g, base = baseC, lambda = l))

Idat <- data.frame(
    y = as.vector(y),
    i = rep(1:p, each = n),
    r = rep(1:n, p))

ifits <- lapply(cgLambdas, function(cm) {
    cat("fitting ... ")
    o <- inla(
        y ~ 0 + f(i, model = cm, replicate = r),
        data = Idat, 
        control.family = list(
            hyper = list(prec = list(
                             initial = 20, fixed = TRUE))))
    cat("ok!\n")
    return(o)
})

iil <- which(lower.tri(diag(p)))
ires <- lapply(ifits, function(ri) {
    hs <- inla.hyperpar.sample(n = 5000, ri, TRUE)
    cs <- t(sapply(1:nrow(hs), function(i)
        basepcor(hs[i, ], p, baseC$iLtheta,
                 d0 = baseC$d0)$base[iil]))
})

cols <- c(rgb(1,0,0,.7), rgb(0,0,1,.7), 1, 6)
xlbs <- list(bquote(rho[2~","~1]),
             bquote(rho[3~","~1]),
             bquote(rho[3~","~2]))

par(mfrow = c(4, 6), mar = c(3,3,0.1,0.1),
    mgp = c(1.5,0.5,0), bty = "n")
for(i in 1:4) {
    for(j in 1:3) {
        thj <- c(th.true[j], th.base[j])
        h.thj <- hist(stan.th.samples[[i]][, j], 100, plot = FALSE)
        ds.thj <- density(stan.th.samples[[i]][, j])
        sj <- inla.smarginal(ifits[[i]]$internal.marginals.hyperpar[[j]])
        plot(h.thj, freq = FALSE, main = "",
             xlim = range(th.base[j], th.true[j],
                          h.thj$breaks, ds.thj$x, sj$x),
             ylim = range(h.thj$density, ds.thj$y, sj$y),
             xlab = as.expression(bquote(theta[j])),
             border = 'transparent')
        lines(ds.thj, col = cols[1], lwd = 2, lty = 2)        
        lines(sj, col = cols[2], lwd = 2, lty = 2)
        abline(v = thj, lty = 3, col = cols[3:4])
        if(j==1)
            legend("topleft", "", bty = "n",
                   title = as.expression(bquote(lambda==.(lambdas[i]))))
    }
    for(j in 1:3) {
        h.cj <- hist(stan.cor.samples[[i]][, j], 100, plot = FALSE)
        ds <- density(stan.cor.samples[[i]][, j])
        di <- density(ires[[i]][, j])
        ccj <- c(corr[iil[j]], baseC$base[iil[j]])
        plot(h.cj, freq = FALSE, main = "",
             xlim = range(ccj, h.cj$breaks, ds$x, di$x),
             ylim = range(h.cj$density, ds$y, di$y),
             xlab = as.expression(xlbs[[j]]),
             border = "transparent")
        lines(ds, col = cols[1], lwd = 2, lty = 2)
        lines(di, col = cols[2], lwd = 2, lty = 2)
        abline(v = ccj, lty = 3, col = cols[3:4])
    }
    if(i==1)
        legend("topright", bty = "n",
               c("STAN", "INLA", "TRUE", "base"),
               col = c(cols), lty = c(2,2,3,3))
}
