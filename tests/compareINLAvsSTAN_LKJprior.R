
library(INLA)
library(graphpcor)
library(coda)
library(rstan)
options(mc.cores = 4L)

eta0s <- c(0.1, 0.5, 1, 1.1, 3, 5, 20, 50, 1000)
names(eta0s) <- sprintf("%1.1f", eta0s)

## notice that eta<1 may leads to trouble near corners
par(mfrow = c(3, 3), mar = c(4,4,2,1), mgp = c(3,1,0), las = 1, bty = 'n')
for(i in 1:length(eta0s)) {
    plot(function(x) dbeta((1+x)/2, eta0s[i], eta0s[i])/2, -1, 1, n = 1+1e4,
         ylim = c(0, max(dbeta(c(1e-4, 0.5, 1-1e-4), eta0s[i], eta0s[i])/2)),
         xlab = 'Correlation', ylab = "marginal prior density",
         main = as.expression(bquote(eta==.(eta0s[i]))))
    plot(function(x) sapply(x, function(r)
        dLKJ(matrix(c(1, r, r, 1), 2), eta0s[i])),
         -1+1e-5, 1-1e-5, n = 1+1e4, lty = 2, col = 2, add = TRUE)
    legend("center", format(integrate(function(x) sapply(x, function(r)
        dLKJ(matrix(c(1, r, r, 1), 2), eta0s[i])),
        -1+1e-7, 1-1e-7)$value), bty = "n")
}

## STAN model definition for generic 'eta'
model.eta <- "
data {
  int<lower=1> p;
  vector[p] y;
  vector[p] mu;
  real<lower=0> eta;
}
parameters {
  cholesky_factor_corr[p] Lcorr;
}
model {
  Lcorr ~ lkj_corr_cholesky(eta);
  y ~ multi_normal_cholesky(mu, Lcorr);
}
generated quantities {
  corr_matrix[p] corr;
  corr = multiply_lower_tri_self_transpose(Lcorr);
}
"

## it takes time to compile the model
system.time(
    stan_cmpld <- 
        stan_model(
            model_code = model.eta,
            model_name = "LKJ"
        )
)

## choose eta>1 (as inla will only fit one mode)
eta.s <- c(1.1, 5, 50, 500)
names(eta.s) <- paste0("eta", eta.s)
eta.s


stan.samples <- vector("list", length(eta.s))

for(i in 1:length(eta.s)) {
    cat("SAMPLING with eta = ", eta.s[i], "\n")
    stan.samples[[i]] <- sampling(
        stan_cmpld,
        data = list(p = as.integer(3), y = rep(0, 3),
                    mu = rep(0, 3), eta = eta.s[i]),
        iter = 30000,
        warmup = 5000,
        chains = 4
    )
    cat("Finish SAMPLING with eta = ", eta.s[i], "\n")
}

## collect and organize the 4 parallel sample chains
rhonams <- c("corr[2,1]", "corr[3,1]", "corr[3,2]")
stan.cor.samples <- lapply(stan.samples, function(s)
    mcmc(Reduce("cbind", extract(s, rhonams))))

## define the cgeneric LKJ models
cLKJetas <- lapply(eta.s, function(e)
    cgeneric("LKJ", n = 3, eta = e))

ires <- lapply(cLKJetas, function(cm) {
    ri <- inla(y ~ 0 + f(i, model = cm),
               data = data.frame(y = NA, i = 1:3),
               control.family = list(
                   hyper = list(prec = list(
                                    initial = 20, fixed = TRUE))))
    hs <- inla.hyperpar.sample(n = 5000, ri, TRUE)
    cs <- t(sapply(1:nrow(hs), function(i)
        tcrossprod(cholcor(hs[i, ], 3))[c(2,3,6)]))
})

cols <- c(rgb(1,0,0,.7), rgb(0,0,1,.7))
xlbs <- list(bquote(rho[2~","~1]),
             bquote(rho[3~","~1]),
             bquote(rho[3~","~2]))

par(mfrow = c(4, 3), mar = c(3,3,1,0.5), mgp = c(1.5,0.5,0))
for(i in 1:4) {
    for(j in 1:3) {
        ds <- density(stan.cor.samples[[i]][, j])
        di <- density(ires[[i]][, j])
        plot(function(x) dbeta((1+x)/2, eta.s[i], eta.s[i])/2,
             -1, 1, n = 1+1e4,
             ylim = c(0, max(dbeta(c(1e-4, 0.5, 1-1e-4),
                                   eta.s[i], eta.s[i])/2)),
             xlab = as.expression(xlbs[[j]]),
             ylab = "marginal prior density")
        lines(ds, col = cols[1], lwd = 2, lty = 2)
        lines(di, col = cols[2], lwd = 2, lty = 3)
        if(j==1)
            legend("topleft", "", bty = "n",
                   title = as.expression(bquote(eta==.(eta.s[i]))))
    }
}
legend("topright", c("Theory", "STAN", "INLA"), bty = "n",
       col = c(1, cols), lty = c(1,2,3))
