
library(spdep)
library(INLA)
library(corGraphs)

inla.setOption(
    num.threads = 1,
    safe = FALSE
)

nxy <- c(100, 200)
nb <- grid2nb(d = nxy, queen = FALSE)
nnb <- card(nb)
n <- length(nnb)

nb.graph <- sparseMatrix(
    i = rep(1:n, nnb),
    j = unlist(nb[nnb>0]),
    x = 1,
    dims = c(n, n)
)

R0 <- inla.as.sparse(Diagonal(n, nnb) - nb.graph)
R0[1:min(n, 5), 1:min(n, 10)]

m1 <- cgeneric_generic0(
    R = inla.as.sparse(Diagonal(n, nnb) * 0.25 + R0),
    param = c(1, 0.5),
    constr = FALSE,
    scale = FALSE
)

theta1 <- 0
Q1 <- cgeneric_get(m1, "Q", theta = theta1, optimize = FALSE)

Q1[1:min(5,n), 1:min(10,n)]

x1 <- inla.qsample(
    Q = Q1,
    n = 1
)

dtest1 <- list(
    i = 1:n,
    y = rpois(n, exp(3 + x1))
)

cpred <- list(link = 1)
cmode <- list(
    theta = theta1,
    restart = TRUE, 
    fixed = FALSE
)
cinla <- list(
        int.strategy = 'eb'
)
ccomp <- list(config = TRUE)

fit.1 <- inla(
    formula = y ~ 1 +
        f(i, model = m1),
    data = dtest1,
    family = "poisson",
    control.predictor = cpred,
    control.inla = cinla,
    control.mode = cmode,
    control.compute = ccomp
)

fg <- y ~ 1 + 
        f(i, model = 'generic0', Cmatrix = Q1, 
          constr = FALSE,
          diagonal = 0,
          hyper = list(
              theta = list(
                  initial = theta1,
                  prior = 'pc.prec',
                  param = c(1, 0.5)
              )
          ))

fit.i <- inla(
    formula = fg,
    data = dtest1,
    family = "poisson", 
    control.predictor = cpred,
    control.inla = cinla,
    control.mode = cmode,
    control.compute = ccomp
)

rbind(fit.1$cpu.used, fit.i$cpu.used)
rbind(fit.1$misc$nfunc, fit.i$misc$nfunc)

c(fit.1$cpu.used["Total"] / fit.1$misc$nfunc,
  fit.i$cpu.used["Total"] / fit.i$misc$nfunc) 

c(fit.1$mode$theta, fit.i$mode$theta)

args(INLA:::plot.inla)
##plot(fit.1, F, F, F, T, F, plot.prior = TRUE, plot.opt.trace = TRUE)
##plot(fit.i, F, F, F, T, F, plot.prior = TRUE, plot.opt.trace = TRUE)

diag(cor(fit.1$summary.random$i,
         fit.i$summary.random$i))

unlist(inla.zmarginal(inla.tmarginal(
    function(x) exp(-x/2),
    fit.1$internal.marginals.hyperpar[[1]]), TRUE))

unlist(inla.zmarginal(inla.tmarginal(
    function(x) exp(-x/2),
    fit.i$internal.marginals.hyperpar[[1]]), TRUE))

lprec.seq <- seq(-3, 5, 0.1)
prec.seq <- exp(lprec.seq)
cg.lprior <- sapply(lprec.seq, function(x)
    cgeneric_get(m1, "log.prior", theta = x))

post1 <- inla.tmarginal(
    exp, fit.1$internal.marginals.hyperpar[[1]])
post1s <- inla.tmarginal(
    function(x) exp(-x/2),
    fit.1$internal.marginals.hyperpar[[1]])

posti <- inla.tmarginal(
    exp, fit.i$internal.marginals.hyperpar[[1]])
postis <- inla.tmarginal(
    function(x) exp(-x/2),
    fit.i$internal.marginals.hyperpar[[1]])

par(mfrow = c(1, 2), mar = c(3, 3, 1.5, 0.5), mgp = c(1.5, 0.5, 0))
plot(posti, type = "l")
lines(post1, col = 2)
lines(prec.seq, exp(cg.lprior - lprec.seq), col = 4, lty = 2)
plot(postis, type = "l")
lines(post1s, col = 2)
lines(exp(-0.5*lprec.seq),
      exp(cg.lprior - 0.5 * lprec.seq), col = 4, lty = 2)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
