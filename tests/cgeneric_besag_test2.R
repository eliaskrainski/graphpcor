
library(spdep)
library(INLA)
library(corGraphs)

inla.setOption(
    safe = FALSE
)

nxy <- c(40, 50)
nb <- grid2nb(d = nxy, queen = FALSE)
nnb <- card(nb)
(n <- length(nnb))

nb.graph <- sparseMatrix(
    i = rep(1:n, nnb),
    j = unlist(nb[nnb>0]),
    x = 1, 
    dims = c(n, n)
)

theta <- 0

Q0 <- inla.as.sparse(
    exp(theta) * (
        Diagonal(n, nnb) - nb.graph
    )
)

cnstr <- list(
    A = matrix(1, 1, n),
    e = 0
)

Qs <- inla.scale.model(
    Q0,
    constr = cnstr
)

Q0[1:min(5, n), 1:min(12, n)]
Qs[1:min(5, n), 1:min(12, n)]

pparam <- c(1, 0.01)

m1 <- cgeneric_besag(
    graph = nb.graph,
    param = pparam,
    constr = !TRUE,
    scale = TRUE,
    debug = !TRUE
)

str(m1)

dtest1 <- list(
    i = 1:n,
    x = inla.qsample(
        n = 1,
        Q = Qs + Diagonal(n, nnb) * 0.0001,
        constr = cnstr
    )[, 1]
)
dtest1$y <- rpois(n, exp(3 + dtest1$x))

r1 <- inla(
    y ~ 0 + f(i, model = m1, constr = FALSE),
    data = dtest1,
    family = 'poisson',
    control.inla = list(
        int.strategy = 'eb'
    ),
    control.mode = list(
        theta = theta,
        fixed = TRUE,
        restart = !TRUE
    ),
    control.compute = list(
        config = TRUE)
)

ri <- inla(
    y ~ 0 +
        f(i, model = 'besag', graph = nb.graph,
          constr = FALSE,
          scale.model = TRUE, diagonal = 0,
          hyper = list(
              theta = list(
                  prior = 'pc.prec',
                  param = pparam
              )
          )
          ),
    data = dtest1,
    family = 'poisson',
    control.inla = list(
        int.strategy = 'eb'),
    control.mode = list(
        theta = theta,
        fixed = TRUE,
        restart = !TRUE
    ),
    control.compute = list(config = TRUE)
)

grep("iagonal", r1$logfile, value = TRUE)
grep("iagonal", ri$logfile, value = TRUE)

grep("onstraint", r1$logfile, value = TRUE)
grep("onstraint", ri$logfile, value = TRUE)

rbind(r1$cpu.used, ri$cpu.used)

print(r1$cpu.used["Total"] / 
      ri$cpu.used["Total"]) ### ;)

c(r1$misc$nfunc,
  ri$misc$nfunc)

c(r1$mode$theta, ri$mode$theta)

diag(cor(r1$summary.random$i,
         ri$summary.random$i))

unlist(inla.zmarginal(inla.tmarginal(
    function(x) exp(-x/2),
    r1$internal.marginals.hyperpar[[1]]), TRUE))

unlist(inla.zmarginal(inla.tmarginal(
    function(x) exp(-x/2),
    ri$internal.marginals.hyperpar[[1]]), TRUE))

all.equal(r1$misc$configs$config[[1]]$Qprior,
          ri$misc$configs$config[[1]]$Qprior)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
