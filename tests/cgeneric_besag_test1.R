
library(spdep)
library(INLA)
library(corGraphs)

inla.setOption(
    num.threads = 1,
    safe = FALSE
)

nxy <- c(5, 5)
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
    constr = FALSE
)

str(m1)
str(m1$f$cgeneric$data$smatrices,1)

str(cgeneric_get(m1, "initial"))

theta1 <- 0
str(cgeneric_get(m1, "mu", theta = theta1))
str(cgeneric_get(m1, "log_prior", theta = theta1))

str(g1ij <- cgeneric_get(m1, "graph"))
str(g1ij)
table(g1ij[[1]]<=g1ij[[2]])

str(cgeneric_get(m1, "graph", optimize = FALSE))

str(cgeneric_get(m1, "Q", theta = theta1))

Q1 <- cgeneric_get(m1, "Q", theta = theta1, optimize = FALSE)
str(Q1)

Q1[1:min(5,n), 1:min(10,n)]
image(Q1)

dtest1 <- list(
    i = 1:n,
    y = rep(NA, n)
)

cpred <- list(link = 1)
cmode <- list(
    theta = theta1,
    restart = FALSE, 
    fixed = TRUE
)
cinla <- list(
        int.strategy = 'eb'
)
ccomp <- list(config = TRUE)
cfam <- list(
    hyper = list(
        prec = list(
            initial = 10,
            fixed = TRUE
        )
    )
)

fit.1 <- inla(
    formula = y ~ 0 +
        f(i, model = m1),
    data = dtest1,
    verbose = !TRUE,
    control.family = cfam,
    control.predictor = cpred,
    control.inla = cinla,
    control.mode = cmode,
    control.compute = ccomp
)

fg <- y ~ 0 + 
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
    control.family = cfam,
    control.predictor = cpred,
    control.inla = cinla,
    control.mode = cmode,
    control.compute = ccomp
)

rbind(fit.1$cpu.used, fit.i$cpu.used)
rbind(fit.1$misc$nfunc, fit.i$misc$nfunc)

c(fit.1$cpu.used["Total"] / fit.1$misc$nfunc,
  fit.i$cpu.used["Total"] / fit.i$misc$nfunc) ### ;)

c(fit.1$mode$theta, fit.i$mode$theta)

names(fit.1$misc$configs$config[[1]])

all.equal(fit.1$misc$configs$config[[1]]$Qprior,
          fit.i$misc$configs$config[[1]]$Qprior)


detach("package:corGraphs", unload = TRUE)
library(corGraphs)
