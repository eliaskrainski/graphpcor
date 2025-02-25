
library(spdep)
library(INLA)
library(graphpcor)

inla.setOption(
    num.threads = 1,
    safe = FALSE
)

nxy <- c(25, 35)
nb <- spdep::grid2nb(d = nxy, queen = FALSE)
nnb <- spdep::card(nb)
n <- length(nnb)

nb.graph <- sparseMatrix(
    i = rep(1:n, nnb),
    j = unlist(nb[nnb>0]),
    x = 1,
    dims = c(n, n)
)

R0 <- inla.as.sparse(Diagonal(n, nnb) - nb.graph)
R0[1:min(n, 5), 1:min(n, 10)]

m1 <- cgeneric(
    model = "generic0",
    R = inla.as.sparse(Diagonal(n, nnb) * 0.25 + R0),
    param = c(1, 0.5),
    constr = FALSE
)

str(m1)
str(m1$f$cgeneric$data$smatrices,1)

str(initial(m1))

mu(m1)

theta1 <- 0
prior(m1, theta = theta1)

str(g1ij <- graph(m1, optimize = TRUE))
str(g1ij)
table(g1ij[[1]]<=g1ij[[2]])

str(graph(m1))

str(precision(m1, theta = theta1))

Q1 <- precision(m1, theta = theta1)
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
    verbose = TRUE,
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

all.equal(fit.1$misc$configs$config[[1]]$Qprior,
          fit.i$misc$configs$config[[1]]$Qprior)


detach("package:graphpcor", unload = TRUE)
library(graphpcor)
