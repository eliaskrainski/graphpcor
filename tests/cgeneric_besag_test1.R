
library(spdep)
library(INLA)
library(corGraphs)

inla.setOption(
    safe = FALSE
)

nxy <- c(100, 150)
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

idx <- which(R0@i >= R0@j) ## should be lower here and transposed later !!!
ord <- order(R0@i[idx])
Rgraph <- inla.as.sparse(
    sparseMatrix(
        i = R0@i[idx[ord]]+1L,
        j = R0@j[idx[ord]]+1L,
        x = R0@x[idx[ord]]
    )
)
str(Rgraph@i)
str(Rgraph@j)
Rgraph[1:min(n, 5), 1:min(n, 10)]

m1 <- do.call(
    "inla.cgeneric.define",
    list(
        model = "inla_cgeneric_besag",
        shlib = paste0(
            system.file(
                "libs/",
                package = "corGraphs"),
            "corGraphs.so"),
        n = n,
        debug = TRUE,
        Rgraph = t(Rgraph),
        param = c(1, 0.5)
    )
)

str(m1)
str(m1$f$cgeneric$data$smatrices,1)

theta1 <- 0

str(cgeneric_get(m1, "initial"))
str(cgeneric_get(m1, "mu", theta = theta1))
str(cgeneric_get(m1, "log_prior", theta = theta1))
str(g1ij <- cgeneric_get(m1, "graph"))
str(g1ij)
table(g1ij[[1]]<=g1ij[[2]])
str(cgeneric_get(m1, "graph", optimize = FALSE))
str(cgeneric_get(m1, "Q", theta = theta1))

str(Q1 <- cgeneric_get(m1, "Q", theta = theta1, optimize = FALSE))

image(Q1)

dtest1 <- list(
    i = 1:n,
    x = inla.qsample(n = 1, Q = Q1 + Diagonal(n, nnb) * 0.0001,
                     constr = list(A = matrix(1, 1, n), e = 0))[, 1]
)
dtest1$y <- rpois(n, exp(3 + dtest1$x))


r1 <- inla(
    y ~ f(i, model = m1, constr = FALSE) - 1,
    data = dtest1,
    family = 'poisson',
    control.inla = list(int.strategy = 'eb'),
    control.mode = list(
        theta = theta1,
        fixed = TRUE
    ),
    control.compute = list(config = TRUE)
)

ri <- inla(
    y ~ f(i, model = 'besag', graph = nb.graph,
          constr = FALSE,
          scale.model = FALSE, diagonal = 0,
          hyper = list(theta = list(prior = 'pc.prec',
                                    param = c(1, 0.5)))) - 1,
    data = dtest1,
    family = 'poisson',
    control.inla = list(int.strategy = 'eb'),
    control.mode = list(
        theta = theta1,
        fixed = TRUE
    ),
    control.compute = list(config = TRUE)
)

rbind(r1$cpu.used, ri$cpu.used)

print(r1$cpu.used["Total"] / 
      ri$cpu.used["Total"]) ### ;)


c(r1$mode$theta, ri$mode$theta)

diag(cor(r1$summary.random$i,
         ri$summary.random$i))

unlist(inla.zmarginal(inla.tmarginal(
    function(x) exp(-x/2),
    r1$internal.marginals.hyperpar[[1]]), TRUE))

unlist(inla.zmarginal(inla.tmarginal(
    function(x) exp(-x/2),
    ri$internal.marginals.hyperpar[[1]]), TRUE))

names(r1$misc$configs$config[[1]])

str(r1$misc$configs$config[[1]]$Qprior)
str(ri$misc$configs$config[[1]]$Qprior)

summary(r1$misc$configs$config[[1]]$Qprior@x/
        ri$misc$configs$config[[1]]$Qprior@x)

all.equal(r1$misc$configs$config[[1]]$Qprior,
          ri$misc$configs$config[[1]]$Qprior)

