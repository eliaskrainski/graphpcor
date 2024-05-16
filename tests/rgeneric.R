
rcar1 <- function(cmd = c("graph", "Q", "mu", "initial",
                           "log.norm.const",
                           "log.prior", "quit"),
                   theta = NULL, ...) {
    graph <- function(n, theta)
        return(data$graph)
    Q <- function(n, theta) {
        qq <- inla.as.sparse(
        (data$C * exp(theta[1]) + data$G) * 1)
##        exp(theta[1]))
        idx <- which(qq@i <= qq@j)
        return(qq@x[idx])
    }
    mu <- function(n, theta)
        return(numeric(0))
    log.norm.const <- function(n, theta)
        return(numeric(0))
    log.prior <- function(n, theta) {
        return(sum(theta -exp(theta)))
    }
    initial <- function(n, theta)
        rep(0, 1)
    quit <- function(n, theta)
        return(invisible())
    cmd <- match.arg(cmd)
    return(
        do.call(
            cmd, 
            args = list(n = data$n,
                        theta = theta)))
}

library(spdep)
library(INLA)

inla.setOption(
    safe = FALSE
)

nxy <- c(3, 4)
nb <- grid2nb(d = nxy)
nnb <- card(nb)
n <- length(nnb)

m1data <- list(
    n = n,
    graph = Diagonal(n) +
        sparseMatrix(
            i = rep(1:n, nnb),
            j = unlist(nb[nnb>0]),
            x = 1)
)
m1data$C <- Diagonal(n, nnb)
m1data$G <- m1data$C - (m1data$graph-Diagonal(n))

m1 <- inla.rgeneric.define(
    model = rcar1,
    optimize = TRUE,
    data = m1data)

str(m1)

str(inla.rgeneric.q(m1, "initial"))
str(inla.rgeneric.q(m1, "mu", theta = c(0)))
str(inla.rgeneric.q(m1, "graph"))
str(inla.rgeneric.q(m1, "Q", theta = c(0)))

library(corGraphs)

m2 <- dag_model(
    graph = matrix(1,2,2),
    sigma.prior.reference = c(1,1),
    sigma.prior.probability = c(.05,.05),
    lambda = 1)

str(cgeneric_get(m2, "initial"))
str(cgeneric_get(m2, "mu", theta = c(0)))
str(cgeneric_get(m2, "graph"))
str(cgeneric_get(m2, "Q", theta = c(0)))

kmodel12 <- kronecker(
    m1,
    m2,
    data = m1data,
    debug = !TRUE)
str(kmodel12)

kmodel21 <- kronecker(
    m2,
    m1,
    data = m1data,
    debug = !TRUE)

str(kmodel21)

str(inla.rgeneric.q(kmodel12, "initial"))
str(inla.rgeneric.q(kmodel12, "mu", theta = c(0,0,0,0)))
str(inla.rgeneric.q(kmodel12, "graph"))
str(inla.rgeneric.q(kmodel12, "Q", theta = c(0,0,0,0)))

str(inla.rgeneric.q(kmodel21, "initial"))
str(inla.rgeneric.q(kmodel21, "mu", theta = c(0,0,0,0)))
str(inla.rgeneric.q(kmodel21, "graph"))
str(inla.rgeneric.q(kmodel21, "Q", theta = c(0,0,0,0)))

image(inla.rgeneric.q(m1, "graph"))
x11()
image(inla.rgeneric.q(kmodel12, "graph"))
x11()
image(inla.rgeneric.q(kmodel21, "graph"))

dataf <- data.frame(
    i = rep(1:n, 2),
    j = rep(1:2, each = n),
    idx = 1:(2*n),
    y = rpois(n * 2, 10)
)

out12 <- inla(
    y ~ f(idx, model = kmodel12), 
    family = 'poisson',
    data = dataf
)

out21 <- inla(
    y ~ f(idx, model = kmodel21), 
    family = 'poisson',
    data = dataf
)

