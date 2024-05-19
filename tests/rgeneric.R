
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
        return(sum(-theta)) #G(1,1)
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

nxy <- c(10, 15)
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
theta1 <- 0
str(inla.rgeneric.q(m1, "mu", theta = theta1))
str(inla.rgeneric.q(m1, "graph"))
str(Q1 <- inla.rgeneric.q(m1, "Q", theta = theta1))

library(corGraphs)

m1.graph <- list(
    c1~c2,
    c2~c3)

m2 <- dag_model(
    graph = m1.graph, 
    sigma.prior.reference = c(1,1,1),
    sigma.prior.probability = c(.05,.05,0.05),
    lambda = 1)

str(cgeneric_get(m2, "initial"))
str(cgeneric_get(m2, "mu", theta = c(0)))
str(cgeneric_get(m2, "graph"))
cgeneric_get(m2, "graph", optimize = !TRUE)
theta2 <- c(1:-1, -0.5, 0.5)
str(cgeneric_get(m2, "Q", theta = theta2))
Q2 <- cgeneric_get(m2, "Q", theta = theta2, optimize = !TRUE)
Q2
solve(Q2)
cov2cor(solve(Q2))

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
str(inla.rgeneric.q(kmodel12, "mu", theta = c(theta1, theta2)))
str(inla.rgeneric.q(kmodel12, "graph"))
str(inla.rgeneric.q(kmodel12, "Q", theta = c(theta1, theta2)))

str(inla.rgeneric.q(kmodel21, "initial"))
str(inla.rgeneric.q(kmodel21, "mu", theta = c(theta2, theta1)))
str(inla.rgeneric.q(kmodel21, "graph"))
str(inla.rgeneric.q(kmodel21, "Q", theta = c(theta2, theta1)))

image(inla.rgeneric.q(m1, "graph"))
x11()
image(inla.rgeneric.q(kmodel12, "graph"))
x11()
image(inla.rgeneric.q(kmodel21, "graph"))

qq12 <- kronecker(Q1, Q2)
Q12 <- inla.rgeneric.q(kmodel12, "Q", theta = c(theta1, theta2))
all.equal(qq12, Q12)

qq21 <- kronecker(Q2, Q1)
Q21 <- inla.rgeneric.q(kmodel21, "Q", theta = c(theta2, theta1))
all.equal(qq21, Q21)

## using Q2 (x) Q1 to sample
x1 <- inla.qsample(n = 1, Q = Q21)[, 1]

dataf <- data.frame(
    idx = 1:(2 * n),
    y2 = rpois(2 * n, exp(3 + x1))
)
## reorder y2 -> y1
dataf$y1 <- as.integer(t(matrix(dataf$y2, n)))

str(dataf)

out12 <- inla(
    y1 ~ f(idx, model = kmodel12), 
    family = 'poisson',
    data = dataf
)

out21 <- inla(
    y2 ~ f(idx, model = kmodel21), 
    family = 'poisson',
    data = dataf
)

rbind(out12$cpu.used,
      out21$cpu.used)

rbind(out12$summary.fix,
      out21$summary.fix)


cbind(true=c(theta1, theta2), out12$summary.hy[, c(1,2,3,5)])
cbind(true=c(theta2, theta1), out21$summary.hy[, c(1,2,3,5)])

tail(out12$logfile, 12)
tail(out21$logfile, 12)

grep("nnz", out12$logfile, value = TRUE)
grep("nnz", out21$logfile, value = TRUE)
