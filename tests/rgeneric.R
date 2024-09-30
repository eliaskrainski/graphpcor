
rcar1 <- function(cmd = c("graph", "Q", "mu", "initial",
                           "log.norm.const",
                           "log.prior", "quit"),
                   theta = NULL, ...) {
    graph <- function(n, theta)
        return(Graph)
    Q <- function(n, theta) {
        qq <- INLA::inla.as.sparse(
        (Dmat * exp(theta[1]) + Gmat) * 1)
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
            args = 
                list(n = n,
                     theta = theta)
        )
    )
}

library(graphpcor)
library(spdep)
library(INLA)

inla.setOption(
    safe = FALSE
)

nxy <- c(30, 50)
nb <- grid2nb(d = nxy)
nnb <- card(nb)
n <- length(nnb)

graph <- sparseMatrix(
    i = rep(1:n, nnb),
    j = unlist(nb[nnb>0]),
    )

m1 <- rgeneric(
    model = rcar1,
    n = n,
    Graph = graph,
    Dmat = Diagonal(n, nnb),
    Gmat = Diagonal(n, nnb) - graph
)

initial(m1)

theta1 <- 0
mu(m1, theta = theta1)

image(graph(m1))

Q1 <- precision(m1, theta = c(0))
image(Q1)

## model 2
m2g <- dtg(
    p1 ~ c1 + c2 + p2,
    p2 ~ c3
)

m2 <- cgeneric(
    model = m2g,
    sigma.prior.reference = c(1,1,1),
    sigma.prior.probability = c(.05,.05,0.05),
    lambda = 1)
(n2 <- m2$f$n)

initial(m2)
mu(m2)
graph(m2, optimize = TRUE)
graph(m2)
precision(m2)

theta2 <- c(seq(1, -1, length = n2), -0.5, 0.0)
Q2 <- precision(m2, theta = theta2)
Q2
solve(Q2)
cov2cor(solve(Q2))

k12 <- kronecker(
    m1,
    m2,
    debug = !TRUE)
str(k12)

k21 <- kronecker(
    m2,
    m1,
    debug = !TRUE)

str(k21)

initial(k12)
initial(k21)

mu(k12)

str(graph(k12))
str(graph(k21))

str(precision(k12))
str(precision(k21))

image(graph(m1))
#x11()
image(graph(k12))
#x11()
image(graph(k21))

qq12 <- as(kronecker(Q1, Q2), 'symmetricMatrix')
Q12 <- as(precision(k12, theta = c(theta1, theta2)), 'symmetricMatrix')
all.equal(qq12, Q12)

qq21 <- as(kronecker(Q2, Q1), 'symmetricMatrix')
Q21 <- as(precision(k21, theta = c(theta2, theta1)), 'symmetricMatrix')
all.equal(qq21, Q21)

initial(k12)

## using Q2 (x) Q1 to sample
x2 <- inla.qsample(n = 1, Q = Q21)[, 1]

cor(t(matrix(x2, n2)))

dataf <- data.frame(
    idx = 1:(n2 * n),
    y2 = x2
)
## reorder y2 -> y1
dataf$y1 <- as.numeric(t(matrix(dataf$y2, n)))

str(dataf)

cfam <- list(hyper = list(prec = list(initial = 20, fixed = TRUE)))

out12 <- inla(
    y1 ~ 0 + f(idx, model = k12), 
    data = dataf,
    control.family = cfam,
    control.mode = list(
        theta = c(theta1, theta2),
        fixed = TRUE
    )
)

out21 <- inla(
    y2 ~ 0 + f(idx, model = k21), 
    data = dataf,
    control.family = cfam,
    control.mode = list(
        theta = c(theta2, theta1),
        fixed = TRUE)
)

all.equal(Q12, precision(out12))
all.equal(Q21, precision(out21))

out12 <- inla(
    y1 ~ 0 + f(idx, model = k12), 
    data = dataf,
    control.family = cfam,
    control.mode = list(
        theta = c(theta1, theta2),
        fixed = !TRUE
    )
)

out21 <- inla(
    y2 ~ 0 + f(idx, model = k21), 
    data = dataf,
    control.family = cfam,
    control.mode = list(
        theta = c(theta2, theta1),
        fixed = !TRUE)
)


rbind(out12$cpu.used,
      out21$cpu.used)

nth1 <- length(theta1)
nth2 <- length(theta2)

rbind(true=c(theta1, theta2), out12$mode$theta)
rbind(true=c(theta2, theta1), out21$mode$theta)

tail(out12$logfile, 12)
tail(out21$logfile, 12)

grep("nnz", out12$logfile, value = TRUE)
grep("nnz", out21$logfile, value = TRUE)


detach("package:graphpcor", unload = TRUE)
library(graphpcor)
