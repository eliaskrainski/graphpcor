## Order of the kronecker product
##  with Q = Q1 (x) Q2, or Q = Q2 (x) Q1.
## We use Q1 in 'generic0' model and
##  Q2 in the group model as 'rw1'

library(INLA)

inla.setOption(
    num.threads = 1L,
    safe = FALSE
)

## Model 1: Besag over a grid
graph <- sparseMatrix(
    i = c(2, 3, 1, 4, 1, 4, 5, 2, 3, 3),
    j = c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5)
)
graph
(n <- nrow(graph))

Q1 <- inla.as.sparse(
    Diagonal(n, 1 + rowSums(graph)) - graph)
Q1[1:min(10, n), 1:min(20, n)]

n2 <- 3
Q2 <- INLA:::inla.rw1(n2)
Q2

Q21 <- as(kronecker(Q2, Q1), "symmetricMatrix")
Q21

(theta0 <- 0)

cfam <- list(
    hyper = list(
        prec = list(
            initial = 10,
            fixed = TRUE)))
cmode <- list(
    theta = theta0,
    fixed = TRUE,
    restart = FALSE)
ccpt <- list(config = TRUE)

dataf <- data.frame(
    y = NA,
    i = rep(1:n, each = n2),
    j = rep(1:n2, n))

ires1 <- inla(
    y ~ 0 + f(i, model = "generic0", Cmatrix = Q1, group = j,
              control.group = list(
                  model = 'rw1',
                  scale.model = FALSE)),
    data = dataf,
    control.family = cfam,
    control.mode = cmode,
    control.compute = ccpt
)
Qinla <- corGraphs::precision(ires1)

all.equal(Q21, Qinla)

