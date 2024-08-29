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

rho <- 0.9
Q2 <- inla.as.sparse(
    sparseMatrix(
        i = c(1,1, 2,2,2, 3,3),
        j = c(1,2, 1,2,3, 2,3),
        x = c(1,  -rho,
              -rho,  1 + rho^2,  -rho,
              -rho,  1)
    )
)
Q2

(n2 <- nrow(Q2))
cov2cor(solve(Q2))

Q21 <- inla.as.sparse(kronecker(Q2, Q1))
Q21

(theta0 <- c(0, -1+2/(1+exp(-rho))))

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
              control.group = list(model = 'ar1')),
    data = dataf,
    control.family = cfam,
    control.mode = cmode,
    control.compute = ccpt
)
Qinla <- inla.as.sparse(ires1$misc$config$config[[1]]$Qprior)

idx.u <- which(Q21@i <= Q21@j)
Q21u <- inla.as.sparse(
    sparseMatrix(
        i = Q21@i[idx.u] + 1L,
        j = Q21@j[idx.u] + 1L,
        x = Q21@x[idx.u]
    )
)

all.equal(Q21u, Qinla)

round(Q21, 2)
round(Q21u, 2)
round(Qinla, 2)

