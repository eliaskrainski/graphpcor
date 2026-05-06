
## Define a model with precision is Q = Q1 (x) Q2
## Use two ways to define it and fit using inla
## Test the cgeneric for each one testing it with
##  the "same" model already implemented in inla
## We set m1 as the 'generic0', with Q1 = 1 * R
##  and Q2 from the 'ar1', group model in inla.
## Note: the m2 model with theta = c(0, 0, x) is
##  equal the (inla) group ar1 model with
##    x = log((1+rho)/(1-rho)) 

library(INLA)
library(graphpcor)

inla.setOption(
    num.threads = 1L,
    safe = FALSE
)

## Model 1 graph
graph <- sparseMatrix(
    i = c(2, 3, 1, 4, 1, 4, 5, 2, 3, 3),
    j = c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5)
)
graph 
(n <- nrow(graph))

Q1 <- as(inla.as.sparse(
    Diagonal(n, 1 + rowSums(graph)) - graph),
    'symmetricMatrix')
Q1[1:min(10, n), 1:min(20, n)]

cfam <- list(
    hyper = list(
        prec = list(
            initial = 10,
            fixed = TRUE)))

out1 <- inla(
    y ~ 0 + f(i, model = 'generic0', Cmatrix = Q1),
    data = list(y = rep(NA, n),
                i = 1:n),
    control.family = cfam,
    control.mode = list(
        theta = c(0),
        fixed = TRUE)
)

all.equal(Q1, cgeneric_Q(out1))

## model 2 definition
m2 <- treepcor(
    p1 ~ c1 - c2
)
m2
dim(m2)
summary(m2)

cgeneric_Q(m2, theta = c(0))
solve(cgeneric_Q(m2, theta = c(0)))
vcov(m2, theta = c(0))

theta.p <- c(0.33)
(V2 <- vcov(m2, theta = theta.p))

C2 <- cov2cor(V2)
C2

m2.cg <- cgeneric(
    m2,
    sigma.prior.reference = c(1, 1),
    sigma.prior.probability = c(0.5, 0.5),
    lambda = 1
)

(n2 <- m2.cg$f$n)

theta.test <- c(seq(-1, 1, length = n2), theta.p)
Q2test <- cgeneric_Q(m2.cg, theta = theta.test)

solve(Q2test)
crossprod(C2 %*% diag(exp(theta.test[1:n2])),
          diag(exp(theta.test[1:n2])))

cov2cor(as.matrix(solve(Q2test)))
C2

theta2 <- c(0.0, 0.0, theta.p) ## unit variance
Q2 <- cgeneric_Q(m2.cg, theta = theta2)
Q2
solve(Q2)

out2 <- inla(
    y ~ 0 + f(i, model = m2.cg),
    data = list(y = rep(NA, n2),
                i = 1:n2),
    control.family = cfam,
    control.mode = list(
        theta = theta2,
        fixed = TRUE)
)
all.equal(Q2, cgeneric_Q(out2))

rho = C2[1, 2]
Q2 * (1 - rho^2)
solve(Q2 * (1 - rho^2))
cov2cor(as.matrix(solve(Q2 * (1 - rho^2))))

Q21 <- as(inla.as.sparse(kronecker(Q2, Q1)),
          'symmetricMatrix')
Q21r <- as(inla.as.sparse(kronecker(Q2 * (1 - rho^2), Q1)),
           'symmetricMatrix')
Q21[1:min(5, n*n2), 1:min(10, n*n2)]

(theta.fixed <- c(
     0 + log(1 - rho^2),
     log((1+rho)/(1-rho))))

cmode <- list(
    theta = theta.fixed,
    fixed = TRUE,
    restart = FALSE)

dataf <- data.frame(
    y = NA,
    i = rep(1:n, each = n2),
    j = rep(1:n2, n))

ires1 <- inla(
    y ~ 0 + f(i, model = "generic0", Cmatrix = Q1, group = j,
              control.group = list(model = 'ar1')),
    data = dataf,
    control.family = cfam,
    control.mode = cmode
)

Qinla1 <- cgeneric_Q(ires1)

all.equal(Sparse(Q21r), Sparse(Qinla1))

## setup the 'cgeneric' kronecker model
m1 <- cgeneric(
    model = 'generic0',
    R = Q1,
    scale = FALSE,
    param = c(1, 0.0))

Q1c <- Sparse(cgeneric_Q(m1, theta = 0))

all.equal(
    cgeneric_Q(out1),
    Q1c
)

k21 <- kronecker(m2.cg, m1)
q21 <- cgeneric_Q(k21, theta = c(theta2))

all.equal(q21, Q21)

ires2 <- inla(
    y ~ 0 + f(i, model = k21),
    data = data.frame(y = NA, i = 1:(n * n2)),
    control.family = cfam,
    control.mode =
        list(theta = c(theta2),
             fixed = TRUE)
)

Qinla2 <- cgeneric_Q(ires2)

all.equal(Sparse(q21), Qinla2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
