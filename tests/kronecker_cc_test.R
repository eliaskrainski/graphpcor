library(spdep)
library(corGraphs)

## Model 1: Besag over a grid
nxy <- c(40, 50)
nb <- grid2nb(d = nxy, queen = FALSE)
nnb <- card(nb)
n <- length(nnb)

nb.graph <- sparseMatrix(
    i = rep(1:n, nnb),
    j = unlist(nb[nnb>0]),
    dims = c(n, n),
    repr = "T"
)

R0 <- Diagonal(n, 1 + nnb) - nb.graph
R0[1:min(5, n), 1:min(20, n)]

## m1 definition
m1 <- cgeneric_generic0(
    debug = FALSE,
    R = R0,
    scale = FALSE,
    param = c(1, 0.5) 
)

initial(m1)

prior(m1, theta = -1)

prior(m1, theta = 1)

(theta1 <- initial(m1))

Q1 <- precision(m1, theta = theta1)
Q1[1:min(n,5), 1:min(n,10)]

## Model 2
m2.dtg <- dtg(
    p1 ~ p2 + c1 + c2,
    p2 ~ c3 - c4)
d2 <- dim(m2.dtg)
(n2 <- d2[1])

summary(m2.dtg)

## m2 definition
m2 <- cgeneric(
    m2.dtg, 
    debug = FALSE,
    sigma.prior.reference = rep(1, n2), 
    sigma.prior.probability = rep(.5, n2),
    lambda = 1)

(n2 <- m2$f$n)

initial(m2)
theta2 <- c(1,0,0,-1, -0.5,0.5)
Q2 <- precision(m2, theta = theta2)
Q2

solve(Q2)

cov2cor(solve(Q2))

Q12 <- kronecker(Q1, Q2)
Q21 <- kronecker(Q2, Q1)

ijo <- order(rep(1:n2, n))
all.equal(Q12[ijo, ijo], Q21)

## The M1 (x) M2 Kronecker product model definition
m1o <- m1; m1o$old = TRUE
m2o <- m2; m2o$old = TRUE

km12o <- kronecker(m1o, m2o)

km12 <- kronecker(m1, m2)

## The M2 (x) M1 Kronecker product model definition
km21o <- kronecker(m2o, m1o)
km21 <- kronecker(m2, m1)

str(initial(km12))

q12o <- precision(km12o, theta = c(theta1, theta2))
q21o <- precision(km21o, theta = c(theta2, theta1))

q12 <- precision(km12, theta = c(theta1, theta2))
q21 <- precision(km21, theta = c(theta2, theta1))

all.equal(Q12, q12o)
all.equal(Q21, q21o)

all.equal(Q12, q12)
all.equal(Q21, q21)

nscn <- 100
th12 <- matrix(rnorm(nscn * (length(theta1) + length(theta2))), nscn)

t12 <- sapply(1:nscn, function(i) {
    theta1 <- th12[i, 1:length(theta1)]
    theta2 <- th12[i, length(theta1) + 1:length(theta2)]
    Q1 <- precision(m1, theta = theta1)
    Q2 <- precision(m2, theta = theta2)
    Q12 <- kronecker(Q1, Q2)
    Q21 <- kronecker(Q2, Q1)    
    t0 <- Sys.time()
    q12 <- precision(km12o, theta = c(theta1, theta2))
    q21 <- precision(km21o, theta = c(theta2, theta1))
    r1 <- c(all.equal(Q12, q12), all.equal(Q21, q21))
    t1 <- Sys.time()
    q12 <- precision(km12, theta = c(theta1, theta2))
    q21 <- precision(km21, theta = c(theta2, theta1))
    t2 <- Sys.time()
    c(r1,
      all.equal(Q12, q12), 
      all.equal(Q21, q21),
      diff(c(t0, t1, t2)))
}
)

str(t12)

table(t12[1, ], t12[2, ])
table(t12[3, ], t12[4, ])

sum(t12[5,])
sum(t12[6,])


detach("package:corGraphs", unload = TRUE)
library(corGraphs)
