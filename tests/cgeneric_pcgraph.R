library(INLA)
library(graphpcor)

g1 <- corgraph(x ~ y, y ~ v, v ~ z, z ~ x)

g1

dim(g1)

summary(g1)

G1 <- precision(g1)
G1

n <- ncol(G1)

L1 <- t(chol(G1 + diag(n)))
L1

L0 <- L1; L0[G1==0] <- 0
sum((L0==0) & (L1!=0))
L0

ne <- sum(L0!=0)-n
ne

## the R code
all.equal(graphpcor:::fiL(L0,n+4), L1)

## the C code
Lf <- .C("fillL", as.integer(4), as.integer(1),
         as.integer(3), as.integer(1), l=L0)$l
all.equal(Lf, L1)

## build the model
cmodel <- cgeneric(
    model = "pcgraph",
    graph = g1,
    lambda = 1,
    theta.base = rep(-1, ne),
    sigma.prior.reference = rep(1, n),
    sigma.prior.probability = rep(0.5, n),
    debug = 1e9)

all.equal(
    cmodel,
    cgeneric(
        model = "pcgraph",
        graph = G1,
        lambda = 1,
        theta.base = rep(-1, ne),
        sigma.prior.reference = rep(1, n),
        sigma.prior.probability = rep(0.5, n),
        debug = 1e9)
)

graph(cmodel)

theta1 <- c(log(diag(L0)), L0[(L0!=0) & lower.tri(L0, diag = FALSE)])
theta1

Q.a <- precision(g1, theta = theta1)
Q.a

all.equal(Q.a, G1 + diag(n))

Q.c <- as.matrix(precision(cmodel, theta = theta1))
Q.c

all.equal(Q.c,
          matrix(G1, n, n))



