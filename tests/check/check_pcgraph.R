library(graphpcor)

g <- corgraph(x ~ y, y ~ v, v ~ z, z ~ x)

g

ne <- dim(g)
ne

summary(g)

## sometimes we need it
G <- Laplacian(g)
G

## alternatively
all.equal(G,
          Laplacian(corgraph(x~y+z,y~v,v~z)))

## now do the Cholesky of some matrix
Q1 <- G + diag(ne[1])
L1 <- t(chol(Q1))
L1

## take out the fill-in
L1a <- L1; L1a[is.zero(G)] <- 0
L1a

## define theta from L1a
theta1 <- c(
    log(diag(L1a)), ## diagonal
    L1a[lower.tri(L1a) & (!is.zero(L1a))] ## lower & 0
)
theta1

Q.a <- precision(g, theta = theta1)
Q.a

## the fill-in 
idx.fill <- which(is.zero(L1a) & (!is.zero(L1)))
idx.fill
nfi <- length(idx.fill)
nfi

## the R code to fill-in
all.equal(graphpcor:::fiL(L1a, idx.fill), L1)

## the C code to fill.in
ij.fill <- cbind(
    row(Q1)[idx.fill],
    col(Q1)[idx.fill])
ij.fill
Lf <- .C("fillL",
         as.integer(ne[1]),
         as.integer(nfi),
         as.integer(ij.fill[, 1]-1),
         as.integer(ij.fill[, 2]-1),
         l=L1a)$l
all.equal(Lf, L1)

## theta (lower L) for the base model
theta0l <- rep(-1, ne[2])

## what it gives
Q0 <- precision(
    g,
    theta = c(rep(0, ne[1]), ## diag
              theta0l)) ## lower L
C0 <- cov2cor(solve(Q0))
C0

## the 'iid' case would be
cov2cor(solve(precision(g, theta = rep(0, sum(ne)))))

## build the model
cmodel <- cgeneric(
    model = "pcgraph",
    graph = g, ## use the corgraph
    lambda = 1,
    base = theta0l,
    sigma.prior.reference = rep(1, ne[1]),
    sigma.prior.probability = rep(0.5, ne[1]))

## anothe way
all.equal(
    cmodel,
    cgeneric(
        model = "pcgraph",
        graph = G, ## using G
        lambda = 1,
        base = theta0l,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## yet another way
all.equal(
    cmodel,
    cgeneric(
        model = "pcgraph",
        graph = G!=0, ## any binary matrix works
        lambda = 1,
        base = theta0l,
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

## specify the base model from C0
all.equal(
    cmodel,
    cgeneric(
        model = "pcgraph",
        graph = Q1, ## using Q1
        lambda = 1,
        base = C0, 
        sigma.prior.reference = rep(1, ne[1]),
        sigma.prior.probability = rep(0.5, ne[1]))
)

graph(cmodel)

##
all.equal(Q.a, Q1)

Q.c <- precision(cmodel, theta = theta1)
Q.c

all.equal(Q1, as.matrix(Q.c))

prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = rnorm(sum(ne)))
prior(cmodel, theta = theta1)

dataf <- list(
    i = 1:ne[1],
    y = rep(NA, ne[1]))

library(INLA)

fit1 <- inla(
    formula = y ~ 0 + f(i, model = cmodel),
    family = 'poisson',
    data = dataf,
    control.mode = list(theta = theta1, fixed = TRUE)
)

all.equal(Q.c, precision(fit1))

## some data
nrep <- 3000
nd <- nrep * ne[1]

xx <- matrix(rnorm(nd), nrep) %*% chol(C0)
cov(xx)

theta.y <- log(5)
datar <- data.frame(
    r = rep(1:nrep, each = ne[1]),
    i = rep(1:ne[1], nrep),
    y = rnorm(nd, 1 + xx, exp(-2*theta.y))
)

m1 <- y ~ f(i, model = cmodel, replicate = r)
fit1 <- inla(
    formula = m1,
    data = datar,
    control.inla = list(int.strategy = 'eb')
)

## transform back thetaL
Hd <- graph2H(g, theta0l)
th0.til <- fit1$mode$theta[1+ne[1]+1:ne[2]]
th0.fit <- Hd$hneg.5 %*% th0.til + theta0l

th.fit <- c(fit1$mode$theta[2:(1+ne[1])], th0.fit)

Qfit <- precision(cmodel, theta = th.fit)

Q0
Qfit

C0
cov2cor(solve(Qfit))
