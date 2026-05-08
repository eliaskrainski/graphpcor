library(graphpcor)

## graph in Example 2.6 of the GMRF book
g <- graphpcor(x1~x2+x3, x4~x2+x3)
(ne <- dim(g))

## Laplacian
(G <- Laplacian(g))

## define initial L
theta.low <- rep(-3, ne[2])
idx.low <- (!is.zero(G) & lower.tri(G))
L1a <- diag(ne[1]:1); L1a[idx.low] <- theta.low
L1a

## the C code to fill.in
idx.fill <- which(
    is.zero(G) &
    (!is.zero(t(chol(G + diag(ne[1]))))))
nfi <- length(idx.fill)
nfi

if(nfi) {
    ij.fill <- cbind(
        row(G)[idx.fill],
        col(G)[idx.fill])
    print(ij.fill)
    Lf <- .C("fillL",
             as.integer(ne[1]),
             as.integer(nfi),
             as.integer(ij.fill[, 1]-1),
             as.integer(ij.fill[, 2]-1),
             l=L1a)$l
    print(Lf)
    Q0 <- crossprod(Lf)
} else {
    Q0 <- crossprod(L1a)
}

Q0

all.equal(
    Q0,
    crossprod(
        graphpcor:::Lprec0(
                        theta.low,
                        iLtheta = g,
                        d0 = ne[1]:1)))

V0 <- chol2inv(t(Lf))
C0 <- cov2cor(V0)
C0

all.equal(cov2cor(V0),
          vcov(g, theta = theta.low),
          check.attributes = FALSE)

all.equal(V0,
          vcov(g, theta = c(log(diag(V0))/2, theta.low)),
          check.attributes = FALSE)

## define the sigma parameters
logsigmas <- log(c(0.3, 0.7, 1.5, 0.9))

## this is used in the cgeneric C code
sr <- exp(logsigmas)/sqrt(diag(V0))
sr

V1 <- diag(sr) %*% V0 %*% diag(sr)
V1

stopifnot(all.equal(sqrt(diag(V1)), exp(logsigmas)))

## base model (theta for L)
theta0l <- rep(-0.3, ne[2])
basepcor(theta0l, iLtheta = g)

## build the cgeneric model
## Note: here 'model' is a 'graphpcor'
vmodel <- cgeneric(
    model = g, ## use the graphpcor
    lambda = 1,
    base = theta0l,
    sigma.prior.probability = rep(0.5, ne[1]),
    debug = 1e9,
    useINLAprecomp = FALSE)

cgeneric_prior(vmodel, theta = c(logsigmas, theta0l))

theta1 <- c(d = logsigmas,
            l = theta.low)
cgeneric_prior(vmodel, theta = theta1)

cgeneric_mu(vmodel)

cgeneric_initial(vmodel)

cgeneric_graph(vmodel) ## the structure of the precision

all.equal(solve(V1),
          as.matrix(cgeneric_Q(vmodel, theta = theta1)))

## Note: here 'model' is a 'graphpcor'
cmodel <- cgeneric(
    model = g, ## use the graphpcor
    lambda = 1,
    base = theta0l,
    debug = 0,
    useINLAprecomp = FALSE)

c(cgeneric_prior(cmodel, theta = theta0l), 
  cgeneric_prior(cmodel, theta = theta.low))
cgeneric_prior(cmodel, theta = cbind(theta0l, theta.low))

cgeneric_prior(cmodel, theta = rnorm(ne[2]*5, ne[2]))

all.equal(solve(C0),
          as.matrix(cgeneric_Q(cmodel, theta = theta.low)))

all.equal(cgeneric_Q(g, theta = theta.low),
          cgeneric_Q(cmodel, theta = theta.low))

dataf <- list(
    i = 1:ne[1],
    y = rep(NA, ne[1]))

library(INLA)

fit1 <- inla(
    formula = y ~ 0 + f(i, model = cmodel),
    control.family = list(hyper = list(
                              prec = list(intial = 20, fixed = TRUE))),
    data = dataf,
    control.mode = list(theta = theta.low, fixed = TRUE)
)

all.equal(cgeneric_Q(cmodel, theta = theta.low), prec(fit1))

