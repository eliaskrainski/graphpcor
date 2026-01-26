
library(spdep)
library(INLA)
library(graphpcor)

inla.setOption(
    num.threads = 1L,
    safe = FALSE
)

## Model 1: Besag over a grid
nxy <- c(40, 50)
nb <- grid2nb(d = nxy, queen = FALSE)
nnb <- card(nb)
n <- length(nnb)

nb.graph <- sparseMatrix(
    i = rep(1:n, nnb),
    j = unlist(nb[nnb>0]),
    x = 1,
    dims = c(n, n)
)

R0 <- inla.as.sparse(Diagonal(n, nnb) - nb.graph)
R0[1:min(5, n), 1:min(20, n)]

## scale first
cntr <- list(A = matrix(1, 1, n), e = 0)
Rs <- inla.scale.model(
    R0,
    constr = cntr
)

## m1 definition
m1 <- cgeneric(
    model = "generic0",
    R = Rs,
    scale = FALSE, ## scale could have been done here
    param = c(1, 0.0) ## to fix it at this value
)

prior(m1, theta = -1.0)
prior(m1, theta = +1.0)

initial(m1)

theta1 <- 0
Q1 <- prec(m1, theta = theta1)

summary(diag(inla.qinv(Q1 + Diagonal(n, 1e-5), cntr)))

## Model 2: 
m2.graph <- graphpcor(
    x1 ~ x2 + x3
)

d2 <- dim(m2.graph)
d2
n2 <- d2[1]

## m2 definition
m2 <- cgeneric(
    model = m2.graph, 
    sigma.prior.reference = rep(1, n2),
    sigma.prior.probability = rep(.5, n2), 
    lambda = 1)

str(initial(m2))
length(theta2 <- c(0.5, 0.3, 0, -1, 1))

Q2 <- prec(m2, theta = theta2)
Q2

solve(Q2)

cov2cor(as.matrix(solve(Q2)))
cov2cor(vcov(m2.graph, theta = theta2[-(1:n2)]))

## The M1 (x) M2 Kronecker product model definition
k12 <- kronecker(m1, m2)

## The M2 (x) M1 Kronecker product model definition
k21 <- kronecker(m2, m1)

### two ways of getting the precision matrix
Q12 <- kronecker(Q1, Q2)
q12 <- prec(k12, theta = c(theta2))
all.equal(Q12, q12)

Q21 <- kronecker(Q2, Q1)
q21 <- prec(k21, theta = c(theta2))
all.equal(Q21, q21)

## reorder test
ijo <- order(rep(1:n2, n))
all.equal(Q12[ijo, ijo], Q21)

ijo2 <- order(rep(1:n, n2))
all.equal(Q21[ijo2, ijo2], Q12)

## using Q1 (x) Q2 to sample
xx <- inla.qsample(
    n = 1,
    Q = Q12 + Diagonal(n*n2, 1e-9),
    constr = k12$f$extraconst)[, 1]
summary(xx)

diag(solve(Q2))
apply(t(matrix(xx, n2)), 2, sd)^2

cov2cor(as.matrix(solve(Q2)))
cor(t(matrix(xx, n2)))

dataf <- list(
    id1 = rep(1:n, each = n2),
    id2 = rep(1:n2, n),    
    idx = 1:(n2 * n), 
    y1 = new('numeric', xx)
)

## reorder y1 -> y2 (to use M2 (x) M1)
dataf$y2 <- as.vector(t(matrix(
    dataf$y1, n2, n)))

str(dataf)

all.equal(dataf$y1, as.vector(t(matrix(dataf$y2, n))))

cfam <- list(
    hyper = list(
        prec = list(
            initial = 10, fixed = TRUE
        )
    )
)
cmode <- list(theta = theta2)

## fit the model using data augmented
Nna <- rep(NA, n * n2)
dataex <- list(
    yy = cbind(c(dataf$y1, Nna), 
               rep(c(NA, 0), each = n * n2)),
    id1 = c(dataf$id1, Nna),
    id2 = c(dataf$id2, Nna),
    id1c = c(Nna, dataf$id1),
    id2c = c(Nna, dataf$id2),
    idx0 = c(Nna, 1:(n*n2)),
    w0 = rep(c(0, -1), each = n*n2)
)
str(dataex)

out2g <- inla(
    yy ~ 
        f(id2, model = m2, group = id1,
          control.group = list(model = 'besag', graph = nb.graph)) +
        f(id2c, copy = 'id2', group = id1c,
          hyper = list(theta = list(initial = -10, fixed = TRUE)),
          control.group = list(model = "iid")) +
        f(idx0, w0, model = 'iid',
          extraconstr = k12$f$extraconstr,
          hyper = list(theta = list(initial = -10, fixed = TRUE))), 
    data = dataex,
    family = rep("gaussian", 2),
    control.family = list(
        cfam,
        list(hyper = list(prec = list(initial = 10, fixed = TRUE)))),
    control.mode = cmode
)

out12 <- inla(
    y1 ~ f(idx, model = k21), 
    data = dataf,
    control.family = cfam,
    control.mode = cmode
)

out21 <- inla(
    y2 ~ f(idx, model = k12), 
    data = dataf,
    control.family = cfam,
    control.mode = cmode
)

rbind(out2g$cpu.used,
      out12$cpu.used,
      out21$cpu.used)
c(out2g$misc$nfunc,
  out12$misc$nfunc,
  out21$misc$nfunc)

rbind(out2g$summary.fix,
      out12$summary.fix,
      out21$summary.fix)

rbind(c(th = theta2),
      out2g$mode$theta,
      out12$mode$theta,
      out21$mode$theta)


diag(solve(Q2))
diag(solve(prec(m2, theta = out12$mode$theta)))
diag(solve(prec(m2, theta = out21$mode$theta)))

cov2cor(as.matrix(solve(Q2)))

cov2cor(vcov(m2.graph, theta = out12$mode$theta[-(1:n2)]))
cov2cor(vcov(m2.graph, theta = out21$mode$theta[-(1:n2)]))

summary(out12$summary.random$idx$mean)
summary(out21$summary.random$idx$mean)

sd(out12$summary.random$idx$mean)
sd(out21$summary.random$idx$mean)

tail(out12$logfile, 12)
tail(out21$logfile, 12)

grep("Number of constraints", out12$logfile, value = TRUE)
grep("Number of constraints", out21$logfile, value = TRUE)

grep("nnz", out12$logfile, value = TRUE)
grep("nnz", out21$logfile, value = TRUE)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
