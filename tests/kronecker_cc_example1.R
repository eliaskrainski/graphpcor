
library(INLA)
library(graphpcor)

inla.setOption(
    num.threads = 1L,
    safe = FALSE
)

## Model 1: iid
n <- 100

## m1 definition
m1 <- cgeneric(
    model = 'iid',
    n = n, 
    param = c(1, 0.0) ## to fix it at this value
)

prior(m1, theta = -1.0)
prior(m1, theta = +1.0)

initial(m1)

theta1 <- 0
Q1 <- precision(m1, theta = theta1)
if(Q1@uplo == "L")
    Q1 <- t(Q1)

## Model 2
m2.graph <- cortree(
    p1 ~ p2 + c1 - c2,
    p2 ~ c3 + c4
)

## m2 definition
m2 <- cgeneric(
    model = m2.graph,
    sigma.prior.reference = c(1,1,1,1),
    sigma.prior.probability = c(.5,.5,0.5,0.5),
    lambda = 1)
(n2 <- m2$f$n)

initial(m2)

length(theta2 <- c(0.7,0.5,0.2,0.6, 0, 1))
Q2 <- precision(m2, theta = theta2)
Q2

solve(Q2)
cov2cor(solve(Q2))

Q12 <- kronecker(Q1, Q2)
Q21 <- kronecker(Q2, Q1)

## The M1 (x) M2 Kronecker product model definition
k12 <- kronecker(m1, m2)

## The M2 (x) M1 Kronecker product model definition
k21 <- kronecker(m2, m1)

### two ways of getting the precision matrix
q12 <- precision(k12, theta = c(theta2))
all.equal(Q12, q12)

q21 <- precision(k21, theta = c(theta2))
all.equal(Q21, q21)

## reorder test
ijo <- order(rep(1:n2, n))
all.equal(q12[ijo, ijo], q21)

ijo2 <- order(rep(1:n, n2))
all.equal((q21)[ijo2, ijo2], q12)

## using Q1 (x) Q2 to sample
xx <- inla.qsample(n = 1, Q = Q12)[, 1]
summary(xx)

cov(t(matrix(xx, n2)))/solve(Q2)-1

cov2cor(solve(Q2))
cor(t(matrix(xx, n2)))

dataf <- data.frame(
    id1 = rep(1:n, each = n2),
    id2 = rep(1:n2, n),
    idx = 1:(n * n2), 
    y1 = xx
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
cmode <- list(theta = rep(0, length(theta2)))##theta2)

out2r <- inla(
    y1 ~ f(id2, model = m2, replicate = id1),
    data = dataf,
    control.family = cfam,
    control.mode = cmode
)

out12 <- inla(
    y1 ~ f(idx, model = k12), 
    data = dataf,
    control.family = cfam,
    control.mode = cmode
)

out21 <- inla(
    y2 ~ f(idx, model = k21), 
    data = dataf,
    control.family = cfam,
    control.mode = cmode
)

rbind(out2r$cpu.used,
      out12$cpu.used,
      out21$cpu.used)
c(out2r$misc$nfunc,
  out12$misc$nfunc,
  out21$misc$nfunc)

rbind(out2r$summary.fix,
      out12$summary.fix,
      out21$summary.fix)

rbind(c(th=theta2),
      out2r$mode$theta,
      out12$mode$theta,
      out21$mode$theta)

cbind(true=theta2,
      out12$summary.hy[, c(1,2,3,5)])
cbind(true=theta2, 
      out21$summary.hy[, c(1,2,3,5)])

tail(out12$logfile, 12)
tail(out21$logfile, 12)

grep("nnz", out12$logfile, value = TRUE)
grep("nnz", out21$logfile, value = TRUE)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
