
library(spdep)
library(INLA)
library(corGraphs)

inla.setOption(
    num.threads = 1L,
    safe = FALSE
)

## Model 1: iid
n <- 100

## m1 definition
m1 <- cgeneric_iid(
    n = n, 
    param = c(1, 0.0) ## to fix it at this value
)

cgeneric_get(m1, "log.prior", theta = -1.0)
cgeneric_get(m1, "log.prior", theta = +1.0)
cgeneric_get(m1, "initial")

theta1 <- 0
Q1 <- cgeneric_get(m1, "Q", theta = theta1, optimize = FALSE)

## Model 2
m2.graph <- list(
    p1~p2+c1-c2,
    p2~c3+c4
)

## m2 definition
m2 <- dcg_model(
    dcg = m2.graph, 
    sigma.prior.reference = c(1,1,1,1),
    sigma.prior.probability = c(.5,.5,0.5,0.5),
    lambda = 1)
(n2 <- m2$f$n)

str(cgeneric_get(m2, "initial"))
length(theta2 <- c(0.7,0.5,0.2,0.6, 0, 1))
Q2 <- cgeneric_get(m2, "Q", theta = theta2, optimize = FALSE)
Q2

solve(Q2)
cov2cor(solve(Q2))

qq12 <- kronecker(Q1, Q2)
qq21 <- kronecker(Q2, Q1)

## The M1 (x) M2 Kronecker product model definition
kmodel12 <- kronecker(m1, m2)

## The M2 (x) M1 Kronecker product model definition
kmodel21 <- kronecker(m2, m1)

### two ways of getting the precision matrix
Q12 <- cgeneric_get(kmodel12, "Q", theta = c(theta2), optimize = FALSE)
all.equal(qq12, Q12)

Q21 <- cgeneric_get(kmodel21, "Q", theta = c(theta2), optimize = FALSE)
all.equal(qq21, Q21)

## reorder test
ijo <- order(rep(1:n2, n))
all.equal(qq12[ijo, ijo], qq21)

ijo2 <- order(rep(1:n, n2))
all.equal(qq21[ijo2, ijo2], qq12)

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
    y1 ~ f(idx, model = kmodel12), 
    data = dataf,
    control.family = cfam,
    control.mode = cmode
)

out21 <- inla(
    y2 ~ f(idx, model = kmodel21), 
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

summary(out12$summary.random$idx$mean)
summary(out21$summary.random$idx$mean)

tail(out12$logfile, 12)
tail(out21$logfile, 12)

grep("nnz", out12$logfile, value = TRUE)
grep("nnz", out21$logfile, value = TRUE)

detach("package:corGraphs", unload = TRUE)
library(corGraphs)
