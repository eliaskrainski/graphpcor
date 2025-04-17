
library(INLA)
library(graphpcor)

inla.setOption(
    num.threads = 1L,
    safe = FALSE
)

## Model 1 graph
graph <- sparseMatrix(
    i = c(1, 2, 1, 3, 3, 5, 2, 4, 6, 5, 4, 6),
    j = c(2, 1, 3, 1, 5, 3, 4, 2, 5, 6, 6, 4)
)
graph*1
nnb <- colSums(graph)
n <- length(nnb)

R1 <- inla.as.sparse(Diagonal(n, 1+nnb) - graph)
R1[1:min(5, n), 1:min(20, n)]

m1 <- cgeneric(
    model = "generic0",
    R = R1,
    constr = FALSE,
    debug = !TRUE,
    scale = FALSE,
    param = c(1, 0.0),
    useINLAprecomp=FALSE
)

str(m1)
head(m1$f$cgeneric$data$smatrices)

theta1 <- 0

str(initial(m1))

str(graph(m1, optimize = TRUE))
str(graph(m1))

str(prec(m1))
str(prec(m1, optimize = TRUE))

str(prec(m1, optimize = TRUE, theta = -1))
str(prec(m1, optimize = TRUE, theta = 0))

if(FALSE) ### test model1
    inla(y ~ 0 + f(i, model = m1), "poisson",
         data = data.frame(y=rpois(n,1), i = 1:n))$cpu.used

Q1 <- prec(m1, theta = theta1)
Q1[1:min(5, n), 1:min(20, n)]

solve(Q1)
s1 <- exp(mean(log(diag(solve(Q1)))))

## a Second model graph:
graph2 <- sparseMatrix(
    i = c(1, 3, 1, 4, 2, 2, 4, 3),
    j = c(3, 1, 4, 1, 3, 4, 2, 2)
)
graph2*1
nnb2 <- colSums(graph2)
n2 <- length(nnb2)

R2 <- inla.as.sparse(Diagonal(n2, 1+nnb2) - graph2)
R2[1:min(5, n2), 1:min(20, n2)]

m2 <- cgeneric(
    model = "generic0",
    R = R2,
    scale = FALSE,
    constr = FALSE,
    param = c(1, 0.5),
    useINLAprecomp=FALSE)

initial(m2)

str(mu(m2))

str(graph(m2, optimize = TRUE))
str(graph(m2))

theta2 <- c(0)
str(Q2 <- prec(m2, theta = theta2))

Q2
solve(Q2)
s2 <- exp(mean(log(diag(solve(Q2)))))

if(TRUE) ## test model2
    inla(y ~ f(i, model = m2),
         data=list( i= 1:n2, y = rpois(n2,1)))$cpu.used

str(Q1Q2 <- kronecker(Q1, Q2))
str(Q2Q1 <- kronecker(Q2, Q1))

kmodel12 <- kronecker(
    m1,
    m2,
    debug = !TRUE)

Q12 <- prec(kmodel12, theta = c(theta2))
all.equal(Q1Q2, Q12)

str(kmodel12)
kmodel12$f$cgeneric$data$ints$n

head(kmodel12$f$cgeneric$data$smatrices$Kgraph)
length(kmodel12$f$cgeneric$data$smatrices$Kgraph)
(141-35) * (5-3) * 2 + 105

M <- kmodel12$f$cgeneric$data$smatrices$Kgraph[3]
M
summary(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+1:M])
summary(diff(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+1:M]))
summary(diff(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+M+1:M]))
summary(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+M+1:M])
summary(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+M+M+1:M])

gg12 <- inla.as.sparse(sparseMatrix(
    kmodel12$f$cgeneric$data$smatrices$Kgraph[3+1:M] + 1,
    kmodel12$f$cgeneric$data$smatrices$Kgraph[3+M+1:M] + 1,
    symmetric = TRUE
))

str(gg12)

image(gg12)

kmodel12$f$cgeneric$data$ints$idx1u
kmodel12$f$cgeneric$data$ints$idx2u
head(kmodel12$f$cgeneric$data$smatrices$Kgraph)
Me <- kmodel12$f$cgeneric$data$smatrices$Kgraph[3]
Me
range(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+1:Me])
range(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+Me+1:Me])
range(kmodel12$f$cgeneric$data$smatrices$Kgraph[3+Me*2+1:Me])

kmodel21 <- kronecker(
    m2,
    m1)

Q21 <- prec(kmodel21, theta = c(theta2))
all.equal(Q2Q1, Q21)

str(kmodel21)
kmodel21$f$cgeneric$data$ints$n

stopifnot(M == kmodel21$f$cgeneric$data$smatrices$Kgraph[3])
summary(kmodel21$f$cgeneric$data$smatrices$Kgraph[3+1:M])
summary(diff(kmodel21$f$cgeneric$data$smatrices$Kgraph[3+1:M]))
summary(kmodel21$f$cgeneric$data$smatrices$Kgraph[3+M+1:M])
summary(kmodel21$f$cgeneric$data$smatrices$Kgraph[3+M+M+1:M])

initial(kmodel12)
mu(kmodel12)
str(graph(kmodel12, optimize = TRUE))
str(graph(kmodel12))

image(Q12)

initial(kmodel21)
mu(kmodel21)
str(graph(kmodel21, optimize = TRUE))
str(graph(kmodel21))

image(graph(m1))
x11()
image(graph(m2))
x11()
image(graph(kmodel12))
x11()
image(graph(kmodel21))

kmodel12 <- kronecker(
    m1,
    m2,
    debug = FALSE)
kmodel21 <- kronecker(
    m2,
    m1,
    debug = FALSE)

data1 <- data.frame(
    i = 1:(n*n2),
    y = rep(NA, n*n2)
)

lkfixed <- list(
    hyper = list(
        prec = list(
            initial = 10, fixed = TRUE
        )
    )
)

out12 <- inla(
    y ~ 0 + f(i, model = kmodel12),
    data = data1,
    control.family = lkfixed,
    control.mode = list(theta = theta2, fixed = TRUE)
)
out21 <- inla(
    y ~ 0 + f(i, model = kmodel21), 
    data = data1,
    control.family = lkfixed,
    control.mode = list(theta = theta2, fixed = TRUE)
)

all.equal(prec(out12), Q1Q2)
all.equal(prec(out21), Q2Q1)

## with data
nrepl <- 1000

## using Q1 (x) Q2 to sample
xx <- array(inla.qsample(n = nrepl, Q = Q12), c(n2, n, nrepl))
xx1 <- matrix(aperm(xx, c(3, 1, 2)), n2*nrepl)
xx2 <- matrix(aperm(xx, c(3, 2, 1)), n*nrepl)

cov2cor(solve(Q1))
cor(xx1)

cov2cor(solve(Q2))
cor(xx2)


rbind(diag(solve(Q1)) * s2,
           diag(cov(xx1)))

rbind(diag(solve(Q2)) * s1,
           diag(cov(xx2)))


data2 <- data.frame(
    i = rep(1:(n*n2), nrepl),
    r = rep(1:nrepl, each = n*n2),
    y1 = as.vector(xx1),
    y2 = as.vector(xx2)
)

out12 <- inla(
    y2 ~ f(i, model = kmodel12, replicate = r),
    data = data2,
    control.family = lkfixed
)

out21 <- inla(
    y1 ~ f(i, model = kmodel21, replicate = r),
    data = data2,
    control.family = lkfixed
)

rbind(out12$cpu.used,
      out21$cpu.used)

c(out12$misc$nfunc, out21$misc$nfunc)

rbind(out12$summary.fix,
      out21$summary.fix)

cbind(c(th=theta2),
      out12$mode$theta,
      out21$mode$theta)

cbind(true=c(theta2),
      out12$summary.hy[, c(1,2,3,5)])

cbind(true=c(theta2),
      out21$summary.hy[, c(1,2,3,5)])

grep("nnz", out12$logfile, value = TRUE)
grep("nnz", out21$logfile, value = TRUE)

plot(out12$internal.marginals.hyperpar[[1]], type = "l",
     bty = 'n')
lines(out21$internal.marginals.hyperpar[[1]], col = 2, lty = 2)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
