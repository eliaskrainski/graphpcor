
library(graphpcor)

treepcor(p1 ~ p2)
treepcor(p1 ~ c2)

treepcor(p1 ~ c1 + c2,
    p2 ~ c3)

treepcor(p1 ~ c1 + c2,
    p2 ~ p1 + c2 + c3)

treepcor(p1 ~ c1 + c2,
    p2 ~ p3 + c2 + c3)

treepcor(p1 ~ p2 + c1 + c2,
    p2 ~ c2 + c3)

g1 <- treepcor(p1 ~ c1 + c2 - c3)

g1

dim(g1)

summary(g1)

plot(g1)

prec(g1)

(q1 <- prec(g1, theta = c(0)))

v1 <- chol2inv(chol(q1))

v1

cov2cor(v1)

vcov(g1)
vcov(g1, theta = 0)
vcov(g1, theta = -1)
vcov(g1, theta = 1)

cov2cor(vcov(g1))
cov2cor(vcov(g1, theta = -1))
cov2cor(vcov(g1, theta = 1))

g2 <- treepcor(p1 ~ p2 + c1 + c2,
          p2 ~ c3 - c4)
g2
dim(g2)
summary(g2)

plot(g2)

prec(g2)
prec(g2, theta = c(0, 0))
prec(g2, theta = c(-1, 1))

solve(prec(g2))

solve(prec(g2, theta = c(0, 0)))
vcov(g2)

chol2inv(chol(prec(g2, theta = c(0, 0))))[1:4, 1:4]
vcov(g2, theta = c(0, 0))

g2
g3 <- treepcor(p1 ~ -p2 + c1 + c2,
          p2 ~ -c3 + c4)
g3
dim(g3)
summary(g3)

plot(g3)

prec(g3)
prec(g3, theta = c(0, 0))

chol2inv(chol(prec(g3, theta = c(0, 0))))[1:4, 1:4]
vcov(g3, theta = c(0, 0))

summary(g2)

g3
drop1(g3) ## to be fixed (do not remove childrens!)

prec(g3)
prec(drop1(g3))

n3 <- dim(g3)[1]
all.equal(
    solve(prec(g2, theta = c(0, 0)))[1:n3, 1:n3],
    solve(prec(g3, theta = c(0, 0)))[1:n3, 1:n3]
)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
