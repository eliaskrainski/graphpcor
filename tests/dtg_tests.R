
library(graphpcor)

dtg(p1 ~ p2)
dtg(p1 ~ c2)

dtg(p1 ~ c1 + c2,
    p2 ~ c3)

dtg(p1 ~ c1 + c2,
    p2 ~ p1 + c2 + c3)

dtg(p1 ~ c1 + c2,
    p2 ~ p3 + c2 + c3)

dtg(p1 ~ p2 + c1 + c2,
    p2 ~ c2 + c3)

g1 <- dtg(p1 ~ c1 + c2 - c3)

g1

dim(g1)

summary(g1)

plot(g1)

precision(g1)

(q1 <- precision(g1, theta = c(0)))

v1 <- chol2inv(chol(q1))

v1

cov2cor(v1)

variance(g1)
variance(g1, theta = 0)
variance(g1, theta = -1)
variance(g1, theta = 1)

cov2cor(variance(g1))
cov2cor(variance(g1, theta = -1))
cov2cor(variance(g1, theta = 1))

g2 <- dtg(p1 ~ p2 + c1 + c2,
          p2 ~ c3 - c4)
g2
dim(g2)
summary(g2)

plot(g2)

precision(g2)
precision(g2, theta = c(0, 0))
precision(g2, theta = c(-1, 1))

solve(precision(g2))

solve(precision(g2, theta = c(0, 0)))
variance(g2)

chol2inv(chol(precision(g2, theta = c(0, 0))))[1:4, 1:4]
variance(g2, theta = c(0, 0))

g2
g3 <- dtg(p1 ~ -p2 + c1 + c2,
          p2 ~ -c3 + c4)
g3
dim(g3)
summary(g3)

plot(g3)

precision(g3)
precision(g3, theta = c(0, 0))

chol2inv(chol(precision(g3, theta = c(0, 0))))[1:4, 1:4]
variance(g3, theta = c(0, 0))

summary(g2)

g3
drop(g3) ## to be fixed (do not remove childrens!)

precision(g3)
precision(drop(g3))

n3 <- dim(g3)[1]
all.equal(
    solve(precision(g2, theta = c(0, 0)))[1:n3, 1:n3],
    solve(precision(g3, theta = c(0, 0)))[1:n3, 1:n3]
)

detach("package:graphpcor", unload = TRUE)
library(graphpcor)
