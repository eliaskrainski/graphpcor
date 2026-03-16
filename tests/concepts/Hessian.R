
library(graphpcor)

g4 <- graphpcor(x1~x2+x3, x4~x2+x3)

p <- 4
m <- 4
theta <- rnorm(m)

H4 <- hessian(g4, theta)
H4

c4 <- basepcor(theta, p, iLtheta = g4)
all.equal(H4, hessian(c4))

c4$base

c2 <- basepcor(theta, p, itheta = g4,
               iunknown = c(1,3))

hessian(c2)
H4[c(1,3), c(1,3)]

c22 <- basepcor(theta[c(1,3)], p,
                itheta = g4, iparam = c(1,1,2,2))
hessian(c22)

c1 <- basepcor(theta[c(1,3)], p, itheta = g4,
               iparam = c(1,1,2,2), iunknown = 1)
hessian(c1)

graphpcor:::dspd(H4)

###

gs1 <- graphpcor(x1~x2+x3+x4+x5)
G1 <- attr(gs1, "graph")
G1

gs2 <- graphpcor(G1[c(2:5,1), c(2:5,1)])
G2 <- attr(gs2, "graph")
G2

th1 <- c(-2,-1,-0.5,0.5)
th1 <- rep(-1,4)
C1 <- basepcor(th1, iLtheta = gs1)
C1

C2 <- basepcor(C1$base[c(2:5,1), c(2:5,1)], iLtheta = gs2)
C2

H1 <- hessian(gs1, C1$theta)
H2 <- hessian(gs2, C2$theta)

round(H1, 5)
round(H2, 5)

str(graphpcor:::dspd(H1))
str(graphpcor:::dspd(H2))
