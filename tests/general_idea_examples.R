library(INLA)
library(corGraphs)

g1 <- list(c1~c2, c2~c3, c3~c4, c4~c5, c5~c6)
G1 <- graph_Laplacian(g1)
G1

g2 <- c(g1, c1~c6)
G2 <- graph_Laplacian(g2)
G2

g3 <- list(c1~c2, c1~c3, c1~c4, c1~c5, c1~c6)
G3 <- graph_Laplacian(g3)
G3

g3o <- inla.qreordering(G3)
G3.reord <- G3[g3o$reo, g3o$reo]
G3.reord

g4 <- list(c1~c3, c2~c3, c3~c5, c4~c5, c5~c6)
G4 <- graph_Laplacian(g4)
G4

g4o <- inla.qreordering(G4)
G4.reord <- G4[g4o$reo, g4o$reo]
G4.reord

par(mfcol = c(2, 3), mar = c(1, 1, 1, 1))
plot(inla.read.graph(G1))
plot(inla.read.graph(G2))
plot(inla.read.graph(G3))
plot(inla.read.graph(G3.reord))
plot(inla.read.graph(G4))
plot(inla.read.graph(G4.reord))

thetas1 <- rep(c(log(2), -1), c(6, 5))
thetas1b <- rep(c(log(2), 1), c(6, 5))
thetas1c <- c(rep(log(2), 6), -1,1,-1,1,-1)
thetas2 <- rep(c(log(2), -1), c(6, 6))
thetas2b <- rep(c(log(2), 1), c(6, 6))
thetas2c <- c(rep(log(2), 6), -1,1,-1,1,-1,1)

c64 <- rgb(seq(0,1,length=64), 0.3,
           seq(1,0,length=64))

par(mfrow = c(4, 4), mar = c(.5,.5,.5,.5))
image(as.matrix(cov2cor(chol2inv(t(dag_L(G1, thetas1))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G2, thetas2))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G3, thetas1))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G4, thetas1))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G1, thetas1b))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G2, thetas2b))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G3, thetas1b))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G4, thetas1b))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G1, thetas1c))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G2, thetas2c))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G3, thetas1c))))), col=c64, axes = FALSE)
image(as.matrix(cov2cor(chol2inv(t(dag_L(G4, thetas1c))))), col=c64, axes = FALSE)
plot(inla.read.graph(G1))
plot(inla.read.graph(G2))
plot(inla.read.graph(G3))
plot(inla.read.graph(G4))

detach("package:corGraphs", unload = TRUE)
library("corGraphs")
