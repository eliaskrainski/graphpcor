library(INLA)
library(graphpcor)

g1 <- graph(c1~c2, c2~c3, c3~c4, c4~c5, c5~c6)
G1 <- Laplacian(g1)
G1

g2 <- graph(c1~c2+c6, c2~c3, c3~c4, c4~c5, c5~c6)
G2 <- Laplacian(g2)
G2

g3 <- graph(c1~c2+c3+c4+c5+c6)
G3 <- Laplacian(g3)
G3

g3o <- inla.qreordering(G3)
G3.reord <- G3[g3o$reo, g3o$reo]
G3.reord

g4 <- graph(c1~c3, c2~c3, c3~c5, c4~c5, c5~c6)
G4 <- Laplacian(g4)
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
image(cov2cor(variance(g1, theta = thetas1)), col=c64, axes = FALSE)
image(cov2cor(variance(g2, theta = thetas2)), col=c64, axes = FALSE)
image(cov2cor(variance(g3, theta = thetas1)), col=c64, axes = FALSE)
image(cov2cor(variance(g4, theta = thetas1)), col=c64, axes = FALSE)
image(cov2cor(variance(g1, theta = thetas1b)), col=c64, axes = FALSE)
image(cov2cor(variance(g2, theta = thetas2b)), col=c64, axes = FALSE)
image(cov2cor(variance(g3, theta = thetas1b)), col=c64, axes = FALSE)
image(cov2cor(variance(g4, theta = thetas1b)), col=c64, axes = FALSE)
image(cov2cor(variance(g1, theta = thetas1c)), col=c64, axes = FALSE)
image(cov2cor(variance(g2, theta = thetas2c)), col=c64, axes = FALSE)
image(cov2cor(variance(g3, theta = thetas1c)), col=c64, axes = FALSE)
image(cov2cor(variance(g4, theta = thetas1c)), col=c64, axes = FALSE)
plot(inla.read.graph(G1))
plot(inla.read.graph(G2))
plot(inla.read.graph(G3))
plot(inla.read.graph(G4))

detach("package:graphpcor", unload = TRUE)
library("graphpcor")
