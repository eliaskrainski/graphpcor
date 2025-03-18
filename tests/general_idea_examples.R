library(graphpcor)
library(INLA)

g1 <- graphpcor(c1~c2, c2~c3, c3~c4, c4~c5, c5~c6)

G1 <- Laplacian(g1)
G1

g2 <- graphpcor(c1~c2+c6, c2~c3, c3~c4, c4~c5, c5~c6)
G2 <- Laplacian(g2)
G2

g3 <- graphpcor(c1~c2+c3+c4+c5+c6)
G3 <- Laplacian(g3)
G3

g4 <- graphpcor(c1~c3, c2~c3, c3~c5, c4~c5, c5~c6)
G4 <- Laplacian(g4)

thetas1 <- rep(-1, 5)
thetas1b <- rep(1, 5)
thetas1c <- c(-1,1,-1,1,-1)
thetas2 <- rep(-1,6)
thetas2b <- rep(1, 6)
thetas2c <- c(-1,1,-1,1,-1,1)

c64 <- rgb(seq(0, 1, length = 64), 0.3,
           seq(1, 0, length = 64))

par(mfrow = c(4, 4), mar = c(.5,.5,.5,.5))
plot(g1)
plot(g2)
plot(g3)
plot(g4)
image(variance(g1, theta = thetas1), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g2, theta = thetas2), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g3, theta = thetas1), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g4, theta = thetas1), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g1, theta = thetas1b), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g2, theta = thetas2b), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g3, theta = thetas1b), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g4, theta = thetas1b), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g1, theta = thetas1c), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g2, theta = thetas2c), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g3, theta = thetas1c), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))
image(variance(g4, theta = thetas1c), zlim = c(-1, 1), col=c64, axes = FALSE)
text(0:5/5, 0:5/5, paste(1:6))

detach("package:graphpcor", unload = TRUE)
library("graphpcor")
