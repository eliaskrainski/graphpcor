library(INLA)
library(corGraphs)

g1 <- list(c1~c2, c2~c3, c3~c4, c4~c5, c5~c6)
G1 <- graph_Laplacian(g1)

g2 <- c(g1, c1~c6)
G2 <- graph_Laplacian(g2)

g3 <- list(c1~c2, c1~c3, c1~c4, c1~c5, c1~c6)
G3 <- graph_Laplacian(g3)

g3o <- inla.qreordering(G3)
G3.reord <- G3[g3o$reo, g3o$reo]

g4 <- list(c1~c3, c2~c3, c3~c5, c4~c5, c5~c6)
G4 <- graph_Laplacian(g4)
g4o <- inla.qreordering(G4)
G4.reord <- G4[g4o$reo, g4o$reo]

par(mfrow = c(2, 3), mar = c(1, 1, 1, 1))
plot(inla.read.graph(G1))
plot(inla.read.graph(G2))
plot(inla.read.graph(G3))
plot(inla.read.graph(G3.reord))
plot(inla.read.graph(G4))
plot(inla.read.graph(G4.reord))

thetas11 <- rep(c(log(2), -1), c(6, 5))
thetas12 <- rep(c(log(2), -1), c(6, 6))

library(ggpubr)

ggarrange(
    image(cov2cor(chol2inv(t(dag_L(G1, thetas11))))),
    image(cov2cor(chol2inv(t(dag_L(G2, thetas12))))),
    image(cov2cor(chol2inv(t(dag_L(G3, thetas11))))),
    image(cov2cor(chol2inv(t(dag_L(G4, thetas11)))))
)

detach("package:corGraphs", unload = TRUE)
library("corGraphs")
