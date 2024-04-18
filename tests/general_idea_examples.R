library(INLA)

g1 <- inla.as.sparse(
    sparseMatrix(
        i = 1:5,
        j = 2:6,
        x = 1L,
        dims = c(6, 6)
    )
); g1 <- g1 + t(g1)

g2 <- inla.as.sparse(
    sparseMatrix(
        i = c(1:5, 1),
        j = c(2:6, 6),
        x = 1L,
        dims = c(6, 6)
    )
); g2 <- g2 + t(g2)

g3 <- inla.as.sparse(
    sparseMatrix(
        i = rep(1, 5),
        j = 2:6,
        x = 1L,
        dims = c(6, 6)
    )
); g3 <- g3 + t(g3)
g3o <- inla.qreordering(g3)
g3.reord <- g3[g3o$reo, g3o$reo]

g4 <- inla.as.sparse(
    sparseMatrix(
        i = c(1, 2, 3, 4, 5),
        j = c(3, 3, 5, 5, 6),
        x = 1L,
        dims = c(6, 6)
    )
); g4 <- g4 + t(g4)
g4o <- inla.qreordering(g4)
g4.reord <- g4[g4o$reo, g4o$reo]

par(mfrow = c(2, 3), mar = c(1, 1, 1, 1))
plot(inla.read.graph(g1))
plot(inla.read.graph(g2))
plot(inla.read.graph(g3))
plot(inla.read.graph(g3.reord))
plot(inla.read.graph(g4))
plot(inla.read.graph(g4.reord))

library(corGraphs)

thetas11 <- rep(1, 6+5)
thetas12 <- rep(1, 6+6)

dag_L(g1, thetas11)

dag_L(g2, thetas12)

dag_L(g3, thetas11)
dag_L(g3.reord, thetas11)

dag_L(g4, thetas11)
dag_L(g4.reord, thetas11)

library(ggpubr)

ggarrange(
    image(cov2cor(chol2inv(t(dag_L(g1, thetas11))))),
    image(cov2cor(chol2inv(t(dag_L(g2, thetas12))))),
    image(cov2cor(chol2inv(t(dag_L(g3, thetas11))))),
    image(cov2cor(chol2inv(t(dag_L(g4, thetas11)))))
)

detach("package:corGraphs", unload = TRUE)
library("corGraphs")
