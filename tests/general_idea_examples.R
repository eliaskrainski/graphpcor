library(INLA)

g1 <- inla.as.sparse(
    sparseMatrix(
        i = 1:5,
        j = 2:6,
        x = 1L,
        dims = c(6, 6)
    )
)

g2 <- inla.as.sparse(
    sparseMatrix(
        i = c(1:5, 1),
        j = c(2:6, 6),
        x = 1L,
        dims = c(6, 6)
    )
)

g3 <- inla.as.sparse(
    sparseMatrix(
        i = rep(1, 5),
        j = 2:6,
        x = 1L,
        dims = c(6, 6)
    )
)

g4 <- inla.as.sparse(
    sparseMatrix(
        i = c(1, 2, 3, 4, 5),
        j = c(3, 3, 5, 5, 6),
        x = 1L,
        dims = c(6, 6)
    )
)

par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
plot(inla.read.graph(g1))
plot(inla.read.graph(g2))
plot(inla.read.graph(g3))
plot(inla.read.graph(g4))

library(corGraphs)
library(ggpubr)

ggarrange(
    image(cov2cor(chol2inv(chol(dag_L(g1))))),
    image(cov2cor(chol2inv(chol(dag_L(g2))))),
    image(cov2cor(chol2inv(chol(dag_L(g3))))),
    image(cov2cor(chol2inv(chol(dag_L(g4)))))
)


