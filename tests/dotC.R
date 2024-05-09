library(INLA)
library(corGraphs)

n <- 5

g <- matrix(c(3,-1,-1,0,0,
              -1,3,0,0,0,
              -1,0,3,-1,0,
              0,0,-1,3,-1,
              0,0,0,-1,3),
            n)
all(g==t(g))

l <- t(chol(g))
l
tcrossprod(l)

iq <- which((g!=0) & lower.tri(g))
iq
nlq <- length(iq)
nlq

ifi <- which((l!=0) & (g==0))
ifi

th <- c(rep(log(5), n), c(-1, .3, -2, -0.5))
th

l0 <- matrix(0, n, n)
diag(l0) <- th[1:n]
l0[iq] <- th[n + 1:nlq]
l0

l1 <- .C("fillL", as.integer(n), as.integer(length(ifi)),
         row(l0)[ifi]-1L, col(l0)[ifi]-1L, l=l0,
   PACKAGE = "corGraphs")$l
l1

g
tcrossprod(l1)

is.matrix(g)
corGraphs:::matrix2graph(g)
corGraphs:::graph_check(corGraphs:::matrix2graph(g))
corGraphs:::graph_elements(
                corGraphs:::matrix2graph(g))


dm <- dag_model(g, rep(1,n), (1:n-0.5)/n, 1, debug = !TRUE)

str(dm$f$cgeneric$data)

cgeneric_get(dm, cmd = "graph")
cgeneric_get(dm, cmd = "Q")
cgeneric_get(dm, cmd = "initial")
cgeneric_get(dm, cmd = "mu")
cgeneric_get(dm, cmd = "log_prior")

str(
    cgeneric_get(
        dm,
        cmd = c("graph", "Q", "initial", "mu", "log_prior")
    )
)


str(
    cgeneric_get(
        dm,
        cmd = c("graph", "Q", "initial", "mu", "log_prior"),
        theta = c(
            rep(log(5), n),
            rep(0.5, dm$f$cgeneric$data$ints$ne)
        )
    )
)


detach("package:corGraphs", unload = TRUE)
library(corGraphs)
