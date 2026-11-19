
par(mfrow = c(2, 3), mar = c(0,0,0,0))
plot(graphpcor(x~y+v, z~y+v))
plot(graphpcor(x~y,x~v,z~y,z~v))
plot(graphpcor(x~y, v~x, y~z, z~v))
plot(graphpcor(y~x, v~x, z~y, z~v))
plot(graphpcor(y~x+z, v~z+x))
plot(graphpcor(y~x+z, v~x, z~v))

par(mfrow = c(2, 3), mar = c(0,0,0,0))
plot(graphpcor("x~y+v", z~y+v))
plot(graphpcor("x~y",x~v,"z~y","z~v"))
plot(graphpcor(x~y, v~x, y~z, z~v))
plot(graphpcor(paste0("x1~x",2:5)))
plot(graphpcor(list(paste0("x1~x",2:5))))
plot(graphpcor(
    sapply(c(paste0("x1~x",3:5),"x1~x2"), as.character)))
graphpcor(lapply(2:5, function(i) paste0("x1~x",i)))

## the graph in Example 2.6 of the GMRF book
g <- graphpcor(x ~ y + v, z ~ y + v)

class(g)
g

par(mfrow = c(2,2), mar = c(0,0,0,0))
plot(g) ## default 
plot(g, Rgraphviz = TRUE)
plot(g, igraph = TRUE)
plot(g, qgraph = TRUE)

summary(g) ## the graph: nodes and edges (nodes ordered as given)

ne <- dim(g)
ne

## sometimes we need it
G <- Laplacian(g)
G 

## alternatively
all.equal(G,
          Laplacian(graphpcor(x~y, v~x, y~z, z~v)))
all.equal(G,
          Laplacian(graphpcor(x~y, x~v, z~y, v~z)))

plot(graphpcor(G)) ## from a matrix

## base model (theta for lower triangle Cholesky)
theta0l <- rep(-0.5, ne[2])

## vcov() method for graphpcor computes the correlation
##  if only theta for lower of L is provided
C0 <- vcov(g, theta = theta0l)
C0

## the precision for a correlation matrix
Q0 <- gphcQ(g, theta = theta0l)
Q0

all.equal(C0, as.matrix(solve(Q0)))

## the Hessian matrix around a base model
b0 <- basepcor(theta0l, p = ne[1], iLtheta = g)
b0

## a base model can also be a matrix
## however it shall give a precision with
## same sparse pattern as the graph
basepcor(C0, p = ne[1], iLtheta = g)

## the 'iid' case would be
vcov(g, theta = rep(0, ne[2]))
vcov(g, theta = rep(0, sum(ne)))

## marginal variance specified throught standard errors
sigmas <- c(0.3, 0.7, 1.2, 0.5)
## the covariance
vcov(g, theta = c(log(sigmas), rep(0, ne[2]))) ## IID
vcov(g, theta = c(log(sigmas), theta0l))

vcov(g, theta = rep(-3, ne[2])) ## no edge 2~3 but high correlation!!!

## A correlation matrix from Figure 3 in
## https://doi.org/10.1002/sim.6423

nc <- 12
Corr <- diag(0.5,nc,nc);
Corr[upper.tri(Corr)] <- c(
    41,
    36,38,
    27,37,45,
    46,33,35,35,
    46,36,48,46,37,
    69,51,22,20,32,25,
    70,52,48,39,48,49,57,
    20,37,22,20,06,37,28,30,
    61,43,45,59,45,45,51,76,31,
    51,34,45,48,28,48,40,48,18,46,
    39,30,51,42,54,48,34,40,37,49,38
)/100; Corr <- Corr + t(Corr)
caus <- c(
    "Oral", "Stom", "Coln", "Rect", "Livr", "Panc",
    "Larx", "Lung", "Prst", "Blad", "Lymp", "Leuk")
dimnames(Corr) <- list(caus, caus)
round(100*Corr)

baseCorr <- basepcor(
    base = Corr, d0 = (nc:1)^2,
    n = 540, gamma = 1.0, start = "empty")

gg <- graphpcor(baseCorr)
gg

par(mfrow = c(2,2), mar = c(0,0,0,0))
plot(gg) ## default 
plot(gg, Rgraphviz = TRUE)
plot(gg, igraph = TRUE)
plot(gg, qgraph = TRUE)

