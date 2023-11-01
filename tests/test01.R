
library(corGraphs)
library(INLA)

mcor <- matrix(c(1, 0.9, 0.9, 1), 2)
sds <- c(3, 0.5)
mcov <- t(mcor * sds) * sds
mcov

n <- 300
m <- nrow(mcov)

ll <- chol(mcov)
xx <- matrix(rnorm(n * m), n) %*% ll

cov(xx)
cor(xx)

dataf <- data.frame(w = runif(n))
dataf$y <- xx[,1] + xx[,2] * dataf$w

dataf$idx1 <- rep(1L, n)
dataf$idx2 <- rep(2L, n)
dataf$repl <- 1:n

ff <- y ~ 0 +
    f(idx1, model = rGmodel, replicate = repl) +
    f(idx2, w, copy = "idx1", replicate = repl)

S <- list(p1 ~ c1 + c2)

SP_plot <- GraphPlot(S, base=0)

par(mar = c(1, 1, 1, 1))
plot(SP_plot$gr, nodeAttrs = SP_plot$nAttrs)

SP <- GraphDens(S)
names(SP)

ArgsList <- list(
    S = S,
    lambda = 7,
    SP = SP,
    Tdist = Tdist,
    GraphPrior = GraphPrior,
    init = 0)

str(ArgsList)

rGmodel <- inla.rgeneric.define(
    CorGraphs.rmodel,
    args = ArgsList)

inla.setOption(safe = FALSE)

hfix <- list(prec = list(initial = 10, fixed = TRUE))

fit <- inla(
    formula = ff,
    data = dataf, 
    control.family = list(list(hyper = hfix)),
    verbose = !TRUE)

summary(fit)

ESTcor <- ThetaCor(
    SP, fit$summary.hyperpar$mean,
    COV = FALSE)
colnames(ESTcor) <- rownames(ESTcor) <-
    SP$STR[[1]][(SP$NP+1):(SP$NP+SP$NC)]
print(ESTcor)

hessian(SP$JD[[1]], rep(0, length(SP$STR[[1]])), SDev = exp(-0.483))
solve(hessian(SP$JD[[1]], rep(0, length(SP$STR[[1]])), SDev = exp(-0.483)))
cov2cor(solve(hessian(SP$JD[[1]], rep(0, length(SP$STR[[1]])), SDev = exp(-0.483))))
