
## test functionalities
library(graphpcor)

g4 <- graphpcor(x1~x2+x3, x4~x2+x3)

dd <- dim(g4)
n <- dd[1]
m <- dd[2]

ths <- rnorm(n)
thc <- rnorm(m)

basepcor(thc, iLtheta = g4)
basepcor(thc, iLtheta = g4, d0 = 10 + m:1)

basepcor(thc[2:m], iLtheta = g4, iparams = c(1,1:(m-1)))
basepcor(thc[3:m], iLtheta = g4, iparams = c(1,1:(m-2),1))

hessian(basepcor(thc[3:m], iLtheta = g4, iparams = c(1,1:(m-2),1)))

solve(cgeneric_Q(cgeneric(g4, lambda = 1, useINLAprecomp = FALSE),
           theta = thc))

cgeneric_Q(cgeneric(g4, lambda = 1, useINLAprecomp = FALSE,
              sigma.reference.probability = rep(0.1, n)),
              theta = c(ths, thc))

cgeneric_Q(cgeneric(g4, lambda = 1, useINLAprecomp = FALSE),
     theta = rep(thc[1], 4))

cgeneric(g4, lambda = 1, iparam = c(1:4, 5,5,5,5), debug = TRUE)
cgeneric(g4, lambda = 1, iparam = c(1:4, 5,5,5,5), debug = TRUE)$f$cgeneric$data

str(cgeneric(g4, lambda = 1, iparam = c(1:4, 5,5,5,5))$f$cgeneric$data$ints)

cgeneric_initial(cgeneric(g4, lambda = 1, iparam = c(1:4, 5,5,5,5),
                 useINLAprecomp = FALSE, debug = TRUE))
cgeneric_initial(cgeneric(g4, lambda = 1, iparam = c(1,1,2,2,3,3,4,4),
                 useINLAprecomp = FALSE, debug = TRUE))

thc

str(cgeneric(g4, lambda = 1, sigma.prior.p = runif(n))$f$cgeneric$data)

str(cgeneric(g4, lambda = 1, sigma.prior.probability = c(0.5),
             iparam = c(1,1,1,1,2,3,4,5))$f$cgeneric$data$int)

cgeneric_initial(cgeneric(g4, lambda = 1, sigma.prior.p = c(0.5),
                 useINLAprecomp = FALSE,
                 iparam = c(1,1,1,1,2,3,4,5)))

all.equal(
    cgeneric_Q(cgeneric(g4, lambda = 1, sigma.prior.prob = rep(0.5,n), 
                  useINLAprecomp = FALSE),
         theta = c(rep(ths[1], n), thc)),
    cgeneric_Q(cgeneric(g4, lambda = 1, sigma.prior.p = c(0.5),
                  useINLAprecomp = FALSE,
                  iparam = c(1,1,1,1,2,3,4,5)),
         theta = c(ths[1], thc))
)

all.equal(
    cgeneric_Q(cgeneric(g4, lambda = 1, sigma.prior.prob = rep(0.5,n), 
                  useINLAprecomp = FALSE),
         theta = c(ths[c(1,2,1,1)], thc)),
    cgeneric_Q(cgeneric(g4, lambda = 1, sigma.prior.p = c(0.5,0.5),
                  useINLAprecomp = FALSE,
                  iparam = c(1,2,1,1,3,4,5,6)),
         theta = c(ths[1:2], thc))
)

all.equal(
    cgeneric_Q(cgeneric(g4, lambda = 1), theta = rep(thc[1], m)), 
    cgeneric_Q(cgeneric(g4, lambda = 1, iparam = c(1:4, 5,5,5,5)),
         theta = thc[1])
)

all.equal(
    cgeneric_Q(cgeneric(g4, lambda = 1,
                  useINLAprecomp = FALSE),
         theta = thc[c(1,1,2,2)]),
    cgeneric_Q(cgeneric(g4, lambda = 1, iparam = c(1:4, 5,5,6,6),
                  useINLAprecomp = FALSE),
         theta = thc[1:2])
)

all.equal(
    cgeneric_Q(cgeneric(g4, lambda = 1,
                  useINLAprecomp = FALSE),
         theta = thc[c(1,2,2,1)]), 
    cgeneric_Q(cgeneric(g4, lambda = 1, iparam = c(1:4, 5,6,6,5),
                  useINLAprecomp = FALSE),
         theta = thc[1:2])
)
