
library(INLA)
library(graphpcor)

pth <- file.path(.libPaths(), "INLA")[1]

system(
    paste0('gcc -Wall -fpic -g -O -c -o cgeneric-demo.o ',
           pth, '/cgeneric/cgeneric-demo.c')
)
system('gcc -shared -o cgeneric-demo.so cgeneric-demo.o')

model <- inla.cgeneric.define(
    model = 'inla_cgeneric_ar1_model',
    shlib = 'cgeneric-demo.so',
    n = 7L,
    debug = FALSE
)

model$f$cgeneric$data$characters <-
    c(model$f$cgeneric$data$characters,
      list("test1", 'test2', 'last')
      )
model$f$cgeneric$data$characters
class(model) <- c("cgeneric", "inla.cgeneric")

(ini <- initial(model))
prior(model, theta = ini)

graph(model, optimize = TRUE)
prec(model, theta = ini, optimize = TRUE)

graph(model)
prec(model, theta = ini)

