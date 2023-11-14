## check consistency by defining the same model
## with different label, ordering

library(corGraphs)

g1 <- list(p1 ~ c1,
           p2 ~ p1 + c2 + c3,
           p3 ~ p2 + c4)
g2 <- list(p1 ~ p2 + c4,
           p2 ~ p3 + c2 + c3,
           p3 ~ c1)
g3 <- list(p1 ~ c4,
           p2 ~ p1 + c2 + c3,
           p3 ~ p2 + c1)
g4 <- list(p1 ~ p2 + c1,
           p2 ~ p3 + c2 + c3,
           p3 ~ c4)

dag_elements(g1)

g1.plot <- GraphPlot(g1)
g2.plot <- GraphPlot(g2)
g3.plot <- GraphPlot(g3)
g4.plot <- GraphPlot(g4)

par(mfrow=c(2,2), mar = c(0, 0, 0, 0))
plot(g1.plot$gr, nodeAttrs = g1.plot$nAttrs)
plot(g2.plot$gr, nodeAttrs = g2.plot$nAttrs)
plot(g3.plot$gr, nodeAttrs = g3.plot$nAttrs)
plot(g4.plot$gr, nodeAttrs = g4.plot$nAttrs)

q1 <- dag_precision(g1, c(0,0,0))
q2 <- dag_precision(g2, c(0,0,0))
q3 <- dag_precision(g3, c(0,0,0))
q4 <- dag_precision(g4, c(0,0,0))

nc <- 4
cov2cor(solve(q1)[1:nc, 1:nc])
cov2cor(solve(q2)[1:nc, 1:nc])
cov2cor(solve(q3)[1:nc, 1:nc])
cov2cor(solve(q4)[1:nc, 1:nc])

d1 <- GraphDens(g1)
d2 <- GraphDens(g2)
d3 <- GraphDens(g3)
d4 <- GraphDens(g4)

ThetaCor(d1, rep(0, 7))
ThetaCor(d2, rep(0, 7))
ThetaCor(d3, rep(0, 7))
ThetaCor(d4, rep(0, 7))

p1 <- GraphPlotPrior(g1)
p2 <- GraphPlotPrior(g2)
p3 <- GraphPlotPrior(g3)
p4 <- GraphPlotPrior(g4)

par(mfrow = c(3,3), mar = c(0,0,0,0))
for(k in 1:length(p2))
   plot(p2[[k]]$gr)
plot(1, type = "n", axes = FALSE, xlab = "", ylab="")
for(k in 1:length(p4))
   plot(p4[[k]]$gr)

detach("package:corGraphs", unload=TRUE)
library(corGraphs)
