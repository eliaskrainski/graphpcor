## check consistency by defining the same model
## with different label, ordering

library(corGraphs)

s1 <- list(p1 ~ c1,
           p2 ~ p1 + c2 + c3,
           p3 ~ p2 + c4)
s2 <- list(p1 ~ p2 + c4,
           p2 ~ p3 + c2 + c3,
           p3 ~ c1)
s3 <- list(p1 ~ c4,
           p2 ~ p1 + c2 + c3,
           p3 ~ p2 + c1)
s4 <- list(p1 ~ p2 + c1,
           p2 ~ p3 + c2 + c3,
           p3 ~ c4)

s1.plot <- GraphPlot(s1)
s2.plot <- GraphPlot(s2)
s3.plot <- GraphPlot(s3)
s4.plot <- GraphPlot(s4)

par(mfrow=c(2,2), mar = c(0, 0, 0, 0))
plot(s1.plot$gr, nodeAttrs = s1.plot$nAttrs)
plot(s2.plot$gr, nodeAttrs = s2.plot$nAttrs)
plot(s3.plot$gr, nodeAttrs = s3.plot$nAttrs)
plot(s4.plot$gr, nodeAttrs = s4.plot$nAttrs)

q1 <- QS(s1, c(0,0,0))
q2 <- QS(s2, c(0,0,0))
q3 <- QS(s3, c(0,0,0))
q4 <- QS(s4, c(0,0,0))

nc <- 4
cov2cor(solve(q1)[1:nc, 1:nc])
cov2cor(solve(q2)[1:nc, 1:nc])
cov2cor(solve(q3)[1:nc, 1:nc])
cov2cor(solve(q4)[1:nc, 1:nc])

d1 <- GraphDens(s1)
d2 <- GraphDens(s2)
d3 <- GraphDens(s3)
d4 <- GraphDens(s4)

ThetaCor(d1, rep(0, 7))
ThetaCor(d2, rep(0, 7))
ThetaCor(d3, rep(0, 7))
ThetaCor(d4, rep(0, 7))

p1 <- GraphPlotPrior(s1)
p2 <- GraphPlotPrior(s2)
p3 <- GraphPlotPrior(s3)
p4 <- GraphPlotPrior(s4)

par(mfrow = c(3,3), mar = c(0,0,0,0))
for(k in 1:length(p2))
   plot(p2[[k]]$gr)
plot(1, type = "n", axes = FALSE, xlab = "", ylab="")
for(k in 1:length(p4))
   plot(p4[[k]]$gr)

detach("package:corGraphs", unload=TRUE)
library(corGraphs)
