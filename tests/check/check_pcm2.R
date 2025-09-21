library(graphpcor)
library(plot3D)

lseq <- c(0.3, 1, 5)
llab <- lapply(lseq, function(l) 
  bquote(lambda == .(l)))
nl <- length(lseq)
lcol <- c(1,1,1)#2, 4)
b <- 5; a <- -b
h  <- b/100
x0 <- c(seq(a+h/2, 0, h), seq(h/2, b-h/2, h))
nx0 <- length(x0)

dfn <- function(lambda) {
  outer(x0, x0, function(a,b) {
      r <- sqrt(a^2+b^2)
      h <- 0.001
      il <- which(r<h)
      r[il] <- h/2 + (r[il] -h/2)/(h-h/2)
      lambda*exp(-lambda*r)/(2*pi*r)})
}

d2 <- lapply(lseq, function(l) dfn(l))

sapply(d2, function(x)
    c(sum(colSums(x*h)*h),
      sum(rowSums(x*h)*h), sum(x*h*h)))

par(mfcol = c(1, 2), mar = c(4,4,1,1), mgp = c(2, 1, 0), cex.lab = 1)
plot(x0, colSums(d2[[nl]])*h, type = 'n')
for(i in 1:nl) {
    lines(x0, colSums(d2[[i]])*h, type = 'l', lty = i)
}
plot(x0, rowSums(d2[[nl]])*h, type = 'n')
for(i in 1:nl) {
    lines(x0, rowSums(d2[[i]])*h, type = 'l', lty = i)
}

par(mfcol = c(1, 3), mar = c(0,0,0,3), cex.lab = 1)
for(i in 1:nl) {
#  bkz <- exp(seq(log(min(d2[[i]]))-1e-9, 
 #                log(max(d2[[i]]))+1e-9, length = 65))
  surf3D(matrix(x0, nx0, nx0), 
         matrix(rep(x0, each = nx0), nx0), 
         matrix(d2[[i]], nx0), 
  #       breaks = bkz, colkey = !TRUE,
         xlab = "", 
         ylab  = "", 
         zlab = 'Density',
         box = TRUE, bty = "b", phi = 20, theta = 30)
  text(-0.15, -0.43, expression(xi[1]), adj = 1)
  text(0.28, -0.37, expression(xi[2]), adj = 1)
  legend("topleft", bty = "n", as.expression(llab[[i]]), cex = 2)
}

if(FALSE) {
    
    library(rgl)
    
    rglwidget()
    
    plot3d(rep(x0, nx0),
           rep(x0, each = nx0),
           d2[[1]],
           "X", "Y", "Z")

}


