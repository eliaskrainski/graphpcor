## for details see
## https://link.springer.com/article/10.1007/s10260-025-00788-y

if(FALSE) {
    
### examples of what is not allowed
    treepcor(p1 ~ p2)
    treepcor(p1 ~ c2)
    
    treepcor(
        p1 ~ c1 + c2,
        p2 ~ c3)
    
    treepcor(
        p1 ~ c1 + c2,
        p2 ~ p1 + c2 + c3)
    
    treepcor(
        p1 ~ c1 + c2,
        p2 ~ p3 + c2 + c3)
    
    treepcor(
        p1 ~ p2 + c1 + c2,
        p2 ~ c2 + c3)

}

### allowed cases

## 3 children and 1 parent
tree1 <- treepcor(p1 ~ c1 + c2 - c3)

tree1

dim(tree1)

summary(tree1)

plot(tree1)

tpcQ(tree1)

(q1 <- tpcQ(tree1, theta = c(0)))

v1 <- chol2inv(chol(q1))

v1

cov2cor(v1)

vcov(tree1, raw = TRUE)
cov2cor(vcov(tree1, raw = TRUE))
vcov(tree1)

vcov(tree1, theta = 0)
vcov(tree1, theta = -1)
vcov(tree1, theta = 1)

cov2cor(vcov(tree1))
cov2cor(vcov(tree1, theta = -1))
cov2cor(vcov(tree1, theta = 1))

## 4 children and 2 parent
tree2 <- treepcor(
    p1 ~ p2 + c1 + c2,
    p2 ~ -c3 + c4)
tree2
dim(tree2)
summary(tree2)

tpcQ(tree2)
tpcQ(tree2, theta = c(0, 0))
tpcQ(tree2, theta = c(-1, 1))

vcov(tree2, raw = TRUE)
cov2cor(vcov(tree2, raw = TRUE))

cov2cor(solve(tpcQ(tree2, theta = c(0,0))))
vcov(tree2, theta = c(0,0))
vcov(tree2, theta = c(log(4:1), 0,0))

vcov(tree2)

tree2

## 4 children and 2 parent (notice the signs)
tree2b <- treepcor(
    p1 ~ -p2 + c1 + c2,
    p2 ~ -c3 + c4)
tree2b
dim(tree2b)
summary(tree2b)

summary(tree2)
summary(tree2b)

par(mfrow = c(1, 2), mar = c(0,0,0,0))
plot(tree2)
plot(tree2b)

tpcQ(tree2)
tpcQ(tree2b)

## prec is not equal
all.equal(tpcQ(tree2, theta = c(0, 0)),
          tpcQ(tree2b, theta = c(0, 0)))

## vcov is equal
all.equal(vcov(tree2, theta = c(0, 0)),
          vcov(tree2b, theta = c(0, 0)))


