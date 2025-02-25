## see https://en.wikipedia.org/wiki/N-sphere

library(numDeriv)

## check with m = 4
to.x <- function(z) {
    n <- length(z)
    stopifnot(n == 4)
    r <- z[1]
    phi <- z[-1]
    x <- numeric(n)
    x[1] <- r * cos(phi[1])
    x[2] <- r * sin(phi[1]) * cos(phi[2])
    x[3] <- r * sin(phi[1]) * sin(phi[2]) * cos(phi[3])
    x[4] <- r * sin(phi[1]) * sin(phi[2]) * sin(phi[3])
    return (x)
}

to.rphi <- function(x) {
    n <- length(x)
    stopifnot(n == 4)
    r <- sqrt(sum(x^2))
    phi <- numeric(n-1)
    phi[1] <- atan2(sqrt(sum(x[4:2]^2)), x[1])
    phi[2] <- atan2(sqrt(sum(x[4:3]^2)), x[2])
    phi[3] <- atan2(sqrt(sum(x[4:4]^2)), x[3])
    return (c(r, phi))
}

## everything assumes that n=4
n <- 4
x <- rnorm(n)
stopifnot(n == 4)

print(mean(abs(to.x(to.rphi(x)) - x)))

jacobian(to.rphi, x)
xx <- to.rphi(x)

## this is a test that forward-Jacobian is 1/backward-Jacobian
print(abs(det(jacobian(to.rphi, x))))
print(1/abs(det(jacobian(to.x, xx))))

## this is using formula for the 'a closed-form expression for the volume element in spherical
## coordinates' in https://en.wikipedia.org/wiki/N-sphere

print(abs(det(jacobian(to.x, xx))))
r <- xx[1]
phi <- xx[-1]
print(r^(n-1) * sin(phi[1])^(n-2) * sin(phi[2])^(n-3))
