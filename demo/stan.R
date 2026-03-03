
## STAN model code without the prior for L, as LCorr
Scode0 <- "
data {
  int<lower=1> p;
  int<lower=1> n;
  vector[p] y[n];
  vector[p] mu;
}
model {
  y ~ multi_normal_cholesky(mu, LCorr);
}
generated quantities {
  corr_matrix[p] Corr = tcrossprod(LCorr);
}
"

## add the pc_correl code
Scode <- stan_add(Scode0, 'pc_correl', lambda = 1, name = "LCorr")
Scode

Sdata0 <- list(
  p = as.integer(3), 
  n = as.integer(1), 
  y = matrix(0,1,3), 
  mu = rep(0.0, 3)
)

baseC <- basecor(rep(0,3), 3)

Sdata <- stan_add(Sdata0, baseC, lambda = 1, name = "LCorr")
str(Sdata)

