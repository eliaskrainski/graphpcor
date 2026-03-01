#' Build STAN code for the PC-prior for a correlation matrix.
#' @param base `basecor` or `basepcor` object to define the
#' base correlation model. See [basecor()] or [basepcor()].
#' @param lambda the parameter for the exponential prior on
#' the radius of the sphere, see details in the
#' PC-multivariate vignette.
#' @param name character to provide the name for the
#' Cholesky of a correlation matrix or the correlation matrix.
#' See details.
#' @details
#' The parametrization is set as in [basecor()] or [basepcor()].
#' If a `basecor` is provided, the prior would be considered for
#' the Cholesky factor of a correlation matrix.
#' If a `basepcor` is provided, the prior would be considered for
#' a correlation matrix (parametrized from it's inverse).
#' The base is used to define an informative prior, as derived in
#' the pcmultivariate vignette.
#' @references
#' Daniel Simpson, H\\aa vard Rue, Andrea Riebler, Thiago G.
#' Martins and Sigrunn H. S\\o rbye (2017).
#' Penalising Model Component Complexity:
#' A Principled, Practical Approach to Constructing Priors.
#' Statistical Science 2017, Vol. 32, No. 1, 1–28.
#' <doi 10.1214/16-STS576>
#' @return a list of two elements, one as a list of three
#' additional code to be added into a STAN code and
#' the other with the required additional data.
stan <- function(base, lambda, name) {
  if(missing(base)) stop("Please provide 'base'!")
  if(missing(lambda)) stop("Please provide 'lambda'!")
  if(missing(name)) stop("Please provide 'name'!")
  if(inherits(base, "basecor")) {
    return(stan_pc_correl(base, lambda, name))
  }
  if(inherits(base, "basepcor")) {
    return(stan_pc_graphpcor(base, lambda, name))
  }
}

stan_pc_correl <- function(base, lambda, name) {
  if(length(lambda)>1) {
    warning('length(lambda)>1, using lambda[1]!')
  }
  lambda <- as.numeric(lambda[1])
  stopifnot(lambda>0)

  ## build the additional data
  aSdata <- list(
    Lcorrel_lambda = lambda
  )
  I0 <- hessian(base)
  I0dec <- dspd(I0)
  aSdata$Lcorrel_theta_dim <- ncol(I0)
  aSdata$Lcorrel_logDetIhalf <- abs(I0dec$logDeterminant) * 0.5
  aSdata$Lcorrel_theta_base <- base$theta
  aSdata$Lcorrel_Ibase_half <- I0dec$sqrt
  names(aSdata) <- gsub("Lcorrel", name,
                        names(aSdata), fixed = TRUE)

  ## data input definitions
  aScode <- list(
  data = "
  int Lcorrel_dim;
  int Lcorrel_theta_dim;
  real<lower=0> Lcorrel_lambda;
  real<lower=0> Lcorrel_logDetIhalf;
  vector<Lcorrel_theta_dim> Lcorrel_theta_base;
  matrix<Lcorrel_theta_dim, Lcorrel_theta_dim> Lcorrel_Ibase_half;
  ")

  ## parameters imput
  aScode$parameters <- "
  vector[Lcorrel_theta_dim] Lcorrel_theta;
  "

  ## transformed parameters
  aScode$"transformed parameters" <- "
  vector[m] Lcorrel_xi = Lcorrel_Ibase_half * (Lcorrel_theta - Lcorrel_theta_base);
  real<lower=0> Lcorrel_radius = dot_self(Lcorrel_xi);
  matrix[Lcorrel_dim,Lcorrel_dim] Lcorrel_A;
  matrix[Lcorrel_dim,Lcorrel_dim] Lcorrel_B;
  cholesky_factor_corr[Lcorrel_dim] Lcorrel;
  Lcorrel[1,1] = 1.0;
  for(i in 1:Lcorrel_dim) {
    Lcorrel_A[i,i] = 1.0;
    Lcorrel_B[i,i] = 1.0;
  }
  for(i in 1:(Lcorrel_dim-1)) {
    for(j in (i+1):Lcorrel_dim) {
        Lcorrel_A[i,j] = 0.0;
        Lcorrel_B[i,j] = 0.0;
        Lcorrel[i,j] = 0.0;
    }
  }
  {
  int k = 0;
  for(j in 1:(Lcorrel_dim-1)) {
    for(i in (j+1):Lcorrel_dim) {
      k = k+1;
      Lcorrel_A[i,j] = tanh(Lcorrel_theta[k]);
      Lcorrel_B[i,j] = sqrt(1-Lcorrel_A[i,j]^2);
    }
  }
  }
  for(j in 2:(Lcorrel_dim-1)) {
    for(i in j:Lcorrel_dim) {
      Lcorrel_B[i,j] *= Lcorrel_B[i,j-1];
    }
  }
  for(i in 1:Lcorrel_dim) {
     Lcorrel[i, 1] = Lcorrel_A[i,1];
  }
  for(j in 2:Lcorrel_dim) {
    for(i in j:Lcorrel_dim) {
      Lcorrel[i,j] = Lcorrel_A[i,j] * Lcorrel_B[i,j-1];
    }
  }
"
  ## model part
  aSdata$model <- "
  Lcorrel_radius ~ exponential(Lcorrel_lambda);
  target += lgamma(Lcorrel_theta_dim * 0.5) -Lcorrel_theta_dim*0.5 * log(pi());
  target += Lcorrel_logDetIhalf - (Lcorrel_theta_dim-1.0)*log(Lcorrel_radius);
  "
  for(i in seq_along(length(aScode))) {
    aScode[[i]] <- gsub("Lcorrel", name, fixed = TRUE)
  }
  return(aScode)

  return(list(data = aSdata, code = aScode))
}
