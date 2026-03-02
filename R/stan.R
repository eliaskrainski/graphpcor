#' Build/add STAN code/data for the correlation's PC-prior.
#' @param x either a STAN code or list with the data used to fit
#' a STAN model.
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
#' @export
stan_add <- function(x, base, lambda, name) {
  if(missing(x)) stop("Please provide 'x'!")
  if(missing(base)) stop("Please provide 'base'!")
  if(missing(lambda)) stop("Please provide 'lambda'!")
  if(missing(name)) stop("Please provide 'name'!")
  if(inherits(base, "basecor") || (base == 'pc_correl')) {
    return(stan_add_pc_correl(x, base, lambda, name))
  }
  if(inherits(base, "basepcor") || (base == 'graphpcor')) {
    return(stan_add_graphpcor(x, base, lambda, name))
  }
}

stan_add_pc_correl <- function(x, base, lambda, name) {
  if(length(lambda)>1) {
    warning('length(lambda)>1, using lambda[1]!')
  }
  lambda <- as.numeric(lambda[1])
  stopifnot(lambda>0)

  if(inherits(x, "list")) {
  ## build the additional data
    aD <- list(
      Lcorrel_dim = ncol(base$base),
      Lcorrel_lambda = lambda)
    I0 <- hessian(base)
    I0dec <- dspd(I0)
    aD$Lcorrel_theta_dim <- as.integer(ncol(I0))
    aD$Lcorrel_logDetIhalf <- abs(I0dec$logDeterminant) * 0.5
    aD$Lcorrel_theta_base <- base$theta
    aD$Lcorrel_Ibase_half <- I0dec$sqrt
    names(aD) <- gsub("Lcorrel", name, names(aD), fixed = TRUE)
    return(c(x, aD))
  }

  if(inherits(x, "character")) {

##    stopifnot(grep(name, x)) ## check

    ## data input definitions
    aC <- list(data = "
    int Lcorrel_dim;
    int Lcorrel_theta_dim;
    real<lower=0> Lcorrel_lambda;
    real<lower=0> Lcorrel_logDetIhalf;
    vector[Lcorrel_theta_dim] Lcorrel_theta_base;
    matrix[Lcorrel_theta_dim, Lcorrel_theta_dim] Lcorrel_Ibase_half;
    ")

    ## parameters input
    aC$parameters <- "
    vector[Lcorrel_theta_dim] Lcorrel_theta;
    "

    ## transformed parameters
    aC$"transformed parameters" <- "
    vector[Lcorrel_theta_dim] Lcorrel_xi = Lcorrel_Ibase_half * (Lcorrel_theta - Lcorrel_theta_base);
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
    aC$model <- "
    Lcorrel_radius ~ exponential(Lcorrel_lambda);
    target += lgamma(Lcorrel_theta_dim * 0.5) -Lcorrel_theta_dim*0.5 * log(pi());
    target += Lcorrel_logDetIhalf - (Lcorrel_theta_dim-1.0)*log(Lcorrel_radius);
  "
    aC$"generated quantities" = ""

    for(i in seq_along(aC)) {
      aC[[i]] <- gsub("Lcorrel", name, aC[[i]], fixed = TRUE)
    }

    csnams <- c("data", "parameters", "transformed parameters",
                "model", "generated quantities")
    stopifnot(all(names(aC)==csnams)) ## check

    sx <- strsplit(x, "}")[[1]]
    i <- lapply(csnams, grep, sx)
    ad0 <- 0
    for(k in seq_along(i)) {
      if(length(i[[k]])==0) {
        aC[[k]] <- paste0("\n", csnams[[k]], "{\n", aC[[k]])
      } else {
        aC[[k]] <- paste(sx[i[[k]]], aC[[k]])
      }
    }
    return(paste0(paste(unlist(aC), collapse="}\n"), "\n}\n"))
  }

}
