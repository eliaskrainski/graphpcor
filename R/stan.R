#' Build/add STAN code/data for the correlation's PC-prior.
#' @param x either a STAN code or list with the data used to fit
#' a STAN model.
#' @param model either a character ("pc_correl" or "graphcor") or
#' a `basecor` or a `basepcor` object to define the
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
#' @example demo/stan.R
#' @export
stan_add <- function(x, model, lambda, name) {
  if(missing(x)) stop("Please provide 'x'!")
  if(missing(model)) stop("Please provide 'model'!")
  if(missing(lambda)) stop("Please provide 'lambda'!")
  if(missing(name)) stop("Please provide 'name'!")
  if(is.character(model)) {
    if(model == 'pc_correl') {
      return(stan_add_pc_correl(x, model, lambda, name))
    }
    if(model == 'graphpcor') {
      return(stan_add_graphpcor(x, model, lambda, name))
    }
  } else {
    if(inherits(model, "basecor")) {
      return(stan_add_pc_correl(x, model, lambda, name))
    }
    if(inherits(model, "basepcor")) {
      return(stan_add_graphpcor(x, model, lambda, name))
    }
  }
  warning("Nothing done!")
  return(x)
}
#' @describeIn stan_add method for `basecor`
stan_add_pc_correl <- function(x, model, lambda, name) {
  if(length(lambda)>1) {
    warning('length(lambda)>1, using lambda[1]!')
  }
  lambda <- as.numeric(lambda[1])
  stopifnot(lambda>0)

  if(inherits(x, "list")) {
  ## build the additional data
    aD <- list(
      Lcorrel_dim = ncol(model$base),
      Lcorrel_lambda = lambda)
    I0 <- hessian(model)
    I0dec <- dspd(I0)
    aD$Lcorrel_theta_dim <- as.integer(ncol(I0))
    aD$Lcorrel_logDetIhalf <- abs(I0dec$logDeterminant) * 0.5
    aD$Lcorrel_theta_base <- model$theta
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

    for(i in seq_along(aC)) {
      aC[[i]] <- gsub("Lcorrel", name, aC[[i]], fixed = TRUE)
    }

    return(stan_add_code(x, aC))
  }

  warning("Nothing done!")
  return(x)
}
#' @describeIn stan_add method for `basepcor`
stan_add_graphpcor <- function(x, model, lambda, name) {
  if(length(lambda)>1) {
    warning('length(lambda)>1, using lambda[1]!')
  }
  lambda <- as.numeric(lambda[1])
  stopifnot(lambda>0)

  if(inherits(x, "list")) {
    ## build the additional data
    aD <- list(
      grpc_dim = ncol(model$base),
      grpc_lambda = lambda)
    I0 <- hessian(model)
    I0dec <- dspd(I0)
    aD$grpc_theta_dim <- as.integer(ncol(I0))
    aD$grpc_logDetIhalf <- abs(I0dec$logDeterminant) * 0.5
    aD$grpc_theta_base <- model$theta
    aD$grpc_Ibase_half <- I0dec$sqrt
    aD$grpc_d0 <- model$d0
    L <- matrix(0, model$p, model$p)
    aD$grpc_ii <- as.array(row(L)[model$iLtheta])
    aD$grpc_jj <- as.array(col(L)[model$iLtheta])
    L[model$iLtheta] <- -1
    L <- L + t(L)
    diag(L) <- 1-colSums(L)
    Lf <- t(chol(L))
    ifl <- setdiff(which((abs(Lf)>0) & lower.tri(L)),
                   model$iLtheta)
    aD$grpc_nfi <- length(ifl)
    aD$grpc_ifi <- as.array(row(L)[ifl])
    aD$grpc_jfi <- as.array(col(L)[ifl])
    return(c(x, aD))
  }

  if(inherits(x, "character")) {

    ##    stopifnot(grep(name, x)) ## check

    ## data input definitions
    aC <- list(data = "
    int grpc_dim;
    int grpc_theta_dim;
    real<lower=0> grpc_lambda;
    real<lower=0> grpc_logDetIhalf;
    vector[grpc_theta_dim] grpc_theta_base;
    matrix[grpc_theta_dim, grpc_theta_dim] grpc_Ibase_half;
    vector[grpc_dim] grpc_d0;
    array[grpc_theta_dim] int grpc_ii;
    array[grpc_theta_dim] int grpc_jj;
    int grpc_nfi;
    array[grpc_nfi] int grpc_ifi;
    array[grpc_nfi] int grpc_jfi;
    ")

    ## parameters input
    aC$parameters <- "
    vector[grpc_theta_dim] grpc_theta;
    "

    ## transformed parameters
    aC$"transformed parameters" <- "
    vector[grpc_theta_dim] grpc_xi = grpc_Ibase_half * (grpc_theta - grpc_theta_base);
    real<lower=0> grpc_radius = dot_self(grpc_xi);
    matrix[grpc_dim,grpc_dim] grpc_L0;
    for(i in 1:grpc_dim) {
      for(j in 1:grpc_dim) {
        if(i==j) {
          grpc_L0[i,i] = grpc_d0[i];
        } else {
          grpc_L0[i,j] = 0.0;
        }
      }
    }
    {
      int i,j;
      for(k in 1:grpc_theta_dim) {
        i = grpc_ii[k];
        j = grpc_jj[k];
        grpc_L0[i,j] = grpc_theta[k];
      }
      if(grpc_nfi>0) {
        for(l in 1:grpc_nfi) {
          i = grpc_ifi[l];
          j = grpc_jfi[l];
          if(j>0) {
            for(k in 1:(j-1)) {
              grpc_L0[i,j] -= grpc_L0[i,k] * grpc_L0[j,k] / grpc_L0[j,j];
            }
          }
        }
      }
    }
    matrix[grpc_dim,grpc_dim] grpc_V0 = chol2inv(grpc_L0);
    vector[grpc_dim] grpc_s0inv = inv_sqrt(diagonal(grpc_V0));
    graphpcor_cov_mat = quad_form_diag(grpc_V0, grpc_s0inv);
    "

    ## model part
    aC$model <- "
    grpc_radius ~ exponential(grpc_lambda);
    target += lgamma(grpc_theta_dim * 0.5) -grpc_theta_dim*0.5 * log(pi());
    target += grpc_logDetIhalf - (grpc_theta_dim-1.0)*log(grpc_radius);
  "


    for(i in seq_along(aC)) {
      aC[[i]] <- gsub("graphpcor_cov_mat", name, aC[[i]], fixed = TRUE)
    }

    return(stan_add_code(x, aC))
  }

  warning("Nothing done!")
  return(x)
}
#' @describeIn stan_add add code at the end of each section
#' @param to_add named list with the code to be added
#' @importFrom utils tail
stan_add_code <- function(x, to_add) {
  all_sec_names <- c(
    "functions", "data", "transformed data",
    "parameters", "transformed parameters",
    "model", "generated quantities")
  stopifnot(sum(names(to_add) %in% all_sec_names)>0)
  names(all_sec_names) <- c(
    "fnc", "dat", "tdt", "par", "tpr", "mod", "gqt")
  nsec <- length(all_sec_names)
  sec_ini <- lapply(all_sec_names, gregexpr, ## not only the first
                    x, fixed = TRUE)
  if(any(sapply(sec_ini, length)>c(1,2,1,2,1,1,1))) {
    stop("Code fix needed!")
    ## to do: if section names are used in code as variable
  }
  if(any(sec_ini$tdt[[1]]>0)) {
    if(length(sec_ini$dat[[1]])==1) {
      sec_ini$dat[[1]] <- -1
    }
    if(length(sec_ini$dat[[1]])==2) {
      sec_ini$dat[[1]] <- sec_ini$dat[[1]][1]
    }
  }
  if(any(sec_ini$tpr[[1]]>0)) {
    if(length(sec_ini$par[[1]])==1) {
      sec_ini$par[[1]] <- -1
    }
    if(length(sec_ini$par[[1]])==2) {
      sec_ini$par[[1]] <- sec_ini$par[[1]][1]
    }
  }
  sec_ini <- sapply(sec_ini, function(x) x[[1]][1])
  icurl <- gregexpr("}", x, fixed = TRUE)[[1]]
  sec_end <- sec_ini
  for(i in 1:(nsec-1)) {
    if(sec_ini[i]>0) {
      a <- icurl[icurl>sec_ini[i]]
      j <- i + which(sec_ini[(i+1):nsec]>0)
      if(length(j)==0) {
        b <- tail(icurl, 1)
      } else {
        b <- icurl[icurl<sec_ini[j[1]]]
      }
      sec_end[i] <- intersect(a, b)
    }
  }
  if(sec_ini[nsec]>0) {
    sec_end[nsec] <- tail(icurl, 1)
  }

  isx <- pmatch(names(to_add), all_sec_names)

  secs <- lapply(1:nsec, function(x) "")
  names(secs) <- all_sec_names
  for(i in 1:nsec) {
    if(sec_ini[i]>0) {
      secs[[i]] <- substr(x, sec_ini[i], sec_end[i]-1)
    }
  }
  for(i in seq_along(to_add)) {
    ini <- sec_ini[isx[i]]
    if(ini<0) {
      secs[[isx[i]]] <- paste0(
        "\n", names(to_add)[i], " {\n",
        to_add[[i]], "\n")
    } else {
      end <- sec_end[isx[i]]
      secs[[isx[i]]] <- paste0(secs[[isx[i]]], to_add[[i]])
    }
  }
  for(i in 1:nsec) {
    if (nchar(secs[[i]])>0)
      secs[[i]] <- paste0(secs[[i]], "}\n")
  }
  return(Reduce("paste0", secs))
}
