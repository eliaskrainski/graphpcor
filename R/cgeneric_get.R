#' Function to extract cgeneric model
#' @param cgeneric_model an object containing the cgeneric model
#' @param theta a numeric vector with theta.
#' If NULL (default) will use the initial from the cgeneric model.
#' @param cmd an string to specify which model element to get
#' @export
cgeneric_get <- function(cgeneric_model,
                         cmd = c("graph", "Q", "initial", "mu", "log.prior"),
                         theta = NULL
                         ) {

  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)

  cgdata <- cgeneric_model$f$cgeneric$data
  stopifnot(!is.null(cgdata))
  stopifnot(!is.null(cgdata$ints))
  stopifnot(!is.null(cgdata$characters))

  cmds <- c("graph", "Q", "initial", "mu", "log_prior")
  cmd <- match.arg(cmd,
                   cmds,
                   several.ok = TRUE)
  stopifnot(length(cmd)>0)

  itheta <- .Call(
    "cgeneric_element_get",
    "initial",
    NULL,
    cgdata$ints,
    cgdata$doubles,
    cgdata$characters,
    cgdata$matrices,
    cgdata$smatrices,
    PACKAGE = "corGraphs"
  )

  if(length(cmd) == 1){
    if(cmd == "initial") {
      return(itheta)
    } else {
      return(
        .Call(
          "cgeneric_element_get",
          cmd,
          theta,
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "corGraphs"
          )
        )
      }
    } else {
    names(cmd) <- cmd
    if(is.null(theta)) {
      theta = itheta
    }
    return(
      lapply(cmd, function(x) {
        .Call(
          "cgeneric_element_get",
          x,
          theta,
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "corGraphs"
        )
      })
    )
  }

}
