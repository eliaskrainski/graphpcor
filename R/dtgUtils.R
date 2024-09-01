#' Directed Tree Graph - DTG definitions an dmethods
#'
#' @export
#' @examples
#' g <- dtg(p1 ~ c1 + c2 - c3)
#' g
dtg <- setClass('dtg', contains = "formula")

setValidity("dtg",
  function(object) {
    if(length(object)<1)
      stop("Please provide a valid object!")
    if(inherits(object, "formula")) {
      object <- list(object)
    }
    lch <- lapply(object, as.character)
    pp <- sapply(lch, function(x) x[2])
    if (length(pp)>length(unique(pp))) {
      stop("Non-unique parent definition!")
    }
    if(!all(sapply(pp, substr, 1, 1) == "p"))
      stop("Please use the letter 'p' for parents!")
    if(any(is.na(as.integer(
      sapply(lch, function(x)
        substring(x[2], 2)))))) {
      stop("Please use integer after letter 'p' for parents!")
    }
    for(i in 1:length(lch)) {
      x <- gsub(" ", "", lch[[i]])
      s0 <- strsplit(x[3], "-", fixed = TRUE)[[1]]
      s1 <- unlist(strsplit(s0, "+", fixed = TRUE))
      c1 <- all(substr(s1, 1, 1) %in% c("p", "c"))
      if(!c1) stop("Invalid variable labeling in", object[[i]], "!")
      cn <- any(is.na(as.integer(substring(s1, 2))))
      if(cn) stop("Invalid variable numbering in", object[[i]], "!")
    }
    return(TRUE)
  }
)

