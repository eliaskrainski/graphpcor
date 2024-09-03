#' Directed Tree Graph - DTG definitions an dmethods
#' @param ... a list of formula used as relationship
#' to define the Directed Tree Graph - DTG.
#' Parent nodes shall be in the right side while children
#' (or parents with a parent) in the left side.
#' Example 1: dtg(p1 ~ c1 + c2
#' Example 2: dtg(p1 ~ c1 + c2 + p2, p2 ~ c3)
#' @export
#' @examples
#' g1 <- dtg(p1 ~ c1 + c2 - c3)
#' g1
#' g2 <- dtg(p1 ~ c1 + c2 + p2, p2 ~ c3 - c4)
#' g2
dtg <- function(...) {

  fch <- as.character(match.call())[-1]
##  fch <- lapply(fch0, as.formula)
  if(length(fch)<1)
    stop("Please provide an argument!")

  ch <- lapply(fch, function(x)
    as.character(as.formula(x)))

  ## left side check
  pp <- sapply(ch, function(x) x[2])
  if (length(pp)>length(unique(pp))) {
    stop("Non-unique parent definition!")
  }
  if(!all(sapply(pp, substr, 1, 1) == "p"))
    stop("Please use the letter 'p' for parents!")
  if(any(is.na(as.integer(
    sapply(ch, function(x)
      substring(x[2], 2)))))) {
    stop("Please use integer after letter 'p' for parents!")
  }

  ## order
  opp <- order(pp)
  fch <- fch[opp]
  ch <- ch[opp]
  pp <- pp[opp]
  m <- length(pp)
  P <- sort(pp)

  ## right side check, and collect terms
  terms.l <- vector("list", m)
  terms.i <- vector("list", m)
  terms.s <- vector("list", m)
  for(i in 1:m) {
    x <- gsub(" ", "", ch[[i]][3])
##    print(c(x = x))
    if(substr(x, 1, 1) != "-")
      x <- paste0("+", x)
  ##  print(c(x = x))
    i0 <- gregexpr("-", x, fixed = TRUE)[[1]]
    i1 <- gregexpr("+", x, fixed = TRUE)[[1]]
    si <- integer(max(i0, i1))
    if(any(i0>0))
      si[i0] <- -1
    if(any(i1>0))
      si[i1] <- +1
    ##print(si)
    terms.s[[i]] <- si[si != 0]
    schi <- strsplit(x, "-", fixed = TRUE)[[1]]
##    print(schi)
    schi <- unlist(strsplit(schi, "+", fixed = TRUE))
  ##  print(schi)
    if(schi[1]=="") schi <- schi[-1]
    ch.l <- substr(schi, 1, 1)
    ##print(c(ch.l))
    if(!all(ch.l %in% c("p", "c")))
      stop("Invalid variable labeling in ~ ", ch[[i]][3])
    jj.i <- as.integer(substring(schi, 2))
    ##print(jj.i)
    if(any(is.na(jj.i)) | any(jj.i<1))
      stop("Invalid variable numbering in ~ ", ch[[i]][3])
    if(any(ch.l == "p")) {
      if(any(jj.i[ch.l == "p"] <= 1))
        stop("Parent id in ~ ", ch[[i]][3], " must be >", 1)
      if(any(jj.i[ch.l == "p"] > m))
        stop("Parent id in ~ ", ch[[i]][3], " should not be >", m)
    }
    terms.l[[i]] <- ch.l
    terms.i[[i]] <- jj.i
  }

##  print(terms.l)
  ##print(terms.i)
  ##print(terms.s)

  ic <- unlist(terms.l) == "c"
  if(any(ic)) {
    n <- max(unlist(terms.i)[ic])
  } else{
    n <- 0L
  }

  trm <- matrix(0L, n+m, m)
  colnames(trm) <- paste0("p", 1:m)
  if(n>0) {
    rownames(trm) <- c(
      paste0("c", 1:n),
      paste0("p", 1:m)
    )
  } else {
    stop("There is no children variable!")
  }

  for(i in 1:m) {
    icc <- terms.l[[i]] == "c"
    if(any(icc)) {
      trm[terms.i[[i]][icc], i] <- terms.s[[i]][icc]
    }
    npc <- sum(!icc)
    if(npc>0) {
      trm[n + terms.i[[i]][!icc], i] <- terms.s[[i]][!icc]
    }
  }

  trm <- trm[-(n+1), , drop = FALSE]

  nr <- rowSums(trm != 0)
  if(any(nr==0)) {
    stop("Missing definition for ",
         rownames(trm)[nr==0])
  }

  class(fch) <- "dtg"
  attr(fch, "childrens") <- n
  attr(fch, "parents") <- m
  attr(fch, "relations") <- trm

  return(fch)

}

#' @export
print.dtg <- function(x, ...) {
  cat("DTG for", attr(x, "childrens"),
      "childrens with",
      attr(x, "parents"), "parents\n")
  for(i in 1:length(x)) {
    cat(x[[i]], "\n")
  }
}

#' @export
summary.dtg <- function(object, ...) {
  attr(object, "relations")
}


