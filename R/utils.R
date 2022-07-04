#' function to transform observed bottom-level series into dummy series.
#' 
#' @import Matrix
#' @param x dhts object
#' @tag multiple-levels
#' @return dummy matrix
cons_realDummy <- function(x) {
  if (!is.dhts(x)) stop("Argument x should be a dhts object.")
  coherent_domain <- x$domain$coherent_domain
  series <- aggdhts(x)
  r = dim(coherent_domain)[1]
  dummy_mat <- as(matrix(0, dim(series)[1], r), "sparseMatrix")
  for (i in 1:dim(series)[1]) {
    dummy_mat[i, which(apply(coherent_domain, 1, function(x) {
      all(x == series[i, ])
    }))] = 1
  }
  dummy_mat
}

#' cost of moving probabilities
#‘
#' function to calculate distance between coherent point
#' @param incoherent_domain coherent domain
#' @param coherent_domain incoherent domain
#' @tag multiple-levels
#' @return distance matrix
cal_costeMatrix <- function(incoherent_domain, coherent_domain) {
  stopifnot(is.coherent_domain(coherent_domain),
            is.incoherent_domain(incoherent_domain))
  r = dim(coherent_domain)[1]
  q = dim(incoherent_domain)[1]
  distance = matrix(NA, nrow = r, ncol = q)
  for (j in 1:r) {
    for (k in 1:q) {
      distance[j, k] = sum(abs(coherent_domain[j, ] - incoherent_domain[k, ]))
    }
  }
  distance
}


#' convert marginal distributions into joint distribution assuming independence
#' 
#' @param x list of distributions of all series
#' @param domain domain of hierarchy, if coherent, only bottom series are used (bottom-up approach).
#' @tag multiple-levels
#' @return joint distribution matrix
#' @export 
marginal2Joint <- function(x, domain){
  cls <- "jdist-ind"
  m <- attr(domain, 'm')
  n <- dim(domain)[2]
  if (is.coherent_domain(domain)){
    x <- x[(n-m+1):n]
    domain <- domain[, (n-m+1):n]
    cls <- "jdist-bu"
  }
  time_window <- dim(x[[1]])[1]
  res <- NULL
  for (j in 1:dim(domain)[1]){
    tmp <- 1
    for (i in 1:length(x)){
      indx <- domain[j, i]
      tmp <- tmp * x[[i]][,indx+1]
    }
    res <- cbind(res, tmp)
  }
  rownames(res) <- NULL
  colnames(res) <- NULL
  structure(res, class=cls)
}

#' function to compute distribution of sum of bottom series assuming independence
#' 
#' @param x list of distributions of some series
#' @param domain should be coherent domain
#' @param which indicating which upper nodes, default NULL means all upper series.
#' @return distribution of upper series.
#' @tag 
marginal2Sum <- function(x, domain, which = NULL){
  stopifnot(is.coherent_domain(domain))
  time_window <- dim(x[[1]])[1]
  
  m <- attr(domain, 'm')
  n <- dim(domain)[2]
  
  d <- domain[,(n-m+1):n]
  
  output <- list()
  if (is.null(which)){
    which <- 1:(n-m)
  }
  for (i in seq_along(which)){
    ds <- unique(domain[,which[i]])
    output[[i]] <- matrix(0, time_window, length(ds))
    colnames(output[[i]]) <- ds
    for (j in 1:length(domain[,which[i]])){
      current_col <- as.character(domain[j, which[i]])
      tmp <- 1
      for (k in 1:m){
        tmp = tmp * x[[n-m+k]][,as.character(d[j, k])]
      }
      output[[i]][, current_col] = output[[i]][, current_col] + tmp
    }
  }
  if (length(output) == 1) return(output[[1]])
  output
}


#' function to convert Joint distribution to marginal Distribution
#' 
#' @param x joint distribution
#' @param domain domain of hierarchy, incoherent or coherent
#' @param which integer indicating which dimension (which column of domain), if
#' NULL, return marginal distribution of all series.
#' @tag multiple-levels
#' @return marginal distribution
Joint2Marginal <- function(x, domain, which=NULL){
  n <- dim(domain)[2]
  m <- attr(domain, "m")
  time_window <- dim(x)[1]
  if (!("jdist-bu" %in% class(x) | "jdist-rec" %in% class(x) | "jdist-ind" %in% class(x))){
    stop("only bottom-up joint distritbution or reconciled distribution is supported!")
  }
  
  if (is.null(which)){
    output <- list()
    for (i in 1:n){
      if ("jdist-bu" %in% class(x)){
        if (i <= n-m){
          output[[i]] <- NULL
          next
        }
      }
      output[[i]] <- Joint2Marginal(x, domain, i)
    }
  } else {
    if ("jdist-bu" %in% class(x)){
      if (which <= n-m) {
        warning("can not obtain marginal distribution of upper levels given bottom up joint distribution")
        return(NULL)
      }
    }
    output <- NULL
    for (i in 1:time_window){
      output <- rbind(output, 
                      sapply(split(x[i,], domain[,which]), sum))
    }
  }
  output
}

