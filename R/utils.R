#' function to transform observed bottom-level series into dummy series.
#' 
#' @param x time series
#' @param hier dhier object
#' @return dummy matrix
cons_realDummy <- function(x, hier) {
  series <- x %*% t(hier$s_mat) 
  dummy_mat <- spMatrix(NROW(series), hier$quantities['r'])
  for (i in 1:dim(series)[1]) {
    dummy_mat[i, which(apply(hier$coherent_domain, 1, function(x) {
      all(x == series[i, ])
    }))] = 1
  }
  dummy_mat
}

#' cost of moving probabilities
#‘
#' function to calculate distance between coherent point
#' @param hier dhier object
#' @return distance matrix
cal_costeMatrix <- function(hier) {
  cd <- hier$coherent_domain
  id <- hier$incoherent_domain

  distance = matrix(NA, nrow = NROW(cd), ncol = NROW(id))
  for (j in 1:NROW(cd)) {
    for (k in 1:NROW(id)) {
      distance[j, k] = sum(abs(cd[j, ] - id[k, ]))
    }
  }
  distance
}


#' convert marginal distributions into joint distribution assuming independence
#' 
#' @param x list of distributions of all series
#' @param hier dhts object
#' @return joint distribution matrix
marginal2Joint <- function(x, hier, method = c("ind", "bu")){
  
  domain <- hier$incoherent_domain
  method <- match.arg(method)
  
  if (method == "ind") {
    domain <- hier$incoherent_domain
  } else {
    n <- NROW(hier$s_mat)
    m <- NCOL(hier$s_mat)
    x <- x[(n-m+1):n]
    domain <- hier$coherent_domain[, (n-m+1):n]
  }
  
  time_window <- NROW(x[[1]])
  res <- NULL
  for (j in 1:NROW(domain)){
    tmp <- 1
    for (i in 1:NCOL(domain)){
      indx <- domain[j, i]
      tmp <- tmp * x[[i]][,paste0(indx)]
    }
    res <- cbind(res, tmp)
  }

  unname(res)
}

#' function to compute distribution of sum of bottom series assuming independence
#' 
#' @param x list of distributions of some series
#' @param hier dhier obj
#' @param which indicating which upper nodes, default NULL means all upper series.
#' @return distribution of upper series.
marginal2Sum <- function(x, hier, which = NULL){
  all_ts <- marginal2Joint(x, hier, method = "bu")
  Joint2Marginal(all_ts, hier, which)
}


#' function to convert Joint distribution to marginal Distribution
#' 
#' @param x joint distribution
#' @param hier dhier object
#' @param which integer indicating which dimension (which column of domain), if
#' NULL, return marginal distribution of all series.
#' @param coherent If the passed joint distribution coherent
#' @return marginal distribution
Joint2Marginal <- function(x, hier, which=NULL, coherent=TRUE){
  n <- NROW(hier$s_mat)
  m <- NCOL(hier$s_mat)
  time_window <- dim(x)[1]
  if (coherent) domain <- hier$coherent_domain
  else domain <- hier$incoherent_domain
  
  if (is.null(which)) which <- 1:n
  
  output <- list()
  
  output <- lapply(which, function(i){
    nums <- sort(unique(domain[,i]))
    
    tmp <- sapply(nums, function(y){
      idx <- domain[,i] == y
      if(sum(idx) == 1) return(x[, idx])
      else return(apply(x[, idx, drop=FALSE], 1, sum))
    }, simplify = "array")
    if (is.null(dim(tmp))) tmp <- matrix(tmp, 1)
    colnames(tmp) <- nums
    tmp
  })
  
  if (length(which) == 1) return(output[[1]])
  names(output) <- colnames(hier$coherent_domain)[which]
  output
}





