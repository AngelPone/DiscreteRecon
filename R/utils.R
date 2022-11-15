#' function to transform observed bottom-level series into dummy series.
#' 
#' @import Matrix
#' @param x dhts object
#' @tag multiple-levels
#' @return dummy matrix
cons_realDummy <- function(x) {
  if (!is.dhts(x)) stop("Argument x should be a dhts object.")
  coherent_domain <- x$meta$coherent_domain
  series <- aggdhts(x)
  r = NROW(coherent_domain)
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
#' @param meta metadata of dhts object
#' @tag multiple-levels
#' @return distance matrix
cal_costeMatrix <- function(meta) {
  cd <- meta$coherent_domain
  id <- meta$incoherent_domain

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
#' @param obj dhts object
#' @param method bu for bottom-up, produce coherent joint distribution; 
#' ind for independent, producing incoherent base joint distribution.
#' @tag multiple-levels
#' @return joint distribution matrix
#' @export 
marginal2Joint <- function(x, meta, method){
  stopifnot(is.list(x), length(x) == NROW(meta$s_mat))
  stopifnot(method %in% c("bu", "ind"))
  
  n <- NROW(meta$s_mat)
  m <- NCOL(meta$s_mat)
  
  
  domain <- meta$incoherent_domain
  cls <- c("incoherent", "jdist")
  if (method == "bu"){
    x <- x[(n-m+1):n]
    domain <- meta$coherent_domain[, (n-m+1):n]
    cls <- c("coherent", "jdist", "bu")
  }
  
  time_window <- dim(x[[1]])[1]
  res <- NULL
  for (j in 1:NROW(domain)){
    tmp <- 1
    for (i in 1:NCOL(domain)){
      indx <- domain[j, i]
      tmp <- tmp * x[[i]][,paste0(indx)]
    }
    res <- cbind(res, tmp)
  }

  structure(unname(res), class=cls)
}

#' function to compute distribution of sum of bottom series assuming independence
#' 
#' @param x list of distributions of some series
#' @param obj dhts obj
#' @param which indicating which upper nodes, default NULL means all upper series.
#' @return distribution of upper series.
#' @tag 
marginal2Sum <- function(x, meta, which = 1){
  all_ts <- marginal2Joint(x, meta, method = "bu")
  Joint2Marginal(all_ts, meta, which)
}


#' function to convert Joint distribution to marginal Distribution
#' 
#' @param x joint distribution
#' @param obj dhts object
#' @param which integer indicating which dimension (which column of domain), if
#' NULL, return marginal distribution of all series.
#' @tag multiple-levels
#' @return marginal distribution
Joint2Marginal <- function(x, meta, which=NULL){
  n <- NROW(meta$s_mat)
  m <- NCOL(meta$s_mat)
  time_window <- dim(x)[1]
  if (!is_jdist(x)){
    stop("x shoule be one kind of joint distribution")
  }
  
  if (is_coherent(x)) domain <- meta$coherent_domain
  else domain <- meta$incoherent_domain
  
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
  names(output) <- colnames(meta$coherent_domain)[which]
  output
}

is_jdist <- function(x){
  "jdist" %in% class(x)
}
is_coherent <- function(x){
  "coherent" %in% class(x)
}

is_coherentJdist <- function(x){
  is_jdist(x) & is_jdist(x)
}
is_incoherentJdist <- function(x){
  ("incoherent" %in% class(x)) & is_jdist(x)
}

