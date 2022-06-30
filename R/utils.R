#' function to transform observed bottom-level series into dummy series.
#' 
#' @import Matrix
#' @param x dhts object
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
#' @param incoherent_domain
#' @param coherent_domain
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
#' @return joint distribution matrix
#' @export 
marginal2Joint <- function(x, domain){
  cls <- "jdist-ind"
  if (is.coherent_domain(domain)){
    domain <- domain[,2:dim(domain)[2]]
    cls <- "jdist-bu"
  }
  x <- ifelse(is.coherent_domain(domain), x[,2:dim(domain)], x)
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

marginal2Sum <- function(x, domain){
  if (length(x) == 1){
    return(x[[1]])
  }
  time_window <- dim(x[[1]])[1]
  
  # expand the domain
  domains <- expand.grid(data.frame(domain))
  
  tprobf <- lapply(x, function(f){ f[1,] })
  basef <- apply(expand.grid(tprobf), 1, prod)
  ans <- sapply(split(basef, apply(domains, 1, sum)), sum)
  for (i in 2:time_window){
    tprobf <- lapply(x, function(f){ f[i,] })
    basef <- apply(expand.grid(tprobf), 1, prod)
    ans <- rbind(ans, sapply(split(basef, apply(domains, 1, sum)), sum))
  }
  rownames(ans) <- NULL
  ans[,sort(colnames(ans))]
}


#' function to convert Joint distribution to marginal Distribtion
#' 
#' @param x joint distribution
#' @param domain domain of hierarchy, incoherent or cohernet
#' @param which integer indicating which dimension (which column of domain), if
#' NULL, return marginal distribution of all series.
#' @return marginal distribution
Joint2Marginal <- function(x, domain, which=NULL){
  time_window <- dim(x)[1]
  marginal <- NULL
  if (is.null(which)){
    res = list()
    for (i in 1:dim(domain)[2]){
      res[[i]] <- Joint2Marginal(x, domain, which = i)
    }
    return(res)
  } else {
    for (i in 1:time_window){
      marginal <- rbind(marginal, 
                        sapply(split(x[i,], domain[,which]), sum))
    }
    return(marginal)
  }
}



#' function to construct summing matrix of a simple two-level hierarchy
#' 
#' @param m number of children nodes
#' @return summing matrix
twolevelhierarchy <- function(m){
  rbind(rep(1, m), diag(rep(1, m)))
}
