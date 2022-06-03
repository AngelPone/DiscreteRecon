#' function to transform observed bottom-level series into dummy series.
#' 
#' @import Matrix
#' @param x dhts object
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

#' construct domain
#' 
#' Method for constructing coherent domain vector given domain matrix of all 
#' bottom-level series and the summing matrix.
#' @param domain_bts domain matrix used in creating dhts object.
#' @param s_mat summing matrix
#' @param coherent logical indicating if produced domain is coherent or not.
#' @return domain object
cons_domain <- function(domain_bts, s_mat, coherent = TRUE, node_names=NULL) {
  if (coherent){
    # df of coherent domain
    df <- expand.grid(
      # domain of each variable
      apply(domain_bts, 2, function(x){x[1]:x[2]}, simplify = FALSE)
    )
    allDomain <- t(s_mat %*% t(df))
  } else {
    df <- t(s_mat %*% t(domain_bts))
    allDomain <- expand.grid(
      # domain of each variable
      apply(df, 2, function(x){x[1]:x[2]}, simplify = FALSE)
    )
  }
  if (is.null(node_names)){
    node_names <- 1:dim(allDomain)[2]
  }
  colnames(allDomain) <- node_names
  return(allDomain)
}


#' cost of moving probabilities
#' 
#' function to calculate distance between coherent point
cal_costeMatrix <- function(incoherent_domain, coherent_domain) {
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


#' method for incoherent probabilistic forecasts
#' @export 
#' @note This method assume the list containing base forecasts is ordered by the 
#' order of summing matrix and probabilities are arranged in order of domain.
marginal2Joint <- function(x){
  time_window <- dim(x[[1]])[1]
  domains <- expand.grid(lapply(x, function(x){0:(dim(x)[2] - 1)}))
  res <- NULL
  for (j in 1:dim(domains)[1]){
    tmp <- 1
    for (i in 1:length(x)){
      indx <- domains[j, i]
      tmp <- tmp * x[[i]][,indx+1]
    }
    res <- cbind(res, tmp)
  }
  # tprobf <- lapply(x, function(f){ f[1,] })
  # basef <- apply(expand.grid(tprobf), 1, prod)
  # for (i in 2:time_window){
  #   tprobf <- lapply(x, function(f){ f[i,] })
  #   basef <- rbind(basef, apply(expand.grid(tprobf), 1, prod))
  # }
  rownames(res) <- NULL
  colnames(res) <- NULL
  res
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
#' @param domain domain of hierarchy
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
