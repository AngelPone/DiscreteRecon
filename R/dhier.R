#' discrete hierarchy constructor.
#' 
#' Create a 'dhier' object 
#' @rdname dhier-class
#' @param domain the domain of bottom-level series, 2 * m matrix.
#' Each column corresponds to a series, and the first row corresponds to
#' the minimums, and the second row corresponds the maximums.
#' @param s_mat the summing matrix
#' @return a dhier object
#' @examples 
#' # create a 3-nodes `dhier` object.
#' s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
#' domain <- matrix(c(0, 1, 0, 1), 2)
#' dhier(s_mat, domain)
#' @export
dhier <- function(s_mat, domain, node_names=NULL){
  stopifnot(NCOL(domain) == NCOL(s_mat),
            NROW(domain) == 2)
  stopifnot(all(domain[1,] <= domain[2,]))
  
  if (is.null(node_names)){
    node_names <- paste0('s', 1:NROW(s_mat))  
  }
  
  idomain <- cons_domain(domain, s_mat, FALSE, node_names)
  cdomain <- cons_domain(domain, s_mat, TRUE, node_names)
  
  
  structure(
    list(domain = domain,
         s_mat = s_mat,
         incoherent_domain = idomain,
         coherent_domain = cdomain,
         coherent_flags = getCoherentFlag(idomain, cdomain, s_mat),
         quantities = c(r=NROW(cdomain), q=NROW(idomain),
                        n=NROW(s_mat), m=NCOL(s_mat))),
    class = c("dhier")
  )
}


#' @export
print.dhier <- function(x){
  cat("Hierarchy with series taking discrete values:\n")
  cat(paste0("n=", unname(x$quantities['n']),
      ", m=", unname(x$quantities['m']),
      ", r=", unname(x$quantities['r']),
      ", q=", unname(x$quantities['q']), "."))
}


#' utility functions
getCoherentFlag <- function(idomain, cdomain, smat){
  output <- list()
  coherent <- NULL
  
  n <- NROW(smat)
  m <- NCOL(smat)
  comp <- t(smat[1:(n-m),] %*% t(idomain[,(n-m+1):n]))
  coherent <- which(rowSums(comp == idomain[,1:(n-m)]) == n-m)
  list(coherent=coherent, incoherent=setdiff(1:dim(idomain)[1], coherent))
}


#' @export
is.dhier <- function(dhier){
  is.element("dhier", class(dhier))
}

#' construct domain
#' 
#' Method for constructing coherent domain vector given domain matrix of all 
#' bottom-level series and the summing matrix.
#' @param domain domain matrix used in creating dhier object.
#' @param s_mat summing matrix
#' @param coherent logical indicating if produced domain is coherent or not.
#' @return domain object
cons_domain <- function(domain, s_mat, coherent, node_names) {
  if (coherent){
    # df of coherent domain
    df <- expand.grid(
      # domain of each variable
      apply(domain, 2, function(x){x[1]:x[2]}, simplify = FALSE)
    )
    allDomain <- t(s_mat %*% t(df))
  } else {
    df <- t(s_mat %*% t(domain))
    allDomain <- expand.grid(
      # domain of each variable
      apply(df, 2, function(x){x[1]:x[2]}, simplify = FALSE)
    )
  }
  allDomain <- unname(as.matrix(allDomain))
  colnames(allDomain) <- node_names
  structure(allDomain, class=ifelse(coherent, "coherent_domain", "incoherent_domain"))
}



