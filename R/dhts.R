#' discrete hierarchical time series constructor.
#' 
#' Create a 'dhts' object 
#' @rdname dhts-class
#' @param domain the domain of bottom-level series, 2 * m matrix.
#' Each column corresponds to a series, and the first row corresponds to
#' the minimums, and the second row corresponds the maximums.
#' @param s_mat the summing matrix
#' @param bts bottom-level time series, T * m
#' @tag multiple-levels
#' @return a dhts object
#' @examples 
#' # create a 3-nodes discrete hts.
#' bts <- matrix(sample(0:1, 100, replace=TRUE), ncol=2)
#' s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
#' domain <- matrix(c(1, 0, 1, 0), 2)
#' dhts(bts, s_mat, domain)
#' @export
dhts <- function(bts, s_mat, domain_bts, node_names=NULL){
  stopifnot(dim(bts)[2] == dim(s_mat)[2],
            dim(domain_bts)[2] == dim(s_mat)[2],
            dim(domain_bts)[1] == 2)
  for (i in 1:dim(s_mat)[2]){
    stopifnot(bts[,i] >= domain_bts[1,i],
              bts[,i] <= domain_bts[2,i])
  }
  
  if (!is.ts(bts)){
    bts <- as.ts(bts)
  }
  
  domain <- list(
    domain_bts = domain_bts,
    incoherent_domain = cons_domain(domain_bts, s_mat, FALSE, node_names = node_names),
    coherent_domain = cons_domain(domain_bts, s_mat, node_names = node_names))
  structure(
    list(bts = bts, 
         domain = domain, 
         s_mat = s_mat),
    class = c("dhts")
  )
}

#' aggregate discrete hierarchical time series
#' 
#' function for aggregting discrete hierarchical time series.
#' 
#' @param dhts a dhts object
aggdhts <- function(x){
  if (!is.dhts(x)) {
    stop("Argument must be a dhts object", call. = FALSE)
  }
  t(x$s_mat %*% t(x$bts))
}


#' A function to check if a object is \code{dhts}.
#' @rdname dhts-class
#' @param dhts object
#' @export
is.dhts <- function(dhts){
  is.element("dhts", class(dhts))
}

#' construct domain
#' 
#' Method for constructing coherent domain vector given domain matrix of all 
#' bottom-level series and the summing matrix.
#' @param domain_bts domain matrix used in creating dhts object.
#' @param s_mat summing matrix
#' @param coherent logical indicating if produced domain is coherent or not.
#' @tag multiple-levels
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
    node_names <- paste0('s', 1:dim(s_mat)[1])
  }
  allDomain <- as.matrix(allDomain)
  colnames(allDomain) <- node_names
  class(allDomain) <- 
  structure(allDomain, class=ifelse(coherent, "coherent_domain", "incoherent_domain"),
            m = dim(s_mat)[2])
}

#' @export
is.coherent_domain <- function(x){
  "coherent_domain" %in% class(x)
}

#' @export
is.incoherent_domain <- function(x){
  "incoherent_domain" %in% class(x)
}



