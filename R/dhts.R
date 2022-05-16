#' discrete hierarchical time series constructor.
#' 
#' Create a 'dhts' object 
#' @rdname dhts-class
#' @param domain the domain of bottom-level series, 2 * m matrix.
#' Each column corresponds to a series, and the first row corresponds to
#' the minimums, and the second row corresponds the maximums.
#' @param s_mat the summing matrix
#' @param bts bottom-level time series, T * m
#' @return a dhts object
#' @examples 
#' # create a 3-nodes discrete hts.
#' bts <- matrix(sample(0:1, 100, replace=TRUE), ncol=2)
#' s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
#' domain <- matrix(c(1, 0, 1, 0), 2)
#' dhts(bts, s_mat, domain)
#' @export
dhts <- function(bts, s_mat, domain_bts){
  stopifnot(dim(bts)[2] == dim(s_mat)[2],
            dim(domain_bts)[2] == dim(s_mat)[2],
            dim(domain_bts)[1] == 2)
  
  if (!is.ts(bts)){
    bts <- as.ts(bts)
  }
  
  domain <- list(
    domain_bts = domain_bts,
    incoherent_domain = cons_domain(domain_bts, s_mat, FALSE),
    coherent_domain = cons_domain(domain_bts, s_mat))
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




