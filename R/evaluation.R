#' Brier Score
#' 
#' function to compute brier score given observations and probabilistic forecasts
#' @param probf joint probabilistic forecasts, T * n
#' @param y dhts
#' @export
brier_score <- function(probf, y) {
  time_window <- dim(probf)[1]
  stopifnot(is(y, "dhts"))
  if (is(probf, "jdist-bu") | is(probf, "jdist-rec") | is(probf, "jdist-td")){
    bs <- sum((probf - cons_realDummy(y))^2)/time_window
  } else if (is(probf, "jdist-ind")){
    flags <- getCoherentFlag(y$domain$incoherent_domain, y$domain$coherent_domain)
    cps <- probf[,flags$coherent]
    ips <- probf[,flags$incoherent]
    cps <- sum((cps - cons_realDummy(y))^2)/time_window
    return(c(coherency=cps, incoherency=sum(ips^2)/time_window))
  }
}

#' utility functions
getCoherentFlag <- function(idomain, cdomain){
  output <- list()
  coherent <- NULL
  for (i in 1:dim(idomain)[1]){
    for (j in 1:dim(cdomain)[1]){
      if (all(idomain[i,] == cdomain[j,])){
        coherent <- c(coherent, i)
      }
    }
  }
  list(coherent=coherent, incoherent=setdiff(1:dim(idomain)[1], coherent))
}


#' Point forecast
#' 
#' function to generate point forecasts of the hierarchy
#' @param dist coherent distribution
#' @param domain coherent domain
#' @return point forecasts of all series
#' @export 
point_forecast <- function(dist, domain){
  fms <- matrix(0, nrow = dim(dist)[1], ncol = dim(domain)[2])
  for (i in 1:dim(domain)[2]){
    fi <- Joint2Marginal(dist, domain, i)
    ds <- as.numeric(colnames(fi))
    fm <- apply(fi, 1, function(x){sum(ds * x)})
    fms[,i] <- fm
  }
  colnames(fms) <- colnames(domain)
  return(fms)
}

#' Point metrics of each series
#' @param dist
#' @param y dhts
#' @param f function to calculate metric
#' @return vector containing metric of each series
#' @export
point_metric <- function(dist, y, f){
  pointf <- point_forecast(dist, y$domain$coherent_domain)
  all_y <- aggdhts(y)
  sapply(1:dim(all_y)[2], function(x){
    f(pointf[,x], all_y[,x])
  })
}

