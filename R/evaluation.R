#' Brier Score
#' 
#' function to compute brier score given observations and probabilistic forecasts
#' @param probf joint probabilistic forecasts, T * n
#' @param y dhts
#' @export
brier_score <- function(probf, y) {
  time_window <- dim(probf)[1]
  stopifnot(is(y, "dhts"))
  sum((probf - cons_realDummy(y))^2)/time_window
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
    fm <- apply(fi, 1, function(x){sum((0:(dim(fi)[2]-1)) * x)})
    fms[,i] <- fm
  }
  return(fms)
}

