#' Brier Score
#' 
#' function to compute brier score given observations and probabilistic forecasts
#' @param probf joint probabilistic forecasts, T * n
#' @param y observations
#' @param hier dhier object
#' @export
brier_score <- function(probf, y, hier) {
  
  if (!is.list(probf)){
    time_window <- dim(probf)[1]
    bs <- sum((probf - cons_realDummy(y, hier))^2)/time_window
    probf <- Joint2Marginal(probf, hier)
  } else {
    probf1 <- marginal2Joint(probf, hier, method = "ind")
    time_window <- dim(probf1)[1]
    flags <- hier$coherent_flags
    cps <- probf1[,flags$coherent]
    ips <- probf1[,flags$incoherent]
    cps <- sum((cps - cons_realDummy(y, hier))^2)/time_window
    bs <- c(coherency=cps, incoherency=sum(ips^2)/time_window)
  }
  all_ts <- y %*% t(hier$s_mat)
  bs_level <- sapply(seq_along(probf), function(x){
    tmp_y <- all_ts[,x]
    sum(sapply(unique(tmp_y), function(g){
      idx <- which(tmp_y == g)
      col_idx <- which(colnames(probf[[x]]) == paste0(g))
      sum((1 - probf[[x]][idx, col_idx])^2) + sum(probf[[x]][idx, -col_idx]^2)
    })) / time_window})
  list(series=bs_level, hierarchy=bs)
}


#' Point forecast
#' 
#' function to generate point forecasts of the hierarchy
#' @param dist coherent distribution or basef
#' @param hier metadata
#' @return point forecasts of all series
#' @export 
point_forecast <- function(dist, hier){
  if (!is.list(dist)) series_dist <- Joint2Marginal(dist, hier)
  else series_dist <- dist
  
  sapply(series_dist, function(x){
    d <- as.integer(colnames(x))
    apply(x, 1, function(y){ sum(y * d) })
  }, simplify = "array")
}

#' Point metrics of each series
#' @param dist predictive distribution
#' @param y observation
#' @param f function to calculate metric
#' @return vector containing metric of each series
#' @export
point_metric <- function(dist, y, hier, f){
  pointf <- point_forecast(dist, hier)
  all_y <- y %*% t(hier$s_mat)
  sapply(1:dim(all_y)[2], function(x){
    f(pointf[,x], all_y[,x])
  })
}

