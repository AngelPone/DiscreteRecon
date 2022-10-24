#' Brier Score
#' 
#' function to compute brier score given observations and probabilistic forecasts
#' @param probf joint probabilistic forecasts, T * n
#' @param y dhts
#' @export
brier_score <- function(probf, y) {
  
  stopifnot(is(y, "dhts"))
  if (is_coherentJdist(probf)){
    time_window <- dim(probf)[1]
    bs <- sum((probf - cons_realDummy(y))^2)/time_window
  } else if (is.list(probf)){
    probf <- marginal2Joint(probf, y$meta, method = "ind")
    time_window <- dim(probf)[1]
    flags <- y$meta$coherent_flags
    cps <- probf[,flags$coherent]
    ips <- probf[,flags$incoherent]
    cps <- sum((cps - cons_realDummy(y))^2)/time_window
    bs <- c(coherency=cps, incoherency=sum(ips^2)/time_window)
  }
  probf <- Joint2Marginal(probf, y$meta)
  all_ts <- y$bts %*% t(y$meta$s_mat)
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
#' @param meta metadata
#' @return point forecasts of all series
#' @export 
point_forecast <- function(dist, meta){
  
  if (is_coherentJdist(dist)) series_dist <- Joint2Marginal(dist, meta)
  else if (is.list(dist)) series_dist <- dist
  
  sapply(series_dist, function(x){
    d <- as.integer(colnames(x))
    apply(x, 1, function(y){ sum(y * d) })
  }, simplify = "array")
  # domain <- meta$coherent_domain
  # fms <- matrix(0, nrow = dim(dist)[1], ncol = dim(domain)[2])
  # for (i in 1:dim(domain)[2]){
  #   fi <- Joint2Marginal(dist, domain, i)
  #   ds <- as.numeric(colnames(fi))
  #   fm <- apply(fi, 1, function(x){sum(ds * x)})
  #   fms[,i] <- fm
  # }
  # colnames(fms) <- colnames(domain)
  # return(fms)
}

#' Point metrics of each series
#' @param dist
#' @param y dhts
#' @param f function to calculate metric
#' @return vector containing metric of each series
#' @export
point_metric <- function(dist, y, f){
  pointf <- point_forecast(dist, y$meta)
  all_y <- aggdhts(y)
  sapply(1:dim(all_y)[2], function(x){
    f(pointf[,x], all_y[,x])
  })
}

