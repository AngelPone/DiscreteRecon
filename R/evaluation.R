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

