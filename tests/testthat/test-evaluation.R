test.brier_score <- function(){
  probf <- matrix(c(0.1, 0.2, 0.3, 0.9, 0.8, 0.7), 3)
  real <- matrix(c(1, 0, 0, 1, 1, 0), 3)
  brier_score(probf, real)
}
