library(dplyr)

prepare_test_data1 <- function(){
  bts <- matrix(sample(0:1, 200, replace=TRUE), ncol=2)
  domain <- matrix(c(0, 1, 0, 1), 2)
  s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
  list(bts = bts, hier = dhier(s_mat, domain))
}

prepare_test_data2 <- function(x, size=2000){
  upper <- rep(0, x)
  domain <- rbind(upper, upper + 2)
  s_mat <- rbind(rep(1,x), diag(x))
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=size, replace = TRUE)
  })
  list(bts = bts, hier = dhier(s_mat, domain))
}

prepare_test_data3 <- function(){
  domain <- rbind(rep(0, 4), rep(1, 4))
  s_mat <- rbind(c(1,1,1,1), c(1,1,0,0), c(0,0,1,1), diag(4))
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=100, replace = TRUE)
  })
  list(bts=bts, hier=dhier(s_mat, domain))
}

prepare_basef <- function(dhts){
  dsupper <- (dhts$hier$s_mat %*% t(dhts$hier$domain))[,2]
  dslower <- (dhts$hier$s_mat %*% t(dhts$hier$domain))[,1]
  ds <- dsupper - dslower + 1
  probf <-
    lapply(1:length(ds), function(x){
      a <- abs(matrix(rnorm(NROW(dhts$bts)*ds[x]), NROW(dhts$bts)))
      a <- apply(a, 1, function(x){x/sum(x)})
      a <- t(a)
      colnames(a) <- dslower[x] : dsupper[x]
      a
    })
  probf
}

prepare_recdist <- function(dhts){
  domain <- dhts$hier$coherent_domain
  r <- dim(domain)[1]
  t(apply(matrix(rnorm(100*r), 100), 1, function(x){abs(x)/sum(abs(x))}))
}

