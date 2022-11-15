library(dplyr)

prepare_test_data1 <- function(){
  bts <- matrix(sample(0:1, 200, replace=TRUE), ncol=2)
  domain <- matrix(c(0, 1, 0, 1), 2)
  s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
  dhts(bts, s_mat, domain)
}

prepare_test_data2 <- function(x, size=2000){
  upper <- rep(0, x)
  domain <- rbind(upper, upper + 2)
  s_mat <- rbind(rep(1,x), diag(x))
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=size, replace = TRUE)
  })
  dhts(bts, s_mat, domain)
}

prepare_test_data3 <- function(){
  domain <- rbind(rep(0, 4), rep(1, 4))
  s_mat <- rbind(c(1,1,1,1), c(1,1,0,0), c(0,0,1,1), diag(4))
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=100, replace = TRUE)
  })
  dhts(bts, s_mat, domain)
}

prepare_basef <- function(dhts){
  dsupper <- (dhts$meta$s_mat %*% t(dhts$meta$domain_bts))[,2]
  dslower <- (dhts$meta$s_mat %*% t(dhts$meta$domain_bts))[,1]
  ds <- dsupper - dslower + 1
  probf <- 1:length(ds) %>%
    lapply(function(x){
      a <- matrix(rnorm(NROW(dhts$bts)*ds[x]), NROW(dhts$bts)) %>% abs() %>% 
        apply(1, function(x){x/sum(x)}) %>% t()
      colnames(a) <- dslower[x] : dsupper[x]
      a
    })
  probf
}

prepare_recdist <- function(dhts){
  domain <- dhts$meta$coherent_domain
  r <- dim(domain)[1]
  structure(t(apply(matrix(rnorm(100*r), 100), 1, function(x){abs(x)/sum(abs(x))})),
            class = c("coherent", "jdist", "rec"))
}

