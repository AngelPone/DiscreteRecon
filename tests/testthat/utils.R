library(dplyr)

prepare_test_data1 <- function(){
  bts <- matrix(sample(0:1, 200, replace=TRUE), ncol=2)
  domain <- matrix(c(0, 1, 0, 1), 2)
  s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
  dhts(bts, s_mat, domain)
}

prepare_test_data2 <- function(x){
  upper <- rep(0, x)
  domain <- rbind(upper, upper + sample(1:2, x, replace = TRUE))
  s_mat <- rbind(rep(1,x), diag(x))
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=100, replace = TRUE)
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
  dsupper <- (dhts$s_mat %*% t(dhts$domain$domain_bts))[,2]
  dslower <- (dhts$s_mat %*% t(dhts$domain$domain_bts))[,1]
  ds <- dsupper - dslower + 1
  probf <- 1:length(ds) %>%
    lapply(function(x){
      a <- matrix(rnorm(100*ds[x]), 100) %>% abs() %>% 
        apply(1, function(x){x/sum(x)}) %>% t()
      colnames(a) <- dslower[x] : dsupper[x]
      a
    })
  probf
}

prepare_recdist <- function(dhts){
  domain <- dhts$domain$coherent_domain
  r <- dim(domain)[1]
  structure(t(apply(matrix(rnorm(100*r), 100), 1, function(x){abs(x)/sum(abs(x))})),
            class = "jdist-rec")
}

test_that("prepare base forecasts", {
  dts <- prepare_test_data2(4)
  dsupper <- (dts$s_mat %*% t(dts$domain$domain_bts))[,2]
  dslower <- (dts$s_mat %*% t(dts$domain$domain_bts))[,1]
  ds <- as.numeric(dsupper - dslower + 1)
  basef <- prepare_basef(dts)
  expect_equal(length(basef), 5)
  for (i in 1:5){
    expect_equal(dim(basef[[i]]), as.numeric(c(100, ds[i])))
    expect_vector(rowSums(basef[[i]]), 1)
    expect_equal(as.numeric(colnames(basef[[i]])), dslower[i]:dsupper[i])
  }
})

test_that("prepare base forecasts for multiple levels", {
  dts <- prepare_test_data3()
  ds <- (dts$s_mat %*% t(dts$domain$domain_bts))[,2] - (dts$s_mat %*% t(dts$domain$domain_bts))[,1] + 1
  basef <- prepare_basef(dts)
  expect_equal(length(basef), 7)
  for (i in 1:7){
    expect_equal(dim(basef[[i]]), as.numeric(c(100, ds[i])))
    expect_vector(rowSums(basef[[i]]), 1)
  }
})