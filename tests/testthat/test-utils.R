library(dplyr)

preparce_test_data1 <- function(){
  bts <- matrix(sample(0:1, 100, replace=TRUE), ncol=2)
  domain <- matrix(c(0, 1, 0, 1), 2)
  s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
  dhts(bts, s_mat, domain)
}

preparce_test_data2 <- function(x){
  upper <- rep(0, x)
  domain <- rbind(upper, upper + sample(1:2, x, replace = TRUE))
  s_mat <- twolevelhierarchy(x)
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=100, replace = TRUE)
  })
  dhts(bts, s_mat, domain)
}

preparce_basef <- function(dhts){
  ds <- apply(dhts$domain$domain_bts, 2, function(x){x[2] - x[1] + 1})
  tds <- rowSums(dhts$domain$domain_bts)[2] - rowSums(dhts$domain$domain_bts)[1] + 1
  probf <- c(tds, ds) %>%
    lapply(function(x){
      matrix(rnorm(100*x), 100) %>% abs() %>% 
        apply(1, function(x){x/sum(x)}) %>% t()
    })
  probf
}

test_that("prepare base forecasts", {
  dts <- preparce_test_data2(4)
  ds <- apply(dts$domain$domain_bts, 2, function(x){x[2] - x[1] + 1})
  tds <- rowSums(dts$domain$domain_bts)[2] - rowSums(dts$domain$domain_bts)[1] + 1
  ds <- as.numeric(c(tds, ds))
  basef <- preparce_basef(dts)
  expect_equal(length(basef), 5)
  for (i in 1:5){
    expect_equal(dim(basef[[i]]), as.numeric(c(100, ds[i])))
    expect_vector(rowSums(basef[[i]]), 1)
  }
})

test_that("convert dhts to dummies", {
  sample_dhts <- preparce_test_data2()
  domain <- sample_dhts$domain$coherent_domain
  
  a <- cons_realDummy(sample_dhts)
  for (i in 1:1:dim(a)[1]){
    expect_equal(domain[which.max(a[i,]), 2:dim(domain)[2]], 
                 as.numeric(sample_dhts$bts[i,]))
  }
})

test_that("calculate cost matrix", {
  sample_dhts <- preparce_test_data2(3)
  d1 <- sample_dhts$domain$incoherent_domain
  d2 <- sample_dhts$domain$coherent_domain
  distance <- cal_costeMatrix(d1, d2)
  expect_equal(dim(distance), c(dim(d2)[1], dim(d1)[1]))
  for (i in 1:dim(d2)[1]){
    for (j in 1:dim(d1)[1]){
      expect_equal(distance[i, j], sum(abs(d2[i,] - d1[j,])))
    }
  }
})

test_that("Marginal distributions to Joint distributions", {
  sample_dhts <- preparce_test_data2(3)
  basef <- preparce_basef(sample_dhts)
  r <- dim(sample_dhts$domain$coherent_domain)[1]
  q <- dim(sample_dhts$domain$incoherent_domain)[1]
  m <- dim(sample_dhts$domain$coherent_domain)[2]
  f1 <- marginal2Joint(basef, sample_dhts$domain$incoherent_domain)
  expect_is(f1, "jdist-ind")
  expect_equal(dim(f1), c(100, q))
  for (i in 1:q){
    a <- sample_dhts$domain$incoherent_domain[i,] %>% as.numeric()
    b <- NULL
    for (i in 1:length(basef)){
      b <- cbind(b, basef[[i]][, a[i]+1])
    }
    expect_vector(
      f1[,i],
      apply(b, 2, prod)
    )
  }
  
  f1 <- marginal2Joint(basef, sample_dhts$domain$coherent_domain)
  for (i in 1:r){
    a <- sample_dhts$domain$coherent_domain[i,2:m] %>% as.numeric()
    b <- NULL
    for (i in 2:length(basef)){
      b <- cbind(b, basef[[i]][, a[i]+1])
    }
    expect_vector(
      f1[,i],
      apply(b, 2, prod)
    )
  }
  expect_is(f1, "jdist-bu")
  expect_equal(dim(f1), c(100, r))
})



test.marginal2Sum <- function(){
  basef <- list(
    rbind(c(0.1, 0.9), c(0.1, 0.9), c(0.1, 0.9)),
    rbind(c(0.2, 0.8), c(0.2, 0.8), c(0.2, 0.8)),
    rbind(c(0.4, 0.6), c(0.3, 0.7), c(0.3, 0.7))
  )
  domain <- rbind(c(0, 0, 0), c(1, 1, 1))
  res <- marginal2Sum(basef, domain)
  expected <- rbind(c(0.008, 0.116, 0.444, 0.432),
                    c(0.006, 0.092, 0.398, 0.504),
                    c(0.006, 0.092, 0.398, 0.504))
  return(all.equal(res, expected, check.attributes=FALSE))
}

test.marginal2Joint <- function(){
  basef <- list(
    rbind(c(0.1, 0.8, 0.1), c(0.1, 0.8, 0.1), c(0.1, 0.8, 0.1)),
    rbind(c(0.2, 0.8), c(0.2, 0.8), c(0.2, 0.8)),
    rbind(c(0.4, 0.6), c(0.4, 0.6), c(0.4, 0.6))
  )
  res <- marginal2Joint(basef)
  expected <- rbind(c(0.008, 0.064, 0.008, 0.032, 0.256, 0.032, 0.012, 0.096, 0.012, 0.048, 0.384, 0.048),
                    c(0.008, 0.064, 0.008, 0.032, 0.256, 0.032, 0.012, 0.096, 0.012, 0.048, 0.384, 0.048),
                    c(0.008, 0.064, 0.008, 0.032, 0.256, 0.032, 0.012, 0.096, 0.012, 0.048, 0.384, 0.048))
  return(all.equal(res, expected, check.attributes=FALSE))
}

test.updateJoint <- function(){
  left <- t(matrix(c(0.02, 0.18, 0.08, 0.4, 0.3, 0.02)))
  right <- matrix(c(0.2, 0.1, 0.38, 0.32), nrow = 1)
  left_domain <- cons_domain(rbind(c(0,0), c(1,2)), twolevelhierarchy(2))
  right_domain <- cons_domain(rbind(c(0,0), c(1,1)), twolevelhierarchy(2))
  update_jointDist(left, right, left_domain, right_domain)
  
}

test.topdowntrain <- function(){
  y <- preparce_test_data1()
  library(dplyr)
  library(Matrix)
  # debug(topdown.train)
  A <- topdown.train(y)
  # undebug(topdown.train)
  print(A)
  basef <- matrix(rnorm(7 * 3), 7, 3) %>%
    apply(1, function(x){abs(x)/sum(abs(x))}) %>%
    t() %>% list()
  print(basef[[1]][1,])
  reconcile.topdown(A, basef)[1,]
}


time.topdowntrain <- function(){
  y <- test.dhts()
  library(dplyr)
  library(Matrix)
  start = Sys.time()
  for (i in 1:1000) topdown.train(y)
  end = Sys.time()
  print(end - start)
}


test.Joint2Marginal <- function(){
  library(dplyr)
  
  
  s_mat <- twolevelhierarchy(7)
  domain <- rbind(rep(0, 7), rep(1, 7))
  
  test_f <- c(8, 2, 2, 2, 2, 2, 2, 2) %>% 
    lapply(function(x){
      rnorm(100*x, x) %>% 
        matrix(100) %>%
        apply(1, function(x){abs(x)/sum(abs(x))}) %>%
        t()
    })
  test_bts <- sample(0:1, 7*100, replace = TRUE) %>% 
    matrix(100) %>%
    dhts(s_mat, domain)
  tmp <- marginal2Joint(test_f)
  tmp2 <- Joint2Marginal(tmp, test_bts$domain$incoherent_domain)
  stopifnot(all(tmp2[[1]] == test_f[[1]]))
  
  tmp <- marginal2Joint(test_f[2:8])
  tmp2 <- Joint2Marginal(tmp, test_bts$domain$coherent_domain)
  stopifnot(all(tmp2[[1]] == test_f[[1]]))
}
