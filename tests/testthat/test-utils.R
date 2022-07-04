source('./utils.R')

test_that("convert dhts to dummies", {
  sample_dhts <- prepare_test_data3()
  domain <- sample_dhts$domain$coherent_domain
  m <- dim(sample_dhts$domain$domain_bts)[2]
  n <- dim(domain)[2]
  a <- cons_realDummy(sample_dhts)
  for (i in 1:1:dim(a)[1]){
    expect_equal(as.numeric(domain[which.max(a[i,]), (n-m+1):n]), 
                 as.numeric(sample_dhts$bts[i,]))
  }
})

test_that("calculate cost matrix", {
  sample_dhts <- prepare_test_data3()
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
  sample_dhts <- prepare_test_data3()
  basef <- prepare_basef(sample_dhts)
  r <- dim(sample_dhts$domain$coherent_domain)[1]
  q <- dim(sample_dhts$domain$incoherent_domain)[1]
  n <- dim(sample_dhts$domain$coherent_domain)[2]
  m <- dim(sample_dhts$domain$domain_bts)[2]
  f1 <- marginal2Joint(basef, sample_dhts$domain$incoherent_domain)
  expect_is(f1, "jdist-ind")
  expect_equal(dim(f1), c(100, q))
  for (i in 1:q){
    a <- sample_dhts$domain$incoherent_domain[i,] %>% as.numeric()
    b <- NULL
    for (j in 1:length(basef)){
      b <- cbind(b, basef[[j]][, as.character(a[j])])
    }
    expect_equal(
      f1[,i],
      apply(b, 1, prod)
    )
  }
  
  f1 <- marginal2Joint(basef, sample_dhts$domain$coherent_domain)
  for (i in 1:r){
    a <- sample_dhts$domain$coherent_domain[i,(n-m+1):n] %>% as.numeric()
    b <- NULL
    for (j in (n-m+1):(n)){
      b <- cbind(b, basef[[j]][, as.character(a[j - (n-m)])])
    }
    expect_equal(
      f1[,i],
      apply(b, 1, prod)
    )
  }
  expect_is(f1, "jdist-bu")
  expect_equal(dim(f1), c(100, r))
})

test_that("Marginal distributions to sum assuming independence", {
  dts <- prepare_test_data3()
  basef <- prepare_basef(dts)
  domain <- dts$domain$coherent_domain
  res <- marginal2Sum(basef, domain)
  
  expect_equal(length(res), 3)
  expect_equal(dim(res[[2]]), c(100, 3))
  expect_equal(dim(res[[3]]), c(100, 3))
  expect_equal(dim(res[[1]]), c(100, 5))
  
  expect_equal(res[[1]][,"0"], apply(
    do.call(cbind, lapply(basef[4:7], function(x){x[,"0"]})), 1, prod
  ))
  expect_equal(res[[1]][,"4"], apply(
    do.call(cbind, lapply(basef[4:7], function(x){x[,"1"]})), 1, prod
  ))
  
  res <- marginal2Sum(basef, domain, 1)
  expect_equal(dim(res), c(100, 5))
  
  res <- marginal2Sum(basef, domain, 2)
  expect_equal(dim(res), c(100, 3))
  
  res <- marginal2Sum(basef, domain, 3)
  expect_equal(dim(res), c(100, 3))
  expect_equal(res[,"2"], apply(
    do.call(cbind, lapply(basef[6:7], function(x){x[,"1"]})), 1, prod
  ))
  expect_equal(res[,"0"], apply(
    do.call(cbind, lapply(basef[6:7], function(x){x[,"0"]})), 1, prod
  ))
  expect_equal(res[,"1"], 
                basef[[6]][,"0"] * basef[[7]][,"1"] + 
                basef[[6]][,"1"] * basef[[7]][,"0"])
})


test_that("Joint Distributions to Marginal distribution of specific dimension for bottom up distribution", {
  dts <- prepare_test_data3()
  basef <- prepare_basef(dts)
  domain <- dts$domain$coherent_domain
  
  jdist <- marginal2Joint(basef, domain)
  mdist <- Joint2Marginal(jdist, domain)
  for (i in 1:length(mdist)){
    if (i <= 3){
      expect(is.null(mdist[[i]]), "shoule be null")
      next
    }
    expect_equal(mdist[[i]], basef[[i]])
    expect_equal(colnames(mdist[[i]]), colnames(basef[[i]]))
  }
  
  
})

test_that("Joint Distributions to Marginal distribution of specific dimension for reconciled distribution", {
  dts <- prepare_test_data3()
  recdist <- prepare_recdist(dts)
  domain <- dts$domain$coherent_domain
  mdist <- Joint2Marginal(recdist, domain)
  
  domain <- dts$domain$coherent_domain
  for(i in seq_along(mdist)){
    cdist <- mdist[[i]]
    expect(!is.null(cdist), "shoule not be null")
    
    for (col in colnames(cdist)){
      expect_equal(rowSums(recdist[,which(domain[,i] == col),drop=FALSE]),
                    cdist[,col])
    }
  }
})


