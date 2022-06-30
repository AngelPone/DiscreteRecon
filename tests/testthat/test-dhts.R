# test dhts object and domain

bts <- matrix(sample(0:1, 100, replace=TRUE), ncol=2)
domain <- matrix(c(0, 1, 0, 1), 2)
s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
a <- dhts(bts, s_mat, domain)
t1 <- a$domain$coherent_domain
t2 <- a$domain$incoherent_domain

test_that("domain", {
  expect_is(t1, "coherent_domain")
  expect_is(t2, "incoherent_domain")
  expect_equal(dim(t1), c(4,3))
  expect_equal(dim(t2), c(12, 3))
  expect_equal(t1[,1], t1[,2] + t1[,3])
})

bts <- matrix(sample(0:1, 700, replace=TRUE), ncol=7)
domain <- rbind(rep(0, 7), rep(1, 7))
s_mat <- twolevelhierarchy(7)
a <- dhts(bts, s_mat, domain)
t1 <- a$domain$coherent_domain
t2 <- a$domain$incoherent_domain

test_that("domain", {
  expect_is(t1, "coherent_domain")
  expect_is(t2, "incoherent_domain")
  expect_equal(dim(t1), c(128, 8))
  expect_equal(dim(t2), c(128*8, 8))
  expect_equal(t1[,1], rowSums(t1[,2:8]))
})


# different domain for bottom level time series
upper <- sample(3:8, 7, replace = TRUE)
domain <- rbind(upper, upper + sample(1:2, 7, replace = TRUE))
s_mat <- twolevelhierarchy(7)
bts <- apply(domain, 2, function(x){
  sample(x[1]:x[2], size=100, replace = TRUE)
})
a <- dhts(bts, s_mat, domain)
t1 <- a$domain$coherent_domain
t2 <- a$domain$incoherent_domain

dslower <- apply(domain, 2, function(x){(x[2] - x[1] + 1)})
dsupper <- (rowSums(domain)[2] - rowSums(domain)[1]) + 1

test_that("domain", {
  expect_is(t1, "coherent_domain")
  expect_is(t2, "incoherent_domain")
  expect_equal(dim(t1), c(prod(dslower), 8))
  expect_equal(dim(t2), as.numeric(c(prod(dslower) * dsupper, 8)))
  expect_equal(t1[,1], rowSums(t1[,2:8]))
})


