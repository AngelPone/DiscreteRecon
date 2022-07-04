#' reconcile base probabilistic forecasts
#' 
#' Method for reconciling base probabilistic forecasts 
#' @import Matrix
#' @param basef joint base probabilistic forecasts assuming independence
#' @param real dhts object holding observed series
#' @param lambda penalty terms of cost, default is 1.
#' @param optimized indicating whether optimized. i.e., only moving probabilities to nearest coherent points.
#' @tag multiple-levels
#' @return r*q transformation matrix with class `rec_mat` 
allreconcile.train <- function(basef, 
                               real,
                               lambda = 1,
                               optimized=TRUE) {
  stopifnot(is(real, "dhts"))
  stopifnot(is(basef, "jdist-ind"))
  real_dummy <- cons_realDummy(real)
  stopifnot(dim(basef)[1] == dim(real_dummy)[1])
  
  cdomain <- real$domain$coherent_domain
  idomain <- real$domain$incoherent_domain
  
  r <- dim(cdomain)[1]
  q <- dim(idomain)[1]
  time_window <- dim(real_dummy)[1]
  costs <- cal_costeMatrix(idomain, cdomain)
  
  vec_c <- matrix(t(costs), nrow = 1)
  
  if (optimized){
    indicators <- getVIndexs(costs)
    Dmat <- 2 * construct_Q(r, basef, indicators)
    D <- construct_D(basef, indicators, real_dummy)
    E <- construct_E(costs, indicators)
    n_eq <- dim(E)[1]
    n_var <- dim(E)[2]
    Amat <- rbind(E, 
               diag(rep(1, n_var)),
               -diag(rep(1, n_var)))
    Amat <- t(Amat)
    bvec <- c(rep(1, n_eq),
              rep(0, n_var),
              rep(-1, n_var))
    dvec <- t(as.matrix(2 * D))
  } else {
    Dmat <- 2 * construct_Q(r, basef)  
    D <- 0
    for (i in 1:time_window) {
      D <- D + real_dummy[i,,drop=FALSE] %*% Matrix::bdiag(replicate(r, basef[i,,drop=FALSE], simplify = FALSE))
    }
    A1 <- diag(q)
    for (i in 1:(r - 1)) {
      A1 <- cbind(A1, diag(q))
    }
    b1 <- rep(1, q)
    
    A2 <- diag(r * q)
    b2 <- rep(0, r * q)
    
    A3 <- -diag(r * q)
    b3 <- rep(-1, r * q)
    Amat <- t(rbind(A1, A2, A3))
    bvec <- c(b1, b2, b3)
    dvec <- t(- lambda * matrix(costs, nrow = 1) + as.matrix(2/time_window * D)) * time_window
    n_eq <- q
  }

  solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = n_eq, factorized = FALSE)$solution
  solution[solution<1e-10] = 0
  if (optimized){
    Ahat <- construct_res(costs, solution, indicators)
  }else{
    Ahat <- t(matrix(solution, q, r))
  }
  structure(list(A=apply(Ahat, 2, function(x){x/sum(x)}),
                 domain=real$domain), class="rec_mat")
}


#' function to train reconciliation model
#' 
#' @param basef list containing base probabilistic forecasts of each series
#' @param real dhts object containing corresponding observations of basef
#' @param lambda penality
#' @param optimized boolean indicating if optimization 
#' @param step_wise boolean indicating if reconciliation step wisely.
#' @tag single-level
#' @return sw_res object or rec_mat object
#' @export
reconcile_train <- function(basef, real, 
                            lambda=1, optimized=TRUE, step_wise=TRUE){
  if (!step_wise) {
    basef <- marginal2Joint(basef, real$domain$incoherent_domain)
    return(allreconcile.train(basef, real, lambda=lambda, 
                              optimized = optimized))
  } else {
    return(stepwise_reconcile.train(basef, real, lambda = lambda, optimized = optimized))
  }
}

#' function to step-sisely reconcile base forecasts
#' 
#' @param basef basef forecasts
#' @param real dhts
#' @param lambda penalty
#' @param optimized if optimization
#' @return sw_res object holding reconciliation matrix of each step
#' @export
#' @tag single-level
stepwise_reconcile.train <- function(basef, real, 
                               lambda=1,
                               optimized=TRUE){
  stopifnot(is(real, "dhts"))
  n = NROW(real$s_mat)
  m = NCOL(real$s_mat)
  domain = real$domain$domain_bts 
  if (n - m > 1){
    stop("Only two-level hierarchy is supported!")
  }
  # list storing all A_hat
  As <- list()
  bf_total <- NULL
  
  for (i in 1:(m-1)){
    bf_left <- basef[[i+1]]
    if (i == m-1){
      rdomain <- domain[, (i+1):m, drop=FALSE]
      bf_right <- basef[[n]]
    } else {
      rdomain <- cons_domain(domain[, (i+1):m], rbind(rep(1, m-i), diag(m-i)), 
                             node_names = c(paste0('s', i+2, '-', m+1), colnames(real$domain$coherent_domain)[(i+2):(m+1)]))
      bf_right <- marginal2Sum(basef[(i+1):n], domain = rdomain)
    }
    if (is.null(bf_total)) bf_total <- basef[[i]]
    hierarchy_domain <- cbind(domain[,i], rowSums(domain[,(i+1):m, drop=FALSE]))
    # input for reconciliation of this step
    
    real_i <- dhts(cbind(real$bts[,i], rowSums(real$bts[,(i+1):m,drop=FALSE])), 
                   rbind(rep(1,2), diag(2)),
                   hierarchy_domain)
    basef_i <- marginal2Joint(list(bf_total, bf_left, bf_right), real_i$domain$incoherent_domain)
    modeli <- allreconcile.train(basef_i, real_i, lambda=lambda, optimized = optimized)
    Ahat <- list(A=modeli,
                 domain=real_i$domain,
                 rdomain=rdomain)
    As[[i]] <- Ahat
    bf_total <- Joint2Marginal(reconcile(modeli, list(bf_total, bf_left, bf_right)), real_i$domain$coherent_domain, 3)
  }
  structure(As, class = 'sw_res')
}

#' generic function to reconcile base forecasts
#'
#' @param x trained reconciliation model
#' @export
reconcile <- function(x, ...){
  UseMethod("reconcile", x)
}

#' function to warn irregular object
#' @param x output of reconcile.train
#' @export
reconcile.default <- function(x, ...){
  warning(paste("reconcile does not know how to handle object of class ", 
                class(x)))
}

#' step wise reconciliation 
#' 
#' @param x sw_res object
#' @param basef listing containing base probabilistic forecasts for each series
#' @export
reconcile.sw_res <- function(x, basef){
  stopifnot(is(x, 'sw_res'))
  n <- length(x)
  bf_total <- NULL
  ans <- list()
  for (i in 1:n){
    if (is.null(bf_total)) {
      bf_total <- basef[[i]]
    } else {
      bf_total <- Joint2Marginal(structure(reconciledi, class="jdist-rec"), subtreedomain$coherent_domain, 3) 
    }
    subtreedomain = x[[i]]$domain
    bf_left <- basef[[i+1]]
    if (i == n){
      bf_right <- basef[[length(basef)]]
    }else{
      bf_right <- marginal2Sum(basef[(i+1):length(basef)], x[[i]]$rdomain)
    }
    
    reconciledi <- reconcile(x[[i]]$A, list(bf_total, bf_left, bf_right))
    ans[[i]] <- list(dist=reconciledi, domain=subtreedomain)
  }
  structure(Step2Joint(ans)$dist, class="jdist-rec")
}


#' regular reconciliation 
#' 
#' @param x rec_mat object
#' @param basef listing containing base probabilistic forecasts for each series
#' @param domain domain of 
#' @return reconciled 
#' @export
reconcile.rec_mat <- function(x, basef){
  basef <- marginal2Joint(basef, x$domain$incoherent_domain)
  structure(t(x$A %*% t(basef)), class="jdist-rec")
}

# function to construct Q matrix used in optimization
# r: cardinality of coherent domain
# pi_hat: T*q base probabilistic forecasts
construct_Q <- function(r, basef, indicators=NULL) {
  if (rankMatrix(basef) != dim(basef)[2]){
    basef <- addProbJitter(basef)
  }
  Q <- t(basef) %*% basef
  if(is.null(indicators)){
    Matrix::bdiag(replicate(r, Q, simplify = FALSE))
  }else {
    Qs <- lapply(indicators, function(indicator){
      Q[indicator,][,indicator]
    })
    bdiag(Qs)
  }
}

#' function to train the naive top-down method based on concurrence.
#' 
#' @param real dhts object
#' @return topdown object 
#' @export
#' @import dplyr
#' @import Matrix
topdown.train <- function(real){
  stopifnot(is(real, "dhts"))
  coherent_domain <- data.frame(unclass(real$domain$coherent_domain))
  q <- dim(coherent_domain)[1]
  a <- dim(coherent_domain)[2]
  y <- real$bts
  concurrence <- vector("numeric", q)
  for (i in 1:dim(y)[1]){
    for (j in 1:q){
      if (all(y[i,] == coherent_domain[j, 2:a])){
        concurrence[j] = concurrence[j] + 1
        break
      }
    }
  }
  data.frame(total=coherent_domain$s1, concurrence) %>%
    group_by(total) %>%
    mutate(freq = sum(concurrence), prob = concurrence / freq) %>%
    pull(prob) %>%
    split(coherent_domain$s1) %>%
    bdiag() %>% as.matrix() %>%
    structure(class = "topdown")
}

#' topdown reconcilation function
#' 
#' @param x topdown object
#' @param basef base forecasts of total level
#' @export
reconcile.topdown <- function(x, basef){
  total_f <- basef[[1]]
  structure(t(x %*% t(total_f)), class='jdist-td')
}

#' function to stepwise update the joint distribution and joint domain.
#'
#' @param ans list containing joint distribution of all sub trees.
#' @return reconciled joint distribution  
Step2Joint <- function(ans){
  library(dplyr)
  # function to update domain
  update_domain <- function(l, r){
    l <- data.frame(unclass(l))
    r <- data.frame(unclass(r))
    colnames(l)[dim(l)[2]] <- "key"
    colnames(r)[1] <- "key"
    ans <- right_join(r, l, by="key")
    ans <- ans %>% select(-c("key")) %>% 
      select(c(3:dim(.)[2], 1:2))
    colnames(ans) <- 1:dim(ans)[2]
    structure(as.matrix(ans), class="coherent_domain")
  }
  
  # function to update joint distribution
  update_jointDist <- function(ldist, rdist, ld, rd){
    ld <- data.frame(unclass(ld))
    rd <- data.frame(unclass(rd))
    colnames(ld)[dim(ld)[2]] <- "key"
    colnames(rd)[1] <- "key"
    time_window <- dim(ldist)[1]
    new_dist <- NULL
    for (i in 1:time_window){
      ld$lprob <- ldist[i,]
      new_dist <- mutate(rd, prob=rdist[i,]) %>% group_by(key) %>%
        mutate(ratio=prob/sum(prob)) %>%
        right_join(ld, by="key") %>%
        mutate(prob = lprob * ratio) %>%
        pull(prob) %>%
        rbind(new_dist, .)
    }
    structure(as.matrix(new_dist), class='jdist-rec')
  }
  
  joint_dist <- ans[[1]]$dist
  joint_domain <- ans[[1]]$domain$coherent_domain
  for (i in 1:(length(ans) - 1)){
    # fl: reconciled forecasts in the left hierarchy (as right node)
    fl <- Joint2Marginal(joint_dist, joint_domain, dim(joint_domain)[2])
    # fr: reconciled forecasts in the right hierarchy (as total)
    rd <- ans[[i+1]]$domain$coherent_domain
    fr <- Joint2Marginal(ans[[i+1]]$dist, rd, 1)
    adj_dist <- (fl + fr)/2
    joint_l <- coherentadj(joint_dist, adj_dist, joint_domain, dim(joint_domain)[2])
    joint_r <- coherentadj(ans[[i+1]]$dist, adj_dist, rd, 1)
    joint_dist <- update_jointDist(joint_dist, joint_r, joint_domain, rd)
    joint_domain <- update_domain(joint_domain, rd)
  }
  rownames(joint_dist) <- NULL
  colnames(joint_dist) <- 1:dim(joint_dist)[2]
  list(domain=joint_domain, dist=joint_dist)
}

#' function to adjust joint distribution given new marginal distribution
#' 
#' @param x original joint distribution
#' @param target new marginal distribution
#' @param domain domain of the hierarchy
#' @param which integer indicating which column in the domain 
#' the marginal distribution refers to.
#' @return adjusted joint distribution
#' @import dplyr
coherentadj <- function(x, target, domain, which){
  
  time_window <- NROW(x)
  new_x <- NULL
  for (i in 1:time_window){
    new_x <- data.frame(domain=domain[,which], prob=x[i,]) %>%
      group_by(domain) %>% 
      mutate(prob = prob * target[i, as.character(pull(cur_group(), domain))] / sum(prob) ) %>%
      pull(prob) %>%
      rbind(new_x, .)
  }
  new_x
}


addProbJitter <- function(x){
  add_x <- matrix(runif(length(x), -min(abs(x))/100, min(abs(x))/100), NROW(x), NCOL(x))
  t(apply((x + add_x), 1, function(x){abs(x)/sum(abs(x))}))
}



#' utility functions
#' @tag multiple-levels
getVIndexs <- function(distance){
  r <- dim(distance)[1]
  q <- dim(distance)[2]
  Avec <- matrix(0, r, q)
  indicators <- list()
  for (i in 1:q){
    if (min(distance[,i]) > 0){
      Avec[which(distance[,i]==min(distance[,i])), i] = 1
    }
  }
  for (i in 1:r){
    indicators[[i]] <- which(Avec[i, ] == 1)
  }
  indicators
}

#' utility functions
construct_D <- function(jdist, indicators, y){
  D_mat <- 0
  time_window <- dim(y)[1]
  for (time in 1:time_window){
    Ds <- lapply(indicators, function(indicator){
      matrix(jdist[time, indicator], nrow=1)
    })
    D_mat <- D_mat + y[time,,drop=FALSE] %*% bdiag(Ds)
  }
  D_mat
}

#' utility functions
construct_E <- function(distance, indicators){
  n <- sum(sapply(indicators, function(x){length(x)}))
  E <- NULL
  for (j in 1:dim(distance)[2]){
    if (min(distance[, j]) == 0) next;
    E_jRow <- matrix(0, ncol=n)
    currentVar = 0
    for (i in 1:length(indicators)) {
      E_jRow[which(indicators[[i]] == j) + currentVar] = 1
      currentVar = currentVar + length(indicators[[i]])
    }
    E <- rbind(E, E_jRow)
  }
  E
}

#' utility functions
construct_res <- function(distance, solution, indicators){
  res <- replace(distance, distance==0, 1)
  res <- replace(res, distance>0, 0)
  currentVar = 1
  for (i in 1:length(indicators)){
    res[i, indicators[[i]]] = solution[currentVar: (currentVar + length(indicators[[i]]) - 1)]
    currentVar = currentVar + length(indicators[[i]])
  }
  res
}

