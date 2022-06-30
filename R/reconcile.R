#' reconcile base probabilistic forecasts
#' 
#' Method for reconciling base probabilistic forecasts 
#' @param basef joint base probabilistic forecasts
#' @param real dhts object holding observed series
#' @param lambda penalty terms of cost, default is 1.
#' @param optimized 
#' @return r*q transformation matrix
allreconcile.train <- function(basef, 
                               real,
                               lambda = 1,
                               optimized=TRUE) {
  
  # basef: T * q
  library(Matrix)
  stopifnot(is(real, "dhts"))
  real_dummy <- cons_realDummy(real)
  stopifnot(dim(basef)[1] == dim(real_dummy)[1])
  
  getVIndexs <- function(distance){
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
  construct_D <- function(pi_hat, indicators){
    D_mat <- 0
    for (time in 1:time_window){
      Ds <- lapply(indicators, function(indicator){
        matrix(pi_hat[time, indicator], nrow=1)
      })
      D_mat <- D_mat + real_dummy[time,,drop=FALSE] %*% bdiag(Ds)
    }
    D_mat
  }
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
  
  r <- dim(real_dummy)[2]
  q <- dim(basef)[2]
  time_window <- dim(real_dummy)[1]
  costs <- cal_costeMatrix(real$domain$incoherent_domain, 
                           real$domain$coherent_domain)
  vec_c <- matrix(t(costs), nrow = 1)
  if (optimized){
    indicators <- getVIndexs(costs)
    Dmat <- 2 * construct_Q(r, basef, indicators)
    D <- construct_D(basef, indicators)
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
  # check the result in the optimization function equals to the result of matrix form.
  test.equalres <- function(){
    getVIndexs2 <- function(distance){
      Avec <- matrix(0, r, q)
      indicators <- list()
      for (i in 1:q){
        Avec[which(distance[,i]==min(distance[,i])), i] = 1
      }
      for (i in 1:r){
        indicators[[i]] <- which(Avec[i, ] == 1)
      }
      indicators
    }
    indicators2 <- getVIndexs2(costs)
    Dmat2 <-2 * construct_Q(r, basef, indicators2)
    dvec2 <- t(as.matrix(2 * construct_D(basef, indicators2)))
    tmp_sol <- solution + 2
    Amat <- construct_res(costs, tmp_sol, indicators)
    vec_A <- matrix(t(Amat), nrow=1)[1,]
    vec_A <- vec_A[vec_A > 0]
    vec_A[vec_A > 1.5] = vec_A[vec_A>1.5] - 2
    total_opt = (-t(dvec2) %*% vec_A + 1/2 * t(vec_A) %*% Dmat2 %*% vec_A + time_window)/time_window
    Ahat <- construct_res(costs, solution, indicators)
    res <- sum((t(Ahat %*% t(basef)) - real_dummy)^2)/time_window
    cat("Optimized matrix form:", as.numeric(total_opt), "\n")
    cat("scalar form:", as.numeric(res), "\n")
  }
  # test.equalres()
  solution[solution<1e-10] = 0
  if (optimized){
    Ahat <- construct_res(costs, solution, indicators)
  }else{
    Ahat <- t(matrix(solution, q, r))
  }
  # test if best
  Ahat2 <- construct_res(costs, rnorm(length(solution)), indicators)
  Ahat2 <- apply(Ahat2, 2, function(x){abs(x)/sum(abs(x))})
  sc1 <- brier_score(t(Ahat2 %*% t(basef)), real)
  sc2 <- brier_score(t(Ahat %*% t(basef)), real)
  cat(length(solution), "," , sc1, ",", sc2, ",", sep = "")
  structure(apply(Ahat, 2, function(x){x/sum(x)}), class="rec_mat")
}


#' function to train reconciliation model
#' 
#' @param basef list containing base probabilistic forecasts of each series
#' @param real dhts object containing corresponding observations of basef
#' @param lambda penality
#' @param optimized boolean indicating if optimization 
#' @param step_wise boolean indicating if reconciliation step wisely.
#' @return sw_res object or rec_mat object
#' @export
reconcile_train <- function(basef, real, 
                            lambda=1, optimized=TRUE, step_wise=TRUE){
  if (!step_wise) {
    basef <- marginal2Joint(basef)
    return(allreconcile.train(basef, real, lambda=lambda, 
                              optimized = optimized))
  } else {
    return(stepwise_reconcile.train(basef, real, lambda = lambda, optimized = optimized))
  }
}

#' function to step-sisely reconcile base forecasts
#' 
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
    bf_right <- marginal2Sum(basef[(i+2):n], domain = domain[, (i+1):m])
    if (is.null(bf_total)) bf_total <- basef[[i]]
    hierarchy_domain <- cbind(domain[,i], rowSums(domain[,(i+1):m, drop=FALSE]))
    # input for reconciliation of this step
    basef_i <- marginal2Joint(list(bf_total, bf_left, bf_right))
    real_i <- dhts(cbind(real$bts[,i], rowSums(real$bts[,(i+1):m,drop=FALSE])), 
                   twolevelhierarchy(2),
                   hierarchy_domain)
    Ahat <- list(A=allreconcile.train(basef_i, real_i, lambda=lambda, optimized = optimized),
                 leftdomain=hierarchy_domain[,1],
                 rightdomain=domain[,(i+1):m])
    cat(brier_score(marginal2Joint(list(bf_left, bf_right)), real_i), '\n', sep = "")
    As[[i]] <- Ahat
    bf_total <- Joint2Marginal(t(Ahat$A %*% t(basef_i)), real_i$domain$coherent_domain, 3)
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
  n <- length(x)
  bf_total <- NULL
  ans <- list()
  for (i in 1:n){
    if (is.null(bf_total)) {
      bf_total <- basef[[i]]
    } else {
      bf_total <- Joint2Marginal(reconciledi, subtreedomain, 3) 
    }
    if(!is.null(dim(x[[i]]$rightdomain))){
      right_domain <- x[[i]]$rightdomain
    } else {
      right_domain <- matrix(x[[i]]$rightdomain, ncol=1)
    }
    subtreedomain = cons_domain(cbind(x[[i]]$leftdomain, rowSums(right_domain)),
                                twolevelhierarchy(2))
    bf_left <- basef[[i+1]]
    bf_right <- marginal2Sum(basef[(i+2):length(basef)], x[[i]]$rightdomain)
    reconciledi <- t(x[[i]]$A %*% t(marginal2Joint(list(bf_total, bf_left, bf_right))))
    ans[[i]] <- list(dist=reconciledi, domain=subtreedomain)
  }
  Step2Joint(ans)$dist
}


#' regular reconciliation 
#' 
#' @param x rec_mat object
#' @param basef listing containing base probabilistic forecasts for each series
#' @export
reconcile.rec_mat <- function(x, basef){
  
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
  coherent_domain <- data.frame(real$domain$coherent_domain)
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
  # tmp <- data.frame(y) %>% group_by(across()) %>% count()
  # for (i in 1:q){
  #   for (j in 1:q){
  #     if (all(tmp[i, 1:(a-1)] == coherent_domain[j, 2:a])){
  #       concurrence[i] <- tmp$n[i]
  #     }
  #   }
  # }
  data.frame(total=coherent_domain$X1, concurrence) %>%
    group_by(total) %>%
    mutate(freq = sum(concurrence), prob = concurrence / freq) %>%
    pull(prob) %>%
    split(coherent_domain$X1) %>%
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
  t(x %*% t(total_f))
}

#' function to stepwise update the joint distribution and joint domain.
#'
#' @param ans list containing joint distribution of all sub trees. 
Step2Joint <- function(ans){
  library(dplyr)
  # function to update domain
  update_domain <- function(l, r){
    l <- data.frame(l)
    r <- data.frame(r)
    colnames(l)[dim(l)[2]] <- "key"
    colnames(r)[1] <- "key"
    ans <- right_join(r, l, by="key")
    ans <- ans %>% select(-c("key")) %>% 
      select(c(3:dim(.)[2], 1:2))
    colnames(ans) <- 1:dim(ans)[2]
    ans
  }
  
  # function to update joint distribution
  update_jointDist <- function(ldist, rdist, ld, rd){
    ld <- data.frame(ld)
    rd <- data.frame(rd)
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
    new_dist
  }
  
  joint_dist <- ans[[1]]$dist
  joint_domain <- ans[[1]]$domain
  for (i in 1:(length(ans) - 1)){
    # fl: reconciled forecasts in the left hierarchy (as right node)
    fl <- Joint2Marginal(joint_dist, joint_domain, dim(joint_domain)[2])
    # fr: reconciled forecasts in the right hierarchy (as total)
    fr <- Joint2Marginal(ans[[i+1]]$dist, ans[[i+1]]$domain, 1)
    adj_dist <- (fl + fr)/2
    joint_l <- coherentadj(joint_dist, adj_dist, joint_domain, dim(joint_domain)[2])
    joint_r <- coherentadj(ans[[i+1]]$dist, adj_dist, ans[[i+1]]$domain, 1)
    joint_dist <- update_jointDist(joint_dist, joint_r, joint_domain, ans[[i+1]]$domain)
    joint_domain <- update_domain(joint_domain, ans[[i+1]]$domain)
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

