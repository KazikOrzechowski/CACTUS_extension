logSumExp <- matrixStats::logSumExp

try_split <- function(){
  #pick a cluster to split proportional to its size. Can't split a cluster < 2
  p_clusters <- sapply(clusters_, function(cluster){
                                      if(length(cluster)>1){length(cluster)/M}else{0}
                                      })
  if(!any(p_clusters>0)){
    #print('No clusters of 2 or more!')
    return()
  }
  chosen_ind <- sample(length(clusters_), 1, prob=p_clusters)
  chosen_clust <- clusters_[[chosen_ind]]
  
  #chose two cells i and j from the cluster uniformly
  chosen_cells <- sample(chosen_clust, 2)
  i <- chosen_cells[1]
  j <- chosen_cells[2]
  
  #environment setting
  environment(accept_MH_split)<- environment()
  environment(launch_split)<- environment()
  environment(launch_merge)<- environment()
  environment(propose_split)<- environment()
  
  #create the launch split and launch merge states
  #launch split
  split_ <- launch_split(i, j, chosen_clust)
  
  if(split_$inf){
    return()
  }
  
  t_split <- split_$t
  B_split <- split_$B
  I_split <- split_$I
  n_split <- split_$n
  
  #launch merge
  merge_ <- launch_merge(chosen_clust)
  B_merge <- merge_$B
  I_merge <- merge_$I
  
  #define proposed state
  propose <- propose_split(chosen_clust, i, j,
                            t_split, n_split, B_split, 
                            I_split)
  
  if(propose$inf){
    return()
  }
  
  t_propose <- propose$t
  B_propose <- propose$B
  I_propose <- propose$I
  log_p <- propose$log_p
  n_proposed <- propose$n
  
  #calc M-H acceptance prob
  
  accept <- accept_MH_split(t_propose,
                            chosen_clust, 
                            n_proposed, 
                            log_p)
  
  #if accepted, assign first empty cluster to proposed cluster with cell i
  if(accept){
    #print('Split accepted')
    first_empty <- which(sapply(clusters_, length)==0, arr.ind=TRUE)[1]
    cluster_i <- chosen_clust[t_propose$clust == 0]
    cluster_j <- chosen_clust[t_propose$clust == 1]
    t_[cluster_i] <<- ifelse(is.na(first_empty), length(clusters_)+1, first_empty)
    B[t_[i],,] <<- B_propose[1,,]
    B[t_[j],,] <<- B_propose[2,,]
    I[t_[i]] <<- I_propose[1]
    I[t_[j]] <<- I_propose[2]
    clusters_[[t_[j]]] <<- cluster_j
    if(is.na(first_empty)){
      clusters_[[length(clusters_)+1]] <<- cluster_i
    }else{
      clusters_[[first_empty]] <<- cluster_i
    }
  }
  
  return()
}


try_merge <- function(){
  #pick two clusters to merge inversely proportional to their size
  p_clusters <- sapply(clusters_, function(cluster){
    ifelse(length(cluster)>0, M/length(cluster), 0)
    })
  
  #if only one cluster exists it has p=1 and we cannot merge
  if(any(p_clusters==1)){
    #print('Only one clust!')
    return()
  }
  chosen <- sample(length(clusters_), 2, prob=p_clusters)
  
  cluster_i <- clusters_[[chosen[1]]]
  cluster_j <- clusters_[[chosen[2]]]
  
  #draw a cell from each chosen cluster uniformly
  i <- resample(cluster_i)
  j <- resample(cluster_j)
  
  chosen_cells <- union(cluster_i, cluster_j)
  
  #environment setting
  environment(accept_MH_merge)<- environment()
  environment(launch_split)<- environment()
  environment(launch_merge)<- environment()
  
  #create the launch split and launch merge states
  #launch split
  split_ <- launch_split(i, j, chosen_cells)
  
  if(split_$inf){
    return()
  }
  
  t_split <- split_$t
  B_split <- split_$B
  I_split <- split_$I
  n_split <- split_$n
  
  #we do not need to launch merge in this case, as the probabilities of parameters in proposed state
  #are independent from the launch state
  proposed <- launch_merge(chosen_cells)
  B_proposed <- proposed$B
  I_proposed <- proposed$I
  
  accept <- accept_MH_merge(t_split, i, j, 
                            cluster_i, cluster_j,
                            chosen_cells,
                            n_split, B_split, I_split)
  
  if(accept){
    #print('Merge accepted')
    B[t_[j],,] <<- B_proposed
    I[t_[j]] <<- I_proposed
    clusters_[[t_[j]]] <<- chosen_cells
    clusters_[[t_[i]]] <<- list()
    t_[cluster_i] <<- t_[j]
  }
  
  return()
}


launch_split <- function(i, j, chosen_clust){
  #randomly assign each cell from cluster to i or j (these stay the same)
  t_split <- data.frame(clust=rbinom(chosen_clust, 1, .5),
                        row.names = str(chosen_clust))
  t_split[str(i),] <- 0
  t_split[str(j),] <- 1
  
  #n_split are the populations of split clusters
  n_split <- 1:2
  n_split[2] <- sum(t_split$clust)
  n_split[1] <- length(chosen_clust) - n_split[2]
  
  #we will try drawing parameters from conditional distributions, not priors
  upi <- as.matrix(g) + BCR[i,,]
  upj <- as.matrix(g) + BCR[j,,]
  
  #B_split are the nuc freqs of split clusters
  B_split <- array(0, dim=c(2, L, 4))
  B_split[1,,] <- t(sapply(1:L, function(l){igraph::sample_dirichlet(1, upi[l,])}))
  B_split[2,,] <- t(sapply(1:L, function(l){igraph::sample_dirichlet(1, upj[l,])}))
  
  #I_split are the cluster-clone assignments of split clusters
  I_split <- 1:2
  
  for(i_ in 1:2){
    log_like <- sapply(1:K, function(k){
      loglike_RNA(ifelse(i_==1, i, j), k, C, theta1, theta0, A, D)
    })
    #here we need to amplify, to not get zeros
    log_like <- log_like - max(log_like)
    I_split[i_] <- extraDistr::rcat(1, exp(log_like)) 
  }
  
  #do i_G intermediate Gibbs sampling steps
  for(i_it in 1:i_G){
    for(cell in setdiff(chosen_clust,c(i,j))){
      old_c <- t_split[str(cell),]
      #here we incorporate the model scRNA likelihood
      c_log_p <- sapply(1:2, function(i_){
                         c(log(n_split[i_]-(old_c==(i_-1))),
                         loglike_BCR(B_split[i_,,], BCR[cell,,]),
                         loglike_RNA(cell, I_split[i_], C, theta1, theta0, A, D))
                         })
      
      for(j_ in 1:3){
        if(is.infinite(max(c_log_p[j_,]))){
          #c_log_p[j_,] <- c(0,0)
          return(list('inf'=TRUE))
        }else{
          c_log_p[j_,] <- c_log_p[j_,] - max(c_log_p[j_,])
        }
      }
      
      log_p <- apply(c_log_p, 2, sum)
      
      p <- exp(log_p - max(log_p))
      
      min_maxed <- max(min(p[2]/sum(p) ,1) ,0)
      new_c <- rbinom(1, 1, prob=min_maxed)
      
      if(new_c != old_c){
         n_split[old_c+1] <- n_split[old_c+1] - 1
         n_split[new_c+1] <- n_split[new_c+1] + 1
         t_split[str(cell),] <- new_c
      }
    }
    for(i_ in 1:2){
      clust_i <- chosen_clust[t_split$clust==(i_-1)]
      
      if(length(clust_i) > 1){
        up_prior <- as.matrix(g) + apply(BCR[clust_i,,], 
                                         c(2,3), 
                                         sum)
      }else{
        up_prior <- as.matrix(g) + BCR[clust_i,,]
      }
      
      B_split[i_,,] <- t(sapply(1:L, function(l){igraph::sample_dirichlet(1, up_prior[l,])}))
    }
    for(i_ in 1:2){
      log_like <- sapply(1:K, function(k){
        loglike_RNA(chosen_clust[t_split$clust==(i_-1)], k, C, theta1, theta0, A, D)
        })
      #here we need to amplify, to not get zeros
      log_like <- log_like - max(log_like)
      I_split[i_] <- extraDistr::rcat(1, exp(log_like)) 
    }
  }
  
  list('B'=B_split, 't'=t_split, 'I'=I_split, 'n'=n_split, 'inf'=FALSE)
}


launch_merge <- function(chosen_clust){
  up_prior <- as.matrix(g) + apply(BCR[chosen_clust,,], 
                                   c(2,3), 
                                   sum)
  B_merge <- t(sapply(1:L, function(l){igraph::sample_dirichlet(1, up_prior[l,])}))
  
  log_like <- sapply(1:K, function(k){loglike_RNA(chosen_clust, k, C, theta1, theta0, A, D)})
  log_like <- log_like - max(log_like)
  I_merge <- extraDistr::rcat(1, exp(log_like))
  
  list('B'=B_merge, 'I'=I_merge)
}


propose_split <- function(chosen_clust, i, j, 
                          t_split, n_split, 
                          B_split,
                          I_split){
  log_p_propose <- 0
  
  for(cell in setdiff(chosen_clust,c(i,j))){
    old_c <- t_split[str(cell),]
    
    c_log_p <- sapply(1:2, function(i_){
      c(log(n_split[i_]-(old_c==(i_-1))),
        loglike_BCR(B_split[i_,,], BCR[cell,,]),
        loglike_RNA(cell, I_split[i_], C, theta1, theta0, A, D))
    })
    
    for(j_ in 1:3){
      if(is.infinite(max(c_log_p[j_,]))){
        #c_log_p[j_,] <- c(0,0)
        return('inf'=TRUE)
      }else{
        c_log_p[j_,] <- c_log_p[j_,] - max(c_log_p[j_,])
      }
    }
    
    log_p <- apply(c_log_p, 2, sum)
    
    p <- exp(log_p - max(log_p))
    
    min_maxed <- max(min(p[2]/sum(p) ,1) ,0)
    
    new_c <- rbinom(1, 1, prob=min_maxed)
    
    if(new_c != old_c){
      n_split[old_c+1] <- n_split[old_c+1] - 1
      n_split[new_c+1] <- n_split[new_c+1] + 1
      t_split[str(cell),] <- new_c
    }
    
    log_p_propose <- log_p_propose + log(p[new_c+1]/sum(p))
  }
  
  for(i_ in 1:2){
    clust_i <- chosen_clust[t_split$clust==(i_-1)]
    
    if(length(clust_i) > 1){
      up_prior <- as.matrix(g) + apply(BCR[clust_i,,], 
                                       c(2,3), 
                                       sum)
    }else{
      up_prior <- as.matrix(g) + BCR[clust_i,,]
    }
    
    B_split[i_,,] <- t(sapply(1:L, function(l){igraph::sample_dirichlet(1, up_prior[l,])}))
  }
  
  for(i_ in 1:2){
    log_like <- sapply(1:K, function(k){
                              loglike_RNA(chosen_clust[t_split$clust==(i_-1)], k, C, theta1, theta0, A, D)
      })
    
    #here we need to amplify, to not get zeros
    like <- exp(log_like - max(log_like))
    I_split[i_] <- extraDistr::rcat(1, like)
  }
  
  list('t'=t_split,
       'B'=B_split,
       'I'=I_split,
       'log_p'=log_p_propose,
       'n'=n_split,
       'inf'=FALSE)
}


accept_MH_split <- function(t_propose,
                            chosen_clust, 
                            n_proposed, log_p){
  
  #calc log of fraction p_bcr
  clust_i <- chosen_clust[t_propose$clust==0]
  clust_j <- chosen_clust[t_propose$clust==1]
  
  nm <- apply(BCR[chosen_cells,,], c(2,3), sum)
  
  if(length(clust_i)>1){
    ni <- apply(BCR[clust_i,,], c(2,3), sum)
  }else{
    ni <- BCR[clust_i,,]
  }
  
  if(length(clust_j)>1){
    nj <- apply(BCR[clust_j,,], c(2,3), sum)
  }else{
    nj <- BCR[clust_j,,]
  }
  
  upm <- as.matrix(g) + nm
  upi <- as.matrix(g) + ni
  upj <- as.matrix(g) + nj
  
  log_pbcr <- sum(lgamma(upi)) + 
              sum(lgamma(upj))  -
              sum(lgamma(upm)) - 
              sum(lgamma(g)) +
              sum(lgamma(apply(upm,1,sum))) +
              sum(lgamma(apply(g,1,sum))) -
              sum(lgamma(apply(upi,1,sum))) - 
              sum(lgamma(apply(upj,1,sum))) 
  #+sum(lfactorial(nm))-sum(lfactorial(ni)) -sum(lfactorial(nj)) +sum(lfactorial(apply(ni,1,sum))) +sum(lfactorial(apply(nj,1,sum)))-sum(lfactorial(apply(nm,1,sum)))
  
  #calc log of fraction p_scRNA
  
  m_log <- m_loglike_RNA(chosen_clust, C, theta1, theta0, K, A, D)
  
  m_log_i <- m_loglike_RNA(clust_i, C, theta1, theta0, K, A, D)
  m_log_j <- m_loglike_RNA(clust_j, C, theta1, theta0, K, A, D)
  
  l_prna <- logSumExp(m_log_i) +logSumExp(m_log_j) -logSumExp(m_log)
  
  #calc log of fraction p_T
  log_pt <- -lchoose(sum(n_proposed)-2, n_proposed[1]-1) -log(sum(n_proposed)-1) -log_p

  #calc the factor
  R <- alpha_0 /K *M *(sum(n_proposed)-1) /n_proposed[1] /n_proposed[2] /(sum(sapply(clusters_, inv_len))**2)
  
  #calc log of acceptance probability
  log_acc <- log_pt +log(R) +log_pbcr +l_prna 
  
  log_acc <- max(-50, log_acc)
  
  runif(1) < exp(log_acc)
}


accept_MH_merge <- function(t_split, i, j,
                            cluster_i, cluster_j,
                            chosen_cells,
                            n_split, B_split, I_split){
  
  log_p_propose <- 0
  
  for(cell in setdiff(chosen_cells, c(i,j))){
    old_c <- t_split[str(cell),]
    
    c_log_p <- sapply(1:2, function(i_){
      c(log(n_split[i_]-(old_c==(i_-1))),
        loglike_BCR(B_split[i_,,], BCR[cell,,]),
        loglike_RNA(cell, I_split[i_], C, theta1, theta0, A, D))
    })
    
    for(j_ in 1:3){
      if(is.infinite(max(c_log_p[j_,]))){
        c_log_p[j_,] <- c(0, 0)
      }else{
        c_log_p[j_,] <- c_log_p[j_,] - max(c_log_p[j_,])
      }
    }
    
    log_p <- apply(c_log_p, 2, sum)
    
    p <- exp(log_p - max(log_p))
    
    new_c <- ifelse(t_[cell]==t_[i], 0, 1)
    
    if(new_c != old_c){
      n_split[old_c+1] <- n_split[old_c+1] - 1
      n_split[new_c+1] <- n_split[new_c+1] + 1
    }
    
    log_p_propose <- log_p_propose + log(p[new_c+1]/sum(p))
  }
  
  #calc log of fraction p_bcr
  nm <- apply(BCR[chosen_cells,,], c(2,3), sum)
  
  if(length(cluster_i)>1){
    ni <- apply(BCR[cluster_i,,], c(2,3), sum)
  }else{
    ni <- BCR[cluster_i,,]
  }
  
  if(length(cluster_j)>1){
    nj <- apply(BCR[cluster_j,,], c(2,3), sum)
  }else{
    nj <- BCR[cluster_j,,]
  }
  
  upm <- as.matrix(g) + nm
  upi <- as.matrix(g) + ni
  upj <- as.matrix(g) + nj
  
  log_pbcr <- sum(lgamma(upm)) +
              sum(lgamma(g))  -
              sum(lgamma(upi)) - 
              sum(lgamma(upj)) +
              sum(lgamma(apply(upi,1,sum))) +
              sum(lgamma(apply(upj,1,sum))) -
              sum(lgamma(apply(upm,1,sum))) -
              sum(lgamma(apply(g,1,sum))) 
  #+sum(lfactorial(ni)) +sum(lfactorial(nj))-sum(lfactorial(nm)) +sum(lfactorial(apply(nm,1,sum))) -sum(lfactorial(apply(ni,1,sum))) -sum(lfactorial(apply(nj,1,sum)))
                         
  #calc log of fraction p_scRNA
  
  m_log <- m_loglike_RNA(chosen_cells, C, theta1, theta0, K, A, D)
  
  m_log_i <- m_loglike_RNA(cluster_i, C, theta1, theta0, K, A, D)
  m_log_j <- m_loglike_RNA(cluster_j, C, theta1, theta0, K, A, D)
  
  l_prna <- logSumExp(m_log) - logSumExp(m_log_i) - logSumExp(m_log_j) 
  
  #calc log of fraction p_T
  log_pt <- lchoose(sum(n_split)-2, n_split[1]-1) + log(sum(n_split)-1) + log_p_propose
  
  #calc the factor
  R <- K /alpha_0 /(sum(n_split)-1) *n_split[1]**2 *n_split[2]**2 /M *(sum(sapply(clusters_, inv_len))**2)
  
  #calc log of acceptance probability
  log_acc <- log_pt +log(R) +log_pbcr +l_prna
  
  log_acc <- max(-50, log_acc)
  
  runif(1) < exp(log_acc)
}


#helpful funcs
str <- function(x) as.character(x)


resample <- function(x) x[sample.int(length(x), 1)]


inv_len <- function(x) ifelse(length(x)==0, 0, 1/length(x))


loglike_RNA <- function(cells, clone, C, theta1, theta0, A, D){
  sum(dbinom(A[,cells],
             D[,cells], 
             prob=C[,clone]*theta1+(1-C[,clone])*theta0,
             log=TRUE))
}


m_loglike_RNA <- function(cells, C, theta1, theta0, K, A, D){
  sapply(1:K, function(k){sum(
                 dbinom(A[,cells], D[,cells],
                 prob=C[,k]*theta1 + (1-C[,k])*theta0, log=TRUE)
                 )})
}


loglike_BCR <- function(B, bcr){
  sum(log(B[bcr==1]))
}
