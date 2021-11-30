source('scRNA_sim.R')


#number of intermediate Gibbs sampling steps
i_G <- 3


try_split <- function(){
  #pick a cluster to split proportional to its size. Can't split a cluster < 2
  p_clusters <- sapply(clusters_, function(cluster){
                                      if(length(cluster)>1){length(cluster)/M}else{0}
                                      })
  if(!any(p_clusters>0)) return()
  chosen_ind <- sample(length(clusters_), 1, prob=p_clusters)
  chosen_clust <- clusters_[[chosen_ind]]
  
  #chose two cells i and j from the cluster uniformly
  chosen_cells <- sample(chosen_clust, 2)
  i <- chosen_cells[1]
  j <- chosen_cells[2]
  
  #X_bcr is a list of one-hot encoded bcr seqs of cells from in chosen_clust, for faster computation
  X_bcr <- list()
  for(cell in chosen_clust){
    temp <- array(0,dim=c(L,4))
    colnames(temp) <- nuc
    for(n in nuc){
      temp[,n][BCR[[n]][BCR[[n]]$cell == cell,2]] <- 1
    }
    X_bcr[[str(cell)]] <- temp
  }
  
  #create the launch split and launch merge states
  #launch split
  split_ <- launch_split(i, j, chosen_clust, X_bcr)
  t_split <- split_$t
  B_split <- split_$B
  I_split <- split_$I
  n_split <- split_$n
  
  #launch merge
  merge_ <- launch_merge(X_bcr, chosen_clust)
  B_merge <- merge_$B
  I_merge <- merge_$I
  
  #define proposed state
  propose <- propose_split(chosen_clust, i, j,
                            t_split, n_split, B_split, 
                            X_bcr, I_split)
  t_propose <- propose$t
  B_propose <- propose$B
  I_propose <- propose$I
  log_p_propose <- propose$log_p
  log_d_propose <- propose$log_d
  n_proposed <- propose$n
  
  #calc M-H acceptance prob
  MH <- calc_MH_split(chosen_ind, chosen_clust,
                log_d_propose, log_p_propose,
                n_proposed, t_propose,
                B_propose, I_propose,
                X_bcr, B_split)
  
  
  #if accepted, assign first empty cluster to proposed cluster with cell i
  print(MH)
  print(MH$acc * MH$ext)
  if(runif(1) < min(1, MH$acc * MH$ext)){
    print('Accepted')
    first_empty <- which(sapply(clusters_, length)==0, arr.ind=TRUE)[1]
    cluster_i <- chosen_clust[t_propose$clust == 0]
    cluster_j <- chosen_clust[t_propose$clust == 1]
    t_[cluster_i] <- ifelse(is.na(first_empty), length(clusters_)+1, first_empty)
    B[t_[i],,] <- B_propose[1,,]
    B[t_[j],,] <- B_propose[2,,]
    I[t_[i]] <- I_propose[1]
    I[t_[j]] <- I_propose[2]
    clusters_[[t_[j]]] <- cluster_j
    if(is.na(first_empty)){
      clusters_ <- append(clusters_, cluster_i)
    }else{
      clusters_[[first_empty]] <- cluster_i
    }
  }
  
  list('cells'=chosen_cells, 
       't_split'=t_split, 
       'B_split'=B_split, 
       'B_merge'=B_merge, 
       'I_split'=I_split,
       'I_merge'=I_merge,
       'BCRs'=X_bcr,
       't_prop'= t_propose,
       'B_prop'= B_propose,
       'I_prop'= I_propose,
       'log_p'=log_p_propose,
       'log_d'=log_d_propose,
       'MH'=MH)
}

try_merge <- function(){
  #pick two clusters to merge inversely proportional to their size
  p_clusters <- sapply(clusters_, function(cluster){
    if(length(cluster)>0){M/length(cluster)}else{0}
    })
  
  #if only one cluster exists it has p=1 and we cannot merge
  if(any(p_clusters==1)) return()
  chosen <- sample(length(clusters_), 2, prob=p_clusters)
  
  cluster_i <- clusters_[[chosen[1]]]
  cluster_j <- clusters_[[chosen[2]]]
  
  #draw a cell from each chosen cluster uniformly
  i <- resample(cluster_i)
  j <- resample(cluster_j)

  #X_bcr is a list of one-hot encoded bcr seqs of cells from in chosen_clust, for faster computation
  X_bcr <- list()
  for(cell in union(cluster_i, cluster_j)){
    temp <- array(0,dim=c(L,4))
    colnames(temp) <- nuc
    for(n in nuc){
      temp[,n][BCR[[n]][BCR[[n]]$cell == cell,2]] <- 1
    }
    X_bcr[[str(cell)]] <- temp
  }
  
  #create the launch split and launch merge states
  #launch split
  split_ <- launch_split(i, j, union(cluster_i, cluster_j), X_bcr)
  t_split <- split_$t
  B_split <- split_$B
  I_split <- split_$I
  n_split <- split_$n
  
  #we do not need to launch merge in this case, as the probabilities of parameters in proposed state
  #are independent from the launch state
  proposed <- launch_merge(X_bcr, union(cluster_i, cluster_j))
  B_proposed <- proposed_$B
  I_proposed <- proposed_$I
  
  list(
    'cells'=c(i,j),
    'X'=X_bcr,
    't_split'=t_split,
    'B_split'=B_split,
    'I_split'=I_split,
    'B_merge'=B_proposed,
    'I_merge'=I_proposed
  )
}

launch_split <- function(i,j,chosen_clust, X_bcr){
  #randomly assign each cell from cluster to i or j (these stay the same)
  t_split <- data.frame(clust=rbinom(chosen_clust, 1, .5), 
                        row.names = str(chosen_clust))
  t_split[str(i),] <- 0
  t_split[str(j),] <- 1
  
  #n_split are the populations of split clusters
  n_split <- 1:2
  n_split[2] <- sum(t_split)
  n_split[1] <- length(chosen_clust) - n_split[2]
  
  #B_split are the nuc probabilities of split clusters
  B_split <- array(0, dim=c(2, L, 4))
  B_split[1,,] <- t(sapply(1:L, function(l){sample_dirichlet(1, g[l,])}))
  B_split[2,,] <- t(sapply(1:L, function(l){sample_dirichlet(1, g[l,])}))
  
  #I_split are the cluster-clone assignments of split clusters
  I_split <- rcat(2, Psi)
  
  #do i_G intermediate Gibbs sampling steps
  for(i_it in 1:i_G){
    for(cell in setdiff(chosen_clust,c(i,j))){
      old_c <- t_split[str(cell),]
      #here we will need to incorporate the model scRNA likelihood
      log_p <- sapply(1:2, function(i_){
                         log(n_split[i_]-(old_c==i_-1))+
                         sum(log(B_split[i_,,][X_bcr[[str(cell)]]==1]))+
                         sum(dbinom(A[cell,], D[cell,], C[,I_split[i_]]*theta1+(1-C[,I_split[i_]])*theta0, log=TRUE))
                  })
      p <- exp(log_p - max(log_p))
      new_c <- rbinom(1, 1, prob=p[2]/sum(p))
      if(new_c != old_c){
         n_split[old_c+1] <- n_split[old_c+1] - 1
         n_split[new_c+1] <- n_split[new_c+1] + 1
         t_split[str(cell),] <- new_c
      }
    }
    for(i_ in 1:2){
      up_prior <- as.matrix(g) + Reduce('+', X_bcr[t_split$clust==i_-1])
      B_split[i_,,] <- t(sapply(1:L, function(l){sample_dirichlet(1, up_prior[l,])}))
    }
    for(i_ in 1:2){
      log_like <- sapply(1:K, function(k){
        sum(dbinom(t(A[chosen_clust[t_split$clust==i_-1],]),
                   t(D[chosen_clust[t_split$clust==i_-1],]), 
                   prob=C[,k]*theta1+(1-C[,k])*theta0,
                   log=TRUE))
        })
      #here we need to amplify, to not get zeros
      log_like <- log_like - max(log_like)
      I_split[i_] <- rcat(1, exp(log_like)) 
    }
  }
  
  list('B'=B_split, 't'=t_split, 'I'=I_split, 'n'=n_split)
}

launch_merge <- function(X_bcr, chosen_clust){
  up_prior <- as.matrix(g) + Reduce('+', X_bcr)
  B_merge <- t(sapply(1:L, function(l){sample_dirichlet(1, up_prior[l,])}))
  
  log_like <- sapply(1:K, function(k){
    sum(dbinom(t(A[chosen_clust,]),
               t(D[chosen_clust,]), 
               prob=C[,k]*theta1+(1-C[,k])*theta0,
               log=TRUE))
  })
  log_like <- log_like - max(log_like)
  I_merge <- rcat(1, exp(log_like))
  
  list('B'=B_merge, 'I'=I_merge)
}


propose_split <- function(chosen_clust, i, j, 
                          t_split, n_split, 
                          B_split, X_bcr,
                          I_split){
  log_p_propose <- 0
  log_d_propose <- 0
  
  for(cell in setdiff(chosen_clust,c(i,j))){
    old_c <- t_split[str(cell),]
    #here we will need to incorporate the model scRNA likelihood
    log_p <- sapply(1:2, function(i_){
      log(n_split[i_]-(old_c==i_-1))+
        sum(log(B_split[i_,,][X_bcr[[str(cell)]]==1]))+
        sum(dbinom(A[cell,], D[cell,], C[,I_split[i_]]*theta1+(1-C[,I_split[i_]])*theta0, log=TRUE))
    })
    p <- exp(log_p - max(log_p))
    new_c <- rbinom(1, 1, prob=p[2]/sum(p))
    if(new_c != old_c){
      n_split[old_c+1] <- n_split[old_c+1] - 1
      n_split[new_c+1] <- n_split[new_c+1] + 1
      t_split[str(cell),] <- new_c
    }
    log_p_propose <- log_p_propose + log(p[new_c+1]/sum(p))
  }
  
  for(i_ in 1:2){
    up_prior <- as.matrix(g) + Reduce('+', X_bcr[t_split$clust==i_-1])
    B_split[i_,,] <- t(sapply(1:L, function(l){sample_dirichlet(1, up_prior[l,])}))
    log_d_propose <- log_d_propose + sum(sapply(1:L, function(l){ddirichlet(
                                                                    B_split[i_,l,],
                                                                    up_prior[l,],
                                                                    log=TRUE)}))
  }
  
  for(i_ in 1:2){
    log_like <- sapply(1:K, function(k){
      sum(dbinom(t(A[chosen_clust[t_split$clust==i_-1],]),
                 t(D[chosen_clust[t_split$clust==i_-1],]), 
                 prob=C[,k]*theta1+(1-C[,k])*theta0,
                 log=TRUE))
    })
    #here we need to amplify, to not get zeros
    like <- exp(log_like - max(log_like))
    I_split[i_] <- rcat(1, like)
    log_p_propose <- log_p_propose + log(like[I_split[i_]] / sum(like))
  }
  
  list('t'=t_split,
       'B'=B_split,
       'I'=I_split,
       'log_p'=log_p_propose,
       'log_d'=log_d_propose,
       'n'=n_split)
}

calc_MH_split <- function(chosen_ind, chosen_clust,
                          log_d_propose, log_p_propose,
                          n_proposed, t_propose,
                          B_propose, I_propose,
                          X_bcr, B_split){
  #calculate q_ratio, the fraction of proposal densities
  up_prior <- as.matrix(g) + Reduce('+', X_bcr)
  log_d <- sum(sapply(1:L, function(l){ddirichlet(B[chosen_ind,l,],
                                                  up_prior[l,],
                                                  log=TRUE)}))
  
  log_like <- sapply(1:K, function(k){sum(dbinom(t(A[chosen_clust,]),
                                                 t(D[chosen_clust,]), 
                                                 prob=C[,k]*theta1+(1-C[,k])*theta0,
                                                 log=TRUE))
  })
  like <- exp(log_like - max(log_like))
  log_p <- log(like[I[chosen_ind]] / sum(like))
  
  q_ratio <- exp(log_d - log_d_propose) * exp(log_p - log_p_propose) 
  
  #calculate p_ratio, the ratio of prior densities
  log_d_prior <- 0
  for(i_ in 1:2){
    log_d_prior <- log_d_prior + sum(sapply(1:L, function(l){ddirichlet(
                                                                        B_split[i_,l,],
                                                                        g[l,],
                                                                        log=TRUE)}))
  }
  
  log_d_prior <- log_d_prior - sum(sapply(1:L, function(l){ddirichlet(
                                                                      B[chosen_ind,l,],
                                                                      g[l,],
                                                                      log=TRUE
                                                                      )}))
  
  p_ratio <- alpha_0 / K / choose(sum(n_proposed)-2, n_proposed[1]-1) / (sum(n_proposed)-1) * exp(log_d_prior)
  
  #calculate l_ratio, the ratio of likelihoods
  log_l_ratio <- 0
  for(cell in chosen_clust){
    log_l_ratio <- log_l_ratio + sum(log(B_propose[t_propose[str(cell),]+1,,][X_bcr[[str(cell)]]==1])) - sum(log(B[t_[cell],,][X_bcr[[str(cell)]]==1]))
    old_k <- I[t_[cell]]
    new_k <- I_propose[t_propose[str(cell),]+1]
    if(old_k != new_k){
      log_l_ratio <- log_l_ratio + sum(dbinom(A[cell,], D[cell,], prob=C[,new_k]*theta1+(1-C[,new_k])*theta0, log=TRUE)) - sum(dbinom(A[cell,], D[cell,], prob=C[,old_k]*theta1+(1-C[,old_k])*theta0, log=TRUE))
    }
  }
  l_ratio <- exp(log_l_ratio)
  
  acc <- q_ratio * l_ratio * p_ratio
  
  ext <- M * (sum(n_proposed)- 1)  / product(n_proposed[1]) / (sum(sapply(clusters_, inv_len)) ** 2)

  list('q'=q_ratio,
       'p'=p_ratio,
       'l'=l_ratio, 
       'acc'=acc,
       'ext'=ext)
}

#helpful funcs
str <- function(x) as.character(x)
resample <- function(x) x[sample.int(length(x), 1)]
inv_len <- function(x) ifelse(length(x)==0, 0, 1/length(x))

unknown <- simulate_t_and_cl()

t_true <- t_
clusters_true <- clusters_

t_ <- unknown$t
clusters_ <- unknown$clust

a <- try_split()

b <- try_merge()

I[t_true[897]]
