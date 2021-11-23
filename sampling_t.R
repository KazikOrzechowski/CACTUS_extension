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
  t_and_B <- launch_split(i, j, chosen_clust, X_bcr)
  t_split <- t_and_B$t
  B_split <- t_and_B$B
  I_split <- t_and_B$I
  
  #launch merge
  B_merge <- launch_merge(X_bcr)
  
  #define proposed state
  
  #assign first empty cluster to cell i and
  first_empty <- which(sapply(clusters_, length)==0, arr.ind=TRUE)[1]
  t_[i] <- if(is.na(first_empty)){length(clusters_)+1}else{first_empty}
 
  
  return(list('cells'=chosen_cells, 
              't_split'=t_split, 
              'B_split'=B_split, 
              'B_merge'=B_merge, 
              'I_split'=I_split,
              'BCRs'=X_bcr))
}

try_merge <- function(){
  #pick two clusters to merge inversely proportional to their size
  p_clusters <- sapply(clusters_, function(cluster){
    if(length(cluster)>0){M/length(cluster)}else{0}
    })
  
  #if only one cluster exists it has p=1 and we cannot merge
  if(any(p_clusters==1)) return()
  chosen <- sample(length(clusters_), 2, prob=p_clusters)
  
  #draw a cell from each chosen cluster uniformly
  i <- resample(clusters_[[chosen[1]]])
  j <- resample(clusters_[[chosen[2]]])
  return(c(i,j))
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
  
  list('B'=B_split, 't'=t_split, 'I'=I_split)
}

launch_merge <- function(X_bcr){
  up_prior <- as.matrix(g) + Reduce('+', X_bcr)
  B_merge <- t(sapply(1:L, function(l){sample_dirichlet(1, up_prior[l,])}))
}

#helpful funcs
str <- function(x) as.character(x)
resample <- function(x) x[sample.int(length(x), 1)]


a <- try_split()

b <- a[[2]]

rcatlp(n=1, c(-3029, -5576, 0))
