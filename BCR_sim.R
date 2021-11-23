library('igraph')
library('extraDistr')
library('matrixStats')

#bcr len
L <- 50

#num of cells
M <- 1000

#vec of nuc names
nuc <- c('A', 'C', 'G', 'T')

#g are the prior parameters of nuc probabilities. dim(g) is L*4
g <- data.frame('A'=rpois(L,2)+1, 
                'C'=rpois(L,3)+1,
                'G'=rpois(L,1)+1,
                'T'=rpois(L,1)+1)

#alpha_0 is the concentration parameter of CRP
alpha_0 <- 10


#simulate the clustering
simulate_t_and_cl <- function(){
  t_ <- 1:M
  clusters_ <- list()
  clusters_[[1]] <- 1
  for(i in 2:M){
    p <- sapply(clusters_, function(j){length(j)/(i-1+alpha_0)})
    p <- append(p, alpha_0/(i-1+alpha_0))
    t_[i] <- rcat(1,p)
    if(t_[i]>length(clusters_)){
      clusters_[[t_[i]]] <- i
    }else{
      clusters_[[t_[i]]] <- append(x=clusters_[[t_[i]]], values=i)
    }
  }
  list('clust'=clusters_, 't'=t_)
}



#B is a list len Q, each entry is L*4 probablities of each nuc at spot l of cluster q
simulate_B <- function(){
  B <- array(dim=c(M, L, 4))
  for(q in 1:length(clusters_)){
    temp <- t(sapply(1:L, function(l){sample_dirichlet(1, g[l,])}))
    colnames(temp) <- nuc
    B[q, ,] <- temp
  }
  B
}


#cell x is mat L*4, each row is one-hot of the nuc at spot l
#we simulate it with
#x <- t(sapply(1:L, function(l){rmultinom(1, 1, B[[t_[i]]][l,])}))
#colnames(x) <- nuc


#simulate the BCR seqs of cells
simulate_BCR <- function(){
  As <- Cs <- Gs <- Ts <- data.frame()
  BCR <- list('A'=As, 'C'=Cs, 'G'=Gs, 'T'=Ts)
  for(i in 1:M){
    temp <- t(sapply(1:L, function(l){rmultinom(1, 1, B[t_[i],l,])}))
    colnames(temp) <- nuc
    for(n in nuc){
      if(length(which(temp[,n]==1,arr.ind=TRUE))){
        BCR[[n]] <- rbind(BCR[[n]], 
                          data.frame(cell=i, position=which(temp[,n]>0, arr.ind=TRUE)))
      }
    }
  }
  BCR
}


#calculate BCR likelihood
calc_bcr_model_like <- function(){
  log_like <- 0
  for(n in 1:4){
    clust_vec <- t_[BCR[[n]]$cell]
    mat_inds <- clust_vec + BCR[[n]]$position*10 - 10
    log_like <- sum(log(B[,, n][mat_inds])) + log_like
  }
  exp(log_like)
}

calc_bcr_cell_like <- function(cell, B){
  log_like <- 0
  for(n in 1:4){
    pos <- BCR[[n]][BCR[[n]]$cell == cell,2]
    log_like <- sum(log(B[,n][pos])) + log_like
  }
  exp(log_like)
}

##################################################################################
t_and_clust <- simulate_t_and_cl()

#t_ is cell to cluster assign
t_ <- t_and_clust$t

#clusters_ is the list of clusters containing their members
clusters_ <- t_and_clust$clust

B <- simulate_B()
BCR <- simulate_BCR()

calc_bcr_model_like()

