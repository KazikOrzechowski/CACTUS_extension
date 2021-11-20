library('igraph')
library('extraDistr')
library('matrixStats')


#bcr len
L <- 5

#num of cells
M <- 10

#vec of nuc names
nuc <- c('A', 'C', 'G', 'T')

#g are the prior parameters of nuc probabilities. dim(g) is L*4
g <- data.frame('A'=c(5,5,1,1,1), 
                'C'=c(1,1,5,5,1),
                'G'=c(5,1,5,1,1),
                'T'=c(1,5,1,5,1))

#alpha_0 is the concentration parameter of CRP
alpha_0 <- 2


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



#B is an array dim M*L*4, each entry is L*4 probablities of each nuc at spot l of cluster m
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
calc_bcr_like <- function(){
  like <- 1
  for(n in 1:4){
    clust_vec <- t_[BCR[[n]]$cell]
    mat_inds <- clust_vec + BCR[[n]]$position*10 - 10
    like <- product(B[,, n][mat_inds]) * like
  }
  like
}


##################################################################################
t_and_clust <- simulate_t_and_cl()

#t_ is cell to cluster assign
t_ <- t_and_clust$t

#clusters_ is the list of clusters containing their members
clusters_ <- t_and_clust$clust

B <- simulate_B()
BCR <- simulate_BCR()

calc_bcr_like()

