library('igraph')
library('extraDistr')
library('matrixStats')


resample_few <- function(x, how_many) x[sample.int(length(x), how_many)]


#simulate the clustering
simulate_t_and_cl <- function(alpha_0, M){
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
simulate_B <- function(M, L, K, I, clusters_, g){
  nuc <- c('A', 'C', 'G', 'T')
  
  B <- array(dim=c(M, L, 4))
  BCR_mut <- list()
  up_prior <- list()
  for(k in 1:K){
    BCR_mut[[k]] <- (1:(L%/%K))+ (k-1)*(L%/%K)
  }
  for(q in 1:length(clusters_)){
    n_possible <- length(BCR_mut[[I[q]]])
    exhibited <- resample_few(BCR_mut[[I[q]]], n_possible%/%2)
    n_possible <- length(BCR_mut[[(I[q]+1)%%K+1]])
    exhibited_2 <- resample_few(BCR_mut[[(I[q]+1)%%K+1]], n_possible%/%2)
    temp <- g
    temp$C[exhibited] <- g$A[exhibited]
    temp$A[exhibited] <- g$C[exhibited]
    temp$G[exhibited_2] <- g$C[exhibited_2]
    temp$C[exhibited_2] <- g$G[exhibited_2]
    freq <- t(sapply(1:L, function(l){sample_dirichlet(1, temp[l,])}))
    colnames(freq) <- nuc
    B[q, ,] <- freq
  }
  list('B'=B, 'BCR_mut'=BCR_mut)
}


#cell x is mat L*4, each row is one-hot of the nuc at spot l
#we simulate it with
#x <- t(sapply(1:L, function(l){rmultinom(1, 1, B[[t_[i]]][l,])}))
#colnames(x) <- nuc


#simulate the BCR seqs of cells
simulate_BCR <- function(M, L, B, t_){
  BCR <- array(0, dim=c(M, L, 4))
  for(i in 1:M){
    BCR[i,,] <- t(sapply(1:L, function(l){rmultinom(1, 1, B[t_[i],l,])}))
  }
  BCR
}


data_simulation <- function(L= 300,
                            M= 200,
                            K= 3,
                            N= 100,
                            g= data.frame('A'=rpois(L,1/100)+1, 
                                'C'=rpois(L,1000)+1,
                                'G'=rpois(L,1/100)+1,
                                'T'=rpois(L,1/100)+1),
                            alpha_0= 4,
                            relax_rate_prior= c(0.5, 9.5),
                            prior0= c(0.2, 99.8),
                            prior1= c(4.5, 5.5),
                            mut_freq= 0.3,
                            av_reads= 1
){
  nuc <- c('A', 'C', 'G', 'T')
  
  t_and_clust <- simulate_t_and_cl(alpha_0, M)
  
  #t_ is cell to cluster assign
  t_ <- t_and_clust$t
  
  #clusters_ is the list of clusters containing their members
  clusters_ <- t_and_clust$clust
  
  #relax_rate imulation
  relax_rate <- rbeta(1, relax_rate_prior[1], relax_rate_prior[2])
  
  #theta priors and simulation
  theta0 <- rbeta(1, prior0[1], prior0[2])
  theta1 <- rbeta(N, prior1[1], prior1[2])
  
  #Psi is the prior distribution of cluster-clone assignment
  Psi <- replicate(K, 1/K)
  
  #Omega and C simulation
  Omega <- array(rbinom(N*K, 1, mut_freq), dim=c(N, K))
  p <- abs(Omega-relax_rate)
  C <- array(rbinom(N*K, 1, prob=p), dim=c(N,K))
  
  #cluster-clone assignment
  I <- 1:M
  I[1:length(clusters_)] <- rcat(length(clusters_), Psi)
  
  #All read counts and mutation reads simulation
  D <- array(rpois(N*M, av_reads), dim=c(N, M))
  p <- C[,I[t_]]*theta1 + (1-C[,I[t_]])*theta0
  A <- array(rbinom(N*M, D, prob=p), dim=c(N,M))
  
  B_and_mut <- simulate_B(M, L, K, I, clusters_, g)
  B <- B_and_mut$B
  B_mut <- B_and_mut$BCR_mut
  
  BCR <- simulate_BCR(M, L, B, t_)
  
  list('clustering'=clusters_,
       't'=t_,
       'rr'=relax_rate,
       'theta0'=theta0,
       'theta1'=theta1,
       'Omega'=Omega,
       'C'=C,
       'I'=I,
       'D'=D,
       'A'=A,
       'B'=B,
       'BCR'=BCR)
}



