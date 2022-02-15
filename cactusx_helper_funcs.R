#Variables:
#C: The observed Omega matrix of clonal profiles
#C: The true C matrix of clonal profiles
#prob_mat: the likelihood of the A and D reads
#S1_list: for each clone: a number of reads w/ mutations in cells from clones w/o the mutation
#S2_list: w/o -||- w/o
#S3_list: w/ -||- w/
#S4_list: w/o -||- w/
#Psi: prior probability of cluster -> clone assignment
#A1: numbers of observed mutations
#B1: numbers of observed reference
#K: number of clones
#Q: number of clusters
#t: cluster -> clone assignment
#it: iteration number

Geweke_Z <- function(X, first=0.1, last=0.5) {
  
  if (is.null(dim(X))){
    X <- as.matrix(X, ncol = 1)}
  N <- nrow(X)
  A <- X[1:floor(first*N), , drop = FALSE]
  B <- X[ceiling(last*N):N, , drop = FALSE]
  
  A_col_mean <- colMeans(A)
  B_col_mean <- colMeans(B)
  A_col_var <- rowSums((t(A) - A_col_mean)^2) / (nrow(A) - 1)#^2
  B_col_var <- rowSums((t(B) - B_col_mean)^2) / (nrow(A) - 1)#^2
  
  min_var <- 10^(-50)
  Z <- (A_col_mean - B_col_mean) / sqrt(A_col_var + B_col_var + min_var)
  
  Z
}


devianceIC <- function(logLik_all, logLik_post) {
  
  logLik_mean = mean(logLik_all)
  logLik_var = var(logLik_all)
  
  p_D_Spiegelhalter = -2 * logLik_mean -  (-2 * logLik_post)
  DIC_Spiegelhalter = -2 * logLik_post + 2 * p_D_Spiegelhalter
  
  p_D_Gelman = 2 * logLik_var
  DIC_Gelman = -2 * logLik_post + 2 * p_D_Gelman
  
  DIC = DIC_Gelman
  
  cat(paste("DIC:", round(DIC, 2), 
            "D_mean:", round(-2 * logLik_mean, 2), 
            "D_post:", round(-2 * logLik_post, 2), 
            "logLik_var:", round(logLik_var, 2), "\n"))
  
  list("DIC" = DIC, 
       "logLik_var" = logLik_var, 
       "D_mean" = -2 * logLik_mean, 
       "D_post" = -2 * logLik_post, 
       "DIC_Gelman" = DIC_Gelman, 
       "DIC_Spiegelhalter" = DIC_Spiegelhalter)
}


get_logLik <- function(A1, B1, C, Assign, theta0, theta1) {
  if (is.null(dim(Assign)) || length(dim(Assign)) == 1) {
    Assign_prob <- matrix(0, length(Assign), ncol(C))
    for (i in seq_len(length(Assign))) {
      Assign_prob[i, Assign[i]] = 1
    }
  } else {
    Assign_prob <- Assign
  }
  
  prob_mat <- C %*% t(Assign_prob)
  
  Lik_mat <- (exp(log(theta1) * A1 + log(1 - theta1) * B1) * prob_mat + 
                exp(log(theta0) * A1 + log(1 - theta0) * B1) * (1 - prob_mat))
  
  logLik <- (sum(log(Lik_mat), na.rm = TRUE) + 
               sum(lchoose(A1 + B1, A1), na.rm = TRUE))
  logLik
}


get_logLik_BCR <- function(BCR, B, clusters_, t_){
  log_lik <- 0
  
  for(cl in clusters_){
    if(length(cl)>1){
      nm <- apply(BCR[cl,,], c(2,3), sum)
        
      log_lik <- log_lik + 
                 sum(sapply(cl, function(cell){
                                sum(log(B[t_[cell],,][BCR[cell,,]==1]))
                                })) +
                 sum(lfactorial(apply(nm, 1, sum))) -
                 sum(lfactorial(nm))
    }
    if(length(cl)==1){
      log_lik <- log_lik + sum(log(B[t_[cl],,][BCR[cl,,]==1]))
    }
  }
  
  log_lik
}

sample_I <- function(){
  I_it <- 1:length(clusters_)
  for(i in 1:length(clusters_)){
    if(length(clusters_[[i]])>0){
      log_like <- sapply(1:K, function(k){
        loglike_RNA(clusters_[[i]], k, C, theta1, theta0, A, D)
      })
      
      #here we need to amplify, to not get zeros
      like <- exp(log_like - max(log_like))
      
      I_it[i] <- rcat(1, like)
    }
  }
  
  I <<- I_it
}


sample_c_and_rr <- function(){
  if (it > 0.1 * min_iter){
    diff0 <- sum(Omega == C)
    diff1 <- sum(Omega != C)
    relax_rate <<- r_rate <- stats::rbeta(1, relax_rate_prior[1] + diff1,
                                          relax_rate_prior[2] + diff0)
    
    #C_prior is the prior probability of C=1 given C:=Omega
    C_prior <- Omega
    C_prior[Omega == 1] <- 1 - r_rate
    C_prior[Omega == 0] <- r_rate
    C_prior_oddlog <- log(C_prior) - log(1 - C_prior)
  }
  
  #Iden_mat is the assignment of cells to clones (through clusters)
  Iden_mat[,] <- 0
  for (j in seq_len(M)) {
    Iden_mat[j, I[t_[j]]] <- 1 
  }
  
  # calculate log_probability matrix with genotype 0 and 1
  P0_mat <- A1 * log(theta0) + B1 * log(1 - theta0) + W_log
  P1_mat <- A1 * log(theta1) + B1 * log(1 - theta1) + W_log
  
  oddR_log <- P1_mat %*% Iden_mat  - P0_mat %*% Iden_mat 
  oddR_log <- oddR_log + C_prior_oddlog
  oddR_log[which(oddR_log > 50)] <- 50
  oddR_log[which(oddR_log < -50)] <- -50
  C_prob_tmp <- exp(oddR_log) / (exp(oddR_log) + 1)
  
  C <<- C_n <- matrix(stats::rbinom(N*K, size = 1, C_prob_tmp), nrow=N)
  
  for (k in seq_len(K)) {
    S1_list[[k]] <<- A1 * (1 - C_n[,k])
    S2_list[[k]] <<- B1 * (1 - C_n[,k])
    S3_list[[k]] <<- A1 * C_n[, k]
    S4_list[[k]] <<- B1 * C_n[, k]
  }
}

sample_theta <- function(){
  S1_wgt <- S2_wgt <- 0 # weighted S1
  S3_wgt <- S4_wgt <- matrix(0, nrow = N, ncol = M)
  for (k in seq_len(K)) {
    idx <- which(I[t_] == k)
    S1_wgt <- S1_wgt + sum(S1_list[[k]][,idx], na.rm = TRUE)
    S2_wgt <- S2_wgt + sum(S2_list[[k]][,idx], na.rm = TRUE)
    S3_wgt[,idx] <- S3_wgt[,idx] + S3_list[[k]][,idx]
    S4_wgt[,idx] <- S4_wgt[,idx] + S4_list[[k]][,idx]
  }
  
  S3_wgt[,] <- rowSums(S3_wgt, na.rm = TRUE)
  S4_wgt[,] <- rowSums(S4_wgt, na.rm = TRUE)
  
  theta0_all[it, 1] <<- stats::rbeta(1, prior0[1] + S1_wgt,
                                     prior0[2] + S2_wgt)
  theta1_all[it,  ] <<- stats::rbeta(rep(1, n_element),
                                     prior1[,1] + S3_wgt[idx_vec],
                                     prior1[,2] + S4_wgt[idx_vec])
}

sample_alpha_0 <- function(alpha_prior=c(1,1)){
  # k is the number of non-empty clusters
  k <- sum(sapply(clusters_, length)>0)
  
  #sample new value of sigma
  xSigma <- rbeta(1, alpha_0+1, M)
  
  y <- alpha_prior[1]
  x <- alpha_prior[2]
  
  from_first <- (runif(1, max=y+k+M*(x-log(xSigma))) < y+k)
  if(from_first){
    alpha_0 <<- rgamma(1, shape=y+k, rate=x-log(xSigma))
  }else{
    alpha_0 <<- rgamma(1, shape=y+k-1, rate=x-log(xSigma))
  }
}


ini_clust <- function(BCR, M, L){
  dfBCR <- data.frame('a'=1:M)
  
  for(i in 1:(4*L)){
    dfBCR[,i] <- 0 
  }
  
  for(i in 1:M){
    dfBCR[i,] <-  BCR[i,,][1:(4*L)]
    
  }
  
  ii <- do.call(grouping, dfBCR)
  
  same_clusts <- list()
  same_t <- 1:M
  
  ends_ <- attr(ii,'ends')
  
  for(i in 1:length(ends_)){
    if(i==1){
      same_clusts[[i]] <- ii[1:ends_[i]]
      same_t[same_clusts[[i]]] <- i
    }else{
      same_clusts[[i]] <- ii[(ends_[i-1]+1):ends_[i]]
      same_t[same_clusts[[i]]] <- i
    }
  }
  
  list('clust'=same_clusts, 't'=same_t)
}