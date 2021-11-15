#Variables:
#Config: The observed Omega matrix of clonal profiles
#Config_new: The true C matrix of clonal profiles
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
#t_: cluster -> clone assignment list
#it: iteration number
#alpha_0: concentration parameter of the clustering


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


get_logLik <- function(A1, B1, Config, Assign, theta0, theta1) {
  if (is.null(dim(Assign)) || length(dim(Assign)) == 1) {
    Assign_prob <- matrix(0, length(Assign), ncol(Config))
    for (i in seq_len(length(Assign))) {
      Assign_prob[i, Assign[i]] = 1
    }
  } else {
    Assign_prob <- Assign
  }
  
  prob_mat <- Config %*% t(Assign_prob)
  
  Lik_mat <- (exp(log(theta1) * A1 + log(1 - theta1) * B1) * prob_mat + 
                exp(log(theta0) * A1 + log(1 - theta0) * B1) * (1 - prob_mat))
  
  logLik <- (sum(log(Lik_mat), na.rm = TRUE) + 
               sum(lchoose(A1 + B1, A1), na.rm = TRUE))
  logLik
}

update_prob_mat <- function(){
  for (k in seq_len(K)) {
    logLik_mat[,k] <- (colSums(S1_list[[k]] * log(theta0) , na.rm = TRUE) +
                         colSums(S2_list[[k]] * log(1 - theta0) , na.rm = TRUE) +
                         colSums(S3_list[[k]] * log(theta1) ,     na.rm = TRUE) +
                         colSums(S4_list[[k]] * log(1 - theta1) , na.rm = TRUE))
    logLik_mat[,k] <- logLik_mat[,k] + log(Psi[k])
  }
  for(q in 1:Q){
    if(sum(colnames(A1) %in% t_[[q]])!=1){
      logLik_mat_q[q,] <- colSums(logLik_mat[colnames(A1) %in% t_[[q]],])
    }else{
      logLik_mat_q[q,] <- logLik_mat[colnames(A1) %in% t_[[q]],]
    }
  }
  
  logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
  logLik_mat_amplify_q <- logLik_mat_q - matrixStats::rowMaxs(logLik_mat_q)
  prob_mat <<- exp(logLik_mat_amplify_q) / rowSums(exp(logLik_mat_amplify_q))
}

sample_cluster_to_clone_assign <- function(){
  for (q in seq_len(Q)) {
    assign_all[it,q] <<- temp <- sample(seq_len(K), 1, replace = TRUE, prob = prob_mat[q,])
    assign_all_j[it,match(t_[[q]], colnames(A1))] <<- temp 
  }
}

sample_c_and_rr <- function(){
  if (it > (0.1 * min_iter + 5)){
    diff0 <- sum((Config == Config_new)[, 2:ncol(Config)])
    diff1 <- sum((Config != Config_new)[, 2:ncol(Config)])
    relax_rate <<- r_rate <- stats::rbeta(1, relax_rate_prior[1] + diff1,
                                          relax_rate_prior[2] + diff0)
    
    #Config_prior is the prior probability of C=1 given Config:=Omega
    Config_prior <- Config
    Config_prior[Config == 1] <- 1 - r_rate
    Config_prior[Config == 0] <- r_rate
    if (keep_base_clone) {
      Config_prior[, 1] <- Config[, 1]}
    Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
  }
  
  #Iden_mat is the assignment of cells to clones (through clusters)
  Iden_mat[,] <- 0
  for (j in seq_len(M)) {
    Iden_mat[j, assign_all[it,match(BCR[BCR$cell == colnames(A1)[j],]$cluster,clusters)]] <- 1 
  }
  
  # calculate log_probability matrix with genotype 0 and 1
  P0_mat <- A1 * log(theta0) + B1 * log(1 - theta0) + W_log
  P1_mat <- A1 * log(theta1) + B1 * log(1 - theta1) + W_log
  
  oddR_log <- P1_mat %*% Iden_mat  - P0_mat %*% Iden_mat 
  oddR_log <- oddR_log + Config_prior_oddlog
  oddR_log[which(oddR_log > 50)] <- 50
  oddR_log[which(oddR_log < -50)] <- -50
  Config_prob_tmp <- exp(oddR_log) / (exp(oddR_log) + 1)
  
  Config_new <<- C_n <- matrix(stats::rbinom(N*K, size = 1, Config_prob_tmp), nrow=N)
  
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
    idx <- which(assign_all_j[it,] == k)
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

sample_alpha_0 <- function(){
  # k is the number of non-empty clusters
  k <- 0
  for (q in seq_len(Q)) {
    if(length(t_[[q]]) > 0){
      k <- k+1
    }
  }
  #sample new value of sigma
  xSigma <- rbeta(1, alpha_0+1, M)
  p <- exp(log(k) - log(M) - log1p(-log(xSigma)) - log1p(k/(M*(1-log(xSigma)))))
  if(rbinom(1,1,p)){
    alpha_0 <<- rgamma(1, shape=k+1, rate=1-log(xSigma))
  }else{
    alpha_0 <<- rgamma(1, shape=k, rate=1-log(xSigma))
  }
}
