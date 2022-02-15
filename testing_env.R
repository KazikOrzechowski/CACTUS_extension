source('BCR_sim.R')
source('cactusx_clone_assignment.R')

library('ggplot2')

###############################################################################
#set hyperparameters

save_path <- '~/'

#number of intermittent Gibbs sampling steps
i_G <- 5

max_iter <- 50000

new_sim <- FALSE
new_runs <- FALSE

n_runs <- 1

################################################################################
#data simulation
if(new_sim){    
  simulated_data <- data_simulation(av_reads=0.5)
}else{
  simulated_data <- readRDS(file='sim_data.RData')
}

###########################################################################
#observed variables
sim_A <- simulated_data$A
sim_D <- simulated_data$D
sim_Omega <- simulated_data$Omega
sim_BCR <- simulated_data$BCR

num_cells <- dim(sim_A)[2]
num_clones <- dim(sim_Omega)[2]
num_variants <- dim(sim_Omega)[1]

###########################################################################
#hidden variables
true_t <- simulated_data$t
true_B <- simulated_data$B[true_t,,]
true_I <- simulated_data$I
true_C <- simulated_data$C
true_t0 <- simulated_data$theta0
true_t1 <- simulated_data$theta1

##############################################################################
#run the model
if(new_runs){
  assignments <- cactus_clone_assignment(A = sim_A,
                                         D = sim_D,
                                         Omega = sim_Omega,
                                         BCR = sim_BCR,
                                         n_chain = n_runs,
                                         relax_rate_prior = c(0.5, 9.5),
                                         alpha_0=NULL,
                                         max_iter = max_iter,
                                         i_G= i_G)
}else{
  assignments <- readRDS(file='runs.RData')
}

###############################################################################
#save results and plots
saveRDS(simulated_data, file='sim_data.RData')
saveRDS(assignments, file="runs.RData")

for(run_ in 1:n_runs){
  error_rate <- sapply(1:max_iter, function(it){sum(true_I[true_t] != assignments[[run_]]$assign_all_j[it,]) / num_cells})
  
  pdf(file=paste0(save_path, str(run_), 'err_rate.pdf'), width=4, height=4)
  print(ggplot() + geom_line(aes(x=1:max_iter, y=error_rate)))
  dev.off()
  
  
  pdf(file=paste0(save_path, str(run_), 'C_prob.pdf'), width=4, height=4)
  print(ggplot() + geom_point(aes(x=1:(num_clones*num_variants), y=dbinom(true_C, 1, assignments[[run_]]$C_prob))))
  dev.off()
  
  
  x <- paste0('clone ', 1:num_clones)
  y <- 1:num_cells
  data <- expand.grid(X=x, Y=y)
  data$prob <- as.vector(t(assignments[[run_]]$prob_mat_j[order(assignments[[run_]]$t),]))
  
  pdf(file=paste0(save_path, str(run_), 'ass_conf.pdf'), width=4, height=4)
  print(ggplot(data) + geom_tile(aes(X, Y, fill=prob)))
  dev.off()
  
  
  inferred_B <- assignments[[run_]]$B[assignments[[run_]]$t,,]
  
  av_err <- sapply(1:num_cells, function(cell){mean(abs(inferred_B[cell,,]-true_B[cell,,])/true_B[cell,,])})
  av_abs_err <- sapply(1:num_cells, function(cell){mean(abs(inferred_B[cell,,]-true_B[cell,,]))})
  
  pdf(file=paste0(save_path, str(run_), 'B_err.pdf'), width=4, height=4)
  print(ggplot() + geom_point(aes(x=1:num_cells, y=av_err)))
  print(ggplot() + geom_point(aes(x=1:num_cells, y=av_abs_err)))
  dev.off()
  
  
  log_Lik <- assignments[[run_]]$logLik[2:max_iter]
  #log_Lik <- assignments[[run_]]$logLik[100*(1:(max_iter/100))]
  
  pdf(file=paste0(save_path, str(run_), 'logLik.pdf'), width=4, height=4)
  print(ggplot() + geom_line(aes(x=2:max_iter, y=log_Lik)))
  #ggplot() + geom_line(aes(x=1:(max_iter/100), y=log_Lik))
  dev.off()
  
  
  tst_df <- visualise_clusts(assignments$clusters, assignments$I, num_clones)
}


ggplot() + geom_line(aes(x=(.5*max_iter):max_iter, y=log_Lik[(.5*max_iter):max_iter]))
##################################################################################

Heatmap(as.matrix(gClonePortionS12118),name = paste0("Proportion"),col = f1,show_column_dend = FALSE,show_row_dend = FALSE,show_row_names = TRUE,show_column_names = TRUE,row_split =colnames(gClonePortionS12118)[gAssignedCloneS12118],row_order = order(rowSums(gCloneTableS12118),decreasing = TRUE),row_title =NULL,cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(if(gCloneTableS12118[i, j]!=0)sprintf("%i", as.matrix(gCloneTableS12118)[i, j]), x, y, gp = gpar(fontsize = 9))
},row_gap = unit(1, "mm"), border = TRUE ,row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),column_title_gp = gpar(fontsize = 9,fontface = "bold"),column_title = "c       CACTUS        ",column_order = colnames(gClonePortionS12118), width = unit(2.5, "cm"), height = unit(12, "cm"),show_heatmap_legend = FALSE)

visualise_clusts <- function(clusters_, I, K){
  non_empty_ind <- sapply(clusters_, function(cl){length(cl)>1})
  
  non_empty <- clusters_[non_empty_ind]
  
  assgn_df <- data.frame(row.names=1:length(non_empty))
  
  for(clone in 1:K){
    new_col <- paste('clone', clone)
    populations <- sapply(non_empty, length)
    assgn_df[,clone] <- ifelse(I[non_empty_ind]==clone, populations, 0)
  }
  
  assgn_df
}

clusts <- assignments$clusters

non_empty_ind <- sapply(clusts, function(cl){length(cl)>1})
non_empty <- clusts[non_empty_ind]

sums <- sapply(non_empty, function(cl){apply(sim_BCR[cl,,], c(2,3), sum)})

tst_df <- visualise_clusts(assignments[[1]]$clusters, assignments[[1]]$I, num_clones)

x <- paste0('clone ', 1:num_clones)
y <- 1:(dim(tst_df)[1])
data <- expand.grid(X=x, Y=y)
data$num <- t(as.matrix(tst_df))[1:length(rownames(data))]

ggplot(data) + geom_tile(aes(X,Y, fill=num))
