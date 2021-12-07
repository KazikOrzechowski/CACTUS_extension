source('scRNA_sim.R')
source('cactusx_clone_assignment.R')

library('ggplot2')

#number of intermittent Gibbs sampling steps
i_G <- 3

########################################################################################
unknown <- simulate_t_and_cl()

t_true <- t_
clusters_true <- clusters_

t_ <- unknown$t
clusters_ <- unknown$clust

#########################################################################################
#only clustering test
true_ass <- I[t_true]

I <- rcat(length(clusters_), Psi)

old <- list()

for(it in 1:300){
  old[[str(it)]] <- I[t_]
  a <- try_split()
  b <- try_merge()
}

error_rate <- sapply(old, function(i){sum(i != true_ass)}) / M

print(error_rate)
###########################################################################3
#whole model test

assignments <- cactus_clone_assignment(A = A,
                                       D = D,
                                       Omega = Omega,
                                       BCR = BCR,
                                       n_chain = 1,n_proc = 30,
                                       relax_rate_prior = c(0.5, 9.5))     

error_rate <- sapply(1:20000, function(it){sum(I[t_] != assignments$assign_all_j[it,]) / M})

pdf(file='~/err_rate.pdf', width=4, height=4)
ggplot() + geom_line(aes(x=1:20000, y=error_rate))
dev.off()

pdf(file='~/C_prob.pdf', width=4, height=4)
ggplot() + geom_point(aes(x=1:(K*N), y=dbinom(C, 1, assignments$C_prob)))
dev.off()

View(assignments$C_prob)

x <- paste0('clone ', 1:K)
y <- 1:M
data <- expand.grid(X=x, Y=y)
data$prob <- as.vector(t(assignments$prob_mat_j[order(assignments$t),]))

pdf(file='~/ass_conf.pdf', width=4, height=4)
ggplot(data) + geom_tile(aes(X, Y, fill=prob))
dev.off()

inferred_B <- assignments$B[assignments$t,,]
true_B <- B[t_,,]

av_err <- sapply(1:M, function(cell){mean(abs(inferred_B[cell,,]-true_B[cell,,]))})

pdf(file='~/B_err.pdf', width=4, height=4)
ggplot() + geom_point(aes(x=1:M, y=av_err))
dev.off()

av_av_err <- mean(av_err)

View(B[1,,])
