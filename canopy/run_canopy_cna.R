library(Canopy)
library(magrittr)
library(doParallel)

subject <- "K4B"
#subject <- "K5B"

reads_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K4B/canopy_input.csv"
#reads_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/canopy_input.csv"

falconx_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K4B/falconx_cna.rds"
#falconx_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K5B/falconx_cna.rds"

falconx_input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K4B/falconx_germline_input.csv"
#falconx_input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K5B/falconx_germline_input.csv"

canopy_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K4B/canopy_output.rds"
#canopy_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/canopy_output.rds"

tree_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K4B/output_tree.rds"
#tree_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/output_tree.rds"


sample <- read.csv( reads_save )
cna <- readRDS( falconx_save )
positions <- read.csv( falconx_input_save )[, c(1,2)]

# >>> map assumed copynumbers >>>

sample$MC <- sample$Mc <- 1

for(chr in names(cna)){
  ascn <- cna[[chr]]$ascn
  tauhat <- cna[[chr]]$tauhat
  
  n <- length( tauhat )
  
  ids <- sample$CHROM == chr
  pos_ids <- positions$CHROM == chr
  
  if( n < 1 ) { next }
  
  tauhat_pos <- positions[pos_ids, 2][tauhat]
  
  pos <- sample$POS < tauhat_pos[1]
  
  sample$Mc[ ids * pos ] <- ascn[1, 1]
  sample$MC[ ids * pos ] <- ascn[2, 1]
  
  if( n < 2 ) { next }
  
  for(i in 2:n){
    pos <- and(sample$POS < tauhat_pos[i], sample$POS > tauhat_pos[i-1])
    
    sample$Mc[ ids * pos ] <- ascn[1, i]
    sample$MC[ ids * pos ] <- ascn[2, i]
  }
  
  pos <- sample$POS > tauhat_pos[ n ]
  
  sample$Mc[ ids * pos ] <- ascn[1, n]
  sample$MC[ ids * pos ] <- ascn[2, n]
}

# <<< map assumed copynumbers <<<

R <- sample$BT %>% as.matrix
X <- sample$AT + sample$BT %>% as.matrix

rownames(X) <- rownames(R) <- 1:nrow(R)
colnames(X) <- colnames(R) <- subject

K <- 3:5

WM = as.matrix( sample$MC ) ## observed major copy number (for CNA regions)
Wm = as.matrix( sample$Mc ) ## observed minor copy number (for CNA regions)

rownames(WM) <- rownames(Wm) <- 1:nrow(R)
colnames(WM) <- colnames(Wm) <- subject

epsilonM <- epsilonm <- 0.01 ## standard deviation of WM and Wm, pre-fixed here

Y <- cbind( rep(0,nrow(X)) , diag(nrow(X)) )
rownames(Y) <- 1:nrow(R)
colnames(Y) <- c("non-cna_region", 1:nrow(R))


# >>> run canopy >>>
registerDoParallel(cores=5)

canopy_list <- foreach::foreach(i=1:5) %dopar% {
          library(Canopy)
          canopy.sample(R, X, 
                        WM, Wm, 
                        epsilonM, epsilonm, Y=Y,
                        K=K, numchain=1,
                        max.simrun = 100000,
                        min.simrun = 20000,
                        writeskip = 1000,
                        projectname = subject)
}

saveRDS(canopy_list, canopy_save)

bic = canopy.BIC(sampchain = canopy, projectname = subject, K = K,
                 numchain = 5, burnin = 0, thin = 1, pdf = TRUE)
optK = K[which.max(bic)]


post = canopy.post(sampchain = canopy, projectname = subject, K = K,
                   numchain = 5, burnin = 20, thin = 1, 
                   optK = optK)
samptreethin = post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
config = post[[3]]
config.summary = post[[4]]
print(config.summary)
# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space
# note: if modes of posterior probabilities aren't obvious, run sampling longer.


#######################################################
#######################################################
#######                                         #######
#######          Tree output and plot           #######
#######                                         #######
#######################################################
#######################################################
# choose the configuration with the highest posterior likelihood
config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree = canopy.output(post, config.i)
pdf.name = paste(subject, '_config_highest_likelihood.pdf', sep='')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)

C <- output.tree$Z

heatmap(C, Colv=NA, 
        labCol= c("Reference\n19.9%", "Clone 1\n17.8%", "Clone 2\n27.5%",
                  "Clone 3\n19.9%", "Clone 4\n14.9%")
)

saveRDS(output.tree, tree_save)
