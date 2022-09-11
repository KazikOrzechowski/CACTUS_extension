library(falconx)
library(magrittr)

#input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K4B/falconx_germline_input.csv"
input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K5B/falconx_germline_input.csv"
#input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K8B/falconx_germline_input.csv"

#output_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K4B/falconx_cna.rds"
output_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K5B/falconx_cna.rds"
#output_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K8B/falconx_cna.rds"

#################################################################################
# Apply FALCON-X to generate allele-specific copy number profiles
#################################################################################

ascn.input <- read.csv( input_save, head=T )
ascn.input[, c("mC", "MC")] <- NA

tauhat.list <- list()
cn.list <- list()

for(chr in unique(ascn.input$CHROM)){
  ids <- ascn.input$CHROM == chr
  
  chr.input <- ascn.input[ids, ]
  
  if( nrow(chr.input) < 2){ next }
  
  readMatrix <- chr.input[,c('AN','BN','AT','BT')]
  biasMatrix <- chr.input[,c('sN','sT')]
  
  try({
  tauhat <- getChangepoints.x( readMatrix, 
                               biasMatrix, 
                               pos=chr.input$POS)
  cn <- getASCN.x(readMatrix, 
                  biasMatrix, 
                  tauhat=tauhat, 
                  pos=chr.input$POS, 
                  threshold = 0.1)
  
  cn.list[[chr]] <- cn
  tauhat.list[[chr]] <- tauhat
  })
}


# cn$tauhat would give the indices of change-points.
# cn$ascn would give the estimated allele-specific copy numbers for each segment.
# cn$Haplotype[[i]] would give the estimated haplotype for the major chromosome in segment i
# if this segment has different copy numbers on the two homologous chromosomes.

# add log2ratio
ascn.input$log2ratio <- log2( ascn.input$MC / ascn.input$mC )

# save results
saveRDS( cn.list , output_save )
