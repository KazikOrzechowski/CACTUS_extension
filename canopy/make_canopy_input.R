library(CODEX)
library(falconx)
library(magrittr)


#path_to_wes <- "C:/Users/kazik/outs/K4B/Strelka_S8934_66747_FL_vs_S8934_FB_somatic_snvs.anno.vcf"
path_to_wes <- "C:/Users/kazik/outs/K5B/Strelka_S8934_1336260_FL_vs_S8934_FB_somatic_snvs.anno.vcf"
#path_to_wes <- "C:/Users/kazik/outs/K8B/Strelka_GS104770_LN_vs_GS104770_FB_somatic_snvs.anno.vcf"

#path_to_pos <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/variants/K4B/intersection_pvalue_001.rds"
path_to_pos <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/variants/K5B/intersection_pvalue_001.rds"
#path_to_pos <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/variants/K8B/intersection_pvalue_001.rds"

#input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K4B/canopy_input.csv"
input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/canopy_input.csv"
#input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K8B/canopy_input.csv"


################# >>> helper functions >>>
extract_info <- function( row ){
  row_split <- strsplit( as.character(row), split = ":")
  
  #return only tier 1 reads
  tier_1 <- strsplit( row_split[[1]] , split = ",")
  
  #reads of each nucleotide
  AU <- tier_1[[5]][1]
  CU <- tier_1[[6]][1]
  GU <- tier_1[[7]][1]
  TU <- tier_1[[8]][1]
  
  result <- c( AU, CU, GU, TU )
  names( result ) <- c("A", "C", "G", "T")
  
  result
}


create_reads <- function( wes_mutations, normal, tumor ){
  n_positions <- nrow( wes_mutations )
  
  reads <- matrix( nrow= n_positions,
                   ncol= 4 )
  
  colnames( reads ) <- c( "AN", "BN", "AT", "BT" )
  
  for( i in 1:n_positions ){
    ref <- wes_mutations$REF[i]
    alt <- wes_mutations$ALT[i]
    
    reads[i, "AN"] <- normal[i, ref] %>% as.numeric
    reads[i, "BN"] <- normal[i, alt] %>% as.numeric
    reads[i, "AT"] <- tumor[i, ref] %>% as.numeric
    reads[i, "BT"] <- tumor[i, alt] %>% as.numeric
    }
  
  reads
}

################# <<< helper functions <<<


################# read data
selected_positions <- readRDS( path_to_pos )


wes_mutations <- read.delim(path_to_wes, 
                            header= FALSE, 
                            na.strings= ".", 
                            comment.char= "#")

wes_mutations <- wes_mutations[, c(1,2,4,5,10,11)]

colnames(wes_mutations) <- c("refContig",	"refPos",
                             "REF",	"ALT",
                             "NORMAL", "TUMOR")


################# take only selected positions
wes_mutations <- merge( wes_mutations, selected_positions, by=c("refContig", "refPos") )


################# extract tumor and normal allele counts
normal_split <- wes_mutations$NORMAL %>% sapply(extract_info) %>% t()
tumor_split <- wes_mutations$TUMOR %>% sapply(extract_info) %>% t()


################# create reads matrix
reads <- create_reads( wes_mutations, normal_split, tumor_split )

my_matrix <- wes_mutations[, 1:4] %>% as.matrix()
colnames( my_matrix ) <- c( "CHROM", "POS", "REF", "ALT" )


# total read depth
Y <- matrix( nrow=nrow(reads),
             ncol=2 ) 

Y[, 1] <- reads[,1] + reads[,2]
Y[, 2] <- reads[,3] + reads[,4]


# save falconx input
cbind( my_matrix, reads) %>% write.csv(file=input_save, row.names=FALSE)


