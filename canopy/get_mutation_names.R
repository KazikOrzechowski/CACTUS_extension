library(magrittr)


# >>> parameters >>>

subject <- "K5B"

path_to_wes <- "C:/Users/kazik/outs/K5B/Strelka_S8934_1336260_FL_vs_S8934_FB_somatic_snvs.anno.vcf"

reads_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/canopy_input.csv"

tree_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/output_tree.rds"

save_path <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/tree_mutations.csv"

edges_save_path <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/tree_edges.csv"

mut_set_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/canopy/K5B/mutation_sets.csv"
# <<< parameters <<<


# >>> read in wes and choose only used variants >>>

wes_mutations <- read.delim(path_to_wes, 
                            header= FALSE, 
                            na.strings= ".", 
                            comment.char= "#")[, c(1,2,8)]

colnames( wes_mutations ) <- c("CHROM", "POS", "INFO")

reads <- read.csv( reads_save )


used <- merge( wes_mutations, reads, by = c("CHROM", "POS"))

# <<< read in wes and choose only used variants <<<


# >>> extract gene names >>>

info_list <- used$INFO %>% 
             strsplit( split='|', fixed=TRUE ) %>%
             sapply( function(entry) entry[4] )

used$GENES <- info_list
used$INFO <- NULL

write.csv( used, 
           save_path,
           row.names = FALSE )

# <<< extract gene names <<<


# >>> read canopy output tree >>>

tree <- readRDS( tree_save )

write.csv( tree[["edge"]], edges_save_path, row.names=FALSE)

# <<< read canopy output tree <<<

# >>> assign mutation sets >>>

mut_sets <- unique( tree[["sna"]][, c(2,3)]) %>% as.data.frame()
mut_sets$MUT_ID <- 1:nrow( mut_sets )

with_mut <- merge( tree[["sna"]], mut_sets, 
                   by=c("sna.st.node", "sna.ed.node"),
                   sort=FALSE)

used$SET <- with_mut$MUT_ID[ order(with_mut$sna) ]

write.csv( used, save_path, row.names=FALSE )
write.csv( mut_sets, mut_set_save, row.names=FALSE)
# <<< assign mutation sets <<<
