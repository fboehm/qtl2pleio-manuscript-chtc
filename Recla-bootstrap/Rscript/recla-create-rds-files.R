## GOALS: Prepare rds files for use with CHTC
# Need to prepare: 
# 1. recla.rds: from downloaded data (in qtl2 format)
# 2. recla-aprobs-chr8.rds: allele probabilities object for chromosome 8
# 3. recla-kinship.rds: loco kinship object

# Creating recla.rds:

library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla/recla.zip")
recla <- read_cross2(file)
saveRDS(object = recla, file = "../data/recla.rds")

# Creating recla-aprobs-chr8.rds

insert_pseudomarkers(recla, step = 0.10) -> pseudomap
probs <- calc_genoprob(recla, map = pseudomap)
aprobs <- genoprob_to_alleleprob(probs)
saveRDS(object = aprobs$`8`, file = "../data/recla-aprobs-chr8.rds")

# Creating recla-kinship.rds

kinship <- calc_kinship(aprobs, "loco")
saveRDS(object = kinship, file = "../data/recla-kinship.rds")

