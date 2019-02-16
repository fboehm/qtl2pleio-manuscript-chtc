
# load data
### genotype probabilities ("probs") in form used by R/qtl2
### see ../R/0_DOQTLprobs2qtl2.R for how they were converted
#PATH_TO_DERIVED_DATA <- "~/Box Sync/attie/attiedo"
PATH_TO_DATA <- "data"
#PATH_TO_DERIVED_DATA <- "~/attie"
#load(file.path(PATH_TO_DATA, "GM_Attie_allele_call_haploprobs_4qtl2_wave5.Rdata"))
### clinical phenotypes + phenotype dictionary
### ("pheno_clin" and "pheno_clin_dict")
#load(file.path(PATH_TO_DATA, "pheno_clin.RData"))
### kinship matrices ("loethod) ("K")
load(file.path(PATH_TO_DATA, "kinship_qtl2.RData"))
### covariate matrix ("covar")
#load(file.path(PATH_TO_DERIVED_DATA, "DerivedData/covar.RData"))
### physical map of the markers in the probs array
load(file.path(PATH_TO_DATA, "probs_pmap.RData"))
