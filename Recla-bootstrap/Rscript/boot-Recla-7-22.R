##
pleio_peak_index <- 773
##



library(qtl2pleio)
library(stringr)
library(dplyr)
library(readr)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
# translate command line args
print(args)
print(args$argname)
(proc_num <- as.numeric(args$argname))
(run_num <- as.numeric(args$run_num))
(nboot_per_job <- args$nboot_per_job) # think 1 boot per job
(nsnp <- as.integer(args$nsnp))
(s1 <- as.integer(args$s1))
###############
library(qtl2)
recla <- readRDS("data/recla.rds")
# make sex a covariate for use in pvl_scan
recla[[6]][ , 1, drop = FALSE] -> sex
sex <- sex == "female"
# insert pseudomarkers
insert_pseudomarkers(recla, step = 0.1) -> pseudomap

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
pp <- readRDS("data/recla-aprobs-chr8.rds")

## ------------------------------------------------------------------------
kinship <- readRDS("data/recla-kinship.rds")

## ------------------------------------------------------------------------
recla$pheno -> ph
log(ph) -> lph
apply(FUN = broman::winsorize, X = lph, MARGIN = 2) -> wlph

## ------------------------------------------------------------------------
phe <- wlph[, c(7, 22)]
gm <- pseudomap$`8`
k <- kinship[[8]]
## ------------------------------------------------------------------------
library(qtl2pleio)


###############
# simulate a phenotype
X1 <- pp[ , , pleio_peak_index]
cbind(X1, unlist(sex)) -> Xpre
## remove subjects with missing values of phenotype
is.na(phe[ , 1]) | is.na(phe[ , 2]) -> missing_indic
phe_nona <- phe[!missing_indic, ]
Xpre_nona <- Xpre[!missing_indic, ]
k_nona <- k[!missing_indic, !missing_indic]
##
gemma2::stagger_mats(Xpre_nona, Xpre_nona) -> X
set.seed(proc_num)
calc_covs(pheno = phe_nona, kinship = k_nona) -> cc_out
(cc_out$Vg -> Vg)
(cc_out$Ve -> Ve)
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve, K = k_nona) -> Sigma
solve(Sigma) -> Sigma_inv
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe_nona)
# Start loop
lrt <- numeric()
for (i in 1:nboot_per_job){
  sim1(X = X, B = B, Vg = Vg, Ve = Ve, kinship = k_nona) -> foo
  matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
  rownames(Ysim) <- rownames(phe_nona)
  colnames(Ysim) <- c("t1", "t2")
  scan_pvl(probs = pp[!missing_indic, , ], pheno = Ysim, covariates = sex[!missing_indic, , drop = FALSE], kinship = k_nona, start_snp1 = s1, n_snp = nsnp) -> loglik
# in above call, s1 & nsnp come from command line args
  calc_lrt_tib(loglik) -> lrt[i]
}

fn_out <- paste0("recla-boot-run", run_num, "_", proc_num, ".txt")
write.table(lrt, fn_out)
q("no")


