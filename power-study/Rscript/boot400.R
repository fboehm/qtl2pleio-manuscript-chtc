# translate command line args
##First read in the arguments listed at the command line
args <- commandArgs(TRUE)
print(args)
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(argname)
proc_num <- as.numeric(argname)
print(proc_num)
nsnp <- as.numeric(nsnp)
print(nsnp)
s1 <- as.numeric(s1)
print(s1)
nboot_per_job <- as.numeric(nboot_per_job)
print(nboot_per_job)
run_num <- as.numeric(run_num)
print(run_num)

###############
# first, load setup.R
source("Rscript/setup2.R")
#source("Rscript/setup-chr17-G83.R")
library(dplyr)
library(qtl2pleio)
load("probs_17.RData")
. -> probs

pp <- probs
pm <- pmap$`17`
kinship <- K$`17`

## Determine which markers are shared between pmap & gmap
snp_g <- dimnames(pp)[[3]]
snp_p <- names(pm)
shared_snps <- intersect(snp_g, snp_p)
pp2 <- pp[ , , snp_g %in% shared_snps]
pm2 <- pm[snp_p %in% shared_snps]

phenames <- c("sim1", "sim2")

nboot_per_pheno <- 400

(trait_file_num <- proc_num %% nboot_per_pheno)
 # here is where we assign the trait file name
(fn <- paste0("Ysim-run", run_num, "_", trait_file_num, ".txt"))

PATH_TO_SIM_DATA <- paste0("run", run_num, "-400sims")

as.matrix(read.table(file.path(PATH_TO_SIM_DATA, fn))) -> phe_pre
phe_pre[ , 1] -> phe1
phe_pre[ , 2] -> phe2

indic <- (is.na(phe1) | is.na(phe2))
(phe <- phe_pre[!indic, ])

rownames(pp2) %in% rownames(phe) -> pp_ind
pp2[pp_ind, , ] -> pp3
dim(pp3)

intersect(rownames(phe), rownames(pp3)) -> rn
pp3[rownames(pp3) %in% rn, , ] -> pp4
phe[rownames(phe) %in% rn, ] -> phe3
rownames(phe3) == rownames(pp4)
cbind(rownames(phe3), rownames(pp4))
pp4[order(rownames(pp4)), , ] -> pp5
phe3[order(rownames(phe3)), ] -> phe4
rownames(phe4) == rownames(pp5)

dim(pp5)
dim(phe4)
print(phe4)
rownames(kinship) %in% rownames(phe4) -> kin_ind
k2 <- kinship[kin_ind, kin_ind]
dim(k2)

k3 <- k2[order(rownames(k2)), order(rownames(k2))]

rownames(k3) == rownames(phe4)

(read.csv(file.path(PATH_TO_SIM_DATA, paste0("pleio-peak-indices-table-run", run_num, ".csv"))) -> pleio_peak_indices)

(fn_ind <- which(pleio_peak_indices$Yfns == fn))
(pleio_peak_index <- pleio_peak_indices$pleio_indices[fn_ind])
# simulate a phenotype
X1 <- pp5[ , , pleio_peak_index]
library(gemma2)
stagger_mats(X1, X1) -> X
# gemma2::stagger_mats(X1, X1) -> X
set.seed(proc_num)
calc_covs(pheno = phe4, kinship = k3) -> cc_out
(cc_out$Vg -> Vg)
(cc_out$Ve -> Ve)
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve, K = k3) -> Sigma
(solve(Sigma) -> Sigma_inv)
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe4)
# Start loop
lrt <- numeric()
for (i in 1:nboot_per_job){
  sim1(X = X, B = B, Vg = Vg, Ve = Ve, kinship = k3) -> foo
  matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
  scan_pvl(probs = pp5, pheno = Ysim, kinship = k3, start_snp1 = s1, n_snp = nsnp) -> loglik_mat
# in above call, s1 & nsnp come from command line args
  calc_lrt(loglik_mat) -> lrt[i]
}

fn_out <- paste0("boot400-run", run_num, "_", proc_num, ".txt")
write.table(lrt, fn_out)
q("no")


