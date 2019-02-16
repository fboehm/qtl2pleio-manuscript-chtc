#!/bin/bash

# untar your R installation
tar -xzf R13.tar.gz
tar -xzf SLIBS.tar.gz
# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/SS:$LD_LIBRARY_PATH
# run R, with the name of your  R script
R CMD BATCH '--args argname='$1' nsnp='$2' s1='$3' nboot_per_job='$4' run_num='$5'' Rscript/boot-Recla-7-10.R 'boot400_run'$5'-'$1'.Rout'
