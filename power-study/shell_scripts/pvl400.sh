#!/bin/bash

# untar your R installation
tar -xzf R.tar.gz
tar -xzf SLIBS.tar.gz
# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export LD_LIBRARY_PATH=$(pwd)/SS:$LD_LIBRARY_PATH
# run R, with the name of your  R script
R CMD BATCH '--args argname='$1' nsnp='$2' s1='$3' run_num='$4'' Rscript/pvl400.R 'pvl400_run'$4'-'$1'.Rout'

