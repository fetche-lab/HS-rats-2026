#!/usr/bin/env Rscript 
## hmm based haplotype reconstruction 
setwd("/home/fetche-lab/Desktop/gn_remote/HS_rats/HS_data_v4/genotypes/new_processing/tests/chr1/results/rqlt2_inputs/") 

#install.packages("qtl2") 
library(qtl2) 

cross <- read_cross2("hs_rn8_cross.yaml") 

pr <- calc_genoprob(cross, error_prob = 0.002) #hmm running  

apr <- genoprob_to_alleleprob(pr) 


