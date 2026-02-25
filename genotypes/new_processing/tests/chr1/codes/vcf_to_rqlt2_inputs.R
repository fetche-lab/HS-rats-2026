#!/usr/bin/env Rscript 

#install.packages("vcfR") 
library(vcfR) 
library(dplyr) 
library(stringr) 

#input 
setwd("/home/fetche-lab/Desktop/gn_remote/HS_rats/HS_data_v4/genotypes/new_processing/tests/chr1/results/intersect_output/")

hs_vcf <- "ratgtex_chr1.common.vcf.gz" 
founders_vcf <- "hs_founders.common.vcf.gz" 

pad_length <- 9 #max pos length 

#format marker id function 
format_marker_id <- function(pos) { 
	paste0("hsr", str_pad(pos, width = pad_length, pad= "0")) 
} 

start_time <- Sys.time() 

#load vcfs 
hs <- read.vcfR(hs_vcf, nrows=1000) 
founders <- read.vcfR(founders_vcf, nrows=1000)

#extract positions 
pos <- as.numeric(getFIX(hs)[,"POS"]) 
chr <- getFIX(hs)[,"CHROM"] 

marker_ids <- format_marker_id(pos) 

#create genetic map 
genetic_map <- data.frame(marker = marker_ids, chr=chr, cM = pos/1e6) #cM as approx map 

write.csv(genetic_map, "../rqlt2_inputs/hs_genetic_map.csv", row.names = F) 

#convert genotypes 
convert_gt <- function(gt) {
	gt <- sub(":.*", "", gt) 
	if (gt %in% c("0/0", "0|0")) return("AA") 
	if (gt %in% c("0/1", "1/0", "0|1", "1|0")) return("AB") 
	if (gt %in% c("1/1", "1|1")) return("BB") 
	return(NA) 
} 

##HS genotypes 
hs_gt <- extract.gt(hs) 
hs_df <- apply(hs_gt, c(1,2), convert_gt)
hs_df <- data.frame(marker = marker_ids, hs_df) 
write.csv(hs_df, "../rqlt2_inputs/hs_genotypes.csv", row.names = F) 

##Founder genotypes 
founders_gt <- extract.gt(founders)
founders_df <- apply(founders_gt, c(1,2), convert_gt) 
founders_df <- data.frame(marker = marker_ids, founders_df) 
write.csv(founders_df, "../rqlt2_inputs/hs_founders_genotypes.csv", row.names = F) 


#track runtime 
end_time <- Sys.time() 
runtime <- end_time - start_time 

print(paste("Total runtime:", runtime))
