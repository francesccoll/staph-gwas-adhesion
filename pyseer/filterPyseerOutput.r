#!/usr/bin/env Rscript

library(optparse)

##############################################################################################################
######						                      PARSING ARGUMENTS AND OPTIONS			                           		######
##############################################################################################################

option_list = list(
  make_option(c("-p", "--pyseer_output"), type="character", action="store", default=NULL, help="Pyseer output file (required)"),
  make_option(c("-t", "--p_value_threshold"), type="double", action="store", default=NULL, help="P-value threshold to use (required). Pyseer variants with lrt-pvalue < threshold will be kept."),
  make_option(c("-o", "--output"), type="character", action="store", default=NULL, help="Output file name (required)."),
  make_option(c("-k", "--keep_top_variants"), type="logical", action="store_true", default=FALSE, help="Whether to keep top associated variants, see --top_percentage option too (optional)"),
  make_option(c("-c", "--top_percentile"), type="double", action="store", default=0.05, help="What percentage of top associated variants to keep, see --keep_top_variants option too (optional). [default= %default]")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


######################################################################################################################################
#####					                              READING INPUT FILES AND ARGUMENTS PASSED TO SCRIPT					                      	#####
######################################################################################################################################

pyseer_output = opt$pyseer_output;
pval_thres = as.numeric(opt$p_value_threshold)
output_file = opt$output;
keep_top_variants = opt$keep_top_variants;
top_percentile = opt$top_percentile;

if(!file.exists(pyseer_output)){ stop(paste(pyseer_output," not found",sep="")); }

print(paste("Reading file ", pyseer_output, sep = ""))
lmT = read.delim(pyseer_output,sep="\t",header=T)
print(paste("Info: ", nrow(lmT), " variants extracted from ", pyseer_output, sep = ""))

ldT = as.data.frame(lmT)

if(keep_top_variants){
  print(paste("Info: keeping top variants of percentile : ", as.character(top_percentile), sep = ""))
  pval_thres = quantile(ldT$lrt.pvalue, probs = top_percentile)
}

print(paste("Filtering Pyseer variants using ", as.character(pval_thres), " p-value threshold", sep = ""))
ldT2 = subset(ldT, lrt.pvalue < pval_thres)

print(paste("Info: ", nrow(ldT2), " variants kept in ", output_file, sep = ""))

print(paste("Saving output file ", output_file, sep = ""))
write.table(ldT2, file = output_file, sep = "\t", col.names = T, row.names = F, quote = F)



