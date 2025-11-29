#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(ggrepel)
library(tidyr)


##############################################################################################################
######						                      PARSING ARGUMENTS AND OPTIONS			                           		######
##############################################################################################################

option_list = list(
  make_option(c("-m", "--pyseer_kmer_map"), type="character", action="store", default=NULL, help="Pyseer kmer output file with kmers mapped to reference genomes (output of annotate_hits_pyseer-runner.py)"),
  make_option(c("-s", "--pyseer_kmer_sum"), type="character", action="store", default=NULL, help="Pyseer annotated kmers summarised by locus (output of summarise_kmers_pyseer.py)"),
  make_option(c("-t", "--p_value_threshold"), type="double", action="store", default=NULL, help="P-value threshold to be used (required). Use for horizontal dotted line."),
  make_option(c("-g", "--p_value_threshold_show"), type="double", action="store", default=NULL, help="P-value threshold to be used (required). The names of genes above this threshold are showed."),
  make_option(c("-o", "--output_plot_prefix"), type="character", action="store", default=NULL, help="Prefix used to name output files (required)."),
  make_option(c("-r", "--reference"), type="character", action="store", default=NULL, help="Reference genome to plot (required)."),
  make_option(c("-x", "--show_top_n_loci"), type="character", action="store", default="20", help="Number of top associated loci in reference genome to show.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


######################################################################################################################################
#####					                              READING INPUT FILES AND ARGUMENTS PASSED TO SCRIPT					                      	#####
######################################################################################################################################

pyseer_kmer_map = opt$pyseer_kmer_map;
pyseer_kmer_sum = opt$pyseer_kmer_sum;
p_value_threshold = as.numeric(opt$p_value_threshold)
p_value_threshold_show = as.numeric(opt$p_value_threshold_show)
output_plot_prefix = opt$output_plot_prefix;
reference_chosen = as.character(opt$reference)
length_reference = as.character(opt$length_reference)
show_top_n_loci = as.numeric(opt$show_top_n_loci)

print(paste("p_value_threshold chosen: ", p_value_threshold, sep = ""))
print(paste("p_value_threshold_show chosen: ", p_value_threshold_show, sep = ""))
print(paste("show_top_n_loci chosen: ", as.character(show_top_n_loci),  sep = ""))

if(!file.exists(pyseer_kmer_map)){ stop(paste(pyseer_kmer_map," not found",sep="")); }
if(!file.exists(pyseer_kmer_sum)){ stop(paste(pyseer_kmer_sum," not found",sep="")); }


######################################################################################################################################
#####						                                PREPARING VARIABLES FOR GGPLOT                                                 #####
######################################################################################################################################
# mapped_kmers = read.delim(pyseer_kmer_map, sep = "\t", header = F) # to be used if .map_kmers.csv files, which lack headers
mapped_kmers = read.delim(pyseer_kmer_map, sep = "\t", header = T)
# colnames(mapped_kmers) = c("variant", "af", "filter.pvalue", "lrt.pvalue", "beta", "beta.std.err", "variant_h2", "lineage", "notes", "mappings"); # to be used if .map_kmers.csv files, which lack headers
print(paste("Number of kmers loaded from ", pyseer_kmer_map, ": ", nrow(mapped_kmers)))

# Extracting reference and pos for kmers
# NOTE: mapping information can be expected in column 10: e.g. MRSA252:748155-748254;;;

mapped_kmers$median_pos = NA
mapped_kmers$reference = NA
mapped_kmers$num_mappings = 1

# e.g. MRSA252:748155-748254;;;
extract_kmer_reference = function(x){
  tmp = unlist(strsplit(x, ":")); reference = as.character(tmp[1]); return(reference);
}

extract_median_kmer_pos = function(x){
  tmp = unlist(strsplit(x, ";")); tmp2 = unlist(strsplit(tmp[1], ":")); positions = unlist(strsplit(tmp2[2], "-"));
  median_pos = median(seq(as.numeric(positions[1]), as.numeric(positions[2]), by=1))
  return(median_pos)
}

mapped_kmers$reference = sapply(mapped_kmers[,10], extract_kmer_reference)
mapped_kmers$median_pos = sapply(mapped_kmers[,10], extract_median_kmer_pos)

# mapped_kmers$reference = sapply(mapped_kmers[,9], extract_kmer_reference)
# mapped_kmers$median_pos = sapply(mapped_kmers[,9], extract_median_kmer_pos)

print("Top ten reference/contigs with mapped kmers: ")
print(sort(table(mapped_kmers$reference), decreasing = T)[1:10])

# PROCESSING KMERS mapping to multiple locations, a new line with such kmer information will be added
# Slow: it can be commented out

# kmers_mul_idx = which(grepl(",", mapped_kmers[,10])==TRUE)
# mapped_kmers_mul = mapped_kmers[kmers_mul_idx,]
# 
# for(k in 1:nrow(mapped_kmers_mul)){
#   print(k)
#   mappings = unlist(strsplit(mapped_kmers_mul[k,10], "\\,"))
#   # keeping first kmer mapping in original line/row
#   mapped_kmers[kmers_mul_idx[k],10] = mappings[1]
#   mapped_kmers[kmers_mul_idx[k],"num_mappings"] = as.character(length(mappings))
#   for(m in 2:length(mappings)){
#    new_kmer_row =  mapped_kmers_mul[k,]
#    new_kmer_row[10] = mappings[m]
#    new_kmer_row["reference"] = extract_kmer_reference(mappings[m])
#    new_kmer_row["median_pos"] = extract_median_kmer_pos(mappings[m])
#    new_kmer_row["num_mappings"] = as.character(length(mappings))
#    mapped_kmers = rbind(mapped_kmers, new_kmer_row)
#   }
# }


######################################################################################################################################
#####						                                   PREPARING LOCI LABELS TO PLOT                                               #####
######################################################################################################################################

sum_kmers = read.delim(pyseer_kmer_sum, sep = "\t", header = T)
dim(sum_kmers)
print(paste("Number of loci with associated kmers loaded from ", pyseer_kmer_sum, ": ", nrow(sum_kmers)))

print("Top ten reference/contigs with mapped kmers: ")
print(sort(table(sum_kmers$ref_id), decreasing = T)[1:10])

# Editing locus plot label to include gene name if available, and locus_tag if gene name not available
# Same code of block as that in plot_loci_kmers_pyseer.R

gene_label = rep("", nrow(sum_kmers))
for(r in 1:nrow(sum_kmers)){
  new_label = "not_assigned"
  # if intergenic region
  if(grepl("\\|", sum_kmers$gene[r])){
    locus_tmp = unlist(strsplit(sum_kmers$locus_tag[r],"\\|")); genes_tmp = unlist(strsplit(sum_kmers$gene[r],"\\|"));
    label_left = genes_tmp[1];
    if(label_left==""){ label_left = locus_tmp[1]; }
    if(length(genes_tmp)==2){ label_right = genes_tmp[2]; } else {label_right = ""; }
    if(label_right==""){ label_right = locus_tmp[2]; }
    new_label = paste(label_left, label_right, sep = "-")
  } else
  { if(sum_kmers$gene[r] == "-"){ new_label = sum_kmers$locus_tag[r]; } else { new_label = sum_kmers$gene[r]; } }
  new_label = gsub("\\|", "\\-", new_label)
  gene_label[r] = new_label
}
which(gene_label == "not_assigned")
# integer(0)
sum_kmers$gene_label = gene_label


# Keeping only kmers mapping to chosen reference genome
mapped_kmers = subset(mapped_kmers, reference == reference_chosen)
sum_kmers = subset(sum_kmers, ref_id == reference_chosen)
sum_kmers = sum_kmers[order(sum_kmers$min_pval),]

print(paste("Number of kmers kept for reference ", reference_chosen, ": ", nrow(mapped_kmers)))
print(paste("Number of loci with associated kmers kept for ", reference_chosen, ": ", nrow(sum_kmers)))

# Making both dataframes (kmers and loci) share the same column names: needed for plotting
# mapped_kmers = mapped_kmers[,c("median_pos","lrt.pvalue","beta","af","mappings")]
mapped_kmers = mapped_kmers[,c("median_pos","lrt.pvalue","beta","af","kmer_map")]; # when using .ann_kmers.csv
sum_kmers = sum_kmers[,c("median_pos","min_pval","mean_beta","mean_af","gene_label")]
colnames(sum_kmers) = c("median_pos","lrt.pvalue","beta","af","mappings")

# Keeping loci to label with top associated kmers, based on option show_top_n_loci
# If p_value_threshold_show is chosen, then gene labels below this threshold are chosen
if(!is.null(p_value_threshold_show)){
  show_rows = which(as.numeric(sum_kmers[,"lrt.pvalue"]) < p_value_threshold_show)
  sum_kmers_show = sum_kmers[show_rows,]
  print(paste("Showing ", length(show_rows), " gene labels below p_value_threshold_show ", p_value_threshold_show,sep = ""))
} else {
  last_row = show_top_n_loci
  if(nrow(sum_kmers) < show_top_n_loci){ last_row = nrow(sum_kmers); }
  sum_kmers_show = sum_kmers[1:last_row,]
}

######################################################################################################################################
#####						                                    CREATING MANHATTAN PLOT USING GGPLOT						                          	#####
######################################################################################################################################

plot = ggplot(mapped_kmers, aes(x=median_pos, y=-log10(lrt.pvalue), colour=beta, size=af)) + 
  geom_point(alpha=1) +
  scale_size("Kmer frequency", range=c(1,5)) + 
  scale_colour_gradient('beta') +
  theme_bw(base_size=14) +
  ggtitle(paste("Pyseer associated kmers mapped to reference ", reference_chosen, sep = "")) +
  xlab("Nucleotide position") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept=-log10(p_value_threshold), linetype="dashed") +
  geom_text_repel(data=sum_kmers_show,aes(x=median_pos,y=-log10(lrt.pvalue),label=mappings))

output_plot_name = paste(output_plot_prefix, ".", reference_chosen, ".plot.pdf", sep = "")
ggsave(output_plot_name, plot = plot, device = "pdf", width = 14, height = 7, dpi = 300, units = "in")


