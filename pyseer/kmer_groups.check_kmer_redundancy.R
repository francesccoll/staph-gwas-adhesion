
# This R script is used to check redundancy (i.e. kmers present in the same set of samples) for specific sets of associated kmers

phenotypes_file = "/User/53.marta_zapotoczna_IE/pyseer/phenotypes_tested.adhesion.txt"
phenotypes = as.vector(as.matrix(read.delim(phenotypes_file, sep = "\n", header = F)))

input_dir = "/User/53.marta_zapotoczna_IE/pyseer/kmers/kmers_assoc_ann_picasso/";
output_dir = "/User/53.marta_zapotoczna_IE/pyseer_follow_up/kmer_groups/";

### Load kmers file

dataset="mzIE_phen_samples";
unit="kmers_30-100bp";
assoc="top"; # subset of kmers annotated
assoc="top2"; # subset of kmers annotated

for(p in 1:length(phenotypes))
{
  
kmers_file = paste(input_dir, paste(dataset ,phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"csv", sep = "."), sep = "")
if(file.exists(kmers_file)){
  kmers = read.delim(kmers_file, sep = "\t", header = T)
  dim(kmers)
}

grep_kmers = paste(paste(kmers$variant, " | ", sep = ""))
output_file = paste(output_dir, paste(dataset ,phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"grep_kmers.csv", sep = "."), sep = "")
if(!file.exists(output_file)){
  write.table(grep_kmers, file =  output_file, sep = "\n", col.names = F, row.names = F, quote = F)
}

}

### In bash terminal: grep original fsm kmer file to extract samples

# Run script: extract_kmers_from_fsm_file.sh.txt

### Load slice of fsm file with kmers of interest

dataset="mzIE_phen_samples";
unit="kmers_30-100bp";
assoc="top"; # subset of kmers annotated
assoc="top2"; # subset of kmers annotated

# for(p in 1:length(phenotypes))
# {

p=1; # Adhesion_to_Fg_1_25_ug_mL_absolute
p=2; # Adhesion_to_Fg_5_ug_mL_absolute
p=3; # Adhesion_to_Fg_D ~ 10000 kmers > top file to be used - 2024 kmers
p=4; # Adhesion_to_Fn_5_ug_mL_absolute
p=5; # Adhesion_to_Fn_20_ug_mL_absolute - 6482 kmers > top file to be used - 229 kmers
p=6; # Adhesion_to_Fn_D

print(phenotypes[p])

output_file = paste(output_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"fsm_kmers.30-100bp.txt", sep = "."), sep = "")
kmers2 = read.delim(output_file, sep = "\t", header = F)
print(dim(kmers2))

kmers_file = paste(input_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"csv", sep = "."), sep = "")
kmers = read.delim(kmers_file, sep = "\t", header = T)
print(dim(kmers))

kmers2[,1] = gsub(" ", "", kmers2[,1])
xxx = match(kmers2[,1], kmers$variant)
kmers2 = kmers2[which(!is.na(xxx)),]
print(dim(kmers2))


# simplifying fsm samples' format: sample_id:kmer_num > sample_id

for(k in 1:nrow(kmers2)){
  print(k)
  samples = unlist(strsplit(kmers2[k,2]," "))
  samples = samples[which(samples!="")]
  samples2 = samples
  for(s in 1:length(samples)){
    sample1 = samples[s];
    sample2 = unlist(strsplit(sample1,":"))[1]
    samples2[s] = sample2
  }
  kmers2[k,2] = paste(sort(unique(samples2)), collapse = " ")
}

# To add kmer ids
kmer_ids = paste("kmer_", seq(1, nrow(kmers2), 1), sep = "")
kmers3 = cbind(kmer_ids, kmers2)
colnames(kmers3)=c("kmer_id","kmer","samples")

output_file = paste(output_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"grep_kmers.rf.csv", sep = "."), sep = "")
write.table(kmers3, file =  output_file, sep = "\n", col.names = F, row.names = F, quote = F)


### Calculate proportion of samples' overlap between kmers: kmer similarity matrix

kmers_sim = mat.or.vec(nrow(kmers3), nrow(kmers3))
for(k1 in 1:nrow(kmers3)){
  print(paste(k1,"/",nrow(kmers3)))
  for(k2 in 1:nrow(kmers3)){
    samples1 = unlist(strsplit(kmers3[k1,"samples"]," "))
    samples2 = unlist(strsplit(kmers3[k2,"samples"]," "))
    samplesc = intersect(samples1, samples2)
    # denominator: kmer with higher AF/number of samples?
    denominator = length(samples1)
    if(length(samples2) > length(samples1)){ denominator = length(samples2); }
    numerator = length(samplesc)
    shared_sam_prop = numerator/denominator
    kmers_sim[k1, k2] = shared_sam_prop
  }
}

hist(c(kmers_sim))

### Grouping kmers by samples' overlap proportion

# See https://rdrr.io/github/rformassspectrometry/MsFeatures/man/groupSimilarityMatrix.html
require("MsFeatures")
kmer_groups = groupSimilarityMatrix(kmers_sim, threshold = 0.8, full = TRUE)
length(unique(kmer_groups))
# [1] 52

# Saving kmer group, along with original Pyseer kmer file
xxx = match(kmers3$kmer, kmers$variant)
identical(kmers$variant[xxx], kmers3$kmer)
kmers4 = cbind(kmers[xxx,], kmers3, kmer_groups)
dim(kmers4)
# [1] 586  15

xxx = match(kmers$variant, kmers3$kmer)
identical(kmers$variant, kmers3$kmer[xxx])
kmers4 = cbind(kmers, kmers3[xxx,], kmer_groups[xxx])
dim(kmers4)

output_file = paste(output_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,".kmer_group.csv", sep = "."), sep = "")
write.table(kmers4, file =  output_file, sep = "\t", col.names = T, row.names = F, quote = F)


## Adding original average phenotype values for all samples containing kmers

output_file = paste(output_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,".kmer_group.csv", sep = "."), sep = "")
kmers4 = read.delim(output_file, sep = "\t", header = T)
dim(kmers4)

phen_file = "/User/53.marta_zapotoczna_IE/pyseer/mzIE_phen_samples.kept.phen";
phen = read.delim(phen_file, sep = " ", header = T)
dim(phen)
# [1] 234   6

kmers4$samples_phen_mean = ""
kmers4$samples_phen_iqr = ""
for(k in 1:nrow(kmers4)){
  k_samples = unlist(strsplit(kmers4$samples[k], " "))
  sss = match(k_samples, phen$samples)
  # k_samples_phen = as.numeric(phen[sss[which(!is.na(sss))], phenotype_org])
  k_samples_phen = as.numeric(phen[sss[which(!is.na(sss))], phenotypes[p]])
  kmers4$samples_phen_mean[k] = mean(k_samples_phen)
  kmers4$samples_phen_iqr[k] = paste(quantile(k_samples_phen), collapse = " ")
}

output_file = paste(output_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"kmer_group.phen.csv", sep = "."), sep = "")
write.table(kmers4, file =  output_file, sep = "\t", col.names = T, row.names = F, quote = F)


# creating matrix of samples (as rows) and kmer groups (columns)

output_file = paste(output_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"kmer_group.phen.csv", sep = "."), sep = "")
kmers4 = read.delim(output_file, sep = "\t", header = T)
dim(kmers4)

kmer_groups = sort(unique(kmers4[,15])); # column 15 expected to have "kmer_group" id
length(kmer_groups)

kmers5 = mat.or.vec(nrow(phen)+1, length(kmer_groups) + 1)
dim(kmers5)
colnames(kmers5) = c("sample", paste("kmer_group_", kmer_groups, sep = ""))
kmers5[,1] = c("MRSA252_BX571856.fasta", phen[,1])

for(g in 1:length(kmer_groups)){
  print(g)
  kmer_group = kmer_groups[g]
  ggg = which(kmers4[,15]==kmer_group)
  kmer_group_samples = unique(unlist(strsplit(kmers4$samples[ggg], " ")))
  rrr = match(kmer_group_samples, kmers5[,1])
  kmers5[rrr,g+1] = kmer_group_samples
}

output_file = paste(output_dir, paste(dataset, phenotypes[p], unit, "pyseer_output.ann_kmers",assoc,"kmer_group.samples.csv", sep = "."), sep = "")
write.table(kmers5, file =  output_file, sep = "\t", col.names = T, row.names = F, quote = F)


