library("dndscv")
library("vcfR")
library("deepSNV")
library("stringr")
library("dplyr") 

### collect calls

wgs_calls <- read.table("data/calls/wgs_calls.csv", sep = ",", stringsAsFactors = F, header = T) 

wgs_calls <- wgs_calls %>%
  mutate(biopsy=substr(sampleID,1,8),
         sampleID = paste0(biopsy,"_Cl.",cluster_id),
         mut_id_sample_id = paste0(biopsy,"_",mut_id)) %>%
  dplyr::select(sampleID,Chrom,Pos,Ref,Alt,mut_id_sample_id) %>%
  setNames(.,c("sampleID","chr","pos","ref","mut","mut_id_sample_id")) %>%
  mutate(chr=gsub("chr","",chr)) %>%
  distinct()

exome_calls <- read.table("data/calls/exome_calls.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(sampleID=substr(sampleID,1,8),
         mut_id_sample_id = paste0(sampleID,"_",mut_id)) %>%
  dplyr::select(sampleID,Chrom,Pos,Ref,Alt,mut_id_sample_id) %>%
  setNames(.,c("sampleID","chr","pos","ref","mut","mut_id_sample_id")) %>%
  filter(!mut_id_sample_id %in% wgs_calls$mut_id_sample_id) %>%
  mutate(chr=gsub("chr","",chr)) %>%
  mutate(sampleID=substr(sampleID,1,8)) %>%
  distinct()

final_muts <- rbind(wgs_calls,exome_calls) %>%
  dplyr::select(sampleID,chr,pos,ref,mut) %>%
  distinct()

# collapse dinucleotides
mutations <- final_muts
mutations <- mutations[order(mutations$sampleID,mutations$chr,mutations$pos),]
ind = which(diff(mutations$pos)==1)

dinuc <- mutations[FALSE,]

for(i in 1:length(ind)){
dinuc[i,] <- mutations[ind[i],]
dinuc[i,"ref"] <- paste0(mutations[ind[i],"ref"],mutations[ind[i]+1,"ref"])
dinuc[i,"mut"] <- paste0(mutations[ind[i],"mut"],mutations[ind[i]+1,"mut"])
}

mutations[ind,] <- dinuc
mutations <- mutations[-c(ind+1),]

final_muts <- mutations
rm(mutations)
########

write.table(x = final_muts, file = "combined_collapsed_dnds_input.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

covs = "data/reference_files/covariates_20pc_GRCh37-38.epi_strict_outliers.Rdat"
load(covs) # it loads and object called scores
refcds_38 = "data/reference_files/refcds_GRCh38-GencodeV18+Appris.rda"


## Run dnds separately on the haemochromatosis and A1AD cohort

cohorts=c("A1AD","HH")

for(i in 1:length(cohorts)){
  
  print(cohorts[i])
  
  if(cohorts[i] == "A1AD") {
    cohort_muts <- final_muts[grepl("PD51606|PD60232|PD60233|PD60802|PD60803", final_muts$sampleID), ]
  } else if(cohorts[i] == "HH") {
    cohort_muts <- final_muts[grepl("PD51605|PD51607|PD51608|PD52285|PD52286", final_muts$sampleID), ]
  }

  dndsout_init = dndscv(cohort_muts, refdb=refcds_38, cv=scores, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T)
  
  # Excluding substitution drivers from the indel model
  exc_ind_all <- dndsout_init$sel_cv$gene_name[which(dndsout_init$sel_cv$qallsubs_cv <= 0.05)]
  
  dndsout <- dndscv(cohort_muts, refdb = refcds_38, cv = scores, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, outmats = T, kc = exc_ind_all, onesided=T)
  
  saveRDS(dndsout, file=paste0("outputs/",cohorts[i],"_dndsout.rds"))
  write.table(dndsout$sel_cv,paste0("outputs/",cohorts[i],"_dndsout_selcv.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
  write.table(dndsout$annotmuts,paste0("outputs/",cohorts[i],"_dndsout_annotmuts.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
  
  ## test for recurrent hotspots
  sitednds_liver = sitednds(dndsout, method = "LNP")
  head(sitednds_liver$recursites)
  
  write.table(sitednds_liver$recursites,paste0("outputs/",cohorts[i],"_sitednds_recursites.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
  
}

# get confidence intervals for selected genes

for(i in 1:length(cohorts)){
  dndsout <- readRDS(paste0("outputs/",cohorts[i],"_dndsout.rds"))
  geneci <- geneci(dndsout, gene_list = c("SERPINA1","ACVR2A","CIDEB","GPAM","FOXO1"), level = 0.95)
  write.table(geneci,paste0("outputs/",cohorts[i],"_per_gene_ci.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
}
