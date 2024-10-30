library("dndscv")
library("data.table")
library("vcfR")
library("deepSNV")
library("stringr")
library("dplyr")

basedir <- "~/Mounts/peter/liver/"

### collect wgs calls

project_name <- "haemochromatosis"
project_no = 2792

indir <- paste0(basedir,"data/",project_name,"/",project_no,"/beta_binomial/snv/ndp_out/")

outdir = paste0(basedir,"data/",project_name,"/",project_no,"/beta_binomial/snv/dnds_out/")

samples <- list.files(indir, pattern = "_ndp_assigned_snvs_and_indels.csv$")

for(i in 1:length(samples)){
  
  print(samples[i])
  donorID = substr(samples[i],1,7)
  
  mutations <- read.csv(paste0(indir,samples[i])) %>%
    filter(Mut_Frags > 1,
           !grepl("PDv38is", sampleID))
  
  write.csv(mutations,paste0(basedir,"data/",project_name,"/",project_no,"/beta_binomial/snv/ndp_out/",donorID,"_calls.csv"),row.names = FALSE)
  
  if(i == 1){
    wgs_muts <- mutations
  }else{
    wgs_muts <- rbind(wgs_muts, mutations)
  }
}

wgs_muts$Project <- 2792

write.csv(wgs_muts,"data/calls/wgs_calls_v2.csv",row.names = FALSE)


### collect exome calls
basedir <- "~/Mounts/peter/liver/"

project_nos = c(3131,2814)
project_names = c("A1AT","haemochromatosis")
projects <- data.frame(project_name = project_names, project_no = project_nos)

for(p in 1:nrow(projects)){

  project_no = projects$project_no[p]
  project_name = projects$project_name[p]

  print(paste0("reading project ",project_no))
  
indir <- paste0(basedir,"data/",project_name,"/",project_no,"/beta_binomial/")

samples <- list.files(paste0(indir,"snv/ndp_out"), pattern = "_bb_pass_snvs_all.csv$")
samples <- substr(samples,1,7)

for(i in 1:length(samples)){
  
  print(paste0("reading sample ",samples[i]))
  
  snvs <- read.csv(paste0(indir,"snv/ndp_out/",samples[i],"_bb_pass_snvs_all.csv"))
  indels <- read.csv(paste0(indir,"indel/ndp_out/",samples[i],"_bb_pass_indels_all.csv"))
  
  all_muts <- rbind(snvs,indels)
  write.csv(all_muts,paste0(indir,"snv/ndp_out/",samples[i],"_bb_pass_snvs_and_indels_all.csv"))
  
  ######
  donorID = substr(samples[i],1,7)
  
  mutations <- all_muts %>%
    filter(Mut_Frags > 1,
           sampleID != "PDv38is") %>%
    dplyr::select(-any_of(c("rho", "qval")))
  
  write.csv(mutations,paste0(basedir,"data/",project_name,"/",project_no,"/beta_binomial/snv/ndp_out/",donorID,"_calls.csv"),row.names = FALSE)
  
  if(i == 1){
    project_muts <- mutations
  }else{
    project_muts <- rbind(project_muts, mutations)
  }
}

project_muts$Project <- project_no

write.csv(project_muts,paste0("data/calls/",project_no,"_calls_v3.csv"),row.names = FALSE)

if(p == 1){
  exome_muts <- project_muts
}else{
  exome_muts <- rbind(exome_muts, project_muts)
}

}

write.csv(exome_muts,"data/calls/exome_calls_v3.csv",row.names = FALSE)

#### check if re-running 3131 resulted in other mutations being different than just those in PD60802b_lo0010 or PD60802b_lo0020

v2 <- read.csv("data/calls/exome_calls_v2.csv")
v3 <- read.csv("data/calls/exome_calls_v3.csv")

nrow(v2)
nrow(v3)

nrow(v3) - nrow(v2)

unique(v2[!v2$mutID_sampleID %in% v3$mutID_sampleID,"mut_id"])
unique(v3[!v3$mutID_sampleID %in% v2$mutID_sampleID,"mut_id"])
v3[!v3$mutID_sampleID %in% v2$mutID_sampleID,]
