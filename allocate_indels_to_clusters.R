
library(data.table)
library(dplyr)

source("~/Mounts/peter/liver/code/pc8_Allocate_new_muts.R")

project_no <- 2792
project_name <- "haemochromatosis"
sample <- c("PD51608")
path <- "~/Mounts/peter/liver/"

for(i in 1:length(sample)){
  
  ls.merge.out <- readRDS(paste0(path,"data/",project_name,"/",project_no,"/NDP/",sample[i],"/object_lymph.merge.ls.dat"))
  post.ls.dist <- readRDS(paste0(path,"data/",project_name,"/",project_no,"/NDP/",sample[i],"/object_post.cluster.pos.dat"))
  new.y <- read.table(paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/indel/ndp_out/",sample[i],"_ndp_alt_bb_flt.csv"),sep=",",header=T)
  new.N <- read.table(paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/indel/ndp_out/",sample[i],"_ndp_depth_bb_flt.csv"),sep=",",header=T)
  
  new.y <- as.matrix(new.y[,7:length(colnames(new.y))])
  new.N <- as.matrix(new.N[,7:length(colnames(new.N))])
  
  new.cluster <- allocate.new.muts.to.cluster(ls.merge.out, post.ls.dist, new.y, new.N)
  
  context <- read.table(paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/indel/ndp_out/",sample[i],"_mut_context_GRCh38.txt"),sep="\t",header=T)
  context$clust_assign <- new.cluster
  context$mut_id <- rownames(context)
  clust_assign <- context[,c("clust_assign","mut_id","chrom","pos")]
  
  write.table(clust_assign,paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/indel/ndp_out/",sample[i],"_indels_to_clust_assign.csv"),sep=",",col.names=T,row.names=F,quote=F)
  
  file.indel <- paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/indel/ndp_out/",sample[i],"_bb_pass_indels_all.csv")
  cluster.file <- paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/indel/ndp_out/",sample[i],"_indels_to_clust_assign.csv")
  
  ### assign clusters to annotated mutations
  indel_tbl <- fread(file.indel) %>%
    dplyr::mutate(pos_id = paste(Chrom, Pos, sep = "_"))
  cluster_tbl <- fread(cluster.file) %>%
    dplyr::mutate(pos_id = paste(chrom, pos, sep = "_")) %>%
    dplyr::select(-mut_id, -chrom, -pos)
  assigned_tbl <- left_join(indel_tbl, cluster_tbl, by = "pos_id") %>%
    dplyr::rename("cluster_id" = clust_assign)
  write.csv(assigned_tbl, file = paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/indel/ndp_out/",sample[i],"_ndp_assigned_indels.csv"), quote = F, row.names = F)
  
  ### combine indels and snvs into one file
  assigned_snvs <- read.table(paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/snv/ndp_out/",sample[i],"_ndp_assigned_muts.csv"),sep=",",header=T)
  assigned_all <- rbind(assigned_snvs,assigned_tbl)
  
  sex_file <- fread(paste0(path,"data/",project_name,"/sex_file.txt"))
  sex <- sex_file[sex_file$patientID==sample[i],"sex"]
  
  if(sex=="F"){
    assigned_all <- assigned_all[assigned_all$Chrom != "chrY",]
  }
  
  write.csv(assigned_all,file = paste0(path,"data/",project_name,"/",project_no,"/beta_binomial/snv/ndp_out/",sample[i],"_ndp_assigned_snvs_and_indels.csv"), quote = F, row.names = F)
  
}
