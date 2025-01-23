library(data.table)
library(dplyr)

path <- "/lustre/scratch126/casm/team154pc/nb15/liver/"

source(
  "/lustre/scratch126/casm/team154pc/nb15/liver/code/pc8_Allocate_new_muts.R"
)

project_no <- 2792
project_name <- "HH_A1AD_rerun"
samples <- list.files(paste0(path, "data/", project_name, "/", project_no,
                             "/NDP"), pattern = "^PD")


# Write a function to find and remove and unwanted S0: annotations in
# complex indel instances
read_and_clean_indel_file <- function(file) {
  lines <- readLines(file)
  cleaned_lines <- gsub("(SO:[^,]+,)+", "complex,", lines)
  return(read.csv(text = cleaned_lines, sep = ",", header = TRUE))
}

for (i in seq_along(samples)) {
  ls_merge_out <- readRDS(paste0(path, "data/", project_name, "/", project_no,
                                 "/NDP/", samples[i],
                                 "/object_lymph.merge.ls.dat"))

  post_ls_dist <- readRDS(paste0(path, "data/", project_name, "/", project_no,
                                 "/NDP/", samples[i],
                                 "/object_post.cluster.pos.dat"))

  new_y <- fread(paste0(path, "data/", project_name, "/", project_no,
                        "/beta_binomial/indel/ndp_out/", samples[i],
                        "_ndp_alt_bb_flt.csv"), sep = ",")

  new_n <- fread(paste0(path, "data/", project_name, "/", project_no,
                        "/beta_binomial/indel/ndp_out/", samples[i],
                        "_ndp_depth_bb_flt.csv"), sep = ",")

  new_y <- as.matrix(new_y[, 7:length(colnames(new_y))])
  new_n <- as.matrix(new_n[, 7:length(colnames(new_n))])

  new_cluster <- allocate.new.muts.to.cluster(ls_merge_out, post_ls_dist,
                                              new_y, new_n)

  context <- read.table(paste0(path, "data/", project_name, "/", project_no,
                               "/beta_binomial/indel/ndp_out/", samples[i],
                               "_mut_context_GRCh38.txt"),
                        sep = "\t",
                        header = TRUE)

  context$clust_assign <- new_cluster
  context$mut_id <- rownames(context)
  clust_assign <- context[, c("clust_assign", "mut_id", "chrom", "pos")]

  write.table(clust_assign, paste0(path, "data/", project_name, "/",
                                   project_no, "/beta_binomial/indel/ndp_out/",
                                   samples[i], "_indels_to_clust_assign.csv"),
              sep = ",",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)

  file_indel <- paste0(path, "data/", project_name, "/", project_no,
                       "/beta_binomial/indel/ndp_out/", samples[i],
                       "_bb_pass_indels_all.csv")

  cluster_file <- paste0(path, "data/", project_name, "/", project_no,
                         "/beta_binomial/indel/ndp_out/", samples[i],
                         "_indels_to_clust_assign.csv")

  ### assign clusters to annotated mutations
  indel_tbl <- read_and_clean_indel_file(file_indel) %>%
    dplyr::mutate(pos_id = paste(Chrom, Pos, sep = "_"))

  cluster_tbl <- fread(cluster_file, sep = ",") %>%
    dplyr::mutate(pos_id = paste(chrom, pos, sep = "_")) %>%
    dplyr::select(-mut_id, -chrom, -pos) %>%
    distinct()

  assigned_tbl <- left_join(indel_tbl, cluster_tbl, by = "pos_id") %>%
    dplyr::rename("cluster_id" = clust_assign)

  write.csv(assigned_tbl, file = paste0(path, "data/", project_name, "/",
                                        project_no,
                                        "/beta_binomial/indel/ndp_out/",
                                        samples[i], "_ndp_assigned_indels.csv"),
            quote = FALSE,
            row.names = FALSE)

  ### combine indels and snvs into one file
  assigned_snvs <- read.table(paste0(path, "data/", project_name, "/",
                                     project_no, "/beta_binomial/snv/ndp_out/",
                                     samples[i],
                                     "_ndp_assigned_muts.csv"),
                              sep = ",",
                              header = TRUE)

  assigned_all <- rbind(assigned_snvs, assigned_tbl)

  sex_file <- fread(paste0(path, "data/", project_name, "/sex_file.txt"))
  sex <- sex_file[sex_file$patientID == samples[i], "sex"]

  if (sex == "F") {
    assigned_all <- assigned_all[assigned_all$Chrom != "chrY", ]
  }

  write.csv(assigned_all, file = paste0(path, "data/", project_name, "/",
                                        project_no,
                                        "/beta_binomial/snv/ndp_out/",
                                        samples[i],
                                        "_ndp_assigned_snvs_and_indels.csv"),
            quote = FALSE, row.names = FALSE)
}
