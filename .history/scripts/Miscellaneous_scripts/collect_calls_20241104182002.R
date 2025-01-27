library("dndscv")
library("data.table")
library("vcfR")
library("deepSNV")
library("stringr")
library("dplyr")

basedir <- "/lustre/scratch126/casm/team154pc/nb15/liver/"

### collect wgs calls

project_name <- "HH_A1AD_rerun"
project_no <- 2792

indir <- paste0(basedir, "data/", project_name, "/", project_no,
                "/beta_binomial/snv/ndp_out/")

outdir <- paste0(basedir, "data/", project_name, "/", project_no,
                 "/beta_binomial/snv/dnds_out/")

samples <- list.files(paste0(indir),
                      pattern = "_ndp_assigned_snvs_and_indels.csv$")

for (i in seq_along(samples)) {
  donor_id <- substr(samples[i], 1, 7)
  print(paste0("reading sample ", donor_id))

  mutations <- fread(paste0(indir, samples[i]), sep = ",") %>%
    filter(
      Mut_Frags > 1,
      !grepl("PDv38is", sampleID)
    )

  write.csv(mutations, paste0(basedir, "data/", project_name, "/", project_no,
                              "/beta_binomial/snv/ndp_out/", donor_id,
                              "_calls.csv"), row.names = FALSE)

  if (i == 1) {
    wgs_muts <- mutations
  } else {
    wgs_muts <- rbind(wgs_muts, mutations)
  }
}

wgs_muts$Project <- 2792

write.csv(wgs_muts, "data/calls/wgs_calls_PDv38is_wgs.csv", row.names = FALSE)


### collect exome calls
basedir <- "/lustre/scratch126/casm/team154pc/nb15/liver/"

project_nos <- c(3131, 2814)
project_names <- c("HH_A1AD_rerun", "HH_A1AD_rerun")
projects <- data.frame(project_name = project_names, project_no = project_nos)

for (p in seq_len(nrow(projects))) {
  project_no <- projects$project_no[p]
  project_name <- projects$project_name[p]

  print(paste0("reading project ", project_no))

  indir <- paste0(basedir, "data/", project_name, "/", project_no,
                  "/beta_binomial/")

  samples <- list.files(paste0(indir, "snv/ndp_out"),
                        pattern = "_bb_pass_snvs_all.csv$")

  samples <- substr(samples, 1, 7)

  for (i in seq_along(samples)) {
    print(paste0("reading sample ", samples[i]))

    snvs <- read.csv(paste0(indir, "snv/ndp_out/", samples[i],
                            "_bb_pass_snvs_all.csv"))
    indels <- read.csv(paste0(indir, "indel/ndp_out/", samples[i],
                              "_bb_pass_indels_all.csv"))

    all_muts <- rbind(snvs, indels)
    write.csv(all_muts, paste0(indir, "snv/ndp_out/", samples[i],
                               "_bb_pass_snvs_and_indels_all.csv"))

    ######
    donor_id <- substr(samples[i], 1, 7)

    mutations <- all_muts %>%
      filter(
        Mut_Frags > 1,
        !grepl("PDv38is", sampleID)
      ) %>%
      dplyr::select(-any_of(c("rho", "qval")))

    write.csv(mutations, paste0(basedir, "data/", project_name, "/",
                                project_no, "/beta_binomial/snv/ndp_out/",
                                donor_id, "_calls.csv"), row.names = FALSE)

    if (i == 1) {
      project_muts <- mutations
    } else {
      project_muts <- rbind(project_muts, mutations)
    }
  }

  project_muts$Project <- project_no

  write.csv(project_muts, paste0("data/calls/", project_no,
                                 "_calls_PDv38is_wgs.csv"), row.names = FALSE)

  if (p == 1) {
    exome_muts <- project_muts
  } else {
    exome_muts <- rbind(exome_muts, project_muts)
  }
}

write.csv(exome_muts, "data/calls/exome_calls_PDv38is_wgs.csv",
          row.names = FALSE)

#### check if re-running 3131 resulted in other mutations being different
#### than just those in PD60802b_lo0010 or PD60802b_lo0020

v2 <- read.csv("data/calls/exome_calls.csv")
v3 <- read.csv("data/calls/exome_calls_PDv38is_wgs.csv") %>%
  filter(
    !grepl("PD63757", sampleID),
    !grepl("PD64544", sampleID)
  )


nrow(v2)
nrow(v3)

nrow(v3) - nrow(v2)

unique(v2[!v2$mutID_sampleID %in% v3$mutID_sampleID, "mut_id"])
unique(v3[!v3$mutID_sampleID %in% v2$mutID_sampleID, "mut_id"])
v3[!v3$mutID_sampleID %in% v2$mutID_sampleID, ]

samples <- unique(v3$sampleID)

for (sample in samples) {
  nmut_bef <- length(unique(v2[v2$sampleID == sample, "mut_id"]))
  nmut_aft <- length(unique(v3[v3$sampleID == sample, "mut_id"]))
  print(paste0(sample, " has ", nmut_aft - nmut_bef,
               " extra mutations when called against PDv38is_wgs vs PDv38is"))
}
