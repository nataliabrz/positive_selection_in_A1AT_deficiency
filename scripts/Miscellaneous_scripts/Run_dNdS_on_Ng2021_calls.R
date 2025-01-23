library("dndscv")
library("vcfR")
library("deepSNV")
library("stringr")
library("dplyr")
library("data.table")

### collect calls

load("data/calls/x.snv.indel.matt.foad.RData")
wgs_calls <- x %>%
  rename(chr = chrom,
         lcmID = sample,
         sampleID = id,
         biopsy = donor,
         mut = alt) %>%
  mutate(
    mut_id = paste0(chr,"_",pos,"_",ref,"_",mut),
    mut_id_sample_id = paste0(biopsy, "_", mut_id),
    pos = as.numeric(pos)) %>%
  dplyr::select(sampleID, lcmID, chr, pos, ref, mut, mut_id_sample_id) %>%
  mutate(chr = gsub("chr", "", chr)) %>%
  dplyr::select(sampleID, lcmID, chr, pos, ref, mut) %>%
  distinct()

rm(x)

# collapsing dinucleotides
mutations <- wgs_calls
mutations <- mutations[order(mutations$lcmID,
                             mutations$chr, mutations$pos), ]
ind <- which(diff(mutations$pos) == 1)

dinuc <- data.table(sampleID = character(length(ind)),
                    lcmID = character(length(ind)),
                    chr = character(length(ind)),
                    pos = integer(length(ind)),
                    ref = character(length(ind)),
                    mut = character(length(ind)))

for (i in seq_along(ind)) {
  dinuc[i, ] <- mutations[ind[i], ]
  dinuc[i, "ref"] <- paste0(mutations[ind[i], "ref"],
                            mutations[ind[i] + 1, "ref"])
  dinuc[i, "mut"] <- paste0(mutations[ind[i], "mut"],
                            mutations[ind[i] + 1, "mut"])
}

mutations[ind, ] <- dinuc
mutations <- mutations[-c(ind + 1), ]

print(paste0(nrow(wgs_calls) - nrow(mutations)," SNVs were collapsed into ",
             "dinucleotides"))
rm(wgs_calls)



# Subset to unique mutations per clone (wgs) or donor (exome)
final_muts <- mutations %>%
  dplyr::select(sampleID, chr, pos, ref, mut) %>%
  distinct()


write.table(x = final_muts, file = "outputs/sld_combined_collapsed_dnds_input.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# loading dndscv covariats
covs <- file.path("data", "reference_files",
                  "covariates_20pc_GRCh37-38.epi_strict_outliers.Rdat")

load(covs) # it loads and object called scores
refcds_37 <- "data/reference_files/refcds_GRCh37-GencodeV18+Appris.rda"



  print("initializing...")
  dndsout_init <- dndscv(final_muts, refdb = refcds_37, cv = scores,
                         max_muts_per_gene_per_sample = Inf,
                         max_coding_muts_per_sample = Inf,
                         outmats = TRUE)

  # Excluding substitution drivers from the indel model
  print("running dndscv after exclusion of substitution drivers...")
  exc_ind_all <- dndsout_init$sel_cv$gene_name[which(dndsout_init$sel_cv$qallsubs_cv <= 0.05)]

  dndsout <- dndscv(final_muts, refdb = refcds_37, cv = scores,
                    max_muts_per_gene_per_sample = Inf,
                    max_coding_muts_per_sample = Inf,
                    outmats = TRUE,
                    kc = exc_ind_all,
                    onesided = TRUE)

  saveRDS(dndsout, file = paste0("outputs/SLD_dndsout.rds"))
  write.table(dndsout$sel_cv,
              paste0("outputs/SLD_dndsout_selcv.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(dndsout$annotmuts, paste0("outputs/SLD_dndsout_annotmuts.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

  ## test for recurrent hotspots
  sitednds_liver <- sitednds(dndsout, method = "LNP")
  head(sitednds_liver$recursites)

  write.table(sitednds_liver$recursites, paste0("outputs/SLD_sitednds_recursites.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# get confidence intervals for selected genes

  geneci <- geneci(dndsout, gene_list = c("SERPINA1", "ACVR2A", "CIDEB", "GPAM", "FOXO1"), level = 0.95)
  write.table(geneci, paste0("outputs/SLD_per_gene_ci.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

