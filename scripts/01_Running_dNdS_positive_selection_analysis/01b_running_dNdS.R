library("dndscv")
library("vcfR")
library("deepSNV")
library("stringr")
library("dplyr")
library("data.table")

### collect calls

wgs_calls <- fread("data/calls/wgs_calls.csv", sep = ",") %>%
  mutate(
    biopsy = substr(sampleID, 1, 8),
    lcmID = sampleID,
    sampleID = paste0(biopsy, "_Cl.", cluster_id),
    mut_id_sample_id = paste0(biopsy, "_", mut_id)
  ) %>%
  dplyr::select(sampleID, lcmID, Chrom, Pos, Ref, Alt, mut_id_sample_id) %>%
  setNames(., c("sampleID", "lcmID", "chr", "pos", "ref", "mut", "mut_id_sample_id")) %>%
  mutate(chr = gsub("chr", "", chr)) %>%
  distinct()

exome_calls <- fread("data/calls/exome_calls.csv", sep = ",") %>%
  mutate(
    lcmID = sampleID,
    sampleID = substr(sampleID, 1, 8),
    mut_id_sample_id = paste0(sampleID, "_", mut_id)
  ) %>%
  dplyr::select(sampleID, lcmID, Chrom, Pos, Ref, Alt, mut_id_sample_id) %>%
  setNames(., c("sampleID", "lcmID", "chr", "pos", "ref", "mut", "mut_id_sample_id")) %>%
  filter(!mut_id_sample_id %in% wgs_calls$mut_id_sample_id) %>%
  mutate(chr = gsub("chr", "", chr)) %>%
  mutate(sampleID = substr(sampleID, 1, 8)) %>%
  distinct()


all_calls <- rbind(wgs_calls, exome_calls) %>%
  dplyr::select(sampleID, lcmID, chr, pos, ref, mut) %>%
  distinct()

# collapsing dinucleotides
mutations <- all_calls
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

print(paste0(nrow(all_calls) - nrow(mutations)," SNVs were collapsed into ",
             "dinucleotides"))
rm(all_calls)

#### Filter to only SNVs where both Ref and Alt are a single base
snv_calls <- mutations[
  nchar(ref) == 1 & nchar(mut) == 1
]

# Check for pairs of SNVs within <= 10 bp
close_snvs <- snv_calls[
  order(lcmID, chr, pos), # Ensure sorting
  .(sampleID, lcmID, chr, pos, next_pos = shift(pos, type = "lead")),
  by = .(lcmID, chr)
][
  !is.na(next_pos) & (next_pos - pos <= 10) # Filter for pairs <= 10 bp apart
]

close_snvs_nonadjacent <- close_snvs[
  (next_pos - pos) > 1, # Filter for pairs with distance > 1
  .(sampleID, lcmID, chr, pos, next_pos, distance = next_pos - pos)
]

# Print the results
close_snvs_nonadjacent

# Generate pos_ids by combining Chrom and Pos
close_snvs_nonadjacent[
  , pos_id := paste(chr, pos, sep = "_")
]

close_snvs_nonadjacent[, chr := sub("chr", "", chr)] 

# Count unique pos_ids
unique_pos_ids <- unique(close_snvs_nonadjacent$pos_id)
num_unique_pos_ids <- length(unique_pos_ids)

context <- fread("outputs/dnds/all_calls_mut_context_GRCh38.txt") %>%
  mutate(chr = gsub("chr", "", chr))

mutations <- mutations %>%
left_join(context, by=c("chr","pos","ref","mut"))

######
# Function to collapse MNVs
######

collapse_mnvs <- function(df, close_snvs) {
  for (i in seq_len(nrow(close_snvs))) {
    # Extract information about the current pair
    row <- close_snvs[i]
    sample <- row$lcmID
    chrom <- as.character(row$chr)
    pos1 <- row$pos
    pos2 <- row$next_pos
    distance <- row$distance
    
    # Find the mutations corresponding to the positions
    mut1 <- df[df$lcmID == sample & df$chr == chrom & df$pos == pos1]
    mut2 <- df[df$lcmID == sample & df$chr == chrom & df$pos == pos2]
    
    # Extract the ref and mut sequences
    ref <- paste0(mut1$ref, substr(mut1$CONTEXT, 11 + 1, 11 + distance - 1), mut2$ref)
    mut <- paste0(mut1$mut, substr(mut1$CONTEXT, 11 + 1, 11 + distance - 1), mut2$mut)
    
    # Create a new entry for the collapsed MNV
    new_entry <- data.table(
      chr = chrom,
      pos = pos1,  # Use the position of the first mutation as the new position
      lcmID = sample,
      sampleID = mut1$sampleID,
      ref = ref,
      mut = mut,
      CONTEXT = mut1$CONTEXT  # Keep the context of the first mutation
    )
    
    # Remove the original two SNVs from the dataframe
    df <- df[!(lcmID == sample & chr == chrom & pos %in% c(pos1, pos2))]
    
    # Add the new MNV to the dataframe
    df <- rbind(df, new_entry, fill = TRUE)
    df[lcmID == sample & chr == chrom & pos == pos1]
  }
  
  return(df)
}

# Use above function to collapse nearby snvs into mnvs
mutations_collapsed <- collapse_mnvs(mutations, close_snvs_nonadjacent)

print(paste0(nrow(mutations) - nrow(mutations_collapsed)," SNVs were collapsed into ",
             "multinucleotide variants"))

# Subset to unique mutations per clone (wgs) or donor (exome)
final_muts <- mutations_collapsed %>%
  dplyr::select(sampleID, chr, pos, ref, mut) %>%
  distinct()


write.table(x = final_muts, file = "outputs/dnds/combined_collapsed_dnds_input.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# loading dndscv covariats
covs <- file.path("data", "reference_files",
                  "covariates_20pc_GRCh37-38.epi_strict_outliers.Rdat")

load(covs) # it loads and object called scores
refcds_38 <- "data/reference_files/refcds_GRCh38-GencodeV18+Appris.rda"


## Running dnds separately on the haemochromatosis and A1AD cohort
cohorts <- c("A1AD", "HH")

for (i in seq_along(cohorts)) {
  print(paste0("Running dndscv on ", cohorts[i], " cohort"))

  if (cohorts[i] == "A1AD") {
    cohort_muts <- final_muts[grepl("PD51606|PD60232|PD60233|PD60802|PD60803",
                                    final_muts$sampleID), ]
  } else if (cohorts[i] == "HH") {
    cohort_muts <- final_muts[grepl("PD51605|PD51607|PD51608|PD52285|PD52286",
                                    final_muts$sampleID), ]
  }

  print("initializing...")
  dndsout_init <- dndscv(cohort_muts, refdb = refcds_38, cv = scores,
                         max_muts_per_gene_per_sample = Inf,
                         max_coding_muts_per_sample = Inf,
                         outmats = TRUE)

  # Excluding substitution drivers from the indel model
  print("running dndscv after exclusion of substitution drivers...")
  exc_ind_all <- dndsout_init$sel_cv$gene_name[which(dndsout_init$sel_cv$qallsubs_cv <= 0.05)]

  dndsout <- dndscv(cohort_muts, refdb = refcds_38, cv = scores,
                    max_muts_per_gene_per_sample = Inf,
                    max_coding_muts_per_sample = Inf,
                    outmats = TRUE,
                    kc = exc_ind_all,
                    onesided = TRUE)

  #saveRDS(dndsout, file = paste0("outputs/dnds/", cohorts[i], "_dndsout.rds"))
  write.table(dndsout$sel_cv,
              paste0("outputs/dnds/", cohorts[i], "_dndsout_selcv.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(dndsout$annotmuts, paste0("outputs/dnds/", cohorts[i], "_dndsout_annotmuts.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

  ## test for recurrent hotspots
  sitednds_liver <- sitednds(dndsout, method = "LNP")
  head(sitednds_liver$recursites)

  write.table(sitednds_liver$recursites, paste0("outputs/dnds/", cohorts[i], "_sitednds_recursites.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
}

# get confidence intervals for selected genes

for (i in 1:seq_along(cohorts)) {
  #dndsout <- readRDS(paste0("outputs/dnds/", cohorts[i], "_dndsout.rds"))
  geneci <- geneci(dndsout, gene_list = c("SERPINA1", "ACVR2A", "CIDEB", "GPAM", "FOXO1"), level = 0.95)
  write.table(geneci, paste0("outputs/dnds/", cohorts[i], "_per_gene_ci.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
}
