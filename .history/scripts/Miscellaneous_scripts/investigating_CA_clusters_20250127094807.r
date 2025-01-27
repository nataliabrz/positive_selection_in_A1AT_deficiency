library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)

################################
#### Plot VAF distributions of the unmerged cluster pairs
################################

all_calls <- fread("data/calls/wgs_calls_PDv38is_wgs.csv", sep = ",") %>%
  dplyr::rename(chr = Chrom,
                pos = Pos,
                ref = Ref,
                mut = Alt) %>%
  dplyr::select(sampleID,chr,pos,ref,mut,mut_id,coverage,Mut_Frags,vaf,Project,cluster_id) 

PD51606_calls <- all_calls %>%
  filter(grepl("PD51606",sampleID)) 

cluster_pairs <- list(c(9, 51), c(8, 15), c(35, 7), c(22, 1), c(12, 18))

# Create an empty list to store the plots
plot_list <- list()

# Loop over each pair to create plots
for (pair in cluster_pairs) {
  
  # Filter data for samples containing at least one of the clusters in the pair
  pair_data <- PD51606_calls %>%
    filter(cluster_id %in% pair) %>%
    group_by(sampleID) %>%
    filter(any(cluster_id == pair[1]) | any(cluster_id == pair[2])) %>%
    ungroup()
  
  # Skip if no data meets the criteria
  if (nrow(pair_data) == 0) next
  
  # Further filter samples to include only if at least one cluster has median VAF > 0.20
  samples_to_plot <- pair_data %>%
    group_by(sampleID, cluster_id) %>%
    summarize(median_vaf = median(vaf)) %>%
    filter(median_vaf > 0.20) %>%
    pull(sampleID) %>%
    unique()
  
  # Filter pair_data to include only the selected samples
  pair_data <- pair_data %>%
    filter(sampleID %in% samples_to_plot)
  
  # Skip if no samples meet the criteria
  if (nrow(pair_data) == 0) next
  
  # Calculate the median VAF for each cluster in each sample
  medians <- pair_data %>%
    group_by(sampleID, cluster_id) %>%
    summarize(median_vaf = median(vaf)) %>%
    ungroup()

  # Define custom colors based on the cluster pair
  custom_colors <- setNames(c("blue", "red3"), as.character(pair))
  
  # Plot density with faceting by sampleID and cap x-axis at 0.5
  plot <- ggplot(pair_data, aes(x = vaf, color = as.factor(cluster_id), fill = as.factor(cluster_id))) +
    geom_density(alpha = 0.4) +
    scale_color_manual(values = custom_colors) +  # Apply custom colors
    scale_fill_manual(values = custom_colors) +   # Apply custom colors
    labs(title = paste("VAF Density for Clusters", pair[1], "vs", pair[2]),
         x = "VAF",
         y = "Density",
         color = "Cluster ID",
         fill = "Cluster ID") +
    xlim(0, 0.5) +  # Cap x-axis at 0.5
    theme_minimal() +
    facet_wrap(~sampleID) +  # Facet by sampleID for separate plots
    geom_vline(data = medians, aes(xintercept = median_vaf, color = as.factor(cluster_id)),
               linetype = "dashed", linewidth = 0.5, alpha = 0.5)  # Add dashed vertical lines at medians
  
  # Add the plot to the list
  plot_list[[paste0("Pair_", pair[1], "_vs_", pair[2])]] <- plot

   #pdf_filename <- paste0("outputs/unmerged_cluster_pairs_VAF_density_plots_med_vaf_0.05_", pair[1], "_vs_", pair[2], ".pdf")
   #ggsave(pdf_filename, plot, width = 10, height = 8)
}

# Combine all plots into one grid using gridExtra
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))

# Save the combined grid of plots into a single PDF
pdf("outputs/unmerged_cluster_pairs_VAF_density_plots_med_vaf_0.2.pdf", width = 12, height = 6)
grid.arrange(combined_plot)
dev.off()


################################
#### Check the number of samples that have each mutation, mean vafs etc
################################

cluster_12_muts <- PD51606_calls %>% filter(cluster_id == 12)

# Count how many samples each mutation appears in
mutation_counts <- cluster_12_muts %>%
  group_by(mut_id) %>%
  summarise(samples_count = n_distinct(sampleID)) %>%
  arrange(desc(samples_count))

# Show the top mutations (those appearing in most samples)
head(mutation_counts)

mean(mutation_counts$samples_count)
max(mutation_counts$samples_count)

# on average, a cluster 12 mutation is sean in 1.28 samples and maximum in 4 samples..


# Filter for mutations in PD51606_lo0065
pd51606_muts <- cluster_12_muts %>% filter(sampleID == "PD51606b_lo0065")

# Check VAF distribution of these mutations across other samples
high_vaf_muts <- pd51606_muts %>% filter(vaf > 0.2)  # Adjust VAF threshold as needed

# Find out how often these high VAF mutations appear at lower VAFs in other samples
other_samples_low_vaf <- cluster_12_muts %>%
  filter(mut_id %in% high_vaf_muts$mut_id & sampleID != "PD51606b_lo0065" & vaf < 0.1)

# View how many such mutations exist
table(other_samples_low_vaf$mut_id)
mean(table(other_samples_low_vaf$mut_id))

# View how many fragments these mutations are usually seen in
mean(other_samples_low_vaf$Mut_Frags)

mean(other_samples_low_vaf$vaf)

mean(other_samples_low_vaf$coverage)



# again, high vaf mutations in the main lcm sample are seen at low vaf usually in only 1 other sample
# Filter for mutations in PD51606_lo0065 and other samples in cluster 12
example_mutation <- cluster_12_muts %>%
  filter(mut_id %in% pd51606_muts$mut_id) %>%
  group_by(mut_id) %>%
  filter(n_distinct(sampleID) == 2)  # Find mutations seen in exactly two samples

# View a sample of these mutations
head(example_mutation |>
select(sampleID, chr, pos, ref, mut, mut_id, cluster_id) |>
filter(mut_id == "chr6_8786959_G_T"))


################################
### examine VAF differences
################################

# Define your cluster pairs
cluster_pairs <- list(c(9, 51), c(8, 15), c(35, 7), c(22, 1), c(12, 18))

# Placeholder for storing results
vaf_differences <- list()

# Loop through each pair in the cluster_pairs list
for (pair in cluster_pairs) {

  cluster_1 <- pair[1]
  cluster_2 <- pair[2]

  pair_data <- PD51606_calls %>%
    filter(cluster_id %in% pair) %>%
    group_by(sampleID) %>%
    filter(any(cluster_id == pair[1]) | any(cluster_id == pair[2])) %>%
    ungroup()
  
  # Further filter samples to include only if at least one cluster has median VAF > 0.20
  samples_to_plot <- pair_data %>%
    group_by(sampleID, cluster_id) %>%
    summarize(median_vaf = median(vaf)) %>%
    filter(median_vaf > 0.20) %>%
    pull(sampleID) %>%
    unique()
  
  # Filter pair_data to include only the selected samples
  pair_data <- pair_data %>%
    filter(sampleID %in% samples_to_plot)
  
  # Calculate the mean and median VAF for each cluster
  vaf_summary <- pair_data %>%
    group_by(cluster_id) %>%
    summarise(mean_vaf = mean(vaf, na.rm = TRUE),
              median_vaf = median(vaf, na.rm = TRUE))
  
  # Calculate the difference in mean and median VAFs between the clusters
  vaf_diff <- vaf_summary %>%
    summarise(mean_vaf_diff = diff(mean_vaf),
              median_vaf_diff = diff(median_vaf))
  
  # Store the results
  vaf_differences[[paste(cluster_1, cluster_2, sep = "_")]] <- vaf_diff
}

# View the results for all pairs
vaf_differences

################################
#### Collect bad mutations
################################

# Define your cluster pairs
cluster_pairs <- list(c(9, 51), c(8, 15), c(35, 7), c(22, 1), c(12, 18))

# Initialize an empty vector to store the bad mutations (mut_id with sample_id)
bad_muts_df <- data.frame(sampleID = character(), mut_id = character(), Mut_Frags = numeric(), vaf = numeric(), cluster_id = numeric(), stringsAsFactors = FALSE)

# Loop through each pair in the cluster_pairs list
for (pair in cluster_pairs) {
  # Extract the clusters for the pair
  cluster_1 <- pair[1]
  cluster_2 <- pair[2]
  
  pair_data <- PD51606_calls %>%
    filter(cluster_id %in% pair) %>%
    group_by(sampleID) %>%
    filter(any(cluster_id == pair[1]) | any(cluster_id == pair[2])) %>%
    ungroup()
  
  # Calculate the median VAF for each sample in the first cluster (cluster_1)
  vaf_summary <- pair_data %>%
    filter(cluster_id == as.character(cluster_1)) %>%
    group_by(sampleID) %>%
    summarise(median_vaf = median(vaf, na.rm = TRUE), .groups = "drop")
  
  # Identify the "bad" mutations where the median VAF is below 0.2 (considered artefactual)
  bad_samples <- vaf_summary %>%
    filter(median_vaf < 0.2) %>%
    pull(sampleID)
  
  # Collect the "sampleID" and "mut_id" for those bad mutations
  bad_muts_in_pair <- cluster_1_muts %>%
    filter(sampleID %in% bad_samples) %>%
    select(sampleID, mut_id, Mut_Frags, vaf, cluster_id)
  
  # Append the bad mutations to the bad_muts_df data frame
  bad_muts_df <- bind_rows(bad_muts_df, bad_muts_in_pair)
}
bad_muts_df

write.csv(bad_muts_df, "bad_mutations.csv", row.names = FALSE)

# Remove bad mutations from the dataset
cleaned_muts <- PD51606_calls %>% 
  filter(!paste(sampleID, mut_id, sep = "_") %in% bad_muts)

# View the bad mutations
bad_muts

# View the cleaned data (without the bad mutations)
cleaned_muts

#############################################

#############################################

all_calls <- fread("/lustre/scratch126/casm/team154pc/nb15/liver/data/HH_A1AD_rerun/2792/beta_binomial/snv/ndp_out/PD51606_bb_pass_snvs_all.csv",sep=",")
only_called <- fread("/lustre/scratch126/casm/team154pc/nb15/liver/data/HH_A1AD_rerun/2792/beta_binomial/snv/dnds_out/PD51606_bb_pass_snvs_dnds_input.csv",sep=",")

all_called <- all_called[all_called$Mut_Frags >1,]

not_called <- all_called[all_called$mutID_sampleID %in% setdiff(all_called$mutID_sampleID, only_called$mutID_sampleID), ]

write.csv(not_called, "not_called_but_rescued_mutations.csv", row.names = FALSE)

nrow(bad_muts_df)
nrow(bad_muts_df[paste0(bad_muts_df$sampleID,"_",bad_muts_df$mut_id) %in% not_called$mutID_sampleID,])


##### write a bed file with just cluster 12 mutations to test cgpvaf -bq flag:.groups

# Select relevant columns and rename them to match BED file format
bed_data <- cluster_12_muts[, .(chr, pos, ref, mut)]
setnames(bed_data, c("chr", "pos", "ref", "mut"), c("#chrom", "pos", "ref", "alt"))

bed_data <- bed_data[order(`#chrom`, pos)]

# Specify the output file path
output_file <- "/lustre/scratch126/casm/team154pc/nb15/liver/data/HH_A1AD_rerun/2792/Mutation_Collection_and_Filtering/240823/cgpvaf/SNV/snv_cgpvaf_PD51606_Cl.12.bed"

# Write to file, ensuring no row names, tab-delimited, and include header
fwrite(bed_data, output_file, sep = "\t", col.names = TRUE)


#############################################
#### check if the same sites appear across patients
#############################################

ca_clusters <- c(9, 8, 35, 22, 12)


# Get unique mutation sites in ca_mutations
unique_ca_mutations <- ca_mutations %>%
  distinct(chr, pos)

# Find matching sites in all_calls from other patients
common_mutations <- all_calls %>%
  filter((chr %in% unique_ca_mutations$chr) & (pos %in% unique_ca_mutations$pos)) %>%
  filter(!grepl("PD51606",sampleID))

# View the first few results
head(common_mutations)
nrow(common_mutations)
length(unique(paste0(common_mutations$chr,"_",common_mutations$pos)))