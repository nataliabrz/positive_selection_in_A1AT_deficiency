---
title: "Plot dN/dS result"
author: "Natalia Brzozowska"
date: "30 October, 2024"
output:
  html_document:
    fig_width: 8
    keep_md: true
    toc: true
    toc_float: true
    toc_collapsed: true
    toc_depth: 4
    theme: lumen
---




# Load data

Read the input files


``` r
refcds_38 <- "data/reference_files/refcds_GRCh38-GencodeV18+Appris.rda"
load(refcds_38)

selcvHH <- read.table("outputs/HH_dndsout_selcv.tsv",sep="\t",header=T)
selcvA1AD <- read.table("outputs/A1AD_dndsout_selcv.tsv",sep="\t",header=T)

outpath <- "outputs/"
```


We define a function to extract genomic coordinates from each element (gene) of RefCDS and apply it to RefCDS object. 


``` r
extract_info <- function(gene_info) {
  data.frame(
    gene_name = gene_info$gene_name,
    chr = gene_info$chr,
    pos = gene_info$intervals_cds[1, 1]
  )
}

gene_info_list <- lapply(RefCDS, extract_info)

gene_info_df <- do.call(rbind, gene_info_list)
colnames(gene_info_df) <- c("gene","chr", "pos")
```

Prepare data frame with dN/dS results for genes of interest in A1AD cohort.


``` r
geneLoc <- gene_info_df |> 
  group_by(gene, chr) |> 
  summarize(pos = min(pos), .groups = "drop")

sigEx <- 0.1

manhattan0 <- tibble(selcvHH) |> mutate(target = "HH") |> 
  bind_rows(tibble(selcvA1AD) |> mutate(target = "A1AD")) |> 
  dplyr::rename(gene = gene_name) |> 
  left_join(geneLoc, by = join_by(gene)) |> 
  mutate(pglobalpos_cv = if_else(pglobalpos_cv == 0, 1e-20, pglobalpos_cv)) |> 
  mutate(chr_ordering = str_replace(chr, "X", "23")) |> 
  mutate(chr_ordering = as.numeric(str_replace(chr_ordering, "Y", "24"))) 

data_cum <- manhattan0 %>% 
  group_by(chr_ordering, chr) %>% 
  summarise(max_bp = max(pos), .groups="drop") %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  dplyr::select(chr_ordering, chr, bp_add)

manhattan <- manhattan0 %>% 
  inner_join(data_cum, by = c("chr", "chr_ordering")) %>% 
  mutate(bp_cum = pos + bp_add) |> 
  mutate(label = if_else((qglobalpos_cv < sigEx & target == "HH") | (qglobalpos_cv < sigEx & target == "A1AD"), gene, ""))

axis_set <- manhattan %>% 
  group_by(chr_ordering, chr) %>% 
  summarize(center = mean(bp_cum), .groups = "drop")

ylim <- manhattan %>% 
  filter(pglobalpos_cv == min(pglobalpos_cv)) %>% 
  mutate(ylim = abs(floor(log10(pglobalpos_cv))) + 2) %>% 
  pull(ylim)


# Create a dataframe with genes to label
genes_to_label <- manhattan %>%
  filter((qglobalpos_cv < sigEx & target == "A1AD") | (gene %in% c("ALB","ACVR2A","FOXO1","CIDEB","GPAM","TNRC6B","SERPINA1","HFE") & target == "A1AD")) %>%
  dplyr::select(gene, chr_ordering, bp_cum, label, pglobalpos_cv, qglobalpos_cv)
```


``` r
pManhattanA1AD <- ggplot(manhattan |> filter(target == "A1AD"), aes(x = bp_cum, y = -log10(pglobalpos_cv), label = label,
                                  color = factor(chr_ordering), size = -log10(pglobalpos_cv))) +
  #geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
  geom_point() +
  geom_text_repel(data = genes_to_label, aes(label = gene, fontface = ifelse(qglobalpos_cv < sigEx, "bold", "plain")), 
                  colour = "black", size = 3, min.segment.length = 0) +  # Use the genes_to_label dataframe
  #geom_text_repel(colour = "black", size = 3) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#66c2a5", "#8da0cb"), 
                                  unique(length(axis_set$chr_ordering)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "dN/dS -log<sub>10</sub>(p)") + 
  ggtitle("Alpha-1 antitrypsin deficiency") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )
```

We do the same for the haemochromatosis cohort


``` r
sig <- 1e-05

genes_to_label <- manhattan %>%
  filter((qglobalpos_cv < sigEx & target == "HH") | (gene %in% c("ALB","ACVR2A","FOXO1","CIDEB","GPAM","TNRC6B","SERPINA1","HFE","FASN") & target == "HH")) %>%
  dplyr::select(gene, chr_ordering, bp_cum, label, pglobalpos_cv, qglobalpos_cv)
```


``` r
pManhattanHH <- ggplot(manhattan |> filter(target == "HH"), aes(x = bp_cum, y = -log10(pglobalpos_cv), label = label,
                                  color = factor(chr_ordering), size = -log10(pglobalpos_cv))) +
  #geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
  geom_point() +
  geom_text_repel(data = genes_to_label, aes(label = gene, fontface = ifelse(qglobalpos_cv < sigEx, "bold", "plain")), colour = "black", size = 3, min.segment.length = 0) +  # Use the genes_to_label dataframe
  #geom_text_repel(colour = "black", size = 3, ylim = c(-log10(sigEx),NA)) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#66c2a5", "#8da0cb"), 
                                  unique(length(axis_set$chr_ordering)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "dN/dS -log<sub>10</sub>(p)") + 
  ggtitle("Haemochromatosis") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  )
```

save the plot for both cohorts side by side


``` r
pManhattanA1AD + pManhattanHH + plot_annotation(tag_levels = 'a')
```

![](/Users/nb15/Documents/Nat_Gen_Revision/code/positive_selection_in_A1AT_deficiency/reports/A1AD_HH_dnds_manhattan_plot_Fig_1c_files/figure-html/pManhattan-1.png)<!-- -->
