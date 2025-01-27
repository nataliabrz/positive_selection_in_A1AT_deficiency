library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Biostrings)
library(stringr)


wgs_calls <- fread("/lustre/scratch126/casm/team154pc/nb15/a1ad_paper/data/calls/wgs_calls.csv", sep = ",") %>%
  mutate(
    biopsy = substr(sampleID, 1, 8),
    sampleID = paste0(biopsy, "_Cl.", cluster_id),
    mut_id_sample_id = paste0(biopsy, "_", mut_id)
  ) %>%
  dplyr::select(sampleID, Chrom, Pos, Ref, Alt, mut_id_sample_id) %>%
  setNames(., c("sampleID", "chr", "pos", "ref", "mut", "mut_id_sample_id")) %>%
  distinct()


final_muts <- rbind(wgs_calls, exome_calls) %>%
  dplyr::select(sampleID, chr, pos, ref, mut) %>%
  distinct()

# get list of positions
this.pos.list = final_muts
this.pos.list = this.pos.list[, c("chr", "pos", "ref", "mut")]

# unique mutations
this.pos.list <- this.pos.list[!duplicated(this.pos.list),]

# define the context
context_list <- this.pos.list  %>%
  dplyr::mutate(start = pos - 10,
                end = pos + 10) %>%
  dplyr::select(chr, start,end)

# get the the sequences
this_range <- GenomicRanges::makeGRangesFromDataFrame(context_list)
this_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg37, this_range)
out_seqs <- GenomicRanges::as.data.frame(this_seq)
out.table <- bind_cols(this.pos.list, out_seqs)
colnames(out.table) <- c("chr", "pos", "ref", "mut", "CONTEXT")

# write the table to file
write.table(out.table, file = "/lustre/scratch126/casm/team154pc/nb15/a1ad_paper/data/calls/ng2021_calls_mut_context_GRCh37.txt", quote = F, row.names = F, sep = "\t")


