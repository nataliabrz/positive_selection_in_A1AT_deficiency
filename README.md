# positive_selection_in_A1AT_deficiency

This is the code to accompany the paper 'Selection of somatic escape variants in SERPINA1 in liver tissue of patients with alpha-1 anti-trypsin deficiency'.

This code allows reproduction of all analyses from filtered and processed substitution (Caveman) and indel (Pindel) calls to generation of figures used in the manuscript.

## Repository Overview

.
├── data/                     # Input data files (e.g., mutation calls, metadata)
├── scripts/                  # Analysis scripts for processing and generating plots
├── outputs/                  # Generated analysis results and figures
├── liver_ndp_tree_building/  # Submodule: Pipeline for phylogenetic tree building (Figure 2a)
├── README.md                 # This document

---

## Data Availability 

Due to the sensitive nature of the raw sequencing data, whole genome and exome sequencing data have been deposited in the European Genome–phenome Archive (EGA):

- WGS data: EGAD00001015430 (https://ega-archive.org/studies/EGAD00001015430)
- Exome data: EGAD00001015431 (https://ega-archive.org/studies/EGAD00001015431)
- Additional datasets for steatotic liver disease: EGAD00001006255 (https://ega-archive.org/studies/EGAD00001006255)

Filtered somatic mutation calls used in downstream analyses have been deposited in Mendeley Data and can be freely accessed here: https://data.mendeley.com/datasets/vhybvj2g9p/1. These files are provided as the starting point for all analyses in this repository.


## Data Processing

The bam files were processed to identify somatic mutations using the following pipeline:

1. **Variant Calling**  
   - **SNVs**: Detected with CaVEMan ([Jones *et al.*, 2016](https://doi.org/10.1002/cpbi.20)) using the following parameters:  
     - `-tum-cn_default 10`  
     - `-norm-cn-default 2`  
     - `-normal-contamination 0.1`  
   - **Indels**: Identified with cgpPindel ([Raine *et al.*, 2015](https://doi.org/10.1002/0471250953.bi1507s52)).

2. **VCF Reflagging**  
Both CaVEMan and cgpPindel were originally designed to detect clonal mutations in tumor samples. To accommodate the low-input DNA pipeline (e.g., laser capture microdissection [LCM] samples), the variant calling thresholds and filters were re-adjusted to interpret the data effectively.  

   - **For SNVs**:  
     - The **MNP flag (Matched normal proportion)** is turned off.  
       
   - **For Indels**:  
     - The **FF016** flag is removed to disable the filter requiring a minimum of 5 supporting reads.  
     - The **FF018** flag is removed to disable the filter requiring a minimum depth of 10 in both tumor and normal samples.  

3. **Bulk-convert VCF**  
   - The VCF files are combined and converted into tabular format. 

4. **Custom Filtering for Low-Input Library Prep Artefacts**  
   - A custom filter is applied to remove artifacts specifically associated with the low-input library prep process, targeting those introduced by cruciform DNA structures. The method is described in [SangerLCMFiltering](https://github.com/MathijsSanders/SangerLCMFiltering).

5. **Additional Filtering**  
   Further filtering is applied to refine the variant calls. The following filters are applied to SNVs and indels:

  - **For SNVs**:  
     - **ASRD Filter**: Requires a median (read length adjusted) alignment score of variant-supporting reads >= 0.87.  
     - **CLPM Filter**: Requires a median number of soft-clipped bases in variant-supporting reads of 0.  
     - **MNP Filter**: Requires normal sample mutant allele proportion <= 0.1 and tumor sample mutant allele proportion minus normal sample mutant allele proportion >= 0.05 (allowing lower VAF than the original filter).  

   - **For Indels**:  
     - **FF016 Filter**: Requires at least 2 fragments and bi-directional support.  
     - **FF018 Filter**: Requires at least 6x depth in both tumor and normal samples.  
     - **HOMO Filter**: Removes 1bp indels at homopolymer runs of 6 or more. Any indel in a repeat of greater than 10 is removed.  

7. **Run cgpVaf**  
   - cgpVaf is used to force-call the SNVs and indels across all samples from the same patient, using a cutoff for read mapping quality (30) and base quality (25) to extract VAF values from the BAM files at specific genomic positions defined in a BED file.

8. **Exact Binomial Filter**  
   - An exact binomial filter is applied to aggregated counts of normal and variant reads across all samples. This filter removes sites where the count distributions are consistent with germline SNPs, retaining only somatic mutations.

9. **Beta-Binomial Filter**  
   - A beta-binomial filter is applied to retain only mutations whose count distributions across samples are consistent with an over-dispersed beta-binomial distribution, indicating the mutation is a true somatic variant rather than a low-level sequencing artefact.

For more details on the specific tools and parameters, please contact us or refer to the relevant publications.

You can adapt this stage to your local setup and the mutation-calling algorithm you're most familiar with (e.g., Mutect2, Strelka, VarScan). As long as these tools are run in 'unmatched' mode, their outputs will be compatible with the rest of the analysis pipeline.


## Objects that need to be downloaded from Mendeley Data

These objects should be downloaded and placed in the cloned repository with the following file structure:

data/calls/exome_calls.csv

data/calls/wgs_calls.csv


# Notes on downsteam analysis

## 01 Running dN/dS to test for positive selection

- Script: scripts/01_dnds_analysis.R

This is a script which takes the calls saved within the data/calls/ folder from Mendeley data and performs positive selection analysis using the dndscv R package.

### 02: Generating Plots

The following scripts reproduce the figures in the manuscript:

- Figure 1b: scripts/fig1b_plot.R
- Figure 1c: scripts/fig1c_plot.R
- Figure 2a: liver_ndp_tree_building/ (submodule; see below for details).
- Ext. Figure 2a: scripts/ext_fig2a_plot.R
- Ext. Figure 3: scripts/ext_fig3_plot.R

All scripts save the generated plots in the outputs/ directory.

---

## Using the Liver Tree-Building Submodule (Figure 2a)

The liver_ndp_tree_building submodule is used to generate phylogenetic trees from NDP clustering results.

### Setting Up the Submodule

To clone the repository with the submodule:
git clone --recurse-submodules https://github.com/yourusername/paper-repo.git

If you’ve already cloned the repository, initialize and update the submodule:
git submodule update --init --recursive

Navigate to the submodule:
cd liver_ndp_tree_building

Follow the instructions in the liver_ndp_tree_building/README.md to run the pipeline.

---

## Miscellaneous_scripts

The 'Run_dNdS_on_Ng2021_calls.R' script reproduces dNdS analysis on previously published 'Convergent somatic mutations in metabolism genes in chronic liver disease' article. The data supporting the paper can be downloaded from Mendeley at https://data.mendeley.com/datasets/283gy325fk/1. To run the script, place the downloaded object in the cloned repository with the following file structure:

data/calls/x.snv.indel.matt.foad.RData

---

## Dependencies

The following R packages are required for the analysis:

- dndscv
- ggplot2
- dplyr
- data.table

Install them using CRAN or Bioconductor:
install.packages(c("ggplot2", "dplyr", "data.table"))
BiocManager::install("dndscv")

---

## Notes on Reproducibility

This repository was designed to ensure reproducibility of all analyses in the manuscript. If you encounter any issues or have questions, please contact the authors.