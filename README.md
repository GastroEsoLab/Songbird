![Drawing of a Bird on a musical staff with a artistic rendition of the algorithm](https://github.com/GastroEsoLab/Songbird/blob/master/NAR-graphical%20abstract.jpg)
Wavelet based segmentation and copy number Estimation of DLP+ scDNASeq data

# Quickstart

## Install
Songbird largely uses packages you can get from Bioconductor & CRAN. However, the per-genome qDNASeq annotations require manual install.
```R
library(BiocManager)
library(remotes)

BiocManager::install("QDNAseq.hg19")
remotes::install_github("asntech/QDNAseq.hg38@main")
remotes::install_github("GastroEsoLab/QDNAseq.hs1@main")
remotes::install_github("GastroEsoLab/Songbird@main")
```

## Bedpe Generation
  - Requires [Sinto](https://timoast.github.io/sinto/basic_usage.html), [GNU Parallel](https://www.gnu.org/software/parallel/), [bedtools > 2.30.0](https://github.com/arq5x/bedtools2/releases/tag/v2.30.0), and [samtools > 1.16.1](https://github.com/samtools/samtools/releases/tag/1.16).
  - The Script will create the individualBams subdirectory inside the output folder
```bash
bash preprocess_sample.sh -b /mnt/data/HEK_Controls/DLPPLUS_Ladder/hTERT-1_139566A/alignment/all_cells_bulk_control.bam -o ~/htert-1/
```

## Run Songbird
```R
library(devtools)
setwd('~/Songbird/')
load_all()

folder = '~/hTert-1/individualBams'
bedpes <- list.files(folder, pattern = 'bedpe$', full.names = T)
bams <- list.files(folder, pattern = 'bam$', full.names = T)
sbird_sce <- Songbird::process.batch(bams = bams, bedpes = bedpes, genome = 'hg38', n_cpu = 36)
sbird_sce <- ploidy_correction(sbird_sce, min_reads = 50000)
sbird_sce <- copyCall(sbird_sce)
sbird_sce <- identify_subclones(sbird_sce)
```

## Final Plotting
For your convience & to help with validating the copy number calls, we provide a couple of plotting functions
```R
# Plots a heatmap of the copy number state of the whole sample split by subclone
plot_heatmap(sbird_sce, assay_name = 'copy', row_split = sbird_sce$subclone)

# Plots an individual cell for comparing 
plot_cell(parameters)
```

# Parameters

Banner designed by [Jennifer Jones](https://github.com/jejonesu)
