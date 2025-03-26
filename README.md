# Songbird
Wavelet based segmentation and copy number Estimation of DLP+ scDNASeq data

# To Do List
- [ ] QC Plotting Functions
  - [ ] Heatmap plotter
  - [ ] Cell Plotter (reads, segmentation, cn)
  - [ ] Quality Scorer
  - [ ] Reads per CN
  - [ ] UMAP and PCA plotter
- [ ] minimize dependencies? (Rphenograph is github install only...)
- [ ] Run on both hg19 and hg38 (and T2T?)
- [ ] Allow to skip WGD & Subclone ID
- [ ] Allow to skip any ploidy estimation (Default to CN 2 or a provided value)
- [ ] Build out Vignette and Documentation
- [ ] Allow recalibration of correction factor for ploidy estimation
- [ ] Add Parallel processing(??)
- [ ] Parameter checks for all user facing functions


# To Done (almost) list!
- [-] CN Caller
  - [x] read in bam file using qDNASeq
  - [x] Calculate GC+Map Correction
  - [-] Calculate per cell ploidy estimation
    - [ ] Figure out correction factor (try HG38 scAbsolute)
  - [x] Instantiate scExperiment Class and populate with metadata
  - [x] Make UBH Matrix and add to scExperiment Class
  - [x] Cluster cells using breakpoint matrix + Phenograph
  - [x] Correct Ploidy Estimation & Perform WGD Detection
  - [-] Produce Final Quality Score
    - Our best cells have middling quality scores in the PBMC dataset
  - [x] Make Copy Matrix
  - [x] Minimal Vignette

# Minimal Vignette
```R
library(Songbird)
folder = '~/hTert-1/individualBams'
bams <- list.files(path = folder, pattern = ".bam$", full.names = T)
bedpes <- list.files(path = folder, pattern = ".bedpe$", full.names = T)
res <- mclapply(1:length(bams), function(i) process.cell(bams[i], bedpes[i], bin.size = 500000, min.svSize = 1e6, min_length = 50, max_length = 1000), mc.cores = 32)
sbird_sce <- create_sce(res)
sbird_sce <- identify_subclones(sbird_sce)
sbird_sce <- ploidy_correction(sbird_sce)
sbird_sce <- copyCall(sbird_sce)
```
