# Songbird
Wavelet based segmentation and copy number Estimation of DLP+ scDNASeq data

# To Do List
- [ ] QC Plotting Functions
  - [x] Heatmap plotter
  - [x] Cell Plotter (reads, segmentation, cn)
  - [ ] Quality Scorer
  - [x] Reads per CN
  - [ ] UMAP and PCA plotter
- [x] minimize dependencies? (Rphenograph is github install only...)
- [ ] Allow to skip WGD & Subclone ID
- [ ] Build out Vignette and Documentation
- [x] Allow recalibration of correction factor for ploidy estimation
- [x] Add Parallel processing
- [ ] Parameter checks for all user facing functions
- [ ] Convert all multicore to pbmcapply
- [ ] Write Sonbird creation functions from matrix
- [x] Wrap bedpe/bam loading into a function 

# To Done (almost) list!
- [x] CN Caller
  - [x] read in bam file using qDNASeq
  - [x] Calculate GC+Map Correction
  - [x] Calculate per cell ploidy estimation
  - [x] Figure out correction factor (try HG38 scAbsolute)
  - [x] Instantiate scExperiment Class and populate with metadata
  - [x] Make UBH Matrix and add to scExperiment Class
  - [x] Cluster cells using breakpoint matrix + Phenograph
  - [x] Correct Ploidy Estimation & Perform WGD Detection
  - [x] Produce Final Quality Score
  - [x] Make Copy Matrix
  - [x] Minimal Vignette
  - [x] Run on both hg38 and T2T
  - [x] Expand Bin Exclusion list with T2T?
  - [x] Allow to skip any ploidy estimation (Default mean ploidy range 2-8)
  - [x] Handle empty bams and bedpes gracefully
  - [x] Add flexible tagmentation overlap (for scAbsolute's silliness)
  - [x] Convert bed handling to data.table (read_table single core is the same speed)
# Minimal Vignette
```R
library(Songbird)
folder = '~/hTert-1/individualBams'
bedpes <- list.files(folder, pattern = 'bedpe$', full.names = T)
bams <- list.files(folder, pattern = 'bam$', full.names = T)
sbird_sce <- Songbird::process.batch(bams = bams, bedpes = bedpes, genome = 'hg38', n_cpu = 36)
sbird_sce <- ploidy_correction(sbird_sce, min_reads = 50000)
sbird_sce <- copyCall(sbird_sce)
sbird_sce <- identify_subclones(sbird_sce)

plot_heatmap(sbird_sce, assay_name = 'copy', row_split = sbird_sce$subclone)
```
