# Songbird
Wavelet based segmentation and copy number Estimation of DLP+ scDNASeq data

# To Do List

## By 2/12/2025
- [x] create function files
- [ ] reduce dependencies & add to DESCRIPTION FILE
- [ ] finalize clustering for WGD detection
  - [x] create function to mix cells from 4 samples
  - [x] create testing metric
  - [x] test mclust vs community detection
  - [ ] test feature selection
- [ ] minimal test on hg38 files
- [ ] minimal documentation

## By 2/28/2025
- [ ] set intercompatibility between hg19 and hg38
- [ ] qc plotting functions
- [ ] allow recalibration of correction factor
- [ ] add vignette
- [ ] check inputs

Wow time passed
## By 3/26/2025

### new structure - steps
  - [x] read in bam file using qDNASeq to get counts
  - [x] Calculate GC+Map correction
  - [x] Instantiate scExperiment Class and populate with metadata and Read/GC Corrected matrix
  - [x] Make UBH Matrix and add to scExperiment Class
  - [x] Calculate Ploidy from matching bedpe file and add quality score
  - [x] Cluster cells using the breakpoints + phenograph
  - [x] For each cluster identify WGD cells
  - [x] Make copy matrix
