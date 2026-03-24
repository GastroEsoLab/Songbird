
inspect_cell <- function(sbird_sce, cell) {
  hmm_cn  <- SummarizedExperiment::assay(sbird_sce, 'hmm_cn')[, cell]
  copy_cn <- SummarizedExperiment::assay(sbird_sce, 'copy')[, cell]
  rd      <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))

  dat <- data.frame(
    bin    = seq_along(hmm_cn),
    chr    = rd$chr,
    hmm    = hmm_cn,
    copy   = copy_cn,
    diff   = hmm_cn - copy_cn
  )

  p1 <- ggplot2::ggplot(dat, ggplot2::aes(x = bin, y = hmm,  colour = factor(hmm)))  +
    ggplot2::geom_point(size = 0.4) + ggplot2::theme_classic() +
    ggplot2::labs(title = 'HMM CN', colour = 'CN') + ggplot2::ylim(0, 10)

  p2 <- ggplot2::ggplot(dat, ggplot2::aes(x = bin, y = copy, colour = factor(copy))) +
    ggplot2::geom_point(size = 0.4) + ggplot2::theme_classic() +
    ggplot2::labs(title = 'copyCall CN', colour = 'CN') + ggplot2::ylim(0, 10)

  p3 <- ggplot2::ggplot(dat[dat$diff != 0, ],
                        ggplot2::aes(x = bin, y = diff, fill = factor(sign(diff)))) +
    ggplot2::geom_col(width = 1) + ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(values = c('-1' = '#496bab', '1' = '#b3402e'),
                               labels = c('copyCall higher', 'HMM higher')) +
    ggplot2::labs(title = 'Disagreement (HMM - copyCall)', fill = NULL)

  cowplot::plot_grid(p1, p2, p3, ncol = 1, align = 'v')
}

out.dir <- '/mnt/data00/karol/songbird/HEK-1_139566A/run'

if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = T)
}

setwd(out.dir)
# library(Songbird)
library(ggplot2)
devtools::load_all("~/Dropbox/Postdoc/git/Songbird/")

folder = gsub("run", "individualBams", out.dir)
bedpes <- list.files(folder, pattern = 'bedpe', full.names = T)
bams <- list.files(folder, pattern = 'bam$', full.names = T)



sbird_sce <- process.batch(bams, bedpes, genome = 'hg38',
                           background_method = 'symmetric',
                           symmetric_x       = 50000,
                           ploidy_method     = 'hmm',
                           min_length = 100,
                           n_cpu             = 50,
                           near_start_exclusion = 10)

# sbird_sce <- ploidy_correction_hmm(sbird_sce, n_cpu = 50, cn_states = 1:8, p_stay = 0.95, max_iter = 10)



# Equal weight between doublet and read depth
sbird_sce <- ploidy_correction_hmm_joint(sbird_sce,
                                         cn_states = 1:8,
                                         p_stay    = 1 - 1e-10,
                                         max_iter  = 0,
                                         alpha     = 0.5,
                                         n_cpu     = 50)


sbird_sce <- copyCall(sbird_sce)
sbird_sce <- identify_subclones(sbird_sce)


plot_hmm_doublet(sbird_sce, cell = colnames(sbird_sce)[cell.id], window = 51, chr = "1") #+ ylim(-0.5,2.5)


# plot_cell(sbird_sce, cell = colnames(sbird_sce)[2])

plot_heatmap(sbird_sce, assay_name = 'copy', row_split = sbird_sce$subclone)
plot_heatmap(sbird_sce, assay_name = 'hmm_cn', row_split = sbird_sce$subclone)


cell.id <- 8

plot_overlap_fit(sbird_sce, cell = colnames(sbird_sce)[cell.id])
plot_cell(sbird_sce, cell = colnames(sbird_sce)[cell.id], chr = 1)
plot_cell(sbird_sce, cell = colnames(sbird_sce)[cell.id], assay = "hmm_cn", chr = 1)
plot_hmm_doublet(sbird_sce, cell = colnames(sbird_sce)[cell.id], window = 51, chr = "1")

inspect_cell(sbird_sce, colnames(sbird_sce)[cell.id])

# Return for cowplot composition
p <- plot_overlap_fit(sbird_sce, cell = colnames(sbird_sce)[cell.id], return_plot = TRUE)
cowplot::plot_grid(p, plot_cell(sbird_sce, cell = colnames(sbird_sce)[cell.id], return_plot = TRUE))




# Check how many breakpoints wavelet found for this cell
cell <- colnames(sbird_sce)[cell.id]
wave <- SummarizedExperiment::assay(sbird_sce, 'segmented')[, cell]
chr1_idx <- SummarizedExperiment::rowData(sbird_sce)$chr == 'chr1'
wave_chr1 <- wave[chr1_idx]
sum(diff(wave_chr1) != 0, na.rm = TRUE)  # number of wavelet breakpoints on chr1

# And plot it
plot_cell(sbird_sce, cell = cell, assay = 'segmented', chr = 'chr1')


















# Read-depth dominated — useful if doublet signal is very noisy
sbird_sce <- ploidy_correction_hmm_joint(sbird_sce, alpha = 0.8, ...)

# alpha = 0 exactly reproduces ploidy_correction_hmm
sbird_sce <- ploidy_correction_hmm_joint(sbird_sce, alpha = 0.0, ...)



