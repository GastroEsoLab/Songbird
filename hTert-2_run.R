out.dir <- '/mnt/data00/karol/songbird/hTERT-2_139566A/run'

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
                           n_cpu             = 50)

sbird_sce <- ploidy_correction_hmm(sbird_sce, n_cpu = 50, cn_states = 1:6)
sbird_sce <- copyCall(sbird_sce)
sbird_sce <- identify_subclones(sbird_sce)





plot_overlap_fit(sbird_sce, cell = colnames(sbird_sce)[6])
plot_cell(sbird_sce, cell = colnames(sbird_sce)[2])
# plot_cell(sbird_sce, cell = colnames(sbird_sce)[2])

plot_heatmap(sbird_sce, assay_name = 'copy', row_split = sbird_sce$subclone)
plot_heatmap(sbird_sce, assay_name = 'hmm_cn', row_split = sbird_sce$subclone)



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

inspect_cell(sbird_sce, colnames(sbird_sce)[6])

# Return for cowplot composition
p <- plot_overlap_fit(sbird_sce, cell = colnames(sbird_sce)[2], return_plot = TRUE)
cowplot::plot_grid(p, plot_cell(sbird_sce, cell = colnames(sbird_sce)[2], return_plot = TRUE))
