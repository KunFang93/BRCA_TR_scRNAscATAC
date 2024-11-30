library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(SeuratDisk)
library(tidyverse)
library(plyr)
library(ggVennDiagram)
library(stringr)
library(ComplexHeatmap)
library(reshape2)
library(readxl)
library(writexl)

TTs.cancer.atac.seu <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.addmotif.0620.rds')
hcs_has_region <- read.csv('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/hcs_has_region.csv')
# choosing background peaks
open.peaks <- AccessiblePeaks(TTs.cancer.atac.seu, idents = c("PT_CS1","PT_CS2","PT_CS3","PT_CS4","PT_CS5"))
# match the overall GC content in the peak set
meta.feature <- GetAssayData(TTs.cancer.atac.seu, assay = "peaks", slot = "meta.features")

peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[hcs_has_region$region, ],
  n = 50000
)

hcs_has_enriched.motifs <- FindMotifs(
  object = TTs.cancer.atac.seu,
  features = hcs_has_region$region
)
motif.to.plot <- hcs_has_enriched.motifs[hcs_has_enriched.motifs$motif.name %in% c("FOS","FOSL1::JUND","FOSL1::JUNB","FOXA1","ESR1","SP2"),]
MotifPlot(
  object = TTs.cancer.atac.seu,
  motifs = head(rownames(motif.to.plot)),
  ncol=2,
)

# visualize
rtcs_promoter_links <- read.csv('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/rtcs_links_hcs_promoter.csv')
rtcs_promoter_links_gr <- GRanges(seqnames = rtcs_promoter_links$seqnames,
              ranges = IRanges(start = rtcs_promoter_links$start, end = rtcs_promoter_links$end),
              strand = rtcs_promoter_links$strand,
              score = rtcs_promoter_links$score,
              group = rtcs_promoter_links$group
              )

Links(TTs.cancer.atac.seu) <- rtcs_promoter_links_gr
cs_list <- c("PT_CS1", "PT_CS2", "PT_CS3", "PT_CS4", "PT_CS5", "PRT_CS1", "RT_CS1", 
             "RT_CS2",  "RT_CS3")
annotations <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/annotation.rds')
Annotation(TTs.cancer.atac.seu) <- annotations
levels(TTs.cancer.atac.seu) <- cs_list
cov_plot_TTs <- CoveragePlot(
  object = TTs.cancer.atac.seu,
  assay = 'ATAC',
  annotation = FALSE,
  peaks = F,
  tile=TRUE,
  region='chr20-57100000-57350000'
)
gene_plot <- AnnotationPlot(
  object = TTs.cancer.atac.seu,
  region = "chr20-57100000-57350000"
)
peak_plot <-  PeakPlot(
  object = TTs.cancer.atac.seu,
  region = "chr20-57100000-57350000",
  assay='peaks'
)
link_plot <- LinkPlot(
  object = TTs.cancer.atac.seu,
  region = "chr20-57100000-57350000",
)
CombineTracks(
  plotlist = list( cov_plot_TTs, gene_plot, peak_plot, link_plot),
  heights = c(12,3,3,5),
)


cov_plot_TTs <- CoveragePlot(
  object = TTs.cancer.atac.seu,
  assay = 'ATAC',
  annotation = FALSE,
  peaks = TRUE,
  tile=TRUE,
  region='CA2',
  extend.upstream = 10000,
  extend.downstream = 100000
)
gene_plot <- AnnotationPlot(
  object = TTs.cancer.atac.seu,
  region = "CA2",
  mode='transcript',
  extend.upstream = 10000,
  extend.downstream = 100000
)

# peak_plot <-  PeakPlot(
#   object = TTs.cancer.atac.seu,
#   region = "chr8:85464007-85581493",
#   assay='peaks'
# )
link_plot <- LinkPlot(
  object = TTs.cancer.atac.seu,
  region = "CA2",
  extend.upstream = 10000,
  extend.downstream = 100000,
)
CombineTracks(
  plotlist = list( cov_plot_TTs, gene_plot, link_plot),
  heights = c(12,3,5),
)
