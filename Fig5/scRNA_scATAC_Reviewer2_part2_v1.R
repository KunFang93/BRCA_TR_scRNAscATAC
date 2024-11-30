library(Signac)
library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(cicero)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(RColorBrewer)
library(writexl)

annotations <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/annotation.rds')
TTs.cancer.atac.seu <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/TTs.cancer.atac.0620.rds')
HCS_df <- read.csv('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/hcs_df.csv')
HCS_peaks <- read.csv('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/hcs_peaks.csv')
# subset RTCS
meta.data.TTs <- TTs.cancer.atac.seu@meta.data
meta.data.RTCS <- meta.data.TTs[meta.data.TTs$cell.state %in% c('RT_CS1','RT_CS2','RT_CS3','PRT_CS1'),]
RTCS.cancer.atac.seu <- subset(TTs.cancer.atac.seu, cells = meta.data.RTCS$barcodes) 


################# co accessibility #######################
# convert to CellDataSet format and make the cicero object
RTCS.cancer.atac.cds <- as.cell_data_set(x = RTCS.cancer.atac.seu)
RTCS.cancer.atac.cicero <- make_cicero_cds(RTCS.cancer.atac.cds, reduced_coordinates = reducedDims(RTCS.cancer.atac.cds)$UMAP)
# get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
# use chromosome 1 to save some time
# omit this step to run on the whole genome
genome <- genome[1:24]
# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)
# run cicero
RTCS.cancer.atac.conns <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTCS.cancer.atac.conns.rds')
RTCS.cancer.atac.ccans <- readRDS('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTCS.cancer.atac.ccans.rds')

# check distribution of coaccess
ggplot(RTCS.cancer.atac.conns, aes(x = "", y = coaccess)) +
  geom_violin(fill = "blue") +
  labs(title = "Violin Plot of coaccess", x = "", y = "coaccess") +
  theme_minimal()

# RTCS.cancer.atac.conns.filt <- RTCS.cancer.atac.conns[RTCS.cancer.atac.conns$coaccess > quantile(RTCS.cancer.atac.conns$coaccess,0.75,na.rm=T),]
# RTCS.cancer.atac.ccans.filt <- generate_ccans(RTCS.cancer.atac.conns.filt)
# RTCS.cancer.atac.links.filt <- ConnectionsToLinks(conns = RTCS.cancer.atac.conns.filt, ccans = RTCS.cancer.atac.ccans.filt)
# Links(RTCS.cancer.atac.seu) <- RTCS.cancer.atac.links.filt
RTCS.cancer.atac.conns <- run_cicero(RTCS.cancer.atac.cicero, genomic_coords = genome.df, sample_num = 150)
RTCS.cancer.atac.ccans <- generate_ccans(RTCS.cancer.atac.conns)
# saveRDS(RTCS.cancer.atac.conns, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTCS.cancer.atac.conns.rds')
# saveRDS(RTCS.cancer.atac.ccans, '/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTCS.cancer.atac.ccans.rds')
RTCS.cancer.atac.links <- ConnectionsToLinks(conns = RTCS.cancer.atac.conns, ccans = RTCS.cancer.atac.ccans)
Links(RTCS.cancer.atac.seu) <- RTCS.cancer.atac.links
SaveH5Seurat(RTCS.cancer.atac.seu,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTCS.cancer.atac.seu.rds')

RTCS.cancer.atac.links.df <- as.data.frame(RTCS.cancer.atac.links)
write.csv(RTCS.cancer.atac.links.df,'/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/RTCS.cancer.atac.links.csv')


