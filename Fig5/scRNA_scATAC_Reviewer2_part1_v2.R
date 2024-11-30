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
library(ComplexHeatmap)

### updated heterogeneity-guide expression signatures (Fig.4f, g)
TTs.cancer.seu <-LoadH5Seurat('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_060223.h5seurat')
RT_CS1_markers <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_marker_all.xlsx',sheet='RT_CS1_all')
RT_CS2_markers <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_marker_all.xlsx',sheet='RT_CS2_all')
RT_CS3_markers <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_marker_all.xlsx',sheet='RT_CS3_all')
# get new module signatures
RT_CS1_markers <-RT_CS1_markers[order(RT_CS1_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
RT_CS1_module <- RT_CS1_markers[RT_CS1_markers$avg_log2FC>=0.25,]
RT_CS2_markers <-RT_CS2_markers[order(RT_CS2_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
RT_CS2_module <- RT_CS2_markers[RT_CS2_markers$avg_log2FC>=0.25,]
RT_CS3_markers <-RT_CS3_markers[order(RT_CS3_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
RT_CS3_module <- RT_CS3_markers[RT_CS3_markers$avg_log2FC>=0.25,]
# update fig3f
TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                 features = list(RT_CS1_module$gene),
                                 name="RT_CS1_module",
                                 assay='SCT')
FeaturePlot(TTs.cancer.seu, features="RT_CS1_module1", repel = TRUE,
            min.cutoff = "q50", max.cutoff = "q99", pt.size=0.5,
            order=T) + theme(legend.position = "right")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                 features = list(RT_CS2_module$gene),
                                 name="RT_CS2_module",
                                 assay='SCT')

FeaturePlot(TTs.cancer.seu, features="RT_CS2_module1", repel = TRUE,
            min.cutoff = "q50", max.cutoff = "q99", pt.size=0.5,
            order=TRUE) + theme(legend.position = "right")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

TTs.cancer.seu <- AddModuleScore(TTs.cancer.seu, 
                                 features = list(RT_CS3_module$gene),
                                 name="RT_CS3_module",
                                 assay='SCT')

FeaturePlot(TTs.cancer.seu, features="RT_CS3_module1", repel = TRUE,
            min.cutoff = "q50", max.cutoff = "q99", pt.size=0.5,
            order=TRUE) + theme(legend.position = "right")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

SaveH5Seurat(TTs.cancer.seu, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/TTs_cancer_052124.h5seurat')
PRT_CS1_markers <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_marker_all.xlsx',sheet='RPT_CS_all')
PRT_CS1_markers <-PRT_CS1_markers[order(PRT_CS1_markers$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
PRT_CS1_module <- PRT_CS1_markers[PRT_CS1_markers$avg_log2FC>=0.25,]
sheetlist <- list("RT_CS1_module"=RT_CS1_module, "RT_CS2_module"=RT_CS2_module,
                  "RT_CS3_module"=RT_CS3_module, "RPT_CS_module"=PRT_CS1_module)
write_xlsx(sheetlist, '/scratch/u/kfang/scRNA_scATAC/10x_result/scRNA/RT_RPT_CS_module_all.052124.xlsx')

### update Fig.4c
# Load dp_all from part5
rtcs1_dp_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_all_v2.xlsx',sheet = "RT_CS1")
rtcs1_dp_all$...1<-NULL
rtcs2_dp_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_all_v2.xlsx',sheet = "RT_CS2")
rtcs2_dp_all$...1<-NULL
rtcs3_dp_all <- read_excel('/scratch/u/kfang/scRNA_scATAC/10x_result/scATAC/Dp_all_v2.xlsx',sheet = "RT_CS3")
rtcs3_dp_all$...1<-NULL
# only keep avg_Log2FC > 0.25
rtcs1_dp_all <- rtcs1_dp_all[rtcs1_dp_all$avg_log2FC>=0.25,]
rtcs2_dp_all <- rtcs2_dp_all[rtcs2_dp_all$avg_log2FC>=0.25,]
rtcs3_dp_all <- rtcs3_dp_all[rtcs3_dp_all$avg_log2FC>=0.25,]
cs_upsetplot <- function(rtcs1_df,rtcs2_df,rtcs3_df,labels){
  # dp_de upset
  distal_lt = list(
    RT_CS1 = rtcs1_df[rtcs1_df$genomeLoc_annot=='Distal',]$region,
    RT_CS2 = rtcs2_df[rtcs2_df$genomeLoc_annot=='Distal',]$region,
    RT_CS3 = rtcs3_df[rtcs3_df$genomeLoc_annot=='Distal',]$region
  )
  distal_m = make_comb_mat(distal_lt)
  proximal_lt = list(
    RT_CS1 = rtcs1_df[rtcs1_df$genomeLoc_annot=='Proximal',]$region,
    RT_CS2 = rtcs2_df[rtcs2_df$genomeLoc_annot=='Proximal',]$region,
    RT_CS3 = rtcs3_df[rtcs3_df$genomeLoc_annot=='Proximal',]$region
  )
  proximal_m = make_comb_mat(proximal_lt)
  promoter_lt = list(
    RT_CS1 = rtcs1_df[rtcs1_df$genomeLoc_annot=='Promoter',]$region,
    RT_CS2 = rtcs2_df[rtcs2_df$genomeLoc_annot=='Promoter',]$region,
    RT_CS3 = rtcs3_df[rtcs3_df$genomeLoc_annot=='Promoter',]$region
  )
  promoter_m = make_comb_mat(promoter_lt)
  top_ha = HeatmapAnnotation(
    "Distal" = anno_barplot(comb_size(distal_m), 
                            gp = gpar(fill = "chocolate"), height = unit(3, "cm"),
                            axis_param=list(gp=gpar(fontsize = 14,fontface="bold"))), 
    "Proximal" = anno_barplot(comb_size(proximal_m), 
                              gp = gpar(fill = "darkorange"), height = unit(3, "cm"),
                              axis_param=list(gp=gpar(fontsize = 14,fontface="bold"))), 
    "Promoter" = anno_barplot(comb_size(promoter_m), 
                              gp = gpar(fill = "orange"), height = unit(3, "cm"),
                              axis_param=list(gp=gpar(fontsize = 14,fontface="bold"))), 
    gap = unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0,
    annotation_name_gp= gpar(fontsize = 16,fontface="bold"))
  # the same for using m2 or m3
  ss = c("RT_CS1"=length(rtcs1_df$region),
         "RT_CS2"=length(rtcs2_df$region),
         "RT_CS3"=length(rtcs3_df$region))
  UpSet(distal_m, 
        top_annotation = top_ha,
        set_order = c("RT_CS1","RT_CS2","RT_CS3"),
        left_annotation = rowAnnotation(
          "Total Counts" = anno_barplot(-ss, 
                                        baseline = 0,
                                        axis_param = list(
                                          at = -1 * labels,
                                          labels = labels,
                                          labels_rot = 0),
                                        border = FALSE, 
                                        gp = gpar(fill = c("cornflowerblue","darkorchid1","deeppink")), 
                                        width = unit(4, "cm")
                                        
          ),
          set_name = anno_text(c("RT_CS1","RT_CS2","RT_CS3"), 
                               location = 0.5, 
                               just = "center",
                               gp = gpar(fontsize = 14,
                                         fontface="bold"),
                               width = max_text_width(set_name(distal_m)) + unit(4, "mm")),
          annotation_name_gp= gpar(fontsize = 14, fontface="bold")
        ), 
        right_annotation = NULL,
        show_row_names = FALSE
  )
}

cs_upsetplot(rtcs1_dp_all,rtcs2_dp_all,rtcs3_dp_all,c(0,5000,10000,15000))

