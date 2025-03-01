Sources codes for Fig.1 

### SessionInfo
```
> sessionInfo()
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 8.7 (Green Obsidian)

Matrix products: default
BLAS:   /hpc/apps/R/4.3.3/lib64/R/lib/libRblas.so 
LAPACK: /hpc/apps/R/4.3.3/lib64/R/lib/libRlapack.so;  LAPACK version 3.11.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Chicago
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] readxl_1.4.3                sctransform_0.4.1           celldex_1.12.0             
 [4] SeuratDisk_0.0.0.9021       writexl_1.5.0               scDblFinder_1.16.0         
 [7] scuttle_1.12.0              SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
[10] Biobase_2.62.0              GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
[13] IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1        
[16] MatrixGenerics_1.14.0       matrixStats_1.2.0           lubridate_1.9.3            
[19] forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
[22] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                
[25] tibble_3.2.1                ggplot2_3.5.0               tidyverse_2.0.0            
[28] Seurat_5.0.3                SeuratObject_5.0.1          sp_2.1-3                   

loaded via a namespace (and not attached):
  [1] spatstat.sparse_3.0-3         bitops_1.0-7                  httr_1.4.7                   
  [4] RColorBrewer_1.1-3            tools_4.3.3                   utf8_1.2.4                   
  [7] R6_2.5.1                      lazyeval_0.2.2                uwot_0.1.16                  
 [10] withr_3.0.0                   gridExtra_2.3                 progressr_0.14.0             
 [13] cli_3.6.2                     spatstat.explore_3.2-7        fastDummies_1.7.3            
 [16] spatstat.data_3.0-4           ggridges_0.5.6                pbapply_1.7-2                
 [19] Rsamtools_2.18.0              scater_1.30.1                 parallelly_1.37.1            
 [22] limma_3.58.1                  rstudioapi_0.16.0             RSQLite_2.3.6                
 [25] generics_0.1.3                BiocIO_1.12.0                 ica_1.0-3                    
 [28] spatstat.random_3.2-3         Matrix_1.6-5                  ggbeeswarm_0.7.2             
 [31] fansi_1.0.6                   abind_1.4-5                   lifecycle_1.0.4              
 [34] yaml_2.3.8                    edgeR_4.0.16                  SparseArray_1.2.4            
 [37] BiocFileCache_2.10.2          Rtsne_0.17                    grid_4.3.3                   
 [40] blob_1.2.4                    promises_1.3.0                dqrng_0.3.2                  
 [43] ExperimentHub_2.10.0          crayon_1.5.2                  miniUI_0.1.1.1               
 [46] lattice_0.22-6                beachmat_2.18.1               cowplot_1.1.3                
 [49] KEGGREST_1.42.0               pillar_1.9.0                  metapod_1.10.1               
 [52] rjson_0.2.21                  xgboost_1.7.7.1               future.apply_1.11.2          
 [55] codetools_0.2-20              leiden_0.4.3.1                glue_1.7.0                   
 [58] data.table_1.15.4             vctrs_0.6.5                   png_0.1-8                    
 [61] spam_2.10-0                   cellranger_1.1.0              gtable_0.3.4                 
 [64] cachem_1.0.8                  S4Arrays_1.2.1                mime_0.12                    
 [67] survival_3.5-8                statmod_1.5.0                 bluster_1.12.0               
 [70] interactiveDisplayBase_1.40.0 fitdistrplus_1.1-11           ROCR_1.0-11                  
 [73] nlme_3.1-164                  bit64_4.0.5                   filelock_1.0.3               
 [76] RcppAnnoy_0.0.22              irlba_2.3.5.1                 vipor_0.4.7                  
 [79] KernSmooth_2.23-22            colorspace_2.1-0              DBI_1.2.2                    
 [82] tidyselect_1.2.1              bit_4.0.5                     compiler_4.3.3               
 [85] curl_5.2.1                    BiocNeighbors_1.20.2          hdf5r_1.3.10                 
 [88] DelayedArray_0.28.0           plotly_4.10.4                 rtracklayer_1.62.0           
 [91] scales_1.3.0                  lmtest_0.9-40                 rappdirs_0.3.3               
 [94] digest_0.6.35                 goftest_1.2-3                 spatstat.utils_3.1-1         
 [97] XVector_0.42.0                htmltools_0.5.8.1             pkgconfig_2.0.3              
[100] sparseMatrixStats_1.14.0      dbplyr_2.5.0                  fastmap_1.1.1                
[103] rlang_1.1.3                   htmlwidgets_1.6.4             shiny_1.8.1.1                
[106] DelayedMatrixStats_1.24.0     zoo_1.8-12                    jsonlite_1.8.8               
[109] BiocParallel_1.36.0           BiocSingular_1.18.0           RCurl_1.98-1.14              
[112] magrittr_2.0.3                GenomeInfoDbData_1.2.11       dotCall64_1.1-1              
[115] patchwork_1.2.0               munsell_0.5.1                 Rcpp_1.0.12                  
[118] viridis_0.6.5                 reticulate_1.35.0             stringi_1.8.3                
[121] zlibbioc_1.48.2               MASS_7.3-60.0.1               AnnotationHub_3.10.1         
[124] plyr_1.8.9                    parallel_4.3.3                listenv_0.9.1                
[127] ggrepel_0.9.5                 deldir_2.0-4                  Biostrings_2.70.3            
[130] splines_4.3.3                 tensor_1.5                    hms_1.1.3                    
[133] locfit_1.5-9.9                igraph_2.0.3                  spatstat.geom_3.2-9          
[136] RcppHNSW_0.6.0                reshape2_1.4.4                ScaledMatrix_1.10.0          
[139] BiocVersion_3.18.1            XML_3.99-0.16.1               scran_1.30.2                 
[142] BiocManager_1.30.22           tzdb_0.4.0                    httpuv_1.6.15                
[145] RANN_2.6.1                    polyclip_1.10-6               future_1.33.2                
[148] scattermore_1.2               rsvd_1.0.5                    xtable_1.8-4                 
[151] restfulr_0.0.15               RSpectra_0.16-1               later_1.3.2                  
[154] viridisLite_0.4.2             memoise_2.0.1                 beeswarm_0.4.0               
[157] AnnotationDbi_1.64.1          GenomicAlignments_1.38.2      cluster_2.1.6                
[160] timechange_0.3.0              globals_0.16.3   
```
