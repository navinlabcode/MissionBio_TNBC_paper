# Mission Bio TNBC
This repository contains the metadata, software, and scripts used in the manuscript: Reconstructing mutational lineages in breast cancer by multi-patient-targeted single-cell DNA sequencing.

#### Metadata
includes clinical patient information, total genes submitted for multi-patient targeted panel, and bulk exome / single cell sequencing metrics.

#### Scripts 
R scripts in this directory include all the code required to reproduce the figures in the Mission Bio TNBC paper.

#### Software
includes [*Tapestri Insights*](https://dl.missionbio.io/insights/Tapestri%20Insights%202.2.dmg)(*v2.2*) software, custom panel file, and Tapestri.cnv R Package.

### _Dependencies_
------------
R scripts are dependent on [*CopyKit*](https://github.com/navinlabcode/copykit)(*v0.1.0*) and Tapestri.cnv, which can be installed by:
``` r
devtools::install_github(repo = "navinlabcode/copykit",ref="f709a48")
```
``` r
devtools::install_local(path = "tapestri.cnv_0.1.0.tar.gz", repos='http://cran.us.r-project.org', clean = TRUE)
```
Session info:
```
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Gini_0.1.0           Rtsne_0.15           ggrepel_0.8.2        magrittr_1.5        
 [5] tapestri.cnv_1.0.0   hdf5r_1.3.2          BiocParallel_1.18.1  ineq_0.2-13         
 [9] umap_0.2.6.0         phangorn_2.5.5       ape_5.4-1            circlize_0.4.10     
[13] ComplexHeatmap_2.0.0 flowCore_1.50.0      cowplot_1.1.1        forcats_0.5.0       
[17] stringr_1.4.0        dplyr_1.0.2          purrr_0.3.4          readr_1.3.1         
[21] tidyr_1.1.1          tibble_3.0.3         ggplot2_3.3.2        tidyverse_1.3.0     

loaded via a namespace (and not attached):
 [1] nlme_3.1-148        matrixStats_0.56.0  fs_1.5.0            bit64_4.0.2        
 [5] lubridate_1.7.9     RColorBrewer_1.1-2  httr_1.4.2          tools_3.6.1        
 [9] backports_1.1.8     R6_2.4.1            DBI_1.1.0           BiocGenerics_0.30.0
[13] colorspace_1.4-1    GetoptLong_1.0.2    withr_2.2.0         tidyselect_1.1.0   
[17] bit_4.0.4           compiler_3.6.1      graph_1.62.0        cli_2.0.2          
[21] rvest_0.3.6         Biobase_2.44.0      xml2_1.3.2          scales_1.1.1       
[25] DEoptimR_1.0-8      mvtnorm_1.1-1       robustbase_0.93-6   quadprog_1.5-8     
[29] askpass_1.1         rrcov_1.5-5         pkgconfig_2.0.3     dbplyr_1.4.4       
[33] rlang_0.4.7         GlobalOptions_0.1.2 readxl_1.3.1        rstudioapi_0.11    
[37] shape_1.4.4         generics_0.0.2      jsonlite_1.7.0      Matrix_1.2-18      
[41] Rcpp_1.0.5          munsell_0.5.0       fansi_0.4.1         reticulate_1.16    
[45] lifecycle_0.2.0     stringi_1.4.6       MASS_7.3-52         blob_1.2.1         
[49] parallel_3.6.1      crayon_1.3.4        lattice_0.20-41     haven_2.3.1        
[53] hms_0.5.3           pillar_1.4.6        igraph_1.2.5        rjson_0.2.20       
[57] corpcor_1.6.9       stats4_3.6.1        fastmatch_1.1-0     reprex_0.3.0       
[61] glue_1.4.1          modelr_0.1.8        png_0.1-7           vctrs_0.3.2        
[65] cellranger_1.1.0    gtable_0.3.0        openssl_1.4.2       clue_0.3-57        
[69] assertthat_0.2.1    broom_0.7.0         RSpectra_0.16-0     pcaPP_1.9-73       
[73] cluster_2.1.0       ellipsis_0.3.1                 
```
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggtree_3.2.1                            ape_5.6-2                               dendextend_1.15.2                      
 [4] dbscan_1.1-10                           copykit_0.1.0                           DNAcopy_1.68.0                         
 [7] Rsubread_2.8.1                          SingleCellExperiment_1.16.0             SummarizedExperiment_1.24.0            
[10] MatrixGenerics_1.6.0                    matrixStats_0.61.0                      ggalt_0.4.0                            
[13] DEGreport_1.33.1                        ggpubr_0.4.0                            Homo.sapiens_1.3.1                     
[16] EnsDb.Hsapiens.v86_2.99.0               org.Hs.eg.db_3.14.0                     GO.db_3.14.0                           
[19] OrganismDbi_1.36.0                      GenomicFeatures_1.46.5                  GenomicRanges_1.46.1                   
[22] GenomeInfoDb_1.30.1                     AnnotationDbi_1.56.2                    IRanges_2.28.0                         
[25] S4Vectors_0.32.3                        Biobase_2.54.0                          BiocGenerics_0.40.0                    
[28] RColorBrewer_1.1-2                      ComplexHeatmap_2.10.0                   forcats_0.5.1                          
[31] stringr_1.4.0                           purrr_0.3.4                             readr_2.1.2                            
[34] tidyverse_1.3.1                         tibble_3.1.6                            useful_1.2.6                           
[37] ggplot2_3.3.5                           cowplot_1.1.1                           tidyr_1.2.0                            
[40] dplyr_1.0.8                    
```

### _Data source_
------------
The original sequencing data from this study has been deposited to the Sequence Read Archive (SRA): PRJNA763862, PRJNA629885.

### _Contact_
------------
For any additional information, please [email](mailto:nnavin@mdanderson.org) corresponding author.
