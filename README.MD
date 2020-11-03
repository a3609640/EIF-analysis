# Using this Package

## System Requirements
In general genome analysis carries relatively high demands for
compute, RAM, and disk I/O resources (but not for graphics resources, at
least not with the software used in this project).

This project also makes use of various resource-intensive R packages.

Nonetheless, the necessary hardware is attainable in high-end
consumer-grade systems.

### Description of Development Systems
The following systems have been used to execute the R scripts
in this project:

1. (verified) System76 "Serval" mobile workstation
    * Intel i7-8700k CPU
    * 64GB RAM (DDR4-3000, non-ECC)
    * Samsung NVMe Pro SSD
    * Pop!_OS 18.04 LTS
    * RStudio
    * R 3.6.3

2. (verification in progress) PowerSpec G460 desktop computer
    * Intel i7-8700k CPU
    * 64GB RAM (DDR4-3200, non-ECC)
    * Intel M.2 SATA SSD
    * Samsung NVMe Evo+ SSD
    * Windows 10 Pro
    * RStudio for Windows
    * R 4.0.3

Additional details of these environments are provided in the "Session Information" section below.

## Datasets
Please download all datasets from the following weblinks and store them at the following directory. The script will read all datasets from there.
```
~/Download/
```
This project makes use of data files that must be fetched from remote sources.  The necessary links are provided below.  Note that the .gz files must be unzipped after they are downloaded.

### TCGA and GTEX DATA

#### TCGA CNV dataset (thresholded)
https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz

#### TCGA CNV dataset
https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz

#### TCGA CNV ratio dataset
https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz

#### TCGA RNA-Seq dataset
https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz

#### TCGA sample type annotation
https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz

#### TCGA OS data ##
https://xenabrowser.net/datapages/?dataset=Survival_SupplementalTable_S1_20171025_xena_sp&host=https%3A%2F%2Fpancanatlas.xenahubs.net

#### TCGA and GTEX RNA-Seq dataset
https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz

#### TCGA and GTEX sample type annotation
download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz

### CPTAC DATA

#### CPTAC LUAD Proteomics
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_Proteome_CDAP_Protein_Report.r1/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv

#### CPTAC LUAD Phosproteomics
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_Phosphoproteome_CDAP_Protein_Report.r1/CPTAC3_Lung_Adeno_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv

#### CPTAC LUAD Sample Annotation
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx

#### CPTAC Clinical Data
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx

## Libraries
The work here depends upon many R libraries.

The following command may be a useful way to install them all:
```
BiocManager::install(c("clusterProfiler",
                       "circlize",
                       "ComplexHeatmap",
                       "corrplot",
                       "data.table",
                       "dendextend",
                       "descr",
                       "devtools"
                       "dplyr",
                       "EnvStats",
                       "eulerr",
                       "factoextra",
                       "FactoMineR",
                       "forcats",
                       "forestmodel",
                       "forestplot",
                       "Hmisc",
                       "ggfortify",
                       "ggplot2",
                       "ggpubr",
                       "ggsignif",
                       "ggthemes",
                       "glmnet",
                       "gplots",
                       "grid",
                       "gridExtra",
                       "igraph",
                       "KEGG.db",
                       "lemon",
                       "limma",
                       "MASS",
                       "missMDA",
                       "org.Hs.eg.db",
                       "parallel",
                       "pca3d",
                       "pheatmap",
                       "plotmo",
                       "RColorBrewer",
                       "ReactomePA",
                       "readr",
                       "readxl",
                       "reshape2",
                       "rgl",
                       "scales",
                       "stringr",
                       "survival",
                       "survivalAnalysis",
                       "survMisc",
                       "survminer",
                       "tidyverse",
                       "vcd",
                       "vip"
), update = TRUE, ask = FALSE)
```

This is also necessary:
```
devtools::install_github("zeehio/facetscales")
```

## Output directory
The script will save all the generated plots into the following output folder.
Please generate the output directory.
```
mkdir -p ~/Documents/EIF_output/{CNV,Expression,PCA,PCA/All,PCA/GTEX,PCA/TCGA,PCA/Lung,KM,Cox,Heatmap,CPTAC}
```

## Instructions
1. Install R (a version verified in one of the systems listed above) and RStudio at the recommended Systems.
2. Download all datasets from the listed web links above and store them under the following directory.
```
~/Download/
```
3. Install all the R packages listed above.
4. Create output directories
```
~/Documents/EIF_output/{CNV,Expression,PCA,PCA/All,PCA/GTEX,PCA/TCGA,PCA/Lung,KM,Cox,Heatmap,CPTAC}
```
5. Open the EIFanalysis.R file in RStudio, and run it with the "source" button.
6. The EIFanalysis.R will produce all the figures from the EIF manuscript and store them in the output folders.

If the root directory paths ```~/Download``` and ```~/Documents/EIF_output``` do not suit, they may be adjusted trivially in these lines near the top of the script:

```
data.file.directory <- "~/Downloads"
output.directory <- "~/Documents/EIF_output"
```

Note that you must still pre-create all output subdirectories, under whatever parent output directory you choose.


## Session Information
### System 1
```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Pop!_OS 18.04 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils
 [8] datasets  methods   base

other attached packages:
 [1] vip_0.2.2              vcd_1.4-8              stringr_1.4.0
 [4] purrr_0.3.4            tidyr_1.1.2            tibble_3.0.4
 [7] tidyverse_1.3.0        survminer_0.4.8        survMisc_0.5.5
[10] survivalAnalysis_0.1.3 scales_1.1.1           rgl_0.100.54
[13] reshape2_1.4.4         readxl_1.3.1           readr_1.4.0
[16] ReactomePA_1.28.0      RColorBrewer_1.1-2     plotmo_3.6.0
[19] TeachingDemos_2.12     plotrix_3.7-8          pheatmap_1.0.12
[22] pca3d_0.10.2           org.Hs.eg.db_3.8.2     nortest_1.0-4
[25] missMDA_1.17           limma_3.40.6           lemon_0.4.5
[28] KEGG.db_3.2.3          igraph_1.2.6           gridExtra_2.3
[31] gplots_3.1.0           glmnet_4.0-2           Matrix_1.2-18
[34] ggthemes_4.2.0         ggsignif_0.6.0         ggpubr_0.4.0
[37] ggfortify_0.4.11       Hmisc_4.4-1            Formula_1.2-4
[40] survival_3.2-7         lattice_0.20-41        forestplot_1.10
[43] checkmate_2.0.0        magrittr_1.5           forestmodel_0.6.2
[46] forcats_0.5.0          FactoMineR_2.3         factoextra_1.0.7
[49] ggplot2_3.3.2          facetscales_0.1.0.9000 eulerr_6.1.0
[52] EnvStats_2.3.1         dplyr_1.0.2            descr_1.1.4
[55] dendextend_1.14.0      data.table_1.13.2      corrplot_0.84
[58] ComplexHeatmap_2.0.0   circlize_0.4.10        clusterProfiler_3.12.0
[61] car_3.0-10             carData_3.0-4          AnnotationDbi_1.46.1
[64] IRanges_2.18.3         S4Vectors_0.22.1       Biobase_2.44.0
[67] BiocGenerics_0.30.0

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.1          bit64_4.0.5             knitr_1.30
  [4] rpart_4.1-15            doParallel_1.0.16       generics_0.0.2
  [7] cowplot_1.1.0           RSQLite_2.2.1           mice_3.11.0
 [10] europepmc_0.4           bit_4.0.4               enrichplot_1.4.0
 [13] webshot_0.5.2           xml2_1.3.2              lubridate_1.7.9
 [16] httpuv_1.5.4            assertthat_0.2.1        viridis_0.5.1
 [19] xfun_0.18               hms_0.5.3               promises_1.1.1
 [22] fansi_0.4.1             progress_1.2.2          caTools_1.18.0
 [25] dbplyr_1.4.4            km.ci_0.5-2             DBI_1.1.0
 [28] htmlwidgets_1.5.2       ellipsis_0.3.1          crosstalk_1.1.0.1
 [31] backports_1.1.10        vctrs_0.3.4             abind_1.4-5
 [34] withr_2.3.0             ggforce_0.3.2           triebeard_0.3.0
 [37] prettyunits_1.1.1       cluster_2.1.0           DOSE_3.10.2
 [40] crayon_1.3.4            ellipse_0.4.2           labeling_0.4.2
 [43] pkgconfig_2.0.3         tweenr_1.0.1            nnet_7.3-14
 [46] rlang_0.4.8             lifecycle_0.2.0         miniUI_0.1.1.1
 [49] modelr_0.1.8            cellranger_1.1.0        polyclip_1.10-0
 [52] lmtest_0.9-38           graph_1.62.0            urltools_1.7.3
 [55] KMsurv_0.1-5            zoo_1.8-8               reprex_0.3.0
 [58] base64enc_0.1-3         ggridges_0.5.2          GlobalOptions_0.1.2
 [61] png_0.1-7               viridisLite_0.3.0       rjson_0.2.20
 [64] bitops_1.0-6            KernSmooth_2.23-17      blob_1.2.1
 [67] shape_1.4.5             qvalue_2.16.0           manipulateWidget_0.10.1
 [70] jpeg_0.1-8.1            rstatix_0.6.0           gridGraphics_0.5-0
 [73] reactome.db_1.68.0      leaps_3.1               memoise_1.1.0
 [76] graphite_1.30.0         plyr_1.8.6              compiler_3.6.3
 [79] clue_0.3-57             cli_2.1.0               htmlTable_2.1.0
 [82] MASS_7.3-51.6           tidyselect_1.1.0        stringi_1.5.3
 [85] GOSemSim_2.10.0         latticeExtra_0.6-29     ggrepel_0.8.2
 [88] fastmatch_1.1-0         tools_3.6.3             rio_0.5.16
 [91] rstudioapi_0.11         foreach_1.5.1           foreign_0.8-76
 [94] scatterplot3d_0.3-41    farver_2.0.3            ggraph_2.0.3
 [97] digest_0.6.26           rvcheck_0.1.8           BiocManager_1.30.10
[100] shiny_1.5.0             Rcpp_1.0.5              broom_0.7.2
[103] later_1.1.0.1           httr_1.4.2              colorspace_1.4-1
[106] polylabelr_0.2.0        rvest_0.3.6             fs_1.5.0
[109] splines_3.6.3           graphlayouts_0.7.0      ggplotify_0.0.5
[112] xtable_1.8-4            jsonlite_1.7.1          tidygraph_1.2.0
[115] UpSetR_1.4.0            flashClust_1.01-2       R6_2.4.1
[118] pillar_1.4.6            htmltools_0.5.0         mime_0.9
[121] glue_1.4.2              fastmap_1.0.1           BiocParallel_1.18.1
[124] codetools_0.2-16        fgsea_1.10.1            mvtnorm_1.1-1
[127] curl_4.3                gtools_3.8.2            zip_2.1.1
[130] GO.db_3.8.2             openxlsx_4.2.2          munsell_0.5.0
[133] DO.db_2.9               GetoptLong_1.0.4        iterators_1.0.13
[136] haven_2.3.1             gtable_0.3.0            tidytidbits_0.2.2     
```

### System 2
```
> sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C
[5] LC_TIME=English_United States.1252

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets
 [9] methods   base

other attached packages:
 [1] vip_0.2.2              vcd_1.4-8              stringr_1.4.0
 [4] purrr_0.3.4            tidyr_1.1.2            tibble_3.0.4
 [7] tidyverse_1.3.0        survminer_0.4.8        survMisc_0.5.5
[10] survivalAnalysis_0.1.3 scales_1.1.1           rgl_0.100.54
[13] reshape2_1.4.4         readxl_1.3.1           readr_1.4.0
[16] ReactomePA_1.32.0      RColorBrewer_1.1-2     plotmo_3.6.0
[19] TeachingDemos_2.12     plotrix_3.7-8          pheatmap_1.0.12
[22] pca3d_0.10.2           org.Hs.eg.db_3.11.4    nortest_1.0-4
[25] missMDA_1.17           limma_3.44.3           lemon_0.4.5
[28] KEGG.db_3.2.4          igraph_1.2.6           gridExtra_2.3
[31] gplots_3.1.0           glmnet_4.0-2           Matrix_1.2-18
[34] ggthemes_4.2.0         ggsignif_0.6.0         ggpubr_0.4.0
[37] ggfortify_0.4.11       Hmisc_4.4-1            Formula_1.2-4
[40] survival_3.2-7         lattice_0.20-41        forestplot_1.10
[43] checkmate_2.0.0        magrittr_1.5           forestmodel_0.6.2
[46] forcats_0.5.0          FactoMineR_2.3         factoextra_1.0.7
[49] ggplot2_3.3.2          facetscales_0.1.0.9000 eulerr_6.1.0
[52] EnvStats_2.4.0         dplyr_1.0.2            descr_1.1.4
[55] dendextend_1.14.0      data.table_1.13.2      corrplot_0.84
[58] ComplexHeatmap_2.4.3   circlize_0.4.10        clusterProfiler_3.16.1
[61] car_3.0-10             carData_3.0-4          AnnotationDbi_1.50.3
[64] IRanges_2.22.2         S4Vectors_0.26.1       Biobase_2.48.0
[67] BiocGenerics_0.34.0

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.1          bit64_4.0.5             knitr_1.30
  [4] rpart_4.1-15            doParallel_1.0.16       generics_0.0.2
  [7] cowplot_1.1.0           RSQLite_2.2.1           mice_3.11.0
 [10] europepmc_0.4           bit_4.0.4               enrichplot_1.8.1
 [13] lubridate_1.7.9         webshot_0.5.2           xml2_1.3.2
 [16] httpuv_1.5.4            assertthat_0.2.1        viridis_0.5.1
 [19] xfun_0.18               hms_0.5.3               promises_1.1.1
 [22] fansi_0.4.1             progress_1.2.2          caTools_1.18.0
 [25] dbplyr_1.4.4            km.ci_0.5-2             DBI_1.1.0
 [28] htmlwidgets_1.5.2       ellipsis_0.3.1          crosstalk_1.1.0.1
 [31] backports_1.1.10        vctrs_0.3.4             abind_1.4-5
 [34] withr_2.3.0             ggforce_0.3.2           triebeard_0.3.0
 [37] prettyunits_1.1.1       cluster_2.1.0           DOSE_3.14.0
 [40] crayon_1.3.4            ellipse_0.4.2           labeling_0.4.2
 [43] pkgconfig_2.0.3         tweenr_1.0.1            nnet_7.3-14
 [46] rlang_0.4.8             lifecycle_0.2.0         miniUI_0.1.1.1
 [49] downloader_0.4          modelr_0.1.8            cellranger_1.1.0
 [52] polyclip_1.10-0         lmtest_0.9-38           graph_1.66.0
 [55] urltools_1.7.3          KMsurv_0.1-5            zoo_1.8-8
 [58] reprex_0.3.0            base64enc_0.1-3         ggridges_0.5.2
 [61] GlobalOptions_0.1.2     png_0.1-7               viridisLite_0.3.0
 [64] rjson_0.2.20            bitops_1.0-6            KernSmooth_2.23-17
 [67] blob_1.2.1              shape_1.4.5             qvalue_2.20.0
 [70] manipulateWidget_0.10.1 jpeg_0.1-8.1            rstatix_0.6.0
 [73] gridGraphics_0.5-0      reactome.db_1.70.0      leaps_3.1
 [76] memoise_1.1.0           graphite_1.34.0         plyr_1.8.6
 [79] compiler_4.0.3          scatterpie_0.1.5        clue_0.3-57
 [82] cli_2.1.0               htmlTable_2.1.0         MASS_7.3-53
 [85] tidyselect_1.1.0        stringi_1.5.3           yaml_2.2.1
 [88] GOSemSim_2.14.2         latticeExtra_0.6-29     ggrepel_0.8.2
 [91] fastmatch_1.1-0         tools_4.0.3             rio_0.5.16
 [94] rstudioapi_0.11         foreach_1.5.1           foreign_0.8-80
 [97] scatterplot3d_0.3-41    farver_2.0.3            ggraph_2.0.3
[100] digest_0.6.26           rvcheck_0.1.8           BiocManager_1.30.10
[103] shiny_1.5.0             Rcpp_1.0.5              broom_0.7.2
[106] later_1.1.0.1           httr_1.4.2              colorspace_1.4-1
[109] polylabelr_0.2.0        rvest_0.3.6             fs_1.5.0
[112] splines_4.0.3           graphlayouts_0.7.0      ggplotify_0.0.5
[115] xtable_1.8-4            jsonlite_1.7.1          tidygraph_1.2.0
[118] flashClust_1.01-2       R6_2.4.1                pillar_1.4.6
[121] htmltools_0.5.0         mime_0.9                glue_1.4.2
[124] fastmap_1.0.1           BiocParallel_1.22.0     codetools_0.2-16
[127] fgsea_1.14.0            mvtnorm_1.1-1           curl_4.3
[130] gtools_3.8.2            zip_2.1.1               GO.db_3.11.4
[133] openxlsx_4.2.2          munsell_0.5.0           DO.db_2.9
[136] GetoptLong_1.0.4        iterators_1.0.13        haven_2.3.1
[139] gtable_0.3.0            tidytidbits_0.2.2
```
