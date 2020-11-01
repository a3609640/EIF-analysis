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

## Session Information
### System 1
```
Session info
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 3.6.3 (2020-02-29)
 os       Pop!_OS 18.04 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language en_US:en
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2020-10-29

Packages
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package          * version    date       lib source
 abind              1.4-5      2016-07-21 [1] CRAN (R 3.6.1)
 AnnotationDbi    * 1.46.1     2019-08-20 [1] Bioconductor
 assertthat         0.2.1      2019-03-21 [1] CRAN (R 3.6.0)
 backports          1.1.10     2020-09-15 [1] CRAN (R 3.6.3)
 base             * 3.6.3      2020-02-29 [4] local
 base64enc          0.1-3      2015-07-28 [1] CRAN (R 3.6.0)
 Biobase          * 2.44.0     2019-05-02 [1] Bioconductor
 BiocGenerics     * 0.30.0     2019-05-02 [1] Bioconductor
 BiocManager        1.30.10    2019-11-16 [1] CRAN (R 3.6.1)
 BiocParallel       1.18.1     2019-08-06 [1] Bioconductor
 bit                4.0.4      2020-08-04 [1] CRAN (R 3.6.3)
 bit64              4.0.5      2020-08-30 [1] CRAN (R 3.6.3)
 bitops             1.0-6      2013-08-17 [1] CRAN (R 3.6.0)
 blob               1.2.1      2020-01-20 [1] CRAN (R 3.6.1)
 broom              0.7.2      2020-10-20 [1] CRAN (R 3.6.3)
 car              * 3.0-10     2020-09-29 [1] CRAN (R 3.6.3)
 carData          * 3.0-4      2020-05-22 [1] CRAN (R 3.6.3)
 caTools            1.18.0     2020-01-17 [1] CRAN (R 3.6.1)
 cellranger         1.1.0      2016-07-27 [1] CRAN (R 3.6.1)
 checkmate        * 2.0.0      2020-02-06 [1] CRAN (R 3.6.1)
 circlize         * 0.4.10     2020-06-15 [1] CRAN (R 3.6.3)
 cli                2.1.0      2020-10-12 [1] CRAN (R 3.6.3)
 clue               0.3-57     2019-02-25 [1] CRAN (R 3.6.0)
 cluster            2.1.0      2019-06-19 [4] CRAN (R 3.6.1)
 clusterProfiler  * 3.12.0     2019-05-02 [1] Bioconductor
 codetools          0.2-16     2018-12-24 [4] CRAN (R 3.6.3)
 colorspace         1.4-1      2019-03-18 [1] CRAN (R 3.6.0)
 compiler           3.6.3      2020-02-29 [4] local
 ComplexHeatmap   * 2.0.0      2019-05-02 [1] Bioconductor
 corrplot         * 0.84       2017-10-16 [1] CRAN (R 3.6.0)
 cowplot            1.1.0      2020-09-08 [1] CRAN (R 3.6.3)
 crayon             1.3.4      2017-09-16 [1] CRAN (R 3.6.0)
 crosstalk          1.1.0.1    2020-03-13 [1] CRAN (R 3.6.1)
 curl               4.3        2019-12-02 [1] CRAN (R 3.6.1)
 data.table       * 1.13.2     2020-10-19 [1] CRAN (R 3.6.3)
 datasets         * 3.6.3      2020-02-29 [4] local
 DBI                1.1.0      2019-12-15 [1] CRAN (R 3.6.1)
 dbplyr             1.4.4      2020-05-27 [1] CRAN (R 3.6.3)
 dendextend       * 1.14.0     2020-08-26 [1] CRAN (R 3.6.3)
 descr            * 1.1.4      2018-01-19 [1] CRAN (R 3.6.1)
 digest             0.6.26     2020-10-17 [1] CRAN (R 3.6.3)
 DO.db              2.9        2019-09-08 [1] Bioconductor
 doParallel         1.0.16     2020-10-16 [1] CRAN (R 3.6.3)
 DOSE               3.10.2     2019-06-24 [1] Bioconductor
 dplyr            * 1.0.2      2020-08-18 [1] CRAN (R 3.6.3)
 ellipse            0.4.2      2020-05-27 [1] CRAN (R 3.6.3)
 ellipsis           0.3.1      2020-05-15 [1] CRAN (R 3.6.3)
 enrichplot         1.4.0      2019-05-02 [1] Bioconductor
 EnvStats         * 2.3.1      2018-07-15 [1] CRAN (R 3.6.1)
 eulerr           * 6.1.0      2020-03-09 [1] CRAN (R 3.6.1)
 europepmc          0.4        2020-05-31 [1] CRAN (R 3.6.3)
 facetscales      * 0.1.0.9000 2020-05-28 [1] Github (zeehio/facetscales@f374846)
 factoextra       * 1.0.7      2020-04-01 [1] CRAN (R 3.6.3)
 FactoMineR       * 2.3        2020-02-29 [1] CRAN (R 3.6.1)
 fansi              0.4.1      2020-01-08 [1] CRAN (R 3.6.1)
 farver             2.0.3      2020-01-16 [1] CRAN (R 3.6.1)
 fastmap            1.0.1      2019-10-08 [1] CRAN (R 3.6.1)
 fastmatch          1.1-0      2017-01-28 [1] CRAN (R 3.6.1)
 fgsea              1.10.1     2019-08-21 [1] Bioconductor
 flashClust         1.01-2     2012-08-21 [1] CRAN (R 3.6.1)
 forcats          * 0.5.0      2020-03-01 [1] CRAN (R 3.6.1)
 foreach            1.5.1      2020-10-15 [1] CRAN (R 3.6.3)
 foreign            0.8-76     2020-03-03 [4] CRAN (R 3.6.3)
 forestmodel      * 0.6.2      2020-07-19 [1] CRAN (R 3.6.3)
 forestplot       * 1.10       2020-07-16 [1] CRAN (R 3.6.3)
 Formula          * 1.2-4      2020-10-16 [1] CRAN (R 3.6.3)
 fs                 1.5.0      2020-07-31 [1] CRAN (R 3.6.3)
 generics           0.0.2      2018-11-29 [1] CRAN (R 3.6.0)
 GetoptLong         1.0.4      2020-10-19 [1] CRAN (R 3.6.3)
 ggforce            0.3.2      2020-06-23 [1] CRAN (R 3.6.3)
 ggfortify        * 0.4.11     2020-10-02 [1] CRAN (R 3.6.3)
 ggplot2          * 3.3.2      2020-06-19 [1] CRAN (R 3.6.3)
 ggplotify          0.0.5      2020-03-12 [1] CRAN (R 3.6.1)
 ggpubr           * 0.4.0      2020-06-27 [1] CRAN (R 3.6.3)
 ggraph             2.0.3      2020-05-20 [1] CRAN (R 3.6.3)
 ggrepel            0.8.2      2020-03-08 [1] CRAN (R 3.6.1)
 ggridges           0.5.2      2020-01-12 [1] CRAN (R 3.6.1)
 ggsignif         * 0.6.0      2019-08-08 [1] CRAN (R 3.6.1)
 ggthemes         * 4.2.0      2019-05-13 [1] CRAN (R 3.6.1)
 glmnet           * 4.0-2      2020-06-16 [1] CRAN (R 3.6.3)
 GlobalOptions      0.1.2      2020-06-10 [1] CRAN (R 3.6.3)
 glue               1.4.2      2020-08-27 [1] CRAN (R 3.6.3)
 GO.db              3.8.2      2019-06-09 [1] Bioconductor
 GOSemSim           2.10.0     2019-05-02 [1] Bioconductor
 gplots           * 3.1.0      2020-09-18 [1] CRAN (R 3.6.3)
 graph              1.62.0     2019-05-02 [1] Bioconductor
 graphics         * 3.6.3      2020-02-29 [4] local
 graphite           1.30.0     2019-05-02 [1] Bioconductor
 graphlayouts       0.7.0      2020-04-25 [1] CRAN (R 3.6.3)
 grDevices        * 3.6.3      2020-02-29 [4] local
 grid             * 3.6.3      2020-02-29 [4] local
 gridExtra        * 2.3        2017-09-09 [1] CRAN (R 3.6.0)
 gridGraphics       0.5-0      2020-02-25 [1] CRAN (R 3.6.1)
 gtable             0.3.0      2019-03-25 [1] CRAN (R 3.6.0)
 gtools             3.8.2      2020-03-31 [1] CRAN (R 3.6.3)
 haven              2.3.1      2020-06-01 [1] CRAN (R 3.6.3)
 Hmisc            * 4.4-1      2020-08-10 [1] CRAN (R 3.6.3)
 hms                0.5.3      2020-01-08 [1] CRAN (R 3.6.1)
 htmlTable          2.1.0      2020-09-16 [1] CRAN (R 3.6.3)
 htmltools          0.5.0      2020-06-16 [1] CRAN (R 3.6.3)
 htmlwidgets        1.5.2      2020-10-03 [1] CRAN (R 3.6.3)
 httpuv             1.5.4      2020-06-06 [1] CRAN (R 3.6.3)
 httr               1.4.2      2020-07-20 [1] CRAN (R 3.6.3)
 igraph           * 1.2.6      2020-10-06 [1] CRAN (R 3.6.3)
 IRanges          * 2.18.3     2019-09-24 [1] Bioconductor
 iterators          1.0.13     2020-10-15 [1] CRAN (R 3.6.3)
 jpeg               0.1-8.1    2019-10-24 [1] CRAN (R 3.6.1)
 jsonlite           1.7.1      2020-09-07 [1] CRAN (R 3.6.3)
 KEGG.db          * 3.2.3      2019-09-22 [1] Bioconductor
 KernSmooth         2.23-17    2020-04-26 [4] CRAN (R 3.6.3)
 km.ci              0.5-2      2009-08-30 [1] CRAN (R 3.6.1)
 KMsurv             0.1-5      2012-12-03 [1] CRAN (R 3.6.1)
 knitr              1.30       2020-09-22 [1] CRAN (R 3.6.3)
 labeling           0.4.2      2020-10-20 [1] CRAN (R 3.6.3)
 later              1.1.0.1    2020-06-05 [1] CRAN (R 3.6.3)
 lattice          * 0.20-41    2020-04-02 [4] CRAN (R 3.6.3)
 latticeExtra       0.6-29     2019-12-19 [1] CRAN (R 3.6.1)
 leaps              3.1        2020-01-16 [1] CRAN (R 3.6.1)
 lemon            * 0.4.5      2020-06-08 [1] CRAN (R 3.6.3)
 lifecycle          0.2.0      2020-03-06 [1] CRAN (R 3.6.1)
 limma            * 3.40.6     2019-07-26 [1] Bioconductor
 lmtest             0.9-38     2020-09-09 [1] CRAN (R 3.6.3)
 lubridate          1.7.9      2020-06-08 [1] CRAN (R 3.6.3)
 magrittr         * 1.5        2014-11-22 [1] CRAN (R 3.6.0)
 manipulateWidget   0.10.1     2020-02-24 [1] CRAN (R 3.6.1)
 MASS               7.3-51.6   2020-04-26 [4] CRAN (R 3.6.3)
 Matrix           * 1.2-18     2019-11-27 [4] CRAN (R 3.6.1)
 memoise            1.1.0      2017-04-21 [1] CRAN (R 3.6.0)
 methods          * 3.6.3      2020-02-29 [4] local
 mice               3.11.0     2020-08-05 [1] CRAN (R 3.6.3)
 mime               0.9        2020-02-04 [1] CRAN (R 3.6.1)
 miniUI             0.1.1.1    2018-05-18 [1] CRAN (R 3.6.1)
 missMDA          * 1.17       2020-05-19 [1] CRAN (R 3.6.3)
 modelr             0.1.8      2020-05-19 [1] CRAN (R 3.6.3)
 munsell            0.5.0      2018-06-12 [1] CRAN (R 3.6.0)
 mvtnorm            1.1-1      2020-06-09 [1] CRAN (R 3.6.3)
 nnet               7.3-14     2020-04-26 [4] CRAN (R 3.6.3)
 nortest          * 1.0-4      2015-07-30 [1] CRAN (R 3.6.1)
 openxlsx           4.2.2      2020-09-17 [1] CRAN (R 3.6.3)
 org.Hs.eg.db     * 3.8.2      2019-06-09 [1] Bioconductor
 parallel         * 3.6.3      2020-02-29 [4] local
 pca3d            * 0.10.2     2020-10-02 [1] CRAN (R 3.6.3)
 pheatmap         * 1.0.12     2019-01-04 [1] CRAN (R 3.6.1)
 pillar             1.4.6      2020-07-10 [1] CRAN (R 3.6.3)
 pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 3.6.1)
 plotmo           * 3.6.0      2020-09-13 [1] CRAN (R 3.6.3)
 plotrix          * 3.7-8      2020-04-16 [1] CRAN (R 3.6.3)
 plyr               1.8.6      2020-03-03 [1] CRAN (R 3.6.1)
 png                0.1-7      2013-12-03 [1] CRAN (R 3.6.0)
 polyclip           1.10-0     2019-03-14 [1] CRAN (R 3.6.1)
 polylabelr         0.2.0      2020-04-19 [1] CRAN (R 3.6.3)
 prettyunits        1.1.1      2020-01-24 [1] CRAN (R 3.6.1)
 progress           1.2.2      2019-05-16 [1] CRAN (R 3.6.0)
 promises           1.1.1      2020-06-09 [1] CRAN (R 3.6.3)
 purrr            * 0.3.4      2020-04-17 [1] CRAN (R 3.6.3)
 qvalue             2.16.0     2019-05-02 [1] Bioconductor
 R6                 2.4.1      2019-11-12 [1] CRAN (R 3.6.1)
 rappdirs           0.3.1      2016-03-28 [1] CRAN (R 3.6.1)
 RColorBrewer     * 1.1-2      2014-12-07 [1] CRAN (R 3.6.0)
 Rcpp               1.0.5      2020-07-06 [1] CRAN (R 3.6.3)
 reactome.db        1.68.0     2019-09-08 [1] Bioconductor
 ReactomePA       * 1.28.0     2019-05-02 [1] Bioconductor
 readr            * 1.4.0      2020-10-05 [1] CRAN (R 3.6.3)
 readxl           * 1.3.1      2019-03-13 [1] CRAN (R 3.6.1)
 reprex             0.3.0      2019-05-16 [1] CRAN (R 3.6.1)
 reshape2         * 1.4.4      2020-04-09 [1] CRAN (R 3.6.3)
 rgl              * 0.100.54   2020-04-14 [1] CRAN (R 3.6.3)
 rio                0.5.16     2018-11-26 [1] CRAN (R 3.6.1)
 rjson              0.2.20     2018-06-08 [1] CRAN (R 3.6.0)
 rlang              0.4.8      2020-10-08 [1] CRAN (R 3.6.3)
 rpart              4.1-15     2019-04-12 [4] CRAN (R 3.6.1)
 RSQLite            2.2.1      2020-09-30 [1] CRAN (R 3.6.3)
 rstatix            0.6.0      2020-06-18 [1] CRAN (R 3.6.3)
 rstudioapi         0.11       2020-02-07 [1] CRAN (R 3.6.1)
 rvcheck            0.1.8      2020-03-01 [1] CRAN (R 3.6.1)
 rvest              0.3.6      2020-07-25 [1] CRAN (R 3.6.3)
 S4Vectors        * 0.22.1     2019-09-09 [1] Bioconductor
 scales           * 1.1.1      2020-05-11 [1] CRAN (R 3.6.3)
 scatterplot3d      0.3-41     2018-03-14 [1] CRAN (R 3.6.1)
 sessioninfo        1.1.1      2018-11-05 [1] CRAN (R 3.6.0)
 shape              1.4.5      2020-09-13 [1] CRAN (R 3.6.3)
 shiny              1.5.0      2020-06-23 [1] CRAN (R 3.6.3)
 splines            3.6.3      2020-02-29 [4] local
 stats            * 3.6.3      2020-02-29 [4] local
 stats4           * 3.6.3      2020-02-29 [4] local
 stringi            1.5.3      2020-09-09 [1] CRAN (R 3.6.3)
 stringr          * 1.4.0      2019-02-10 [1] CRAN (R 3.6.0)
 survival         * 3.2-7      2020-09-28 [1] CRAN (R 3.6.3)
 survivalAnalysis * 0.1.3      2020-06-29 [1] CRAN (R 3.6.3)
 survminer        * 0.4.8      2020-07-25 [1] CRAN (R 3.6.3)
 survMisc         * 0.5.5      2018-07-05 [1] CRAN (R 3.6.1)
 TeachingDemos    * 2.12       2020-04-07 [1] CRAN (R 3.6.3)
 tibble           * 3.0.4      2020-10-12 [1] CRAN (R 3.6.3)
 tidygraph          1.2.0      2020-05-12 [1] CRAN (R 3.6.3)
 tidyr            * 1.1.2      2020-08-27 [1] CRAN (R 3.6.3)
 tidyselect         1.1.0      2020-05-11 [1] CRAN (R 3.6.3)
 tidytidbits        0.2.2      2020-04-15 [1] CRAN (R 3.6.3)
 tidyverse        * 1.3.0      2019-11-21 [1] CRAN (R 3.6.3)
 tools              3.6.3      2020-02-29 [4] local
 triebeard          0.3.0      2016-08-04 [1] CRAN (R 3.6.1)
 tweenr             1.0.1      2018-12-14 [1] CRAN (R 3.6.1)
 UpSetR             1.4.0      2019-05-22 [1] CRAN (R 3.6.1)
 urltools           1.7.3      2019-04-14 [1] CRAN (R 3.6.1)
 utils            * 3.6.3      2020-02-29 [4] local
 vcd              * 1.4-8      2020-09-21 [1] CRAN (R 3.6.3)
 vctrs              0.3.4      2020-08-29 [1] CRAN (R 3.6.3)
 vip              * 0.2.2      2020-04-06 [1] CRAN (R 3.6.3)
 viridis            0.5.1      2018-03-29 [1] CRAN (R 3.6.0)
 viridisLite        0.3.0      2018-02-01 [1] CRAN (R 3.6.0)
 webshot            0.5.2      2019-11-22 [1] CRAN (R 3.6.1)
 withr              2.3.0      2020-09-22 [1] CRAN (R 3.6.3)
 xfun               0.18       2020-09-29 [1] CRAN (R 3.6.3)
 xml2               1.3.2      2020-04-23 [1] CRAN (R 3.6.3)
 xtable             1.8-4      2019-04-21 [1] CRAN (R 3.6.0)
 zip                2.1.1      2020-08-27 [1] CRAN (R 3.6.3)
 zoo                1.8-8      2020-05-02 [1] CRAN (R 3.6.3)
```

### System 2
```
```
