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
    * Pop!_OS 17.10

2. (verification in progress) PowerSpec G460 desktop computer
    * Intel i7-8700k CPU
    * 64GB RAM (DDR4-3200, non-ECC)
    * Intel M.2 SATA SSD
    * Samsung NVMe Evo+ SSD
    * Windows 10 Pro
    * RStudio for Windows
    * R 4.0.3


## DATASETS
Please download all datasets from the following weblinks and store them at the following directory. The script will read all datasets from there.
```
~/Download/
```
This project makes use of data files that must be fetched from remote sources.  The necessary links are provided below.  Note that the .gz files must be unzipped after they are downloaded.

### Xena DATA
https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz
https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz
https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
https://toil.xenahubs.net/download/GTEX_phenotype.gz
https://toil.xenahubs.net/download/gtex_RSEM_Hugo_norm_count.gz
https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
https://xenabrowser.net/datapages/?dataset=Survival_SupplementalTable_S1_20171025_xena_sp&host=https%3A%2F%2Fpancanatlas.xenahubs.net

### CPTAC DATA

#### CPTAC LUAD Proteomics
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_Proteome_CDAP_Protein_Report.r1/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv

#### CPTAC LUAD Phosproteomics
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_Phosphoproteome_CDAP_Protein_Report.r1/CPTAC3_Lung_Adeno_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv

#### CPTAC LUAD Sample Annotation
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx

#### CPTAC Clinical Data
https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx


## LIBRARIES
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
                       "Hmisc",
                       "ggfortify",
                       "ggplot2",
                       "ggpubr",
                       "ggsignif",
                       "ggthemes",
                       "glmnet",
                       "gplots",
                       "gridExtra",
                       "igraph",
                       "KEGG.db",
                       "lemon",
                       "limma",
                       "missMDA",
                       "org.Hs.eg.db",
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

## OUTPUT
The script will save all the generated plots into the following output folder.
Please generate the output directory.
```
mkdir -p ~/Documents/EIF_output/{All,CNV,Expression,PCA,KM,Cox,Heatmap,CPTAC}
```

## INSTRUCTIONS
1. Install R (a version verified in one of the systems listed above) and RStudio at the recommended Systems.
2. Download all datasets from the listed web links above and store them under the following directory.
```
~/Download/
```
3. Install all the R packages listed above.
4. Create output directories
```
~/Documents/EIF_output/{All,CNV,Expression,PCA,KM,Cox,Heatmap,CPTAC}
```
5. Open the EIFanalysis.R file in RStudio, and run it with the "source" button.
6. The EIFanalysis.R will produce all the figures from the EIF manuscript and store them in the output folders.
