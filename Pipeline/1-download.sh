#!/bin/sh


## download all datasets from the following weblinks

### create the directory to store all downloaded datasets
data_file_directory="$HOME/Downloads/Test"
mkdir -p $data_file_directory


### TCGA and GTEX DATA
#### TCGA CNV dataset (thresholded)
wget https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz -P $data_file_directory

#### TCGA CNV dataset
wget https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz -P $data_file_directory

#### TCGA CNV ratio dataset
wget https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz -P $data_file_directory

#### TCGA RNA-Seq dataset
wget https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz -P $data_file_directory

#### TCGA sample type annotation
wget https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz -P $data_file_directory

#### TCGA OS data ##
wget https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp -P $data_file_directory

#### TCGA and GTEX RNA-Seq dataset
wget https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz -P $data_file_directory

#### TCGA and GTEX sample type annotation
wget https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz -P $data_file_directory

### CPTAC DATA

#### CPTAC LUAD RNA-Seq data (Gillette et al., 2020)
wget https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/RNA.xlsx -P $data_file_directory

#### CPTAC LUAD Proteomics (Gillette et al., 2020)
wget https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/Protein.xlsx -P $data_file_directory

#### CPTAC LUAD Proteomics
wget https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_Proteome_CDAP_Protein_Report.r1/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv -P $data_file_directory

#### CPTAC LUAD Phosproteomics (Gillette et al., 2020)
wget https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/Phos.xlsx -P $data_file_directory

#### CPTAC LUAD Sample Annotation
wget https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx -P $data_file_directory

#### CPTAC Clinical Data
wget https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx -P $data_file_directory

### CCLE DATA

#### CCLE RNA-Seq data from DepMap Public 20Q4 20Q3
#wget -O CCLE_expression_full.csv https://ndownloader.figshare.com/files/#27902097 -P $data_file_directory #DepMap Public 21Q2

wget https://ndownloader.figshare.com/files/24613349 -O "$data_file_directory/CCLE_expression_full.csv" #DepMap Public 20Q3

#### CCLE annotation data
#wget -O sample_info.csv https://ndownloader.figshare.com/files/27902376 -P #$data_file_directory #DepMap Public 21Q2

wget https://ndownloader.figshare.com/files/24613394 -O "$data_file_directory/sample_info.csv" #DepMap Public 20Q3

#### CCLE proteomics data
wget https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz -P $data_file_directory


gunzip $data_file_directory/*.gz
