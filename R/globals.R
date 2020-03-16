EIF.gene <- c("EIF4E",
              "EIF4G1",
              "EIF4A1",
              "EIF4EBP1")
names(EIF.gene) <- EIF.gene

black_bold_tahoma_7 <- element_text(color = "black",
                                    face  = "bold",
                                    size  = 7)
black_bold_tahoma_12 <- element_text(color = "black",
                                     face  = "bold",
                                     size  = 12)
black_bold_tahoma_16 <- element_text(color = "black",
                                     face  = "bold",
                                     size  = 16)
black_bold_tahoma_16_right <- element_text(color = "black",
                                           face  = "bold",
                                           size  = 16,
                                           angle = 90)
black_bold_tahoma_16_45 <- element_text(color = "black",
                                        face  = "bold",
                                        size  = 16,
                                        angle = 45,
                                        hjust = 1)
black_bold_tahoma_16_90 <- element_text(color = "black",
                                        face  = "bold",
                                        size  = 16,
                                        angle = 90,
                                        hjust = 1,
                                        vjust = 0.5)
black_bold_tahoma_18 <- element_text(color = "black",
                                     face  = "bold",
                                     size  = 18)

# data.table = FALSE gives data.frame
# download:
#  https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
#  https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
#  https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
#  https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz

EbPanCanIlluminaHiseq <- fread("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", data.table = FALSE)
Gistic2CopyNumberAllThresholdedByGenes <- fread("Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", data.table = FALSE)
TcgaPhenotypeDenseData <- readr::read_tsv("TCGA_phenotype_denseDataOnlyDownload.tsv")
TcgaTargetGtexPhenoType <- read_tsv("TcgaTargetGTEX_phenotype.txt")
TcgaTargetGtexRsemHugoNormCount <- fread("TcgaTargetGtex_RSEM_Hugo_norm_count", data.table = FALSE) 
