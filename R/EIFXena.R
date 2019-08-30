## the following script perform PCA on RNA-Seq data of seven EIF genes
## from SKCM amd GTEX dataset with R package "ggfortify".
library(EnvStats)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gridExtra)
library(reshape2)
library(readr)
library(survival)
library(survMisc)
library(survminer)

## read.csv will transform characters into factors
get.EIF.TCGA.GTEX.RNAseq.long <- function () {
  EIF.TCGA.GTEX <- read_csv("project-data/EIFTCGAGTEX.csv")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c(
                                           "sample",
                                           "study",
                                           "sample.type",
                                           "primary.disease",
                                           "variable",
                                           "value"
                                            )
  EIF.TCGA.GTEX.RNAseq.long <- na.omit(EIF.TCGA.GTEX.RNAseq.long)
  tumor.type <- c(
                  "Metastatic",
                  "Primary Tumor",
                  #"Recurrent Tumor",
                  "Solid Tissue Normal",
                  "Normal Tissue"
                  )
  EIF.TCGA.GTEX.RNAseq.long <- EIF.TCGA.GTEX.RNAseq.long[
    EIF.TCGA.GTEX.RNAseq.long$sample.type %in% tumor.type, ]
  EIF.TCGA.GTEX.RNAseq.long <- droplevels(EIF.TCGA.GTEX.RNAseq.long)
  EIF.TCGA.GTEX.RNAseq.long$sample.type <- factor(
    EIF.TCGA.GTEX.RNAseq.long$sample.type, levels = tumor.type)
  return(EIF.TCGA.GTEX.RNAseq.long)
}

##
get.EIF.TCGA.GTEX.RNAseq.tissue <- function (tissue) {
  EIF.TCGA.GTEX <- read_csv("project-data/EIFTCGAGTEX.csv")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c(
    "sample",
    "study",
    "sample.type",
    "primary.disease",
    "variable",
    "value"
  )
  EIF.TCGA.GTEX.tissue <- EIF.TCGA.GTEX.RNAseq.long[grep(tissue,
    EIF.TCGA.GTEX.RNAseq.long$primary.disease), ]
  EIF.TCGA.GTEX.tissue <- na.omit(EIF.TCGA.GTEX.tissue)
  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
  #  "Recurrent Tumor",
    "Solid Tissue Normal",
    "Normal Tissue"
  )
  EIF.TCGA.GTEX.tissue <- EIF.TCGA.GTEX.tissue[
    EIF.TCGA.GTEX.tissue$sample.type %in% tumor.type, ]
  EIF.TCGA.GTEX.tissue <- droplevels(EIF.TCGA.GTEX.tissue)
  EIF.TCGA.GTEX.tissue$sample.type <- factor(
    EIF.TCGA.GTEX.tissue$sample.type, levels = tumor.type)
  return(EIF.TCGA.GTEX.tissue)
}
# x <- get.EIF.TCGA.GTEX.RNAseq.tissue("Skin")

##
get.EIF.TCGA.RNAseq.long <- function (x) {
  EIF.TCGA.GTEX <- read_csv("project-data/EIFTCGAGTEX.csv")
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA', ]
  EIF.TCGA.1 <- EIF.TCGA[!is.na(EIF.TCGA$primary_disease_or_tissue), ]
  EIF.TCGA.RNAseq.long <- melt(EIF.TCGA.1[, 1:10])
  colnames(EIF.TCGA.RNAseq.long) <- c(
                                      "sample",
                                      "study",
                                      "sample.type",
                                      "primary.disease",
                                      "variable",
                                      "value"
                                            )
  EIF.TCGA.RNAseq.long <- EIF.TCGA.RNAseq.long[
    EIF.TCGA.RNAseq.long$variable %in% x, ]
  return(EIF.TCGA.RNAseq.long)
}

##
get.EIF.GTEX.RNAseq.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path(
    "project-data",
    "EIFTCGAGTEX.csv"),
    header = TRUE,
    sep = ",")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c(
    "sample",
    "study",
    "sample.type",
    "primary.disease",
    "variable",
    "value")
  EIF.GTEX.RNAseq.long <- EIF.TCGA.GTEX.RNAseq.long[
    EIF.TCGA.GTEX.RNAseq.long$study == 'GTEX', ]
  EIF.GTEX.RNAseq.long <- na.omit(EIF.GTEX.RNAseq.long)
  #  EIF.GTEX.RNAseq.long <- EIF.GTEX.RNAseq.long[EIF.GTEX.RNAseq.long$sample.type %in% tumor.type,]
  EIF.GTEX.RNAseq.long <- droplevels(EIF.GTEX.RNAseq.long)
  return(EIF.GTEX.RNAseq.long)
}

##
get.EIF.TCGA.GTEX.score.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data",
    "EIFTCGAGTEX.csv"),
    header = TRUE,
    sep = ",")
  EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX
  EIF.TCGA.GTEX.score$EIF4A1 <-
    (EIF.TCGA.GTEX$EIF4A1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <-
    "EIF4A1:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4G1 <-
    (EIF.TCGA.GTEX$EIF4G1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <-
    "EIF4G1:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4EBP1 <-
    (EIF.TCGA.GTEX$EIF4EBP1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <-
    "EIF4EBP1:EIF4E ratio"
  EIF.TCGA.GTEX.score$RPS6KB1 <-
    (EIF.TCGA.GTEX$RPS6KB1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <-
    "RPS6KB1:EIF4E ratio"
  EIF.TCGA.GTEX.score$MYC <-
    (EIF.TCGA.GTEX$MYC - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <-
    "MYC:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4E <-
    (EIF.TCGA.GTEX$EIF4E - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <-
    "EIF4E:EIF4E ratio"
  EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:10])
  colnames(EIF.TCGA.GTEX.score.long) <-
    c("sample",
      "study",
      "sample.type",
      "primary.disease",
      "variable",
      "value")
  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
    "Recurrent Tumor",
    "Solid Tissue Normal",
    "Normal Tissue",
    "Cell Line"
  )
  EIF.TCGA.GTEX.score.long <- EIF.TCGA.GTEX.score.long[EIF.TCGA.GTEX.score.long$sample.type %in% tumor.type, ]
  EIF.TCGA.GTEX.score.long <- droplevels(EIF.TCGA.GTEX.score.long)
  EIF.TCGA.GTEX.score.long$sample.type <-
    factor(EIF.TCGA.GTEX.score.long$sample.type,
      levels = tumor.type)
  return(EIF.TCGA.GTEX.score.long)
}

##
get.EIF.TCGA.GTEX.score.tissue <- function (tissue) {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data",
    "EIFTCGAGTEX.csv"),
    header = TRUE,
    sep = ",")
  EIF.TCGA.GTEX.tissue <- EIF.TCGA.GTEX[grep(tissue,
    EIF.TCGA.GTEX$primary_disease_or_tissue), ]
  EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.tissue
  EIF.TCGA.GTEX.score$EIF4A1 <-
    (EIF.TCGA.GTEX.tissue$EIF4A1 - EIF.TCGA.GTEX.tissue$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <-
    "EIF4A1:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4G1 <-
    (EIF.TCGA.GTEX.tissue$EIF4G1 - EIF.TCGA.GTEX.tissue$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <-
    "EIF4G1:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4EBP1 <-
    (EIF.TCGA.GTEX.tissue$EIF4EBP1 - EIF.TCGA.GTEX.tissue$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <-
    "EIF4EBP1:EIF4E ratio"
  EIF.TCGA.GTEX.score$RPS6KB1 <-
    (EIF.TCGA.GTEX.tissue$RPS6KB1 - EIF.TCGA.GTEX.tissue$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <-
    "RPS6KB1:EIF4E ratio"
  EIF.TCGA.GTEX.score$MYC <-
    (EIF.TCGA.GTEX.tissue$MYC - EIF.TCGA.GTEX.tissue$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <-
    "MYC:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4E <-
    (EIF.TCGA.GTEX.tissue$EIF4E - EIF.TCGA.GTEX.tissue$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <-
    "EIF4E:EIF4E ratio"
  EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:10])
  colnames(EIF.TCGA.GTEX.score.long) <-
    c("sample",
      "study",
      "sample.type",
      "primary.disease",
      "variable",
      "value")
  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
 #   "Recurrent Tumor",
    "Solid Tissue Normal",
    "Normal Tissue"
  )
  EIF.TCGA.GTEX.score.long <- EIF.TCGA.GTEX.score.long[EIF.TCGA.GTEX.score.long$sample.type %in% tumor.type, ]
  EIF.TCGA.GTEX.score.long <- droplevels(EIF.TCGA.GTEX.score.long)
  EIF.TCGA.GTEX.score.long$sample.type <-
    factor(EIF.TCGA.GTEX.score.long$sample.type,
      levels = tumor.type)
  return(EIF.TCGA.GTEX.score.long)
}

##
get.EIF.TCGA.score.long <- function (x) {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data",
    "EIFTCGAGTEX.csv"),
    header = TRUE,
    sep = ",")
  EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX
  EIF.TCGA.GTEX.score$RPS6KB1 <-
    (EIF.TCGA.GTEX$EIF4G1 - EIF.TCGA.GTEX$EIF4EBP1)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <-
    "EIF4G1:EIF4EBP1"
  EIF.TCGA.GTEX.score$EIF4A1 <-
    (EIF.TCGA.GTEX$EIF4A1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <-
    "EIF4A1:EIF4E"
  EIF.TCGA.GTEX.score$EIF4G1 <-
    (EIF.TCGA.GTEX$EIF4G1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <-
    "EIF4G1:EIF4E"
  EIF.TCGA.GTEX.score$EIF4EBP1 <-
    (EIF.TCGA.GTEX$EIF4EBP1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <-
    "EIF4EBP1:EIF4E"
  EIF.TCGA.GTEX.score$MYC <-
    (EIF.TCGA.GTEX$MYC - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <-
    "MYC:EIF4E"
  EIF.TCGA.GTEX.score$EIF4E <-
    (EIF.TCGA.GTEX$EIF4E - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <-
    "EIF4E:EIF4E"
  EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:10])
  colnames(EIF.TCGA.GTEX.score.long) <-
    c("sample",
      "study",
      "sample.type",
      "primary.disease",
      "variable",
      "value")
  EIF.TCGA.score.long <-
    EIF.TCGA.GTEX.score.long[EIF.TCGA.GTEX.score.long$study == 'TCGA', ]
  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
    # "Recurrent Tumor",
    # "Normal Tissue",
    "Solid Tissue Normal"
  )
  EIF.TCGA.score.long <-
    EIF.TCGA.score.long[EIF.TCGA.score.long$sample.type %in% tumor.type, ]
  EIF.TCGA.score.long <- droplevels(EIF.TCGA.score.long)
  EIF.TCGA.score.long$sample.type <- factor(
    EIF.TCGA.score.long$sample.type, 
    levels = tumor.type
  )
  return(EIF.TCGA.score.long)
}

##
get.EIF.GTEX.score.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data",
    "EIFTCGAGTEX.csv"),
    header = TRUE,
    sep = ",")
  EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX
  EIF.TCGA.GTEX.score$EIF4A1 <-
    (EIF.TCGA.GTEX$EIF4A1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <-
    "EIF4A1:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4G1 <-
    (EIF.TCGA.GTEX$EIF4G1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <-
    "EIF4G1:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4EBP1 <-
    (EIF.TCGA.GTEX$EIF4EBP1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <-
    "EIF4EBP1:EIF4E ratio"
  EIF.TCGA.GTEX.score$RPS6KB1 <-
    (EIF.TCGA.GTEX$RPS6KB1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <-
    "RPS6KB1:EIF4E ratio"
  EIF.TCGA.GTEX.score$MYC <-
    (EIF.TCGA.GTEX$MYC - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <-
    "MYC:EIF4E ratio"
  EIF.TCGA.GTEX.score$EIF4E <-
    (EIF.TCGA.GTEX$EIF4E - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <-
    "EIF4E:EIF4E ratio"
  EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:10])
  colnames(EIF.TCGA.GTEX.score.long) <-
    c("sample",
      "study",
      "sample.type",
      "primary.disease",
      "variable",
      "value")
  EIF.GTEX.score.long <-
    EIF.TCGA.GTEX.score.long[EIF.TCGA.GTEX.score.long$study == 'GTEX', ]
  EIF.GTEX.score.long <- na.omit(EIF.GTEX.score.long)
  EIF.GTEX.score.long <- droplevels(EIF.GTEX.score.long)
  return(EIF.GTEX.score.long)
}

##
EIF.TCGA.GTEX <- read_csv("project-data/EIFTCGAGTEX.csv")
EIF.TCGA.GTEX <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$EIF4E != 0, ]
EIF.TCGA.GTEX <-  subset(EIF.TCGA.GTEX, subset = study %in% c('TCGA', "GTEX"))
EIF.TCGA.GTEX$sample_type <- as.factor(EIF.TCGA.GTEX$sample_type)
levels(EIF.TCGA.GTEX$sample_type)
# EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA', ]
tumor.type <- c(
  "Metastatic",
  "Primary Tumor",
  #"Recurrent Tumor",
  "Solid Tissue Normal",
  "Normal Tissue"
)
EIF.TCGA.GTEX <-  subset(EIF.TCGA.GTEX, subset = sample_type %in% tumor.type)
EIF.TCGA.GTEX <- droplevels(EIF.TCGA.GTEX)

cor.test(EIF.TCGA.GTEX$EIF4A1, EIF.TCGA.GTEX$MYC, method = "pearson")
black_bold_tahoma_16 <- element_text(
  color  = "black",
  face   = "bold",
  family = "Tahoma",
  size   = 16
)

ggscatter(EIF.TCGA.GTEX, 
          x        = "EIF4EBP1", 
          y        = "MYC", 
          size     = 0.3,
          color    = "sample_type", 
          palette  = "jco",
          facet.by = "sample_type", #scales = "free_x",
          add      = "reg.line", 
          conf.int = TRUE) +
  stat_cor(aes(color   = "sample_type"), 
               method  = "pearson", 
               label.y = 6)

EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA', ]
ggscatter(EIF.TCGA, 
          x        = "EIF4EBP1", 
          y        = "MYC", 
          size     = 0.3,
          color    = "primary_disease_or_tissue", 
          palette  = "jco",
          facet.by = "primary_disease_or_tissue", #scales = "free_x",
          add      = "reg.line", 
          conf.int = TRUE) +
  stat_cor(aes(color   = "primary_disease_or_tissue"), 
               method  = "pearson", 
               label.y = 6)


    ##
plotEIF.RNAseq.TCGA.GTEX <-  function (x) {
  name <- deparse(substitute(x))
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 9
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 9,
    angle  = 45,
    hjust  = 1
  )
  p1 <- ggplot(data = x,
    aes(x     = sample.type,
        y     = value,
        fill = sample.type,
        color = sample.type)) +
    stat_n_text(size = 6, fontface = "bold") + 
    facet_grid(~ variable,
               scales = "free",
               space  = "free") +
    facet_wrap(~ variable, 
               ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      width    = .1,
      color    = "black",
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
         y = paste("log2(TPM)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    )
  p1 <- p1 +  stat_compare_means(
    comparisons = list(
                      c("Metastatic", "Normal Tissue"),
                      c("Primary Tumor", "Normal Tissue"),
                      c("Metastatic", "Primary Tumor")
                      ),
    method = "t.test", label = "p.signif")
  print(p1)
  p2 <- ggplot(data = x,
    aes(x     = variable,
        y     = value,
        color = variable)) +
    stat_n_text() + 
    facet_grid(~ sample.type,
      scales = "free",
      space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
         y = paste("log2(value)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    )
  p2 <- p2 + stat_compare_means(
    comparisons = list(
                      c("EIF4A1", "EIF4E"),
                      c("EIF4G1", "EIF4E"),
                      c("EIF4EBP1", "EIF4E"),
                      c("RPS6KB1", "EIF4E"),
                      c("MYC", "EIF4E")
                      ),
    method = "t.test", 
    label = "p.signif")
  print(p2)
}

##
plotEIF.RNAseq.TCGA.GTEX.tissue <-  function (tissue) {
  get.EIF.TCGA.GTEX.RNAseq.tissue <- function (tissue) {
    EIF.TCGA.GTEX <- read_csv("project-data/EIFTCGAGTEX.csv")
    EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
    colnames(EIF.TCGA.GTEX.RNAseq.long) <- c(
      "sample",
      "study",
      "sample.type",
      "primary.disease",
      "variable",
      "value"
    )
    EIF.TCGA.GTEX.tissue <- EIF.TCGA.GTEX.RNAseq.long[grep(tissue,
      EIF.TCGA.GTEX.RNAseq.long$primary.disease), ]
    EIF.TCGA.GTEX.tissue <- na.omit(EIF.TCGA.GTEX.tissue)
    tumor.type <- c(
      "Metastatic",
      "Primary Tumor",
      # "Recurrent Tumor",
      "Solid Tissue Normal",
      "Normal Tissue"
    )
    EIF.TCGA.GTEX.tissue <- EIF.TCGA.GTEX.tissue[
      EIF.TCGA.GTEX.tissue$sample.type %in% tumor.type, ]
    EIF.TCGA.GTEX.tissue <- droplevels(EIF.TCGA.GTEX.tissue)
    EIF.TCGA.GTEX.tissue$sample.type <- factor(
      EIF.TCGA.GTEX.tissue$sample.type, levels = tumor.type)
    return(EIF.TCGA.GTEX.tissue)
  }
  x <- get.EIF.TCGA.GTEX.RNAseq.tissue (tissue)
  name <- deparse(substitute(x))
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    size   = 9
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    size   = 9,
    angle  = 45,
    hjust  = 1
  )
  p1 <- ggplot(data = x,
    aes(x     = sample.type,
        y     = value,
        color = sample.type)) +
    stat_n_text() + 
    facet_grid(~ variable,
               scales = "free",
               space  = "free") +
    facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = tissue,
         y = paste("log2(TPM)")) +
    ylim(7, 20)+
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    )
  p1 <- p1 +  stat_compare_means(
    comparisons = list(
      c("Metastatic", "Normal Tissue"),
      c("Primary Tumor", "Normal Tissue"),
      c("Metastatic", "Primary Tumor")
    ),
    method = "t.test", label = "p.signif")
  print(p1)
  p2 <- ggplot(data = x,
    aes(x     = variable,
        y     = value,
        color = variable)) +
    stat_n_text() + 
    facet_grid(~ sample.type,
      scales = "free",
      space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
  #  geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = tissue,
         y = paste("log2(value)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    )
  p2 <- p2 + stat_compare_means(
    comparisons = list(
      c("EIF4A1", "EIF4E"),
      c("EIF4G1", "EIF4E"),
      c("EIF4EBP1", "EIF4E"),
      c("RPS6KB1", "EIF4E"),
      c("MYC", "EIF4E")
    ),
    method = "t.test", 
    label = "p.signif")
  print(p2)
}

##
plot.box.EIF.RNAseq.TCGA <-  function (x) {
  x <- x[x$value != 0, ]
  # name <- deparse(substitute(x))
  # y <- x[x$variable == "EIF4E", ]
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  # reorder bars by explicitly ordering factor levels 
  x$primary.disease <- as.factor(x$primary.disease)
  mean <- within(x[x$variable == "EIF4E", ], # TCGAstudy is one column in df2
    primary.disease <- reorder(primary.disease, value, median))
  mean$primary.disease <- as.factor(mean$primary.disease)
  neworder <- levels(mean$primary.disease)
  x.ordered <- factor(x$primary.disease, levels = neworder)
  
  p1 <- ggplot(data = x,
    aes(x     = x.ordered,
        y     = value,
        color = variable)) +
    stat_n_text(size = 5, fontface = "bold", hjust = 0) + 
    # geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
         size  = .75,
     # width    = .5,
      position = position_dodge(width = .9)
    ) +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "primary disease",
         y = paste("log2(RNA counts)")) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.title    = element_blank(),
      legend.text     = black_bold_tahoma_12,
      legend.position = "top",
      strip.text      = black_bold_tahoma_12
    )
  print(p1)

}

##
plotEIF.RNAseq.TCGA <-  function (x) {
  name <- deparse(substitute(x))
  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
    #  "Recurrent Tumor",
    "Normal Tissue",
    "Solid Tissue Normal")
  x <- x[x$sample.type %in% tumor.type, ]
  x <- droplevels(x)
  x$sample.type <- factor(x$sample.type, levels = tumor.type)
  #y <- x[x$variable == "EIF4E", ]
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_16 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 16
  )
  black_bold_tahoma_16_90 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 16,
    angle  = 90,
    hjust  = 1,
    vjust  = 0.5
  )
  p1 <- ggplot(data = x,
    aes(x     = sample.type,
        y     = value,
        fill  = sample.type,
        color = sample.type)) +
    stat_n_text(size = 6, fontface = "bold",angle = 90, hjust = 0) + 
    facet_grid(. ~ variable,
               scales = "free",
               space  = "free") +
    # facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      width    = .1,
      color    = "black",
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
         y = paste("log2(RNA counts)")) +
    scale_x_discrete(labels = c("Metastatic", 
                                "Primary Tumor", 
                                "Solid Tissue \n (Normal)")
    ) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_16,
      axis.title.x    = element_blank(),
      axis.title.y    = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16_90,
      axis.text.y     = black_bold_tahoma_16_90,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_16
    ) +
    
    stat_compare_means(
      comparisons = list(
        c("Metastatic", "Solid Tissue Normal"),
        c("Primary Tumor", "Solid Tissue Normal"),
        c("Metastatic", "Primary Tumor")
    ),
      method = "t.test", label = "p.signif", size = 6)
  
  p2 <- ggplot(data = x,
    aes(x     = variable,
        y     = value,
        fill  = variable,
        color = variable)) +
    stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) + 
    facet_grid(~ sample.type,
               scales = "free",
               space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      width    = .1,
      color    = "black",
      position = position_dodge(width = .9)
    ) +
    labs(x = "EIF4F complex",
         y = paste("log2(RNA counts)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_16,
      axis.title      = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16_90,
      axis.text.y     = black_bold_tahoma_16_90,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_16
    ) +
    stat_compare_means(comparisons = list(
      c("EIF4A1", "EIF4E"),
      c("EIF4G1", "EIF4E"),
      c("EIF4EBP1", "EIF4E")
    ),
      method = "t.test", label = "p.signif", size = 6, hjust = 0)
  print(p1)
  print(p2)
}
plotEIF.RNAseq.TCGA (get.EIF.TCGA.RNAseq.long(EIF.gene))

  ##
plot.box.EIF.score.TCGA <-  function (x) {
  x <- x[x$value != 0, ]
  name <- deparse(substitute(x))
  ratio <- c("EIF4A1:EIF4E", "EIF4G1:EIF4E", "EIF4EBP1:EIF4E")
  x <- x[x$variable %in% ratio, ]
  # y <- x[x$variable == "EIF4E", ]
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  # reorder bars by explicitly ordering factor levels 
  x$primary.disease <- as.factor(x$primary.disease)
  mean <- within(x[x$variable == "EIF4EBP1:EIF4E", ], # TCGAstudy is one column in df2
    primary.disease <- reorder(primary.disease, value, median))
  mean$primary.disease <- as.factor(mean$primary.disease)
  neworder <- levels(mean$primary.disease)
  x.ordered <- factor(x$primary.disease, levels = neworder)
  
  p1 <- ggplot(data = x,
    aes(x     = x.ordered,
        y     = value, fill  = variable,
        color = variable)) +
    stat_n_text(size = 5, fontface = "bold", hjust = 0) + 
    # geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      position = position_dodge(width = .9)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "primary disease",
         y = paste("log2(ratio)")) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.title    = element_blank(),
      legend.position = "top",
      legend.text     = black_bold_tahoma_12,
      strip.text      = black_bold_tahoma_12
    )
  print(p1)
  
}

##
plotEIF.score.TCGA <-  function (x) {
  name <- deparse(substitute(x))
  ratio <- c("EIF4A1:EIF4E", "EIF4G1:EIF4E", "EIF4EBP1:EIF4E", "EIF4G1:EIF4EBP1")
  x <- x[x$variable %in% ratio, ]
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_16 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 16
  )
  black_bold_tahoma_16_90 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 16,
    angle  = 90,
    hjust  = 1,
    vjust  = 0.5
  )
  p1 <- ggplot(data = x,
    aes(x     = sample.type,
        y     = value,
        fill  = sample.type,
        color = sample.type)) +
    stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) + 
    facet_grid(~ variable,
      scales = "free",
      space  = "free") +
    facet_wrap(~ variable,
      ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      width     = .1,
      color    = "black"
    ) +
    labs(x = "sample type",
      y = paste("log2(ratio)")) +
    scale_x_discrete(labels = c("Metastatic", 
      "Primary Tumor", 
      "Solid Tissue \n (Normal)")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_16,
      axis.title.x    = element_blank(),
      axis.title.y    = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16_90,
      axis.text.y     = black_bold_tahoma_16_90,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(comparisons = list(
      c("Metastatic", "Solid Tissue Normal"),
      c("Primary Tumor", "Solid Tissue Normal"),
      c("Metastatic", "Primary Tumor")
    ),
      method = "t.test", label = "p.signif", size = 6, hjust = 0)
  
  p2 <- ggplot(data = x,
    aes(x     = variable,
      y     = value,
      fill  = variable,
      color = variable)) +
    stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) + 
    facet_grid(~ sample.type,
      scales = "free",
      space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      width     = .1,
      color    = "black"
    ) +
    labs(x = "EIF complex",
      y = paste("log2(ratio)")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_16,
      axis.title      = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16_90,
      axis.text.y     = black_bold_tahoma_16_90,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(comparisons = list(
      c("EIF4A1:EIF4E",
        "EIF4E:EIF4E"),
      c("EIF4G1:EIF4E",
        "EIF4E:EIF4E"),
      c("EIF4EBP1:EIF4E",
        "EIF4E:EIF4E")
    ),
      method = "t.test", label = "p.signif", size  = 6, hjust = 0)
  print(p1)
  print(p2)
}

##
plotEIF.score.TCGA.GTEX.tissue <-  function (tissue) {
  get.EIF.TCGA.GTEX.score.tissue <- function (tissue) {
    EIF.TCGA.GTEX <- read.csv(file.path("project-data",
      "EIFTCGAGTEX.csv"),
      header = TRUE,
      sep = ",")
    EIF.TCGA.GTEX.tissue <- EIF.TCGA.GTEX[grep(tissue,
      EIF.TCGA.GTEX$primary_disease_or_tissue), ]
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.tissue
    EIF.TCGA.GTEX.score <- tibble::add_column(EIF.TCGA.GTEX.score, 
                                              EIF4F = NA, 
                                              .after = "EIF4EBP1")
    EIF.TCGA.GTEX.score$EIF4A1 <-
      (EIF.TCGA.GTEX.tissue$EIF4A1 - EIF.TCGA.GTEX.tissue$EIF4E)
    EIF.TCGA.GTEX.score$EIF4F <-
      (EIF.TCGA.GTEX.tissue$EIF4G1 - EIF.TCGA.GTEX.tissue$EIF4EBP1)
    EIF.TCGA.GTEX.score$EIF4G1 <-
      (EIF.TCGA.GTEX.tissue$EIF4G1 - EIF.TCGA.GTEX.tissue$EIF4E)
    EIF.TCGA.GTEX.score$EIF4EBP1 <-
      (EIF.TCGA.GTEX.tissue$EIF4EBP1 - EIF.TCGA.GTEX.tissue$EIF4E)
    EIF.TCGA.GTEX.score$RPS6KB1 <-
      (EIF.TCGA.GTEX.tissue$RPS6KB1 - EIF.TCGA.GTEX.tissue$EIF4E)
    EIF.TCGA.GTEX.score$MYC <-
      (EIF.TCGA.GTEX.tissue$MYC - EIF.TCGA.GTEX.tissue$EIF4E)
    EIF.TCGA.GTEX.score$EIF4E <-
      (EIF.TCGA.GTEX.tissue$EIF4E - EIF.TCGA.GTEX.tissue$EIF4E)
    colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <-
      "EIF4A1:EIF4E ratio"
    colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <-
      "EIF4G1:EIF4E ratio"
    colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <-
      "EIF4EBP1:EIF4E ratio"
    colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <-
      "RPS6KB1:EIF4E ratio"
    colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <-
      "MYC:EIF4E ratio"
    colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <-
      "EIF4E:EIF4E ratio"
    colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4F")] <-
      "EIF4G1:EIF4EBP1 ratio"
    EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:11])
    colnames(EIF.TCGA.GTEX.score.long) <-
      c("sample",
        "study",
        "sample.type",
        "primary.disease",
        "variable",
        "value")
    tumor.type <- c(
      "Metastatic",
      "Primary Tumor",
      "Recurrent Tumor",
      "Solid Tissue Normal",
      "Normal Tissue",
      "Cell Line"
    )
    EIF.TCGA.GTEX.score.long <- EIF.TCGA.GTEX.score.long[EIF.TCGA.GTEX.score.long$sample.type %in% tumor.type, ]
    EIF.TCGA.GTEX.score.long <- droplevels(EIF.TCGA.GTEX.score.long)
    EIF.TCGA.GTEX.score.long$sample.type <-
      factor(EIF.TCGA.GTEX.score.long$sample.type,
        levels = tumor.type)
    return(EIF.TCGA.GTEX.score.long)
  }
  x <- get.EIF.TCGA.GTEX.score.tissue(tissue)
  name <- deparse(substitute(x))
  gene.number <- nlevels(x$variable)
  metastatic.number <- nrow(x[x$sample.type == "Metastatic", ])/gene.number
  primary.tumor.number <- nrow(x[x$sample.type == "Primary Tumor", ])/gene.number
  recurrent.tumor.number <-
    nrow(x[x$sample.type == "Recurrent Tumor", ])/gene.number
  solid.tissue.normal.number <-
    nrow(x[x$sample.type == "Solid Tissue Normal", ])/gene.number
  normal.tissue.number <-
    nrow(x[x$sample.type == "Normal Tissue", ])/gene.number
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  p1 <- ggplot(data = x,
    aes(x     = sample.type,
        y     = value,
        color = sample.type)) +
    facet_grid(~ variable,
      scales = "free",
      space  = "free") +
    facet_wrap(~ variable,
      ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = tissue,
         y = paste("log2(value)")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_discrete(
      labels = c(
        "Metastatic"          = paste("Metastatic \n n= ",
          metastatic.number),
        "Primary Tumor"       = paste("Primary Tumor \n n= ",
          primary.tumor.number),
        "Solid Tissue Normal" = paste("Solid Tissue Normal \n n= ",
          solid.tissue.normal.number),
        "Normal Tissue"       = paste("Normal Tissue \n n= ",
          normal.tissue.number)
      )
    ) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(comparisons = list(
      c("Metastatic", "Normal Tissue"),
      c("Primary Tumor", "Normal Tissue"),
    #  c("Solid Tissue Normal", "Normal Tissue"),
      c("Metastatic", "Primary Tumor")
    ),
      method = "t.test", label = "p.signif")
  print(p1)
  p2 <- ggplot(data = x,
    aes(x     = variable,
      y     = value,
      color = variable)) +
    facet_grid(~ sample.type,
      scales = "free",
      space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = tissue,
         y = paste("log2(value)")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(comparisons = list(
      c("EIF4A1:EIF4E ratio",
        "EIF4E:EIF4E ratio"),
      c("EIF4G1:EIF4E ratio",
        "EIF4E:EIF4E ratio"),
      c("EIF4EBP1:EIF4E ratio",
        "EIF4E:EIF4E ratio"),
      c("RPS6KB1:EIF4E ratio",
        "EIF4E:EIF4E ratio"),
      c("MYC:EIF4E ratio",
        "EIF4E:EIF4E ratio")
    ),
      method = "t.test", label = "p.signif")
  print(p2)
}

##
plotEIF.GTEX <-  function (x) {
  name <- deparse(substitute(x))
  gene.number <- nlevels(x$variable)
  cell.line.number <- 
    nrow(x[x$sample.type == "Cell Line", ])/gene.number
  normal.tissue.number <- 
    nrow(x[x$sample.type == "Normal Tissue", ])/gene.number
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  ggplot(data = x,
    aes(x     = sample.type,
      y     = value,
      color = sample.type)) +
    facet_grid(~ variable,
      scales = "free",
      space  = "free") +
    facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
      y = paste("log2(TPM)")) +
    scale_x_discrete(labels = c(
      "Cell Line"     = paste("Cell Line \n n= ",
        cell.line.number),
      "Normal Tissue" = paste("Normal Tissue \n n= ",
        normal.tissue.number)
    )) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(method = "anova")
}

##
plot.EIF.seq.all.samples <- function(x) {
  plotEIF.TCGA.GTEX(x) +
    stat_compare_means(method = "anova")
}

##
plot.EIF.seq.each.tumor <- function(x, y) {
  z = y
  print(z)
  m <- x[x$primary.disease == y, ]
  name <- deparse(substitute(m))
  y <- m[m$variable == "EIF4E", ]
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  p1 <- ggplot(data = m,
    aes(x     = sample.type,
        y     = value,
        color = sample.type)) +
    stat_n_text() + 
    facet_grid(~ variable,
               scales = "free",
               space  = "free") +
    facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
         y = paste("log2(RNA counts)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    scale_x_discrete(limit = c("Primary Tumor", "Solid Tissue Normal"),
      labels = c("Tumor","Normal")) +
    stat_compare_means(comparisons = list(
      c("Primary Tumor", 
        "Solid Tissue Normal")),
      method = "t.test", label = "p.signif")
  p2 <- ggplot(data = m,
    aes(x     = variable,
        y     = value,
        color = variable)) +
    stat_n_text() + 
    facet_grid(~ sample.type,
               scales = "free",
               space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "EIF4F complex",
         y = paste("log2(RNA counts)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = element_blank(),
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(comparisons = list(
      c("EIF4A1", "EIF4E"),
      c("EIF4G1", "EIF4E"),
      c("EIF4EBP1", "EIF4E")
    ),
      method = "t.test", label = "p.signif")
  print(p1 +  ggtitle(z))
  print(p2 +  ggtitle(z))
}

##
get.disease.list <- function (x) {
  y <- get.EIF.TCGA.RNAseq.long(x)
  y$primary.disease <- as.factor(y$primary.disease)
  disease.list <- levels(y$primary.disease)
  names(disease.list) <- disease.list
  return(disease.list)
  }

##
plot.EIF.seq.all.normal <- function(x) {
  plotEIF.GTEX(x) +
    stat_compare_means(method = "t.test")
}

#########################################################################
##  Kaplan-Meier curve with clinic and EIF RNASeq data all tumor group ##
#########################################################################
plot.km.EIF.all.tumors <- function(EIF) {
  EIF.TCGA.GTEX <- read_csv(
    "project-data/EIFTCGAGTEX.csv", 
    col_types = cols(OS      = col_number(), 
      OS.time = col_number())
  )
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA', ]
  EIF.TCGA <-
    EIF.TCGA[EIF.TCGA$sample_type != "Solid Tissue Normal",]
  EIF.TCGA <- droplevels(EIF.TCGA)
  df <- na.omit(EIF.TCGA)
  number <- nrow(df)
  sub <- round(number / 5, digits = 0)
  # bottom.label <- paste("Bottom 20%, n = ", sub)
  # top.label <- paste("Top 20%, n = ", sub)
  df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.2)] = "Bottom 20%"
  df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(
    face   = "bold",
    size   = 16,
    family = "Tahoma",
    colour = "black"
  )
  print(
    ggplot2::autoplot(
      km, 
      xlab = "Days",
      ylab = "Survival Probability",
      main = paste0("All TCGA cancer studies (",
        number,
        " cases)"),
      xlim = c(0, 4000)) +
      theme_bw() +
      theme(
        plot.title           = black.bold.12pt,
        axis.title           = black.bold.12pt,
        axis.text            = black.bold.12pt,
        axis.line.x          = element_line(color  = "black"),
        axis.line.y          = element_line(color  = "black"),
        panel.grid           = element_blank(),
        strip.text           = black.bold.12pt,
        legend.text          = black.bold.12pt ,
        legend.title         = black.bold.12pt ,
        legend.position      = c(1, 1),
        legend.justification = c(1, 1)) +
      guides(fill = FALSE) +
      scale_color_manual(
        values = c("red", "blue"),
        name   = paste(EIF, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(paste("Bottom 20%, n = ", sub), 
          paste("Top 20%, n = ", sub))
      ) +
      geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 4000,
        y        = 0.8,
        label    = paste("log-rank test \n p.val = ", p.val),
        size     = 6.5,
        hjust    = 1,
        fontface = "bold"
      )
  )
  # rho = 1 the Gehan-Wilcoxon test
  print(EIF)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

##
plot.km.EIF.each.tumor <- function(EIF, tumor) {
  EIF.TCGA.GTEX <- read_csv(
    "project-data/EIFTCGAGTEX.csv", 
    col_types = cols(OS      = col_number(), 
      OS.time = col_number())
  )
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA', ]
  EIF.TCGA <-
    EIF.TCGA[EIF.TCGA$primary_disease_or_tissue == tumor,]
  EIF.TCGA <-
    EIF.TCGA[EIF.TCGA$sample_type != "Solid Tissue Normal",]
  EIF.TCGA <- droplevels(EIF.TCGA)
  df <- na.omit(EIF.TCGA)
  number <- nrow(df)
  sub <- round(number / 5, digits = 0)
  bottom.label <- paste("Bottom 20%, n = ", sub)
  top.label <- paste("Top 20%, n = ", sub)
  df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.2)] = "Bottom 20%"
  df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <-
    survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <-
    survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(
    face   = "bold",
    size   = 16,
    family = "Tahoma",
    colour = "black"
  )
  print(
    ggplot2::autoplot(
      km,
      xlim = c(0, 4000),
      xlab = "Days",
      ylab = "Survival Probability",
      main = paste0(tumor, " (",
        number, " cases)")
    ) +
      theme_bw() +
      theme(
        plot.title           = black.bold.12pt,
        axis.title           = black.bold.12pt,
        axis.text            = black.bold.12pt,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        strip.text           = black.bold.12pt,
        legend.text          = black.bold.12pt ,
        legend.title         = black.bold.12pt ,
        legend.position      = c(1, 1),
        legend.justification = c(1, 1)
      ) +
      guides(fill = FALSE) +
      scale_color_manual(
        values  = c("red", "blue"),
        name    = paste(EIF, "mRNA expression"),
        breaks  = c("Bottom 20%", "Top 20%"),
        labels  = c(bottom.label, top.label)
      ) +
      #scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+
      geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 4000,
        y        = 0.8,
        label    = paste("log-rank test \n p.val = ", p.val),
        size     = 6.5,
        hjust    = 1,
        fontface = "bold"
      )
  )
  # rho = 1 the Gehan-Wilcoxon test
  print(EIF)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

#################################################################
##  PCA plots on EIF4F RNA-seq data from TCGA and GTEx groups  ##
#################################################################
plot.EIF.TCGA.GTEX.PCA.all <- function () {
  EIF.TCGA.GTEX <- read.csv(file = "project-data/EIFTCGAGTEX2.csv",
                            header = TRUE,
                            sep = "\t")
  EIF.TCGA.GTEX <- as.data.frame(EIF.TCGA.GTEX)
  EIF.TCGA.GTEX <- EIF.TCGA.GTEX[ ,3:10]
  #  EIF.proteomics$gene_id <- as.factor(EIF.proteomics$gene_id)
  EIF.TCGA.GTEX$X_sample_type <-
    as.factor(EIF.TCGA.GTEX$X_sample_type)
  EIF.TCGA.GTEX$X_sample_type <- factor(
    EIF.TCGA.GTEX$X_sample_type,
    levels = c("Normal Tissue",
               "Primary Tumor",
               "Metastatic",
               "Solid Tissue Normal")
    )
  EIF.TCGA.GTEX <- na.omit(EIF.TCGA.GTEX)
  df1 <- EIF.TCGA.GTEX[, -8]
  # df1 <- as.numeric(df1)
  # Do PCA
  PCA <- prcomp(df1)
  # Extract PC axes for plotting
  PCAvalues <- data.frame(Sample.type = EIF.TCGA.GTEX$X_sample_type, PCA$x)
  PCAvalues$Sample.type <- factor(PCAvalues$Sample.type, 
    levels=c("Normal Tissue", "Solid Tissue Normal", "Primary Tumor", "Metastatic"), 
    labels=c("Healthy Tissue (GTEx)", "Solid Normal Tissue (TCGA)", 
      "Primary Tumor (TCGA)", "Metastatic Tumor (TCGA)"))
  # Extract loadings of the variables
  PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
  # Plot
  p <- ggplot(PCAvalues, 
         aes(x      = PC1, 
             y      = PC2, 
             colour = Sample.type)) +
    geom_point(size = 2) +
    geom_segment(data     = PCAloadings, 
                 aes(x    = 0, 
                     y    = 0, 
                     xend = (PC1*5),
                     yend = (PC2*5)), 
                 arrow    = arrow(length = unit(1/3, "picas")),
                 color    = "black") +
    annotate("text",
             size     = 6,
             fontface = "bold",
             x        = (PCAloadings$PC1*5), 
             y        = (PCAloadings$PC2*5),
             label    = PCAloadings$Variables) +
    stat_n_text(geom = "label") +
    ggtitle ("All TCGA and GTEX cases") +
    theme(
      plot.background  = element_blank(),
      plot.title       = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold"),
      panel.background = element_rect(fill    = 'transparent',
                                      color   = 'black',
                                      size    = 1),
      axis.title       = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold"),
      axis.text        = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold"),
      legend.title     = element_blank(),
      legend.position  = c(0.75, 0.9),
      legend.background  = element_blank(),
      legend.text      = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold",
                                      hjust   = 0),
      legend.key       = element_blank()) +
    scale_x_continuous(limits = c(-8, 8)) +
    scale_y_continuous(limits = c(-6, 5))
  print(p)
    }
plot.EIF.TCGA.GTEX.PCA.all()
  
##
plot.EIF.TCGA.GTEX.PCA.each.tissue <- function (x) {
  EIF.TCGA.GTEX <- read.csv(file = "project-data/EIFTCGAGTEX2.csv",
    header = TRUE,
    sep = "\t")
  EIF.TCGA.GTEX <- as.data.frame(EIF.TCGA.GTEX)
  EIF.TCGA.GTEX <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$X_primary_site == x, ]
  # EIF.TCGA.GTEX <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$TCGA_GTEX_main_category != "TCGA Lung Squamous Cell Carcinoma", ]
  EIF.TCGA.GTEX <- EIF.TCGA.GTEX[ ,3:12]
  #  EIF.proteomics$gene_id <- as.factor(EIF.proteomics$gene_id)
  EIF.TCGA.GTEX$X_sample_type <- as.factor(EIF.TCGA.GTEX$X_sample_type)
  EIF.TCGA.GTEX$X_sample_type <- factor(EIF.TCGA.GTEX$X_sample_type,
                                        levels = c("Normal Tissue",
                                                   "Primary Tumor",
                                                   "Metastatic",
                                                   "Solid Tissue Normal"))
  EIF.TCGA.GTEX <- na.omit(EIF.TCGA.GTEX)
  df1 <- EIF.TCGA.GTEX[, -c(8:11)]
  # df1 <- as.numeric(df1)
  # Do PCA
  PCA <- prcomp(df1)
  # Extract PC axes for plotting
  PCAvalues <- data.frame(Sample.type = EIF.TCGA.GTEX$X_sample_type, PCA$x)
  PCAvalues$Sample.type <- factor(PCAvalues$Sample.type, 
    levels=c("Normal Tissue", "Solid Tissue Normal", "Primary Tumor", "Metastatic"), 
    labels=c("Healthy Tissue (GTEx)", "Solid Normal Tissue (TCGA)", 
      "Primary Tumor (TCGA)", "Metastatic Tumor (TCGA)"))
  # Extract loadings of the variables
  PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
  # Plot
  p <- ggplot(PCAvalues, aes(x      = PC1, 
                        y      = PC2, 
                        colour = Sample.type)) +
    geom_segment(data     = PCAloadings, 
                 aes(x    = 0, 
                     y    = 0, 
                     xend = (PC1*3),
                     yend = (PC2*3)), 
                 arrow    = arrow(length = unit(1/3, "picas")),
                 color    = "black") +
    geom_point(size = 2) +
    annotate("text",
              size     = 6,
              fontface = "bold",
              x        = (PCAloadings$PC1*3), 
              y        = (PCAloadings$PC2*3),
              label    = PCAloadings$Variables) +
    ggtitle (x) +
    theme(
      plot.background  = element_blank(),
      plot.title       = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold"),
      panel.background = element_rect(fill    = 'transparent',
                                      color   = 'black',
                                      size    = 1),
      axis.title       = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold"),
      axis.text        = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold"),
      legend.title     = element_blank(),
      legend.position  = c(0.7, 0.9),
      legend.background  = element_blank(),
      legend.text      = element_text(colour  = "black",
                                      size    = 18,
                                      face    = "bold",
                                      hjust   = 0),
      legend.key       = element_blank()) +
    scale_x_continuous(limits = c(-4, 3)) +
    scale_y_continuous(limits = c(-3, 3))
  print(p)
 }

plot.EIF.TCGA.GTEX.PCA.each.tissue ("Lung")
tissue.list <- function() {
  EIF.TCGA.GTEX <- read.csv(file = "project-data/EIFTCGAGTEX2.csv",
    header = TRUE,
    sep = "\t")
  EIF.TCGA.GTEX$X_sample_type <-
    as.factor(EIF.TCGA.GTEX$X_sample_type)
  EIF.TCGA.GTEX$X_sample_type <- factor(
    EIF.TCGA.GTEX$X_sample_type,
    levels = c(
      "Normal Tissue",
      "Primary Tumor",
      "Metastatic",
      "Solid Tissue Normal"
    ))
  EIF.TCGA.GTEX <- na.omit(EIF.TCGA.GTEX)
  EIF.TCGA.GTEX$X_primary_site <- droplevels(EIF.TCGA.GTEX$X_primary_site)
  tissue.list <- levels(EIF.TCGA.GTEX$X_primary_site)
  tissue.list <- tissue.list[ -1]
  return(tissue.list)
}
lapply(tissue.list(), plot.EIF.TCGA.GTEX.PCA.each.tissue)

##
plot.EIF.PCA <- function () {
  EIF.proteomics <- read.csv(file.path("project-data",
    "proteomics.csv"),
    header = TRUE,
    sep = ",")
  EIF.proteomics <- as.data.frame(EIF.proteomics)
  #  EIF.proteomics$gene_id <- as.factor(EIF.proteomics$gene_id)
  EIF.proteomics <- na.omit(EIF.proteomics)
  df1 <- EIF.proteomics[, -1]
  df <- as.numeric(df)
  ggplot2::autoplot(prcomp(df1))
  ggplot2::autoplot(
    prcomp(df),
    data                 = DNFA.TCGA.GTEX.skin,
    colour               = 'sample_type',
    loadings             = TRUE,
    loadings.colour      = 'black',
    loadings.label.vjust = 0,
    loadings.label       = TRUE,
    loadings.label.size  = 6
  ) +
    theme(
      plot.background  = element_blank(),
      panel.background =
        element_rect(
          fill  = 'transparent',
          color = 'black',
          size  = 1
        ),
      axis.title   = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold"
      ),
      axis.text    = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold"
      ),
      legend.title = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold"
      ),
      legend.text  = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold",
        hjust   = 1
      ),
      legend.key   = element_blank()
    )
}

####################################################
####################################################
EIF <- "EIF4A1"
EIF.gene <- c(
              "EIF4E",
              "EIF4G1",
              "EIF4A1",
              "EIF4EBP1"
              )
names(EIF.gene) <- EIF.gene

plotEIF.RNAseq.TCGA.GTEX (get.EIF.TCGA.GTEX.RNAseq.long())
  plotEIF.RNAseq.TCGA.GTEX.tissue ("Colon")

plot.box.EIF.RNAseq.TCGA (get.EIF.TCGA.RNAseq.long(EIF.gene))
plotEIF.RNAseq.TCGA (get.EIF.TCGA.RNAseq.long(EIF.gene))

plot.box.EIF.score.TCGA (get.EIF.TCGA.score.long(EIF.gene))
plotEIF.score.TCGA (get.EIF.TCGA.score.long(EIF.gene))

# plotEIF.score.TCGA (get.EIF.TCGA.GTEX.score.tissue("Lung"))
plotEIF.score.TCGA.GTEX.tissue("Breast")

plot.EIF.seq.all.normal (get.EIF.GTEX.RNAseq.long())
plot.EIF.seq.all.normal (get.EIF.GTEX.score.long())

plot.EIF.seq.each.tumor (
  x = get.EIF.TCGA.RNAseq.long(EIF.gene),
  y = "Skin Cutaneous Melanoma")

lapply(get.disease.list(EIF.gene),
       plot.EIF.seq.each.tumor,
       x = get.EIF.TCGA.RNAseq.long(EIF.gene))

####################################################
plot.km.EIF.all.tumors("EIF4E")
lapply(EIF.gene, plot.km.EIF.all.tumors)

plot.km.EIF.each.tumor("MYC", "Breast Invasive Carcinoma")
lapply(get.disease.list("EIF4E"),
  plot.km.EIF.each.tumor,
  EIF = "MYC")

lapply(EIF.gene,
  plot.km.EIF.each.tumor,
  tumor = "Lung Adenocarcinoma")

####################################################
plot.EIF.TCGA.GTEX.PCA.all()

plot.EIF.TCGA.GTEX.PCA.each.tissue ("Lung")    
lapply(tissue.list(), plot.EIF.TCGA.GTEX.PCA.each.tissue)

