library(readr)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)

TCGA_Breast_BI_Proteome_itraq <- fread(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Proteome.itraq.tsv",
  header = T
)

TCGA_Breast_BI_Phosphopeptide_itraq <- fread(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Phosphoproteome.phosphopeptide.itraq-1.tsv",
  header = T
)

TCGA_Breast_BI_Phosphosite_itraq <- fread(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq-1.tsv", 
  header = T)

TCGA_Breast_BI_Proteome_summary <- read_delim(
  "~/Documents/translation/CPTAC/TCGA_Breast_BI_Proteome.summary.csv", 
  "\t", escape_double = FALSE, trim_ws = TRUE)


TCGA_Breast_BI_Phosphosite_itraq_1 <- TCGA_Breast_BI_Phosphosite_itraq[
  grep("EIF4", TCGA_Breast_BI_Phosphosite_itraq$Gene), ]
TCGA_Breast_BI_Proteome_itraq_2 <- TCGA_Breast_BI_Proteome_itraq_1[ ,-c(224:229)]
x <- TCGA_Breast_BI_Proteome_itraq_2[TCGA_Breast_BI_Proteome_itraq_2$Gene == "Median", ]
x1 <- x[, -1]
x2 <- c(t(x1))
TCGA_Breast_BI_Proteome_itraq_3 <- t(TCGA_Breast_BI_Proteome_itraq_2[-1, -1])
y <- TCGA_Breast_BI_Proteome_itraq_2$Gene
y1 <- y[c(-1)]
colnames(TCGA_Breast_BI_Proteome_itraq_3) <- y1
TCGA_Breast_BI_Proteome_itraq_4 <- t(TCGA_Breast_BI_Proteome_itraq_3) 
TCGA_Breast_BI_Proteome_itraq_5 <- TCGA_Breast_BI_Proteome_itraq_4 + rep(x2, each = nrow(TCGA_Breast_BI_Proteome_itraq_4))
TCGA_Breast_BI_Proteome_itraq_6 <- TCGA_Breast_BI_Proteome_itraq_5[ ,
  grepl("Unshared Log Ratio", colnames(TCGA_Breast_BI_Proteome_itraq_5))]
TCGA_Breast_BI_Proteome_itraq_7 <- melt(as.matrix(TCGA_Breast_BI_Proteome_itraq_6))
p1 <- ggplot(data = TCGA_Breast_BI_Proteome_itraq_7,
  aes(x     = Var1,
    y     = value,
    color = Var1))   +
geom_violin(trim = FALSE) +
  geom_boxplot(
    alpha    = .01,
    size     = .75,
    width    = .5,
    position = position_dodge(width = .9)
  )

  
plot.CPTAC.iTRAQ <- function(status, data) {
  Proteome_itraq_EIF4 <- data[
    grep("EIF4", data$Gene), ]
  Proteome_itraq_EIF4 <- as.data.frame(Proteome_itraq_EIF4)
  Proteome_itraq_EIF4_1 <- Proteome_itraq_EIF4[ ,
    grepl("Log Ratio", colnames(Proteome_itraq_EIF4))]
  Proteome_itraq_EIF4_2 <- Proteome_itraq_EIF4_1[ ,
    grepl("Unshared Log Ratio", colnames(Proteome_itraq_EIF4_1))]
  rownames(Proteome_itraq_EIF4_2) <- Proteome_itraq_EIF4$Gene
  Proteome_itraq_EIF4_3 <- as.data.frame(t(Proteome_itraq_EIF4_2))
  Proteome_itraq_EIF4_4 <- melt(as.matrix(Proteome_itraq_EIF4_3))
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
  p1 <- ggplot(data = Proteome_itraq_EIF4_4,
    aes(x     = Var2,
      y     = value,
      color = Var2)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "protein name",
      y = paste("log2 ratio", status)) +
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
  #  p1 <- p1 + stat_compare_means(method = "anova")
  print(p1)
}

plot.CPTAC.iTRAQ ("Proteome itraq", TCGA_Breast_BI_Proteome_itraq)

plot.CPTAC <- function(status, data) {
  Proteome_EIF4 <- data [grep("EIF4", data$Gene),]
  Proteome_EIF4_2 <- Proteome_EIF4 [,
    grep("Spectral Counts", names(Proteome_EIF4), value = TRUE)]
  rownames(Proteome_EIF4_2) <- Proteome_EIF4$Gene
  Proteome_EIF4_3 <- as.data.frame(t(Proteome_EIF4_2))
  Proteome_EIF4_4 <- melt(as.matrix(Proteome_EIF4_3[c(-1, -38),]))
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
  p1 <- ggplot(data = Proteome_EIF4_4,
    aes(x     = Var2,
        y     = value,
        color = Var2)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "protein name",
      y = paste(status, "Spectral Counts")) +
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
  #  p1 <- p1 + stat_compare_means(method = "anova")
  print(p1)
}

plot.CPTAC ("Peptide", TCGA_Breast_BI_Proteome_summary)
plot.CPTAC ("Phosphoeptide", TCGA_Breast_BI_Phosphoproteome_summary)
