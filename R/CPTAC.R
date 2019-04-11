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
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE
  )

TCGA_Breast_BI_Phosphoproteome_summary <- read_delim(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Phosphoproteome.summary.csv", 
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE
  )

#####################################
## import data from Nature article ##
#####################################

CPTAC_BC_proteomics <- read_excel(
  "Documents/translation/Proteogenomics connects somatic mutations to signalling in breast cancer/nature18003-s2/CPTAC_BC_SupplementaryTable03.xlsx")

CPTAC_BC_proteomics_1 <- CPTAC_BC_proteomics[, c(10, 12:122)]
CPTAC_BC_proteomics_2 <- CPTAC_BC_proteomics_1[
  grep("EIF4", CPTAC_BC_proteomics_1$geneName), ]
## remove duplicated EIF4G1 and EIF4H rows
CPTAC_BC_proteomics_2 <- CPTAC_BC_proteomics_2 [-c(1,3), ]
CPTAC_BC_proteomics_3 <- CPTAC_BC_proteomics_2 [ ,-1]
rownames(CPTAC_BC_proteomics_3) <- CPTAC_BC_proteomics_2$geneName
CPTAC_BC_proteomics_4 <- melt(as.matrix(CPTAC_BC_proteomics_3))
CPTAC_BC_proteomics_5 <- CPTAC_BC_proteomics_4
CPTAC_BC_proteomics_5$type <- "NA"
CPTAC_BC_proteomics_5$type[grep("TCGA", CPTAC_BC_proteomics_5$Var2)] <- "tumor"
CPTAC_BC_proteomics_5$type[grep("CPTAC", CPTAC_BC_proteomics_5$Var2)] <- "normal"
CPTAC_BC_proteomics_5$type <- as.factor(CPTAC_BC_proteomics_5$type)
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
p1 <- ggplot(data = CPTAC_BC_proteomics_5,
  aes(x     = Var1,
      y     = value,
      color = type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(
    alpha    = .01,
    size     = .75,
    width    = .5,
    position = position_dodge(width = .9)
  ) +
  labs(x = "protein name",
    y = paste("log2 ratio")) +
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



#################################
## Use data from CPTAC website ##
#################################

  
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
  Proteome_EIF4 <- data [grep("EIF4", data$Gene), ]
  Proteome_EIF4_2 <- Proteome_EIF4 [ ,
    grep("Spectral Counts", names(Proteome_EIF4), value = TRUE)]
  rownames(Proteome_EIF4_2) <- Proteome_EIF4$Gene
  Proteome_EIF4_3 <- as.data.frame(t(Proteome_EIF4_2))
  Proteome_EIF4_4 <- melt(as.matrix(Proteome_EIF4_3[-38, ]))
  Proteome_EIF4_4$type <- "tumor"
  Proteome_EIF4_4$type[grep("263", Proteome_EIF4_4$Var1)] <- "normal"
  Proteome_EIF4_4$type <- as.factor(Proteome_EIF4_4$type)
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
    aes(x     = type,
        y     = value,
        color = Var2)) +
    facet_grid(~ Var2,
      scales = "free",
      space  = "free") +
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
