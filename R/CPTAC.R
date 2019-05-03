library(readr)
library(readxl)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)

BRCA_Proteome_sample <- read_delim(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Phosphoproteome.sample.csv", 
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE
  )

BRCA_Proteome_summary <- read_delim(
  "~/Documents/translation/CPTAC/TCGA_Breast_BI_Proteome.summary.csv", 
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE
)

BRCA_Proteome_itraq <- fread(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Proteome.itraq.tsv",
  header = T
  )



BRCA_Phosphoproteome_summary <- read_delim(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Phosphoproteome.summary.csv", 
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE
  )

BRCA_Phosphopeptide_itraq <- fread(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Phosphoproteome.phosphopeptide.itraq-1.tsv",
  header = T
  )

BRCA_Phosphosite_itraq <- fread(
  "Documents/translation/CPTAC/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq-1.tsv", 
  header = T
  )





CPTAC_BC_somatic_mutations <- read_excel(
  "Documents/translation/Proteogenomics connects somatic mutations to signalling in breast cancer/nature18003-s2/CPTAC_BC_SupplementaryTable01.xlsx")


##########################################################
## use data from CPTAC Breast Cancer Confirmatory Study ##
##########################################################
## CPTAC_BCprospective_Proteome used different iTRAQ labeling scheme 
## and give different data from the TCGA dataset
CPTAC2_Breast_Prospective_Collection_BI_Proteome <- read_delim(
  "Documents/translation/CPTAC/CPTAC2_Breast_Prospective_Collection_BI_Proteome.summary.csv", 
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE)
Proteome_EIF4 <- CPTAC2_Breast_Prospective_Collection_BI_Proteome [
  grep("EIF4", CPTAC2_Breast_Prospective_Collection_BI_Proteome$Gene), ]
Proteome_EIF4_2 <- Proteome_EIF4 [ ,
  grep("Spectral Counts", names(Proteome_EIF4), value = TRUE)]
rownames(Proteome_EIF4_2) <- Proteome_EIF4$Gene
Proteome_EIF4_3 <- as.data.frame(t(Proteome_EIF4_2))
Proteome_EIF4_4 <- melt(as.matrix(Proteome_EIF4_3[-17, ]))


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
    y     = log2(value),
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
    y = paste("Spectral Counts")) +
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

plot.CPTAC ("Peptide", BRCA_Proteome_summary)
plot.CPTAC ("Phosphoeptide", BRCA_Phosphoproteome_summary)


################################################################
## extract protein quantity from iTRAQ and MS/MS spectra data ##
################################################################
data <- BRCA_Proteome_summary
Proteome_EIF4 <- data[grep("EIF4", data$Gene), ]
Proteome_EIF4_2 <- Proteome_EIF4[ ,
  grep("Spectral Counts", names(Proteome_EIF4), value = TRUE)]
rownames(Proteome_EIF4_2) <- Proteome_EIF4$Gene
colnames(Proteome_EIF4_2) <- str_remove(colnames(Proteome_EIF4_2), 
  "Spectral Counts")
Proteome_EIF4_2 [[38]] <- NULL


data <- BRCA_Proteome_itraq
Proteome_itraq_EIF4 <- data[
  grep("EIF4", data$Gene), ]
Proteome_itraq_EIF4 <- as.data.frame(Proteome_itraq_EIF4)
Proteome_itraq_EIF4_1 <- Proteome_itraq_EIF4[ ,
  grepl("Log Ratio", colnames(Proteome_itraq_EIF4))]
Proteome_itraq_EIF4_2 <- Proteome_itraq_EIF4_1[ ,
  grepl("Unshared Log Ratio", colnames(Proteome_itraq_EIF4_1))]
rownames(Proteome_itraq_EIF4_2) <- Proteome_itraq_EIF4$Gene
# convert all value to non-log transforms
Proteome_itraq_EIF4_3 <- exp(Proteome_itraq_EIF4_2)
colnames(Proteome_itraq_EIF4_3) <- str_remove(colnames(Proteome_itraq_EIF4_2), 
  " Unshared Log Ratio")
colnames(Proteome_itraq_EIF4_3) <- str_remove(colnames(Proteome_itraq_EIF4_3), "\\.1")
colnames(Proteome_itraq_EIF4_3) <- str_remove(colnames(Proteome_itraq_EIF4_3), "\\.2")
ncol(Proteome_itraq_EIF4_3)



## the following function draws the sum of all ratios
BRCA_Proteome_ratiosum <- BRCA_Proteome_sample
EIF4F_list <- rownames(Proteome_itraq_EIF4_3)
for(y in EIF4F_list){
  for(x in 1:37)
    {
    group_item_name <- function(x){
      BRCA_Proteome_sample[x, c("114", "115", "116")]} 
  ## have to use unlist to convert into vector
    v <- as.vector (unlist(group_item_name(x))) 
    x1 <- rowSums(Proteome_itraq_EIF4_3[y, v])
    BRCA_Proteome_ratiosum [x ,y] <- x1
    message("x=", x)
  }
}
BRCA_Proteome_ratiosum <- BRCA_Proteome_ratiosum[ ,EIF4F_list]
BRCA_Proteome_ratiosum_2 <- BRCA_Proteome_ratiosum + 1
rownames(BRCA_Proteome_ratiosum_2) <- colnames(Proteome_EIF4_2)
BRCA_Proteome_ratiosum_3 <- t(BRCA_Proteome_ratiosum_2)
