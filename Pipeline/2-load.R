#### Library Preparation ####
library(AnnotationDbi)
library(car)
library(clusterProfiler)
library(circlize) ## for color options
library(ComplexHeatmap)
library(corrplot)
library(data.table)
library(dendextend)
library(descr)
library(dplyr)
library(EnvStats)
library(eulerr)
library(facetscales)
library(factoextra)
library(FactoMineR)
library(forcats) # change the order of x-axis
library(forestmodel)
library(forestplot)
library(Hmisc)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggthemes) ## color-blind options
library(gplots)
library(gridExtra)
library(igraph)
#library(KEGG.db)
library(lemon) ## coord_capped_cart(bottom='both', left='both')
library(limma)
library(missMDA)
library(nortest) # test for normal distribution
library(org.Hs.eg.db)
library(pca3d)
library(RColorBrewer)
library(ReactomePA)
library(readr)
library(readxl)
library(reshape2)
library(rgl)
library(scales) # Log scaling of the y axis
library(survival)
library(survivalAnalysis)
library(tidyverse)
library(vcd)
library(vip)


#### Directory Preparation ####
data.file.directory <- "~/Downloads/Test"
output.directory <- "~/Documents/EIF_output"


#### Format Preparation ####
black_bold_tahoma_7 <- function() {
  return <- (
    element_text(
      color = "black",
      face = "bold",
      size = 7
    ))
}

black_bold_12 <- function() {
  return(
    element_text(
      color = "black",
      face = "bold",
      size = 12
    )
  )
}

black_bold_12_45 <- function() {
  return(
    element_text(
      color = "black",
      face = "bold",
      size = 12,
      angle = 45,
      hjust = 1
    )
  )
}

black_bold_16 <- function() {
  return(
    element_text(
      color = "black",
      face = "bold",
      size = 16
    )
  )
}

black_bold_16_right <- function() {
  return(
    element_text(
    color = "black",
    face = "bold",
    size = 16,
    angle = 90
  ))
}

black_bold_16_45 <- function() {
  return(
    element_text(
      color = "black",
      face = "bold",
      size = 16,
      angle = 45,
      hjust = 1
    )
  )
}

black_bold_16_90 <- function() {
  return(element_text(
    color = "black",
    face = "bold",
    size = 16,
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ))
}

black_bold_18 <- function() {
  return(
    element_text(
      color = "black",
      face = "bold",
      size = 18
    )
  )
}

color <- function() {
  n <- 32  # TODO: this variable is not used
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
  # TODO: col_vector, as defined in this scope, is not used
  col_vector <- unlist(mapply(
    brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  ))
}
col_vector <- color()
