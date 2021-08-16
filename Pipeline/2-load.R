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
#library(glmnet)
library(gplots)
library(gridExtra)
library(igraph)
library(KEGG.db)
library(lemon) ## coord_capped_cart(bottom='both', left='both')
library(limma)
library(missMDA)
library(nortest) # test for normal distribution
library(org.Hs.eg.db)
library(pca3d)
#library(pheatmap)
library(RColorBrewer)
library(ReactomePA)
library(readr)
library(readxl)
library(reshape2)
library(rgl)
library(scales) # Log scaling of the y axis
library(survival)
library(survivalAnalysis)
#library(survMisc)
#library(survminer)
library(tidyverse)
library(vcd)
library(vip)

data.file.directory <- "~/Downloads/Test"
output.directory <- "~/Documents/EIF_output"

# p <-NCmisc::list.functions.in.file("EIFanalysisv2.R", alphabetic = TRUE)

# TODO: Make all uses of 'stringsAsFactors' explicit.
#
# This code was originally developed with R 3.6.3, then later
# evaluated against R 4.0.x.  From the R 4.0.0 release notes:
#
# https://cran.r-project.org/doc/manuals/r-devel/NEWS.html
# -----
# R now uses a stringsAsFactors = FALSE default, and hence by
# default no longer converts strings to factors in calls to
# data.frame() and read.table().  A large number of packages
# relied on the previous behaviour and so have needed/will
# need updating.
# -----
#
# Therefore, for the time being, here we set a global option
# to restore the old default behavior.  Note that this option
# is deprecated and will no longer work in R 4.1.
options(stringsAsFactors = TRUE)


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
  return(element_text(
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