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
# library(drawProteins)
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
library(glmnet)
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
library(pheatmap)
library(plotmo) # for plot_glmnet
library(RColorBrewer)
library(ReactomePA)
library(readr)
library(readxl)
library(reshape2)
library(rgl)
library(scales) # Log scaling of the y axis
library(survival)
library(survivalAnalysis)
library(survMisc)
library(survminer)
library(tidyverse)
library(vcd)
library(vip)

data.file.directory <- "~/Downloads"
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
##############################################


## Figure 1 ##
################################################################
## stacked bar plots for eIF4F CNV status across tumor groups ##
################################################################
plot.bargraph.EIF.CNV.sum <- function(EIF) {
  pan.TCGA.CNV <- function(EIF) {
    # download https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
    TCGA.CNV <- fread(
      file.path(data.file.directory, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"),
      data.table = FALSE
    )
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    TCGA.CNV <- as.data.frame(TCGA.CNV)
    TCGA.CNV1 <- TCGA.CNV[
      !duplicated(TCGA.CNV$Sample),
      !duplicated(colnames(TCGA.CNV))
    ]
    row.names(TCGA.CNV1) <- TCGA.CNV1$Sample
    TCGA.CNV1$Sample <- NULL
    TCGA.CNV_transpose <- data.table::transpose(TCGA.CNV1)
    rownames(TCGA.CNV_transpose) <- colnames(TCGA.CNV1)
    colnames(TCGA.CNV_transpose) <- rownames(TCGA.CNV1)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    TCGA.CNV.sampletype <- merge(TCGA.CNV_transpose,
      TCGA.sampletype,
      by    = "row.names",
      all.x = TRUE
    )
    TCGA.CNV.anno <- as.data.frame(TCGA.CNV.sampletype)
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    sample.type.list <- levels(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno$primary.disease <- as.factor(TCGA.CNV.anno$primary.disease)
    cancer.type.list <- levels(TCGA.CNV.anno$primary.disease)
    return(TCGA.CNV.anno)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV(EIF)

  TCGA.CNV.anno.EIF <- function(EIF) {
    TCGA.CNV.anno.subset <- TCGA.CNV.anno[
      !TCGA.CNV.anno$sample.type %in% "Solid Tissue Normal",
    ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[
      ,
      colnames(TCGA.CNV.anno.subset) %in% c(
        EIF,
        "sample.type",
        "primary.disease"
      )
    ]
    return(EIF.TCGA.CNV.anno.subset)
  }
  EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.EIF(EIF)

  make.CNV.sum.plot <- function(EIF) {
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.CNV.anno.subset.long$primary.disease
    )
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c(
      "sample.type", "primary.disease",
      "variable", "CNV"
    )
    CNV.sum <- table(EIF.TCGA.CNV.anno.subset.long[, c("CNV", "variable")])
    CNV.sum <- as.data.frame(CNV.sum)
    # CNV.sum$TCGAstudy <- str_remove(CNV.sum$TCGAstudy, regex('_.*\n*.*'))
    CNV.sum$CNV <- factor(CNV.sum$CNV, levels = c("-2", "-1", "0", "1", "2"))
    CNV.sum$variable <- factor(CNV.sum$variable,
      levels = c("TP53","EIF3D", "EIF4E3","EIF4E2","EIF4E","EIF4A1","EIF4G3","EIF4G2","MYC","EIF4EBP1","EIF4H","EIF4B","EIF4A2","EIF4G1")
    )
    # reorder bars by explicitly ordering factor levels
    p1 <- ggplot(CNV.sum, aes(
      fill = CNV,
      y = Freq,
      x = variable
    )) +
      geom_bar(stat = "identity", position = "fill") +
      geom_col() +
      geom_text(aes(label = paste0(Freq / 100, "%")),
        position = position_stack(vjust = 0.5), size = 4
      ) +
      # scale_y_continuous(labels = scales::percent_format())+
      labs(
        x = "Tumor types (TCGA pan cancer atlas 2018)",
        y = "All TCGA tumors combined"
      ) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title.x = black_bold_16(),
        axis.title.y = element_blank(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = black_bold_16(),
        legend.position = "top",
        legend.justification = "left",
        legend.box = "horizontal",
        strip.text = black_bold_16()
      ) +
      guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
      scale_fill_manual(
        name = "Copy number variation",
        breaks = c("-2", "-1", "0", "1", "2"),
        labels = c("Deep del\n 0", "Shallow del\n 1", "Diploid\n 2", "Gain\n 3", "Amp\n 3+"),
        values = c("darkblue", "blue", "lightgreen", "red", "darkred")
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "CNV"),
      filename = "EIFCNVsum.pdf",
      plot = p1,
      width = 9,
      height = 9,
      useDingbats = FALSE
    )
  }
  make.CNV.sum.plot(EIF)
}

plot.bargraph.EIF.CNV.TCGA <- function(EIF) {
  pan.TCGA.CNV <- function() {
    # download https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
    TCGA.pancancer <- fread(
      file.path(data.file.directory, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"),
      data.table = FALSE
    )
    TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.pancancer1 <- TCGA.pancancer[
      !duplicated(TCGA.pancancer$Sample),
      !duplicated(colnames(TCGA.pancancer))
    ]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$Sample
    TCGA.pancancer1$Sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    TCGA.sampletype,
                                    by    = "row.names",
                                    all.x = TRUE
    )
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    sample.type.list <- levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    cancer.type.list <- levels(TCGA.RNAseq.anno$primary.disease)
    return(TCGA.RNAseq.sampletype)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV()
  
  pancancer.TCGA.EIF <- function() {
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno.subset <- TCGA.CNV.anno # [
    #  TCGA.CNV.anno$sample.type %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Recurrent Tumor"), ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[
      ,
      colnames(TCGA.CNV.anno.subset) %in% c(
        EIF,
        "sample.type",
        "primary.disease"
      )
    ]
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.CNV.anno.subset.long$primary.disease
    )
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c(
      "sample.type",
      "primary.disease",
      "variable",
      "CNV"
    )
    
    CNV.sum <- table(EIF.TCGA.CNV.anno.subset.long[, c("CNV", "primary.disease")])
    CNV.sum <- as.data.frame(CNV.sum)
    # CNV.sum$TCGAstudy <- str_remove(CNV.sum$TCGAstudy, regex('_.*\n*.*'))
    CNV.sum$primary.disease <- ordered(CNV.sum$primary.disease, levels = rev(levels(factor(CNV.sum$primary.disease))))
    CNV.sum$CNV <- factor(CNV.sum$CNV, levels = c("-2", "-1", "0", "1", "2"))
    return(CNV.sum)
  }
  CNV.sum <- pancancer.TCGA.EIF()
  
  levels(CNV.sum$CNV)
  # reorder bars by explicitly ordering factor levels
  make.plot <- function(EIF) {
    p1 <- ggplot(
      CNV.sum,
      aes(
        fill = CNV, order = as.numeric(CNV),
        y = Freq,
        x = primary.disease
      )
    ) +
      geom_bar(stat = "identity", position = "fill") +
      labs(
        x = "Tumor types (TCGA pan cancer atlas 2018)",
        y = paste0("Percentages of ", EIF, " CNVs")
      ) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = element_blank(),
        axis.text.x = black_bold_12(),
        axis.text.y = black_bold_12(),
        #axis.line.x = element_line(color = "black"),
        #axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = black_bold_12(),
        legend.position = "top",
        legend.justification = "left",
        legend.box = "horizontal",
        strip.text = black_bold_12()
      ) +
      scale_y_continuous(labels = scales::percent_format()) +
      guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
      scale_fill_manual(
        name = "Copy number variation",
        breaks = c("-2", "-1", "0", "1", "2"),
        labels = c(
          "Deep del\n 0", "Shallow del\n 1",
          "Diploid\n 2", "Gain\n 3", "Amp\n 3+"
        ),
        values = c(
          "darkblue", "blue",
          "lightgreen", "red",
          "darkred"
        )
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "CNV"),
      filename = paste0(EIF, "pancancerCNV.pdf"),
      plot = p1,
      width = 7.5,
      height = 9,
      useDingbats = FALSE
    )
  }
  make.plot(EIF)
}

plot.matrix.EIF.CNV.corr <- function(EIF) {
  pan.TCGA.CNV <- function(EIF) {
    # https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz
    TCGA.CNV <- fread(
      file.path(data.file.directory, "Gistic2_CopyNumber_Gistic2_all_data_by_genes"),
      data.table = FALSE
    )
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    TCGA.CNV <- as.data.frame(TCGA.CNV)
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    TCGA.CNV1 <- TCGA.CNV[
      !duplicated(TCGA.CNV$Sample),
      !duplicated(colnames(TCGA.CNV))
    ]
    row.names(TCGA.CNV1) <- TCGA.CNV1$Sample
    TCGA.CNV1$Sample <- NULL
    TCGA.CNV_transpose <- data.table::transpose(TCGA.CNV1)
    rownames(TCGA.CNV_transpose) <- colnames(TCGA.CNV1)
    colnames(TCGA.CNV_transpose) <- rownames(TCGA.CNV1)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    TCGA.CNV.sampletype <- merge(TCGA.CNV_transpose,
      TCGA.sampletype,
      by    = "row.names",
      all.x = TRUE
    )
    TCGA.CNV.anno <- as.data.frame(TCGA.CNV.sampletype)
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    sample.type.list <- levels(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno$primary.disease <- as.factor(TCGA.CNV.anno$primary.disease)
    cancer.type.list <- levels(TCGA.CNV.anno$primary.disease)
    return(TCGA.CNV.anno)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV(EIF)

  TCGA.CNV.anno.EIF <- function(EIF) {
    TCGA.CNV.anno.subset <- TCGA.CNV.anno[
      !TCGA.CNV.anno$sample.type %in% "Solid Tissue Normal",
    ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[
      ,
      colnames(TCGA.CNV.anno.subset) %in% c(
        EIF,
        "sample.type",
        "primary.disease"
      )
    ]
    return(EIF.TCGA.CNV.anno.subset)
  }
  EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.EIF(EIF)

  plot.EIF.CNV.cor <- function() {
    df1 <- EIF.TCGA.CNV.anno.subset[1:(length(EIF.TCGA.CNV.anno.subset) - 2)]
    # correlation plot
    # res <- cor(df1,  method = "pearson")
    cor_5 <- rcorr(as.matrix(df1), type = "pearson")
    M <- cor_5$r
    p_mat <- cor_5$P
    pdf(file.path(output.directory, "CNV", "EIFCNVcormatrix.pdf"),
    width = 9,
    height = 9,
    useDingbats = FALSE
    )
    corrplot(
      M,
      method      = "color",  
      cl.pos      = "n", # remove color legend
      tl.cex      = 1,
      number.cex  = 1,
      addgrid.col = "gray",
      addCoef.col = "black",
      tl.col      = "black",
      type        = "lower",
      order       = "FPC", tl.srt = 0,
      p.mat       = p_mat,
      sig.level   = 0.05, # insig = "blank"
    )
    dev.off()
  }
  plot.EIF.CNV.cor()
}

plot.boxgraph.EIF.CNVratio.TCGA <- function(EIF) {
  pan.TCGA.CNV <- function() {
    # download https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz
    TCGA.pancancer <- fread(
      file.path(data.file.directory, "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena"),
      data.table = FALSE
    )
    TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.pancancer1 <- TCGA.pancancer[
      !duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))
    ]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    TCGA.sampletype,
                                    by    = "row.names",
                                    all.x = TRUE
    )
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    sample.type.list <- levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    cancer.type.list <- levels(TCGA.RNAseq.anno$primary.disease)
    return(TCGA.RNAseq.sampletype)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV()
  
  pancancer.TCGA.EIF <- function() {
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno.subset <- TCGA.CNV.anno # [
    #  TCGA.CNV.anno$sample.type %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Recurrent Tumor"), ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[
      ,
      colnames(TCGA.CNV.anno.subset) %in% c(
        EIF,
        "sample.type",
        "primary.disease"
      )
    ]
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.CNV.anno.subset.long$primary.disease
    )
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c(
      "sample.type",
      "primary.disease",
      "variable",
      "CNV"
    )
    return(EIF.TCGA.CNV.anno.subset.long)
  }
  EIF.TCGA.CNV.anno.subset.long <- pancancer.TCGA.EIF()
  class(EIF.TCGA.CNV.anno.subset.long$CNV)
  # reorder bars by explicitly ordering factor levels
  
  make.plot <- function(EIF) {
    sts <- boxplot.stats(EIF.TCGA.CNV.anno.subset.long$CNV)$stats
    
    f1 <- factor(EIF.TCGA.CNV.anno.subset.long$primary.disease)
    f.ordered1 <- fct_rev(f1)
    p1 <- ggplot(
      data = EIF.TCGA.CNV.anno.subset.long,
      aes(
        y = 2**CNV,
        x = f.ordered1,
        color = primary.disease
      )
    ) +
      ylim(0, 3) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      stat_n_text(
        size = 5,
        fontface = "bold",
        hjust = 0
      ) +
      geom_boxplot(
        alpha = .01, 
        outlier.colour = NA,
        # size     = .75,
        # width    = 1,
        position = position_dodge(width = .9)
      ) +
      labs(
        x = "primary disease",
        y = paste(EIF, "CNV ratio", "(tumor/normal)")
      ) +
      # scale_color_manual(values = col_vector) +
      coord_cartesian(ylim = c(sts[2] / 2, max(sts) * 1.05)) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = element_blank(),
        axis.text.x = black_bold_12(),
        axis.text.y = black_bold_12(),
        #axis.line.x = element_line(color = "black"),
        #axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = black_bold_12(),
        legend.position = "none",
        legend.justification = "left",
        legend.box = "horizontal",
        strip.text = black_bold_12()
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "CNV"),
      filename = paste0(EIF, "pancancerCNVratio.pdf"),
      plot = p1,
      width = 7,
      height = 9,
      useDingbats = FALSE
    )
  }
  make.plot(EIF)
}

##################################################################
## violin plot for EIF expression in tumors vs adjacent normals ##
##################################################################
plot.boxgraph.EIF.RNAseq.TCGA <- function(EIF.gene) {
  pan.TCGA.gene <- function(EIF.gene) {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.pancancer <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    )
    # use sample column for the rowname
    TCGA.pancancer1 <- TCGA.pancancer[
      !duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))
      ]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.sampletype <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    subset <- TCGA.sampletype[, c(
      "_sample_type",
      "primary disease or tissue",
      "_primary_site",
      "_study"
    ),
    drop = FALSE
    ]
    row.names(subset) <- row.names(TCGA.sampletype)
    colnames(subset) <- c(
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    )
    
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    subset,
                                    by    = "row.names",
                                    all.x = TRUE
    )
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    
    
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[, colnames(TCGA.RNAseq.anno) %in% c(
      EIF.gene,
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    )]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(
      EIF.gene,
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    )]
    # EIF.TCGA.RNAseq.anno.subset$delta <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4G1 - 2**EIF.TCGA.RNAseq.anno.subset$EIF4E - 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!is.na(EIF.TCGA.RNAseq.anno.subset$primary.site), ]
    
    EIF.TCGA.RNAseq.anno.subset.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long <- EIF.TCGA.RNAseq.anno.subset.long[
      !EIF.TCGA.RNAseq.anno.subset.long$value == 0,
      ]
    EIF.TCGA.RNAseq.anno.subset.long$sample.type <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$sample.type
    )
    EIF.TCGA.RNAseq.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.disease
    )
    EIF.TCGA.RNAseq.anno.subset.long$primary.site <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.site
    )
    
    return(EIF.TCGA.RNAseq.anno.subset.long)
  }
  TCGA.RNAseq.anno <- pan.TCGA.gene(EIF.gene)
  
  pancancer.TCGA.EIF <- function() {
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "TCGA", ]
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno.subset [
      !TCGA.RNAseq.anno$sample.type %in% c("Cell Line"),
      ]
    TCGA.RNAseq.anno.subset$sample.type <- as.character(TCGA.RNAseq.anno.subset$sample.type)
    TCGA.RNAseq.anno.subset$sample.type[TCGA.RNAseq.anno.subset$sample.type != "Solid Tissue Normal"] <- "Tumor"
    TCGA.RNAseq.anno.subset$sample.type[TCGA.RNAseq.anno.subset$sample.type == "Solid Tissue Normal"] <- "Normal"
    TCGA.RNAseq.anno.subset$sample.type <- as.factor(TCGA.RNAseq.anno.subset$sample.type)
    return(TCGA.RNAseq.anno.subset)
  }
  pancancer.TCGA.EIF.long <- pancancer.TCGA.EIF()
  
  # reorder bars by explicitly ordering factor levels

    make.plot1 <- function(EIF) {
      pancancer.TCGA.EIF.long1 <- pancancer.TCGA.EIF.long[
        pancancer.TCGA.EIF.long$variable %in% EIF,
        ]
      levels(pancancer.TCGA.EIF.long1$variable)
      
      f1 <- factor(pancancer.TCGA.EIF.long1$primary.disease)
      f.ordered1 <- fct_rev(f1)
      p1 <- ggplot(
        data = pancancer.TCGA.EIF.long1,
        aes(
          x = f.ordered1,
          # x     = x.ordered, # order primary disease
          y = 2**value, 
          color = sample.type, 
          #fill = factor(sample.type)
        )
      ) +
        scale_y_continuous(
          trans = log2_trans(),
          #limits = c(2**11, 2**17),# for 4g
          limits = c(2**7, 2**14),# for eif4E
          labels = label_comma()
        ) +
        #stat_n_text(
        #  size = 5,
        #  fontface = "bold",
        #  hjust = 0
        #) +
        geom_boxplot(
          #alpha = .1,
          #fill = sample.type,
          outlier.shape = NA,
          #size = .75,
          #width = 1,
          position = "dodge"
        ) +
        #scale_color_manual(
        #  values = c(
        #    "Tumor" = "#CC79A7",
        #    "Normal" = "#0072B2")
        #  breaks = c("EIF4EBP1", "EIF4E", "EIF4G1", "EIF4A1")
        #) +
        scale_color_manual(values = c("Tumor" = "#CC79A7", "Normal" = "#0072B2"),
                           breaks = c("Tumor", "Normal"),
                           labels = c("Tumor\n", "Normal\n")
                           ) +
        labs(
          x = "primary disease",
          y = paste(EIF, "expression (RNA-Seq counts)")
        ) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title = black_bold_12(),
          axis.title.x = black_bold_12(),
          axis.title.y = element_blank(),
          axis.text.x = black_bold_12(),
          axis.text.y = black_bold_12(),
          # axis.line.x          = element_line(color = "black"),
          # axis.line.y          = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = black_bold_12(),
          legend.position = "top",
          legend.justification = "left",
          legend.box = "horizontal",
          strip.text = black_bold_12()
        )
        #geom_signif(comparisons=list(c("Tumor", "Normal")))
      print(p1)
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = paste0(EIF, "tumorvsnormal.pdf"),
        plot = p1,
        width = 7.5,
        height = 9,
        useDingbats = FALSE
      )
      }
    #make.plot1("EIF4G1")
  lapply(EIF.gene, make.plot1)
}



# Figure 1
plot.bargraph.EIF.CNV.sum(c("TP53", "EIF4A1","EIF4A2", "EIF4E", "EIF4E2", "EIF4E3", 
                            "MYC", "EIF3D", "EIF4EBP1", "EIF4G1","EIF4G2", "EIF4G3",
                            "EIF4H", "EIF4B"))

lapply(c("TP53", "EIF4A1","EIF4A2", "EIF4E", "EIF4E2", "EIF4E3", 
         "MYC", "EIF3D", "EIF4EBP1", "EIF4G1","EIF4G2", "EIF4G3",
         "EIF4H","EIF4B"), plot.bargraph.EIF.CNV.TCGA)

plot.matrix.EIF.CNV.corr(c("TP53", "EIF4A1","EIF4A2", "EIF4E", "EIF4E2", "EIF4E3", 
                           "MYC", "EIF3D", "EIF4EBP1", "EIF4G1", "EIF4G2", "EIF4G3",
                           "EIF4H", "EIF4B"))

plot.boxgraph.EIF.RNAseq.TCGA (c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
                                 "EIF4G1","EIF4G2","EIF4G3",
                                 "EIF4B","EIF4H",
                                 "EIF4A1","EIF4A2","EIF3D","TP53","MYC"))

lapply(c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
         "EIF4G1","EIF4G2","EIF4G3",
         "EIF4B","EIF4H",
         "EIF4A1","EIF4A2","EIF3D","TP53","MYC"),
       plot.boxgraph.EIF.CNVratio.TCGA)



## Figure 2 ##
###########################
##  KM survival analyses ##
###########################
plot.km.EIF.all.tumors <- function(EIF) {
  pan.TCGA.gene <- function(EIF) {
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      file.path(data.file.directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
      data.table = FALSE
    )
    TCGA.RNAseq1 <- TCGA.RNAseq[
      !duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))
    ]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1 <- as.data.frame(TCGA.RNAseq1)
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    # download https://xenabrowser.net/datapages/?dataset=Survival_SupplementalTable_S1_20171025_xena_sp&host=https%3A%2F%2Fpancanatlas.xenahubs.net
    TCGA.OS <- fread(
      file.path(data.file.directory, "Survival_SupplementalTable_S1_20171025_xena_sp"),
      data.table = FALSE
    )
    TCGA.OS1 <- TCGA.OS[
      !duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))
    ]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1 <- as.data.frame(TCGA.OS1)
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[, c("OS", "OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE
    )
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    TCGA.OS.sampletype$Row.names <- NULL
    TCGA.OS.sampletype$sample.type <- as.factor(TCGA.OS.sampletype$sample.type)
    # remove "solid tissue normal from dataset "
    TCGA.OS.sampletype <-
      TCGA.OS.sampletype[TCGA.OS.sampletype$sample.type != "Solid Tissue Normal", ]
    TCGA.OS.sampletype$sample.type <- droplevels(TCGA.OS.sampletype$sample.type)
    levels(TCGA.OS.sampletype$sample.type)
    
    ## combine OS, sample type and RNAseq data ##
    TCGA.RNAseq.OS.sampletype <- merge(TCGA.RNAseq_transpose,
                                       TCGA.OS.sampletype,
                                       by    = "row.names",
                                       all.x = TRUE
    )
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease
    )
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(EIF)
  # test for normal distribution by Anderson–Darling test
  plot(density(df[[EIF]]), xlab = paste(EIF, "log2(TPM)"), main = NULL)
  test <- ad.test(df[[EIF]])
  shapiro.test(df[[EIF]][0:5000])
  print(paste(EIF, test))
  
  plot.KM <- function(EIF) {
    # df <- subset(df, OS.time <= 4000)
    number <- nrow(df)
    sub <- round(number / 5, digits = 0)
    # bottom.label <- paste("Bottom 20%, n = ", sub)
    # top.label <- paste("Top 20%, n = ", sub)
    df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.2)] <- "Bottom 20%"
    df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.8)] <- "Top 20%"
    df$SurvObj <- with(df, Surv(OS.time, OS == 1))
    df <- na.omit(df)
    km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
    stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
    p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
    p.val <- signif(p.val, 3)
    
    
    KM <- ggplot2::autoplot(
      km,
      censor = FALSE,
      xlab = "Days",
      ylab = "Survival Probability",
      #main = "All TCGA cancer studies",
      main = paste0("All TCGA cancer studies (", number, " cases)"),
      # xlim = c(0, 4100),
      color = strata
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text = black_bold_16(),
        #axis.line.x = element_line(color = "black"),
        #axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        strip.text = black_bold_16(),
        legend.text = black_bold_16(),
        legend.title = black_bold_16(),
        legend.position = c(0.9, 0.98),
        legend.justification = c(1, 1)
      ) +
      guides(fill = FALSE) +
      scale_color_manual(
        values = c("red", "blue"),
        name = paste(EIF, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(
          paste("Bottom 20%, n =", sub),
          paste("Top 20%, n =", sub)
        )
      ) +
      # scale_x_continuous(expand = c(0, 0), limits = c(0, 4100)) +
      # scale_y_continuous(expand = c(0, 0),
      #                   limits = c(0, 1.05),
      #                   labels = scales::percent) +
      # geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 10000,
        y        = 0.75,
        label    = paste("log-rank test \n p.val = ", p.val),
        size     = 6.5,
        hjust    = 1,
        fontface = "bold"
      )
    
    print(KM)
    ggplot2::ggsave(
      path = file.path(output.directory, "KM"),
      filename = paste(EIF, " all tumors KM.pdf"),
      plot = KM,
      width = 6,
      height = 6,
      useDingbats = FALSE
    )
  }
  plot.KM(EIF)
}

plot.km.EIF.each.tumor <- function(EIF, tumor) {
  pan.TCGA.gene <- function(EIF, tumor) {
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      file.path(data.file.directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
      data.table = FALSE
    )
    TCGA.RNAseq1 <- TCGA.RNAseq[
      !duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))
    ]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1 <- as.data.frame(TCGA.RNAseq1)
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    # download https://xenabrowser.net/datapages/?dataset=Survival_SupplementalTable_S1_20171025_xena_sp&host=https%3A%2F%2Fpancanatlas.xenahubs.net
    TCGA.OS <- fread(
      file.path(data.file.directory, "Survival_SupplementalTable_S1_20171025_xena_sp"),
      data.table = FALSE
    )
    TCGA.OS1 <- TCGA.OS[
      !duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))
    ]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1 <- as.data.frame(TCGA.OS1)
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[, c("OS", "OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE
    )
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    TCGA.OS.sampletype$Row.names <- NULL
    TCGA.OS.sampletype$sample.type <- as.factor(TCGA.OS.sampletype$sample.type)
    # remove "solid tissue normal from dataset "
    TCGA.OS.sampletype <-
      TCGA.OS.sampletype[TCGA.OS.sampletype$sample.type != "Solid Tissue Normal", ]
    TCGA.OS.sampletype$sample.type <- droplevels(TCGA.OS.sampletype$sample.type)
    levels(TCGA.OS.sampletype$sample.type)
    
    ## combine OS, sample type and RNAseq data ##
    TCGA.RNAseq.OS.sampletype <- merge(TCGA.RNAseq_transpose,
                                       TCGA.OS.sampletype,
                                       by    = "row.names",
                                       all.x = TRUE
    )
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease
    )
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[TCGA.RNAseq.OS.sampletype$primary.disease %in% tumor, ]
    return(TCGA.RNAseq.OS.sampletype)
  }
  
  plot.KM <- function(EIF, tumor) {
    # df <- subset(df, OS.time <= 2000)
    df <- pan.TCGA.gene(EIF, tumor)
    number <- nrow(df)
    sub <- round(number / 5, digits = 0)
    # bottom.label <- paste("Bottom 20%, n = ", sub)
    # top.label <- paste("Top 20%, n = ", sub)
    df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.2)] <- "Bottom 20%"
    df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.8)] <- "Top 20%"
    df$SurvObj <- with(df, Surv(OS.time, OS == 1))
    df <- na.omit(df)
    km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
    stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
    p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
    p.val <- signif(p.val, 3)
    
    KM <- ggplot2::autoplot(
      km,
      censor = FALSE,
      xlab = "Days",
      ylab = "Survival Probability",
      # xlim = c(0, 2100),
      main = tumor,
      color = strata )+
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text = black_bold_16(),
        #axis.line.x = element_line(color = "black"),
        #axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        strip.text = black_bold_16(),
        legend.text = black_bold_16(),
        legend.title = black_bold_16(),
        legend.position = c(0.9, 0.98),
        legend.justification = c(1, 1)
      ) +
      guides(fill = FALSE) +
      scale_color_manual(
        values = c("red", "blue"),
        name = paste(EIF, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(
          paste("Bottom 20%, n =", sub),
          paste("Top 20%, n =", sub)
        )
      ) +
      # scale_x_continuous(expand = c(0, 0), limits = c(0, 2100)) +
      # scale_y_continuous(expand = c(0, 0),
      #                   limits = c(0, 1.05),
      #                   labels = scales::percent) +
      # geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 7000,
        y        = 0.75,
        label    = paste("log-rank test \n p.val = ", p.val),
        size     = 6.5,
        hjust    = 1,
        fontface = "bold"
      )
    print(KM)
    ggplot2::ggsave(
      path = file.path(output.directory, "KM"),
      filename = paste(EIF, tumor, "KM.pdf"),
      plot = KM,
      width = 6,
      height = 6,
      useDingbats = FALSE
    )
  }
  plot.KM(EIF, tumor)
}

##########################
## Cox regression model ##
##########################
plot.univariate <- function(df, covariate_names, data.np, output.file, plot.title, x.tics, x.range) {
  result <- map(
    vars(
      EIF4E, EIF4E2, EIF4E3, 
      EIF4G1, EIF4G2, EIF4G3, 
      EIF4A1, EIF4A2, EIF3D, 
      EIF3E, EIF4EBP1, EIF4EBP2, #PABPC1, 
      MKNK1, MKNK2, EIF4B, EIF4H, 
      MTOR, #RPS6KB1, 
      MYC
    ),
    function(by) {
      analyse_multivariate(df,
                           vars(OS.time, OS),
                           covariates = list(by), # covariates expects a list
                           covariate_name_dict = covariate_names
      )
    }
  )
  
  HR.table <- function(x) {
    list <- as.data.frame(result[[x]]["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
    return(list)
  }
  a <- lapply(c(1:18), HR.table)
  b <- do.call(rbind.data.frame, a)
  
  # Testing proportional Hazards assumption on univariate analysis
  univ_formulas <- sapply(
    covariate_names,
    function(x) as.formula(paste("Surv(OS.time, OS)~", x))
  )
  univ_models <- lapply(
    univ_formulas,
    function(x) {
      cox.zph(coxph(x, data = df))
    }
  )
  coxassump <- function(x) {
    c <- print(univ_models[[x]])
    return(c)
  }
  univ_results <- lapply(covariate_names, coxassump)
  d <- do.call(rbind.data.frame, univ_results)
  e <- d[-grep("GLOBAL", rownames(d)), ]
  rownames(e) <- gsub("\\..*", "", rownames(e))
  f <- e[, 3, drop = FALSE]
  colnames(f) <- "pinteraction"
  
  data <- merge(b, f, by.x = "factor.id", by.y = "row.names")
  data[, 4:11] <- round(data[, 4:11], digits = 3)
  data[, 4:6] <- round(data[, 4:6], digits = 2)
  data$np <- data.np
  data <- as.data.frame(data)
  data$HRCI <- paste0(data$HR, " (", data$Lower_CI, "-", data$Upper_CI, ")")
  data$p[data$p < 0.001] <- "<0.001"
  data$pinteraction[data$pinteraction < 0.001] <- "<0.001"
  data <- data[order(data$HR, decreasing = TRUE), ]
  tabletext1 <- cbind(
    c("Gene", data$factor.id),
    c("No. of\nPatients", data$np),
    c("Hazard Ratio\n(95% CI)", data$HRCI),
    c("P Value", data$p),
    c("P Value for\nInteraction", data$pinteraction)
  )
  
  pdf(
    file = output.file,
    width = 14,
    height = 12,
    onefile = F
  )
  p <- forestplot(
    labeltext = tabletext1,
    graph.pos = 3, graphwidth = unit(12, "cm"),
    hrzl_lines = list(
      "1" = gpar(lwd = 1, col = "black"),
      "2" = gpar(lwd = 1, col = "black")
    ),
    mean = c(NA, data$HR),
    lower = c(NA, data$Lower_CI),
    upper = c(NA, data$Upper_CI),
    title = plot.title,
    xlab = "<---Good prognosis---    ---Poor prognosis--->",
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = fpColors(box = "black", lines = "black"),
    xticks = x.tics,
    # xlog = 0,
    clip = x.range,
    zero = 1,
    cex = 1.2,
    lineheight = "auto", # height of the graph
    boxsize = 0.2,
    colgap = unit(6, "mm"), # the gap between column
    lwd.ci = 2,
    ci.vertices = FALSE,
    ci.vertices.height = 0.02,
    new_page = getOption("forestplot_new_page", FALSE)
  )
  dev.off()
}

plot.multivariate <- function(df, covariate_names, data.np, output.file, plot.title, x.tics, x.range) {
  df %>%
    analyse_multivariate(vars(OS.time, OS), 
                         covariates = vars(EIF4E, EIF4E2, EIF4E3, 
                                           EIF4G1, EIF4G2, EIF4G3, 
                                           EIF4A1, EIF4A2, EIF3D, 
                                           EIF3E, EIF4EBP1, EIF4EBP2, #PABPC1, 
                                           MKNK1, MKNK2, EIF4B, EIF4H,
                                           MTOR, #RPS6KB1, 
                                           MYC), 
                         covariate_name_dict = covariate_names) -> result1
  data <- as.data.frame(result1["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
  
  # Testing proportional Hazards assumption on univariate analysis
  mv_fit <- coxph(Surv(OS.time, OS) ~  EIF4E + EIF4E2 + EIF4E3 + 
                    EIF4G1 + EIF4G2 + EIF4G3 + 
                    EIF4A1 + EIF4A2 + EIF3D + 
                    EIF3E + EIF4EBP1 + EIF4EBP2 + #PABPC1 + 
                    MKNK1 + MKNK2 + EIF4B + EIF4H +
                    MTOR + #RPS6KB1 + 
                    MYC, data = df)
  test.ph <- cox.zph(mv_fit)
  test <- print(test.ph)
  f <- test[, 3, drop = FALSE]
  colnames(f) <- "pinteraction"
  
  data <- merge(data, f, by.x = "factor.id", by.y = "row.names")
  data[, 4:11] <- round(data[, 4:11], digits = 3)
  data[, 4:6] <- round(data[, 4:6], digits = 2)
  data$np <- data.np
  data <- as.data.frame(data)
  data$HRCI <- paste0(data$HR, " (", data$Lower_CI, "-", data$Upper_CI, ")")
  data$p[data$p < 0.001] <- "<0.001"
  data$pinteraction[data$pinteraction < 0.001] <- "<0.001"
  data <- data[order(data$HR, decreasing = TRUE), ]
  tabletext1 <- cbind(
    c("Gene", data$factor.id),
    c("No. of\nPatients", data$np),
    c("Hazard Ratio\n(95% CI)", data$HRCI),
    c("P Value", data$p),
    c("P Value for\nInteraction", data$pinteraction)
  )
  
  pdf(
    file = output.file,
    width = 14,
    height = 12,
    onefile = F
  )
  p <- forestplot(
    labeltext = tabletext1,
    graph.pos = 3, graphwidth = unit(12, "cm"),
    hrzl_lines = list(
      "1" = gpar(lwd = 1, col = "black"),
      "2" = gpar(lwd = 1, col = "black")
    ),
    #  "3.75" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
    #  "3" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
    #  "7" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
    #  "9" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922")),
    mean = c(NA, data$HR),
    lower = c(NA, data$Lower_CI),
    upper = c(NA, data$Upper_CI),
    title = plot.title,
    xlab = "     <---Good prognosis---    ---Poor prognosis--->",
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = fpColors(box = "black", lines = "black"),
    xticks = x.tics,
    clip = x.range,
    zero = 1,
    cex = 1.2,
    lineheight = "auto", # height of the graph
    boxsize = 0.2,
    colgap = unit(6, "mm"), # the gap between column
    lwd.ci = 2,
    ci.vertices = FALSE,
    ci.vertices.height = 0.02,
    new_page = getOption("forestplot_new_page", FALSE)
  )
  dev.off()
  print(p)
}

plot.coxph.EIF.all.tumors <- function() {
  pan.TCGA.gene <- function(EIF) {
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      file.path(data.file.directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
      data.table = FALSE
    )
    TCGA.RNAseq1 <- TCGA.RNAseq[
      !duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))
    ]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1 <- as.data.frame(TCGA.RNAseq1)
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    TCGA.OS <- fread(
      file.path(data.file.directory, "Survival_SupplementalTable_S1_20171025_xena_sp"),
      data.table = FALSE
    )
    TCGA.OS1 <- TCGA.OS[
      !duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))
    ]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1 <- as.data.frame(TCGA.OS1)
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[, c("OS", "OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE
    )
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    TCGA.OS.sampletype$Row.names <- NULL
    TCGA.OS.sampletype$sample.type <- as.factor(TCGA.OS.sampletype$sample.type)
    # remove "solid tissue normal from dataset "
    TCGA.OS.sampletype <-
      TCGA.OS.sampletype[TCGA.OS.sampletype$sample.type != "Solid Tissue Normal", ]
    TCGA.OS.sampletype$sample.type <- droplevels(TCGA.OS.sampletype$sample.type)
    levels(TCGA.OS.sampletype$sample.type)
    
    ## combine OS, sample type and RNAseq data ##
    TCGA.RNAseq.OS.sampletype <- merge(TCGA.RNAseq_transpose,
                                       TCGA.OS.sampletype,
                                       by    = "row.names",
                                       all.x = TRUE
    )
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease
    )
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(c(
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF4G1", "EIF4G2", "EIF4G3", 
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2", #"PABPC1", 
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", #"RPS6KB1", 
    "MYC"
  ))
  
  covariate_names <- c(
    EIF4E = "EIF4E", EIF4E2 = "EIF4E2", EIF4E3 = "EIF4E3", 
    EIF4G1 = "EIF4G1", EIF4G2 = "EIF4G2", EIF4G3 = "EIF4G3",
    EIF4A1 = "EIF4A1", EIF4A2 = "EIF4A2", EIF3D = "EIF3D", 
    EIF3E = "EIF3E", EIF4EBP1 = "EIF4EBP1", EIF4EBP2 = "EIF4EBP2",#PABPC1 = "PABPC1", 
    MKNK1 = "MKNK1", MKNK2 = "MKNK2", EIF4B = "EIF4B", EIF4H = "EIF4H",
    MTOR = "MTOR", #RPS6KB1 = "RPS6KB1", 
    MYC = "MYC"
  )
  
  plot.univariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 10235,
    output.file = file.path(output.directory, "Cox", "EIFUniCox.pdf"),
    plot.title = "Univariate Cox proportional-hazards regression analysis (all tumor types)",
    x.tics = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8),
    x.range = c(0.6, 1.8))
  
  plot.multivariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 10235,
    output.file = file.path(output.directory, "Cox", "EIFmultiCox.pdf"),
    plot.title = "Multivariate Cox proportional-hazards regression analysis (all tumor types)",
    x.tics = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8),
    x.range = c(0.6, 1.8))
}

plot.coxph.EIF.each.tumor <- function(tumor) {
  pan.TCGA.gene <- function(EIF, tumor) {
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      file.path(data.file.directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
      data.table = FALSE
    )
    TCGA.RNAseq1 <- TCGA.RNAseq[
      !duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))
    ]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1 <- as.data.frame(TCGA.RNAseq1)
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    TCGA.OS <- fread(
      file.path(data.file.directory, "Survival_SupplementalTable_S1_20171025_xena_sp"),
      data.table = FALSE
    )
    TCGA.OS1 <- TCGA.OS[
      !duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))
    ]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1 <- as.data.frame(TCGA.OS1)
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[, c("OS", "OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE
    )
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    TCGA.OS.sampletype$Row.names <- NULL
    TCGA.OS.sampletype$sample.type <- as.factor(TCGA.OS.sampletype$sample.type)
    # remove "solid tissue normal from dataset "
    TCGA.OS.sampletype <-
      TCGA.OS.sampletype[TCGA.OS.sampletype$sample.type != "Solid Tissue Normal", ]
    TCGA.OS.sampletype$sample.type <- droplevels(TCGA.OS.sampletype$sample.type)
    levels(TCGA.OS.sampletype$sample.type)
    
    ## combine OS, sample type and RNAseq data ##
    TCGA.RNAseq.OS.sampletype <- merge(TCGA.RNAseq_transpose,
                                       TCGA.OS.sampletype,
                                       by    = "row.names",
                                       all.x = TRUE
    )
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease
    )
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[TCGA.RNAseq.OS.sampletype$primary.disease %in% tumor, ]
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(c(
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF4G1", "EIF4G2", "EIF4G3", 
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2", #"PABPC1", 
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", #"RPS6KB1", 
    "MYC"
  ), tumor)
  # df$`EIF4E+EIF4EBP1` <- log2(2**df$EIF4E + 2**df$EIF4EBP1 -2 + 1)
  df <- df[c(
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF4G1", "EIF4G2", "EIF4G3", 
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2",#"PABPC1", 
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", #"RPS6KB1", 
    "MYC", # "EIF4E+EIF4EBP1",
    "OS", "OS.time"
  )]
  # Use survivalAnalysis package to draw forest plot of multiple univariate #
  covariate_names <- c(
    EIF4E = "EIF4E", EIF4E2 = "EIF4E2", EIF4E3 = "EIF4E3", 
    EIF4G1 = "EIF4G1", EIF4G2 = "EIF4G2", EIF4G3 = "EIF4G3",
    EIF4A1 = "EIF4A1", EIF4A2 = "EIF4A2", EIF3D = "EIF3D", 
    EIF3E = "EIF3E", EIF4EBP1 = "EIF4EBP1", EIF4EBP2 = "EIF4EBP2", #PABPC1 = "PABPC1", 
    MKNK1 = "MKNK1", MKNK2 = "MKNK2", EIF4B = "EIF4B", EIF4H = "EIF4H",
    MTOR = "MTOR", #RPS6KB1 = "RPS6KB1", 
    MYC = "MYC"
  )
  
  # TODO: Use file.path below (noting that paste0() has been convenient
  #       for filename construction in the current implementation.)
  plot.univariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 517,
    output.file = paste0(output.directory, "/Cox/", tumor, "EIFUniCox.pdf"),
    plot.title = paste("Univariate Cox proportional-hazards regression analysis", tumor),
    x.tics = c(0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8),
    x.range = c(0.4, 2.8))
  
  plot.multivariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 517,
    output.file = paste0(output.directory, "/Cox/", tumor, "EIFmultiCox.pdf"),
    plot.title = paste("Multivariate Cox proportional-hazards regression analysis", tumor),
    x.tics = c(0.4, 0.8, 1.2, 1.6, 2.4, 2.8, 3.2),
    x.range = c(0.4, 3.2))
}


# Figure 2
plot.km.EIF.all.tumors("HNRNPL")
lapply(c("EIF4G1","EIF4G2", "EIF4G3",
         "EIF4A1","EIF4A2", 
         "EIF4E", "EIF4E2", "EIF4E3", 
         "EIF3D", "EIF3E",
         "EIF4EBP1", "EIF4EBP2", 
         "EIF4H", "EIF4B", "MYC",
         "PABPC1", "MKNK1", "MKNK2"), 
       plot.km.EIF.all.tumors)

lapply(
  c(
    "EIF4G1","EIF4G2", "EIF4G3",
    "EIF4A1","EIF4A2", 
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF3D", "EIF3E","EIF4EBP1", "EIF4EBP2", 
    "EIF4H", "EIF4B", "MYC",
    "PABPC1", "MKNK1", "MKNK2"
  ),
  plot.km.EIF.each.tumor,
  tumor = "lung adenocarcinoma"
)

plot.coxph.EIF.all.tumors()

plot.coxph.EIF.each.tumor(c("lung adenocarcinoma"))

## Figure 3 ##
##############################################
## boxplot for EIF expression across tumors ##
##############################################
plot.boxgraph.EIF.RNAseq.TCGA.GTEX <- function(EIF.gene) {
  
  pan.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.pancancer <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    )
    # use sample column for the rowname
    TCGA.pancancer1 <- TCGA.pancancer[
      !duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))
    ]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)

    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.sampletype <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    subset <- TCGA.sampletype[, c(
      "_sample_type",
      "primary disease or tissue",
      "_primary_site",
      "_study"
    ),
    drop = FALSE
    ]
    row.names(subset) <- row.names(TCGA.sampletype)
    colnames(subset) <- c(
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    )


    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
      subset,
      by    = "row.names",
      all.x = TRUE
    )
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)


    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[, colnames(TCGA.RNAseq.anno) %in% c(
      EIF.gene,
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    )]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) + 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) - 1)
    EIF.TCGA.RNAseq.anno.subset$`EIF4E-EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) - 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1))
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(
      EIF.gene,
      "EIF4E+EIF4EBP1",
      "EIF4E-EIF4EBP1",
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    )]
    # EIF.TCGA.RNAseq.anno.subset$delta <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4G1 - 2**EIF.TCGA.RNAseq.anno.subset$EIF4E - 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!is.na(EIF.TCGA.RNAseq.anno.subset$primary.site), ]

    EIF.TCGA.RNAseq.anno.subset.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long <- EIF.TCGA.RNAseq.anno.subset.long[
      !EIF.TCGA.RNAseq.anno.subset.long$value == 0,
    ]
    EIF.TCGA.RNAseq.anno.subset.long$sample.type <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$sample.type
    )
    EIF.TCGA.RNAseq.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.disease
    )
    EIF.TCGA.RNAseq.anno.subset.long$primary.site <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.site
    )

    return(EIF.TCGA.RNAseq.anno.subset.long)
  }
  TCGA.RNAseq.anno <- pan.TCGA.gene()

  pancancer.TCGA.EIF <- function() {
    TCGA.RNAseq.anno <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "TCGA", ]
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno [
      !TCGA.RNAseq.anno$sample.type %in% c("Solid Tissue Normal"),
    ]
    return(TCGA.RNAseq.anno.subset)
  }
  pancancer.TCGA.EIF.long <- pancancer.TCGA.EIF()

  GTEX.EIF <- function() {
    TCGA.RNAseq.anno <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "GTEX", ]
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[
      TCGA.RNAseq.anno$sample.type %in% "Normal Tissue",
    ]
    return(TCGA.RNAseq.anno.subset)
  }
  GTEX.EIF.long <- GTEX.EIF()

  # reorder bars by explicitly ordering factor levels
  make.plot <- function() {
      pancancer.TCGA.EIF.long1 <- pancancer.TCGA.EIF.long[
        pancancer.TCGA.EIF.long$variable %in% EIF.gene,
      ]
      levels(pancancer.TCGA.EIF.long1$variable)
      pancancer.TCGA.EIF.long1$variable <- factor(pancancer.TCGA.EIF.long1$variable,
        levels = rev(c(
          "EIF4EBP1","EIF4E3","EIF4E2","EIF4E", 
          #"EIF4E-EIF4EBP1",
          #"EIF4E+EIF4EBP1",
          "EIF4G3","EIF4G2","EIF4G1","EIF4A2", "EIF4A1"
        ))
      )

      f1 <- factor(pancancer.TCGA.EIF.long1$primary.disease)
      #f.ordered1 <- fct_rev(f1)
      f.ordered1 <- factor(f1)
      p1 <- ggplot(data = pancancer.TCGA.EIF.long1,
                  aes(
                      x = f.ordered1,
                      # x     = x.ordered, # order primary disease
                      y = 2**value)
                ) +
        scale_y_continuous(
          trans = log2_trans(),
          limits = c(2**4, 2**17),
          labels = label_comma()
        ) +
        stat_n_text(
          size = 5,
          fontface = "bold",
          hjust = 0
        ) +
        geom_boxplot(aes(
          colour = factor(variable),
          #fill = factor(variable)
        ),
          # alpha = .01,
          #width    = 0.95,
          outlier.shape = NA,
          position = position_dodge(width = 1)
        ) +
        #scale_color_brewer(palette="Dark2") +
        #scale_fill_brewer(palette = "Set1") +
        scale_color_manual(
          values = c(
            "EIF4EBP1" = "#B997C7",
            "EIF4E3" = "#824D99",
            "EIF4E2" = "#4E78C4",
            "EIF4E" = "#57A2AC", 
            "EIF4G3" = "#7EB875",
            "EIF4G2" = "#D0B541",
            "EIF4G1" = "#E67F33",
            "EIF4A2" = "#CE2220", 
            "EIF4A1" = "#521A13"
          ),
          breaks = rev(c("EIF4EBP1", "EIF4E3","EIF4E2","EIF4E",
                         "EIF4G3", "EIF4G2","EIF4G1","EIF4A2","EIF4A1"))
        ) +
        #scale_fill_manual(
        #  values = c(
        #    "EIF4EBP1" = "#CC79A7",
        #    "EIF4E" = "#0072B2",
        #    "EIF4G1" = "#009E73",
        #    "EIF4A1" = "#D55E00"
        #  ),
        #  breaks = c("EIF4EBP1", "EIF4E","EIF4E2","EIF4E3",
        #             "EIF4G1", "EIF4G2","EIF4G3","EIF4A1","EIF4A2")
        #) +
        # scale_color_manual(
        #  values = c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00"),
        #  breaks = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1"),
        #  labels = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1")) + #for color-blind palettes
        labs(
          x = "primary disease",
          y = paste("Normalized expression (RNA-Seq counts)")
        ) +
        #coord_flip() +
        theme_bw() +
        theme(
          plot.title = black_bold_12(),
          axis.title.x = element_blank(),
          axis.title.y = black_bold_12(),
          axis.text.x = black_bold_12_45(),
          axis.text.y = black_bold_12(),
          # axis.line.x          = element_line(color = "black"),
          # axis.line.y          = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = black_bold_12(),
          legend.position = c(0, 0.15),
          legend.direction = "horizontal",
          legend.justification = c(0, 1),
          legend.box = "horizontal",
          strip.text = black_bold_12()
        )
      print(p1)
      
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorexpression.pdf",
        plot = p1,
        width = 16.5,
        height = 8,
        useDingbats = FALSE)
  }
  make.plot()
}

plot.boxgraph.EIF.RNAseq.TCGA.GTEX (c("EIF4A1","EIF4A2", "EIF4E", "EIF4E2", "EIF4E3", 
                                      "EIF4EBP1", "EIF4G1","EIF4G2", "EIF4G3"))

##########################################
## boxplot for EIF ratios across tumors ##
##########################################
plot.boxgraph.EIF.ratio.TCGA.GTEX.2 <- function(EIF.gene) {
  EIF.gene <- c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
    "EIF4G1","EIF4G2","EIF4G3","EIF3D",
    "EIF4A1","EIF4A2")
  pan.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.pancancer <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    )
    TCGA.pancancer1 <- TCGA.pancancer[
      !duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))
      ]
    
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.sampletype <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    subset <- TCGA.sampletype[, c(
      "_sample_type",
      "primary disease or tissue",
      "_primary_site", 
      "_study"
    ),
    drop = FALSE
    ]
    row.names(subset) <- row.names(TCGA.sampletype)
    colnames(subset) <- c("sample.type", 
                          "primary.disease", 
                          "primary.site", 
                          "study")
    
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    subset,
                                    by    = "row.names",
                                    all.x = TRUE
    )
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    # TCGA.RNAseq.anno <- na.omit(TCGA.RNAseq.anno)
    
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    levels(TCGA.RNAseq.anno$primary.disease)
    TCGA.RNAseq.anno$primary.site <- as.factor(TCGA.RNAseq.anno$primary.site)
    levels(TCGA.RNAseq.anno$primary.site)
    row.names(TCGA.RNAseq.anno) <- TCGA.RNAseq.anno$Row.names
    TCGA.RNAseq.anno$Row.names <- NULL
    
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[
      ,
      colnames(TCGA.RNAseq.anno) %in% c(
        EIF.gene,
        "sample.type",
        "primary.disease",
        "primary.site",
        "study"
      )
      ]
    
    plot.cor <- function(){
      EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[
        ,
        colnames(TCGA.RNAseq.anno) %in% c(
          EIF.gene,
          "sample.type",
          "primary.disease",
          "primary.site",
          "study"
        )
        ]
      
      EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[EIF.TCGA.RNAseq.anno.subset$study == "TCGA", ]
      EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[EIF.TCGA.RNAseq.anno.subset$sample.type != "Solid Tissue Normal", ]
      EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[
        ,
        colnames(TCGA.RNAseq.anno) %in% c(
          EIF.gene)
        ]
      
      cor_5 <- rcorr(as.matrix(EIF.TCGA.RNAseq.anno.subset), type = "pearson")
      M <- cor_5$r
      p_mat <- cor_5$P
      #pdf(file.path(output.directory, "CNV", "EIFCNVcormatrix.pdf"),
      #    width = 9,
      #    height = 9,
      #    useDingbats = FALSE
      #)
      corrplot(
        M,
        method      = "color",
        tl.cex      = 1,
        number.cex  = 1,
        addgrid.col = "gray",
        addCoef.col = "black",
        tl.col      = "black",
        #type        = "lower",
        order       = "FPC", tl.srt = 90,
        p.mat       = p_mat,
        sig.level   = 0.05, # insig = "blank"
    )
    }
    #plot.cor()
    
    
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    EIF.TCGA.RNAseq.anno.subset$sum <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - 1)
    EIF.TCGA.RNAseq.anno.subset$minus <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) - 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) + 1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!is.na(EIF.TCGA.RNAseq.anno.subset$primary.site), ]
    
    EIF.TCGA.GTEX.score <- EIF.TCGA.RNAseq.anno.subset
    # A:E
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4E2)
    EIF.TCGA.GTEX.score$`EIF4A2:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A2 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4A2:\nEIF4E2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A2 - EIF.TCGA.RNAseq.anno.subset$EIF4E2)
    
    # G:E
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E2)
    EIF.TCGA.GTEX.score$`EIF4G3:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G3 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G3:\nEIF4E2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G3 - EIF.TCGA.RNAseq.anno.subset$EIF4E2)
    
    # A:G
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4G1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4A2:\nEIF4G1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A2 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4G2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4G2)
    EIF.TCGA.GTEX.score$`EIF4A2:\nEIF4G2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A2 - EIF.TCGA.RNAseq.anno.subset$EIF4G2)
    
    # EBP:E
    EIF.TCGA.GTEX.score$`EIF4E:\nEIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4E - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4EBP1:\nEIF4E2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - EIF.TCGA.RNAseq.anno.subset$EIF4E2)
    
    # E2:E1
    EIF.TCGA.GTEX.score$`EIF4E2:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4E2 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    
    # G2:G1
    EIF.TCGA.GTEX.score$`EIF4G2:\nEIF4G1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G2 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4G3` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4G3)
    EIF.TCGA.GTEX.score$`EIF4G2:\nEIF3D` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G2 - EIF.TCGA.RNAseq.anno.subset$EIF3D)
    
      # A2:A1
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4A2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4A2)
    
    # EIF.TCGA.GTEX.score$`PABPC1:\nEIF4E` <-
    #  (EIF.TCGA.RNAseq.anno.subset$PABPC1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E+EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E+EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$sum)
    
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E-EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$minus)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E-EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$minus)
    
    EIF.TCGA.GTEX.score$`EIF4E:\nEIF4E+EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4E - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4EBP1:\nEIF4E+EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - EIF.TCGA.RNAseq.anno.subset$sum)

    
    ratio <- c("EIF4A1:\nEIF4E", "EIF4A1:\nEIF4E2", 
               "EIF4A2:\nEIF4E", "EIF4A2:\nEIF4E2",          
               "EIF4G1:\nEIF4E", "EIF4G1:\nEIF4E2", 
               "EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2", 
               "EIF4A1:\nEIF4G1", "EIF4A2:\nEIF4G1", 
               "EIF4A1:\nEIF4G2", "EIF4A2:\nEIF4G2",
               "EIF4E:\nEIF4EBP1", "EIF4EBP1:\nEIF4E2",        
               "EIF4E2:\nEIF4E", "EIF4G2:\nEIF4G1", "EIF4G1:\nEIF4G3", "EIF4G2:\nEIF3D", "EIF4A1:\nEIF4A2",          
               "EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E-EIF4EBP1", "EIF4A1:\nEIF4E-EIF4EBP1" )
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.score[, c(
      ratio, "sample.type", "primary.disease", "primary.site", "study"),
    drop = FALSE
    ]
    EIF.TCGA.GTEX.score <- na.omit(EIF.TCGA.GTEX.score)
    return(EIF.TCGA.GTEX.score)
  }
  TCGA.RNAseq.anno <- pan.TCGA.gene()
  
  pancancer.TCGA.EIF.ratio <- function() {
    
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "TCGA", ]
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno.subset [
      !TCGA.RNAseq.anno.subset$sample.type %in% c("Cell Line"),
      ]
    TCGA.RNAseq.anno.subset$sample.type <- as.character(TCGA.RNAseq.anno.subset$sample.type)
    TCGA.RNAseq.anno.subset$sample.type[TCGA.RNAseq.anno.subset$sample.type != "Solid Tissue Normal"] <- "Tumor"
    TCGA.RNAseq.anno.subset$sample.type[TCGA.RNAseq.anno.subset$sample.type == "Solid Tissue Normal"] <- "NAT"
    TCGA.RNAseq.anno.subset$sample.type <- as.factor(TCGA.RNAseq.anno.subset$sample.type)
    
    
    #TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "TCGA", ]
    # TCGA.RNAseq.anno.subset$study <- NULL
    #TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno.subset[
    #  !TCGA.RNAseq.anno.subset$sample.type %in% "Solid Tissue Normal",
    #  ]
    EIF.TCGA.GTEX.score.long <- melt(TCGA.RNAseq.anno.subset)
    # reorder bars by explicitly ordering factor levels
    EIF.TCGA.GTEX.score.long$primary.disease <- as.factor(
      EIF.TCGA.GTEX.score.long$primary.disease
    )
    EIF.TCGA.GTEX.score.long <- droplevels(EIF.TCGA.GTEX.score.long)
    levels(EIF.TCGA.GTEX.score.long$primary.disease)
    return(EIF.TCGA.GTEX.score.long)
  }
  pancancer.TCGA.EIF.ratio.long <- pancancer.TCGA.EIF.ratio()
  
  make.plot <- function() {
    make.plot1 <- function() {
      #pancancer.TCGA.EIF.ratio.long$label <- sub(
      #  ".*\n", "",
      #  pancancer.TCGA.EIF.ratio.long$variable
      #)
      #pancancer.TCGA.EIF.ratio.long$label <- factor(
      #  pancancer.TCGA.EIF.ratio.long$label,
      #  levels = c("EIF4G1","EIF4G2", "EIF4EBP1", "EIF4E",, "EIF4E2")
      #)
      pancancer.TCGA.EIF.ratio.long1 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4G1:\nEIF4E","EIF4A1:\nEIF4E","EIF4A2:\nEIF4E",
          "EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2" ,"EIF4G1:\nEIF4G3"
        ),
        ]
      pancancer.TCGA.EIF.ratio.long1$variable <- factor(
        pancancer.TCGA.EIF.ratio.long1$variable,
        levels = c(
          "EIF4G1:\nEIF4E","EIF4A1:\nEIF4E","EIF4A2:\nEIF4E",
          #"EIF4G1:\nEIF4E2","EIF4A1:\nEIF4E2","EIF4A2:\nEIF4E2",
          "EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2" , "EIF4G1:\nEIF4G3"
          #"EIF4A1:\nEIF4E2","EIF4A2:\nEIF4E2", "EIF4EBP1:\nEIF4E2",
          #"EIF4E2:\nEIF4E", #"EIF4G2:\nEIF4G1", "EIF4A2:\nEIF4A1"
          #"EIF4A1:\nEIF4G1", "EIF4A2:\nEIF4G1",
          #"EIF4A1:\nEIF4G2", "EIF4A2:\nEIF4G2"
          #"EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1",
          #"EIF4G1:\nEIF4E-EIF4EBP1", "EIF4A1:\nEIF4E-EIF4EBP1"
        )
      )

      
      f1 <- factor(pancancer.TCGA.EIF.ratio.long1$primary.disease)
      f.ordered1 <- fct_rev(f1)
      
      p1 <- ggplot(
        data = pancancer.TCGA.EIF.ratio.long1,
        aes(
          x = f.ordered1,
          y = 2**value,
          # fill  = variable,
          color = sample.type
        )
      ) +
        geom_boxplot(
          #alpha = .1,
          #fill = sample.type,
          outlier.shape = NA,
          #size = .75,
          #width = 1,
          position = "dodge"
        ) +
        #geom_hline(
        #  data = pancancer.TCGA.EIF.ratio.long1,
        #  aes(yintercept = hline),
        #  linetype = "dashed"
        #) +
        # geom_text(data = pancancer.TCGA.EIF.ratio.long1, aes(f.ordered1, hline, label = hline), vjust = 2, hjust = -0.2) +
        scale_color_manual(values = c("Tumor" = "#CC79A7", "NAT" = "#0072B2"),
                           breaks = c("Tumor", "NAT")) +
        # for color-blind palettes
        facet_wrap(~variable, scales = "free_x") +
        ggplot2::facet_grid(~variable, scales = "free_x", space = "free") +
        # ggplot2::facet_grid_sc(cols = vars(label), shrink = TRUE,
        #              scales = list(y = scales_y),
        #              space = "free") +
        guides(colour = guide_legend(nrow = 1)) +
        labs(x = "primary disease",
             y = "Ratio of RNA counts"
             ) +
        coord_flip(ylim = c(0,25)) +
        theme_bw() +
        theme(
          plot.title = black_bold_12(),
          axis.title.x = black_bold_12(),
          axis.title.y = element_blank(),
          axis.text.x = black_bold_12(),
          axis.text.y = black_bold_12(),
          #axis.line.x = element_line(color = "black"),
          #axis.line.y = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.justification = "left",
          legend.box = "horizontal",
          legend.text = black_bold_12(),
          #strip.background = element_blank(),
          strip.text.x = black_bold_12()
        )
      print(p1)
      
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorratio.pdf",
        plot = p1,
        width = 18,
        height = 8,
        useDingbats = FALSE)
      
      
      pancancer.TCGA.EIF.ratio.long2 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4G2:\nEIF4G1", 
          "EIF4E2:\nEIF4E", "EIF4A1:\nEIF4A2",
          "EIF4E:\nEIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1"#,
          #"EIF4G1:\nEIF4E-EIF4EBP1", "EIF4A1:\nEIF4E-EIF4EBP1"
        ),
        ]
      pancancer.TCGA.EIF.ratio.long2$variable <- factor(
        pancancer.TCGA.EIF.ratio.long2$variable,
        levels = c(
          "EIF4G2:\nEIF4G1", 
          "EIF4E2:\nEIF4E", "EIF4A1:\nEIF4A2",
          "EIF4E:\nEIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1"
        )
      )

      
      f2 <- factor(pancancer.TCGA.EIF.ratio.long2$primary.disease)
      f.ordered2 <- fct_rev(f2)
      
      p2 <- ggplot(
        data = pancancer.TCGA.EIF.ratio.long2,
        aes(
          x = f.ordered2,
          y = 2**value,
          # fill  = variable,
          color = sample.type
        )
      ) +
        geom_boxplot(
          #alpha = .1,
          #fill = sample.type,
          outlier.shape = NA,
          #size = .75,
          #width = 1,
          position = "dodge"
        ) +
        geom_hline(
          data = pancancer.TCGA.EIF.ratio.long2,
          aes(yintercept = 4),
          linetype = "dashed"
        ) +
        # geom_text(data = pancancer.TCGA.EIF.ratio.long1, aes(f.ordered1, hline, label = hline), vjust = 2, hjust = -0.2) +
        scale_color_manual(values = c("Tumor" = "#CC79A7", "NAT" = "#0072B2"),
                           breaks = c("Tumor", "NAT")) +
        # for color-blind palettes
        facet_wrap(~variable, scales = "free_x") +
        ggplot2::facet_grid(~variable, scales = "free_x", space = "free") +
        # ggplot2::facet_grid_sc(cols = vars(label), shrink = TRUE,
        #              scales = list(y = scales_y),
        #              space = "free") +
        guides(colour = guide_legend(nrow = 1)) +
        labs(x = "primary disease",
             y = "Ratio of RNA counts"
        ) +
        coord_flip(ylim = c(0,25)) +
        theme_bw() +
        theme(
          plot.title = black_bold_12(),
          axis.title.x = black_bold_12(),
          axis.title.y = element_blank(),
          axis.text.x = black_bold_12(),
          axis.text.y = black_bold_12(),
          #axis.line.x = element_line(color = "black"),
          #axis.line.y = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.justification = "left",
          legend.box = "horizontal",
          legend.text = black_bold_12(),
          #strip.background = element_blank(),
          strip.text.x = black_bold_12()
        )
      print(p2)
      
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorratio2.pdf",
        plot = p2,
        width = 18,
        height = 8,
        useDingbats = FALSE)
      
      
      pancancer.TCGA.EIF.ratio.long3 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4E:\nEIF4EBP1", "EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2", 
          "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", "EIF4A1:\nEIF4A2"
        ),
        ]
      pancancer.TCGA.EIF.ratio.long3$variable <- factor(
        pancancer.TCGA.EIF.ratio.long3$variable,
        levels = c("EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2",
                   "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", 
                   "EIF4A1:\nEIF4A2","EIF4E:\nEIF4EBP1")
        )
      f3 <- factor(pancancer.TCGA.EIF.ratio.long3$primary.disease)
      f.ordered3 <- fct_rev(f3)
      p3 <- ggplot(
        data = pancancer.TCGA.EIF.ratio.long3,
        aes(
          x = f.ordered3,
          y = 2**value,
          # fill  = variable,
          color = sample.type
        )
      ) +
        geom_boxplot(
          #alpha = .1,
          #fill = sample.type,
          outlier.shape = NA,
          #size = .75,
          #width = 1,
          position = "dodge"
        ) +
        geom_hline(
          data = pancancer.TCGA.EIF.ratio.long3,
          aes(yintercept = 1),
          linetype = "dashed"
        ) +
        # geom_text(data = pancancer.TCGA.EIF.ratio.long1, aes(f.ordered1, hline, label = hline), vjust = 2, hjust = -0.2) +
        scale_color_manual(values = c("Tumor" = "#CC79A7", "NAT" = "#0072B2"),
                           breaks = c("Tumor", "NAT")) +
        # for color-blind palettes
        facet_wrap(~variable, scales = "free_x") +
        ggplot2::facet_grid(~variable, scales = "free_x", space = "free") +
        # ggplot2::facet_grid_sc(cols = vars(label), shrink = TRUE,
        #              scales = list(y = scales_y),
        #              space = "free") +
        guides(colour = guide_legend(nrow = 1)) +
        labs(x = "primary disease",
             y = "Ratio of RNA counts"
        ) +
        coord_flip(ylim = c(0,5)) +
        theme_bw() +
        theme(
          plot.title = black_bold_12(),
          axis.title.x = black_bold_12(),
          axis.title.y = element_blank(),
          axis.text.x = black_bold_12(),
          axis.text.y = black_bold_12(),
          #axis.line.x = element_line(color = "black"),
          #axis.line.y = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.justification = "left",
          legend.box = "horizontal",
          legend.text = black_bold_12(),
          #strip.background = element_blank(),
          strip.text.x = black_bold_12()
        )
      print(p3)
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorratio3.pdf",
        plot = p3,
        width = 18,
        height = 8,
        useDingbats = FALSE)
      
      
      pancancer.TCGA.EIF.ratio.long4 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4E:\nEIF4EBP1", "EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2",
          "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E","EIF4G2:\nEIF3D" 
        ),
        ]
      pancancer.TCGA.EIF.ratio.long4$variable <- factor(
        pancancer.TCGA.EIF.ratio.long4$variable,
        levels = c("EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2",
                   "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", 
                    "EIF4E:\nEIF4EBP1", "EIF4G2:\nEIF3D")
      )
      
      f4 <- factor(pancancer.TCGA.EIF.ratio.long4$primary.disease)
      f.ordered4 <- fct_rev(f4)
      
      p4 <- ggplot(
        data = pancancer.TCGA.EIF.ratio.long4,
        aes(
          x = f.ordered4,
          y = 2**value,
          # fill  = variable,
          color = sample.type
        )
      ) +
        geom_boxplot(
          #alpha = .1,
          #fill = sample.type,
          outlier.shape = NA,
          #size = .75,
          #width = 1,
          position = "dodge"
        ) +
        geom_hline(
          data = pancancer.TCGA.EIF.ratio.long4,
          aes(yintercept = 1),
          linetype = "dashed"
        ) +
        # geom_text(data = pancancer.TCGA.EIF.ratio.long1, aes(f.ordered1, hline, label = hline), vjust = 2, hjust = -0.2) +
        scale_color_manual(values = c("Tumor" = "#CC79A7", "NAT" = "#0072B2"),
                           breaks = c("Tumor", "NAT")) +
        # for color-blind palettes
        facet_wrap(~variable, scales = "free_x") +
        ggplot2::facet_grid(~variable, scales = "free_x", space = "free") +
        # ggplot2::facet_grid_sc(cols = vars(label), shrink = TRUE,
        #              scales = list(y = scales_y),
        #              space = "free") +
        guides(colour = guide_legend(nrow = 1)) +
        labs(x = "primary disease",
             y = "Ratio of RNA counts"
        ) +
        coord_flip(ylim = c(0, 25)) +
        theme_bw() +
        theme(
          plot.title = black_bold_12(),
          axis.title.x = black_bold_12(),
          axis.title.y = element_blank(),
          axis.text.x = black_bold_12(),
          axis.text.y = black_bold_12(),
          #axis.line.x = element_line(color = "black"),
          #axis.line.y = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.justification = "left",
          legend.box = "horizontal",
          legend.text = black_bold_12(),
          #strip.background = element_blank(),
          strip.text.x = black_bold_12()
        )
      print(p4)
      
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorratio4.pdf",
        plot = p4,
        width = 18,
        height = 8,
        useDingbats = FALSE)
    }
    make.plot1()
  }
  make.plot()
}
plot.boxgraph.EIF.ratio.TCGA.GTEX.2(c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
                                      "EIF4G1","EIF4G2","EIF4G3","EIF3D",
                                      "EIF4A1","EIF4A2"))
##################################################################
## violin plot for EIF expression in tumors vs adjacent normals ##
##################################################################
plot.violingraph.EIF.RNAseq.TCGA <- function(EIF.gene) {
  tissue.GTEX.TCGA.gene <- function() {
    TCGA.GTEX.anno <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names

    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    subset <- TCGA.GTEX.anno[, c("_sample_type", "_primary_site"),
      drop = FALSE
    ]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]

    # TCGA.GTEX <- as.data.frame(TCGA.GTEX)

    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE
    )
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)]
    TCGA.GTEX.Lung.sampletype <- na.omit(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  EIF.TCGA.RNAseq.anno <- tissue.GTEX.TCGA.gene()

  get.subset.data <- function() {
    EIF.TCGA.RNAseq.anno <- as.data.frame(EIF.TCGA.RNAseq.anno)
    row.names(EIF.TCGA.RNAseq.anno) <- EIF.TCGA.RNAseq.anno$Row.names
    EIF.TCGA.RNAseq.anno$Row.names <- NULL
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno[, c(
      EIF.gene,
      "sample.type",
      "primary.site"
    ),
    drop = FALSE
    ]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0,
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        # "Normal Tissue",
        "Solid Tissue Normal"
      ),
    ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- as.factor(as.character(
      EIF.TCGA.RNAseq.anno.subset$sample.type
    ))
    EIF.TCGA.RNAseq.anno.subset.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long$primary.site <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.site
    )
    return(EIF.TCGA.RNAseq.anno.subset.long)
  }
  EIF.TCGA.RNAseq.anno.subset.long <- get.subset.data()

  make.plot <- function() {
    EIF.TCGA.RNAseq.anno.subset.long1 <- EIF.TCGA.RNAseq.anno.subset.long[
      EIF.TCGA.RNAseq.anno.subset.long$variable %in% c(EIF.gene),]
    EIF.TCGA.RNAseq.anno.subset.long1$variable <- factor(
      EIF.TCGA.RNAseq.anno.subset.long1$variable,
      levels = c("EIF4G1", "EIF4G2","EIF4G3","EIF4A1","EIF4A2", "EIF4B","EIF4H",
                 "EIF4E", "EIF4E2", "EIF4E3", "EIF4EBP1" ,"EIF4EBP2", "EIF3D")
    )
    p1 <- ggplot(
      data = EIF.TCGA.RNAseq.anno.subset.long1,
      aes(
        x = sample.type,
        y = 2**value,
        color = sample.type,
        fill = sample.type
      )
    ) +
      stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
      ggplot2::facet_grid(. ~ variable,
        scales = "free",
        space  = "free"
      ) +
      # facet_wrap(~ variable, ncol = 6) +
      geom_violin(trim = TRUE) +
      # scale_color_manual(values = c("#D55E00","#0072B2","#CC79A7","#009E73")) + #for color-blind palettes
      geom_boxplot(
        alpha = .01,
        width = .25,
        color = "black",
        position = position_dodge(width = .9)
      ) +
      labs(
        x = "sample type",
        y = "normalized RNA counts"
      ) +
      scale_x_discrete(labels = c(
        "Metastatic Tumor",
        "Primary Tumor",
        # "Normal Tissue",
        "Adjacent Normal"
      )) +
      scale_y_continuous(
        trans = log2_trans(),
        labels = label_comma(),
        breaks = c(1, 128, 2048, 32768),
        # labels = c("1","8","64","512"),
        position = "left"
      ) +
      scale_color_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + # for color-blind palettes
      scale_fill_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + # for color-blind palettes
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title.x = element_blank(),
        axis.title.y.right = black_bold_16_right(),
        axis.text.x = black_bold_16_90(),
        axis.text.y = black_bold_16_90(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_16()
      ) +
      stat_compare_means(
        comparisons = list(
          c("Metastatic", "Solid Tissue Normal"),
          c("Primary Tumor", "Solid Tissue Normal"),
          c("Metastatic", "Primary Tumor")
        ),
        method = "t.test",
        label = "p.signif",
        size = 6
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "Expression"),
      filename = "EIFexpressionviolin.pdf",
      plot = p1,
      width = 18,
      height = 9,
      useDingbats = FALSE
    )
  }
  make.plot()
}
plot.violingraph.EIF.RNAseq.TCGA(c("EIF4A1", "EIF4A2", 
                                   "EIF4E", "EIF4E2", "EIF4E3", 
                                   "EIF3D", "EIF4EBP1", #"EIF4EBP2", 
                                   "EIF4G1", "EIF4G2", "EIF4G3",
                                   "EIF4H", "EIF4B"))
#############################################################
## violin plot for EIF ratio in tumors vs adjacent normals ##
#############################################################
plot.violingraph.EIF.ratio.TCGA <- function(EIF.gene) {

  tissue.GTEX.TCGA.gene <- function() {
    TCGA.GTEX.anno <- read_tsv(file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt"))
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names

    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names
    subset <- TCGA.GTEX.anno[, c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE
    )
    # check the name of the last column
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()

  get.EIFratio.anno.data <- function() {
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[, c(
      EIF.gene,
      "sample.type",
      "primary.site"
    ),
    drop = FALSE
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        "Solid Tissue Normal"
      ),
    ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- as.factor(as.character(
      EIF.TCGA.RNAseq.anno.subset$sample.type
    ))
    EIF.TCGA.RNAseq.anno.subset$sum <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score <- EIF.TCGA.RNAseq.anno.subset
    
    # A:E
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4A2:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A2 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    
    # G:E
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G3:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G3 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G3:\nEIF4E2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G3 - EIF.TCGA.RNAseq.anno.subset$EIF4E2)

    # G2:G1
    EIF.TCGA.GTEX.score$`EIF4G2:\nEIF4G1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G2 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4G3` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4G3)
    EIF.TCGA.GTEX.score$`EIF4G2:\nEIF3D` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G2 - EIF.TCGA.RNAseq.anno.subset$EIF3D)
    
    # E2:E1
    EIF.TCGA.GTEX.score$`EIF4E2:\nEIF4E` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4E2 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    
    # A2:A1
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4A2` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4A2)
    # EBP:E
    EIF.TCGA.GTEX.score$`EIF4E:\nEIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4E - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    # EIF.TCGA.GTEX.score$`PABPC1:\nEIF4E` <-
    #  (EIF.TCGA.RNAseq.anno.subset$PABPC1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E+EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E+EIF4EBP1` <-
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$sum)
    
    ratio <- c(
      "EIF4A1:\nEIF4E", "EIF4A2:\nEIF4E", 
      "EIF4G1:\nEIF4E", 
      "EIF4G3:\nEIF4E","EIF4G3:\nEIF4E2",
      "EIF4E:\nEIF4EBP1", 
      "EIF4G2:\nEIF4G1", "EIF4G1:\nEIF4G3","EIF4G2:\nEIF3D",
      "EIF4E2:\nEIF4E", "EIF4A1:\nEIF4A2",
      "EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1"
    )
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.score[,c(ratio, "sample.type", "primary.site")]
    EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score)
    return(EIF.TCGA.GTEX.score.long)
  }
  EIF.TCGA.GTEX.score.long <- get.EIFratio.anno.data()

  make.plot <- function() {
    new.label <- c(
      `EIF4A1:EIF4E` = "EIF4A1:\nEIF4E",
      `EIF4A2:EIF4E` = "EIF4A2:\nEIF4E", 
      `EIF4G1:EIF4E` = "EIF4G1:\nEIF4E",
      `EIF4G3:EIF4E` = "EIF4G3:\nEIF4E",
      `EIF4G3:EIF4E2` = "EIF4G3:\nEIF4E2",
      `EIF4E:EIF4EBP1` = "EIF4E:\nEIF4EBP1",
      `EIF4G2:EIF4G1` = "EIF4G2:\nEIF4G1",
      `EIF4G1:EIF4G3` = "EIF4G1:\nEIF4G3",
      `EIF4E2:EIF4E` = "EIF4E2:\nEIF4E",
      `EIF4A1:EIF4A2` = "EIF4A1:\nEIF4A2",
      `EIF4G1:EIF4E+EIF4EBP1` = "EIF4G1:\nEIF4E+EIF4EBP1",
      `EIF4A1:EIF4E+EIF4EBP1` = "EIF4A1:\nEIF4E+EIF4EBP1"
    )
    EIF.TCGA.GTEX.score.long1 <- EIF.TCGA.GTEX.score.long[
      EIF.TCGA.GTEX.score.long$variable %in% c(
        "EIF4A1:\nEIF4E", "EIF4A2:\nEIF4E", 
        "EIF4G1:\nEIF4E", 
        "EIF4G3:\nEIF4E","EIF4G3:\nEIF4E2",
        "EIF4E:\nEIF4EBP1", 
        "EIF4G2:\nEIF4G1", "EIF4G1:\nEIF4G3",
        "EIF4E2:\nEIF4E", "EIF4A1:\nEIF4A2",
        "EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1"
      ),
    ]
    EIF.TCGA.GTEX.score.long1$variable <- factor(
      EIF.TCGA.GTEX.score.long1$variable,
      levels = c(
        "EIF4G1:\nEIF4E","EIF4A1:\nEIF4E", "EIF4A2:\nEIF4E", 
        "EIF4G3:\nEIF4E","EIF4G3:\nEIF4E2","EIF4G1:\nEIF4G3",
        "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", "EIF4A1:\nEIF4A2",
        "EIF4E:\nEIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1"
      )
    )
    p1 <- ggplot(
      data = EIF.TCGA.GTEX.score.long1,
      aes(
        x = sample.type,
        y = 2**value,
        color = sample.type,
        fill = sample.type
      )
    ) +
      geom_violin(trim = FALSE) +
      geom_boxplot(
        alpha = 0.01,
        width = 0.25,
        color = "black",
        outlier.colour = NA
      ) +
      stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
      ggplot2::facet_grid(~variable,
        scales = "free",
        space  = "free"
      ) +
      #facet_wrap(~variable,
      #  labeller = labeller(variable = as_labeller(new.label))
      #) +
      scale_color_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + # for color-blind palettes
      scale_fill_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + # for color-blind palettes
      labs(
        x = "sample type",
        y = "ratio of RNA counts"
      ) +
      scale_x_discrete(labels = c(
        "Metastatic Tumor",
        "Primary Tumor",
        "Adjacent Normal"
      )) +
      scale_y_continuous(
        trans = log2_trans(),
        # labels = label_comma(),
        breaks = c(0.125, 1, 4, 8, 64, 512),
        labels = c("0.125", "1", "4", "8", "64", "512"),
        position = "left"
      ) +
      geom_hline(yintercept = c(1, 4), linetype = "dashed") +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title.x = element_blank(),
        axis.title.y = black_bold_16(),
        axis.text.x = black_bold_16_90(),
        axis.text.y = black_bold_16_90(),
        #axis.line.x = element_line(color = "black"),
        #axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_16()
      ) +
      stat_compare_means(
        comparisons = list(
          c("Metastatic", "Solid Tissue Normal"),
          c("Primary Tumor", "Solid Tissue Normal"),
          c("Metastatic", "Primary Tumor")
        ),
        method = "t.test",
        label = "p.signif",
        size = 6,
        hjust = 0
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "Expression"),
      filename = "EIFratioviolin.pdf",
      plot = p1,
      width = 18,
      height = 8,
      useDingbats = FALSE
    )

   
  }
  make.plot()
}
plot.violingraph.EIF.ratio.TCGA(c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
                                  "EIF4G1","EIF4G2","EIF4G3", "EIF3D",
                                  "EIF4A1","EIF4A2"))

## Figure 4 ##
#################################################################
##  PCA plots on EIF4F RNA-seq data from TCGA and GTEx groups  ##
#################################################################
plot.EIF.TCGA.PCA.all.tumor <- function(EIF.list) {
  tissue.GTEX.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX.anno <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno, drop = FALSE) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names
    subset <- TCGA.GTEX.anno[, c("_sample_type", "primary disease or tissue"),
                             drop = FALSE
    ]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    TCGA.GTEX <- as.data.frame(TCGA.GTEX, drop = FALSE)
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX <- as.data.frame(TCGA.GTEX, drop = FALSE)
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE
    )
    
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  
  get.EIF.TCGA.GTEX <- function(EIF.list) {
    # EIF.list <- c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","PABPC1","MKNK1","MKNK2", "MYC")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[, c(EIF.list, "sample.type", "primary.site"), drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0,
    ]
    # EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site %in% c("Brain","Blood","Pancreas"), ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        "Normal Tissue",
        "Solid Tissue Normal"
      ),
    ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c(
        "Normal Tissue",
        "Solid Tissue Normal",
        "Primary Tumor",
        "Metastatic"
      ),
      labels = c(
        "Healthy Tissue (GTEx)",
        "Adjacent Normal Tissue (TCGA)",
        "Primary Tumor (TCGA)",
        "Metastatic Tumor (TCGA)"
      )
    )
    EIF.TCGA.RNAseq.anno.subset$primary.site <- as.factor(EIF.TCGA.RNAseq.anno.subset$primary.site)
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)
  
  plot.primary.PCA <- function() {
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Primary Tumor (TCGA)"),
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.list, "sample.type", "primary.site")]
    EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
    ## remove the last two columns
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset) - 2)]
    rownames(df1) <- NULL
    plot.pca.factomineR <- function() {
      res.pca <- PCA(df1,
                     scale.unit = TRUE,
                     ncp = 10,
                     graph = FALSE
      )
      
      biplot <- fviz_pca_biplot(res.pca,
                                axes = c(1, 2),
                                labelsize = 5,
                                col.ind = EIF.TCGA.RNAseq.anno.subset$primary.site,
                                title = "PCA - Biplot (Primary tumors)",
                                palette = col_vector,
                                pointshape = 20,
                                pointsize = 0.75,
                                # addEllipses = TRUE,
                                label = "var",
                                col.var = "black",
                                repel = TRUE
      ) +
        xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16(),
          legend.title = element_blank(),
          legend.position = c(0, .625),
          legend.background = element_blank(),
          legend.text = black_bold_16()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = "EIFPCAprimary.pdf",
        plot = biplot,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      
      eig <- fviz_eig(res.pca,
                      labelsize = 6,
                      geom = "bar",
                      width = 0.7,
                      addlabels = TRUE
      ) +
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16()
        )
      print(eig)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = "EIFeigprimary.pdf",
        plot = eig,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      
      var <- get_pca_var(res.pca)
      pdf(file.path(
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = "EIFcorprimary.pdf"
      ),
      width = 9,
      height = 9,
      useDingbats = FALSE
      )
      corrplot(var$cos2, # cos2 is better than contribute
               title = "PCA (Primary tumors)",
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      dev.off()
      corrplot(var$cos2, # cos2 is better than contribute
               title = "PCA (Primary tumors)",
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      
      contribplot <- function(x) {
        fviz_contrib(res.pca,
                     choice = "var",
                     axes = x,
                     top = 10,
                     fill = "lightblue",
                     color = "black"
        ) +
          theme_minimal() +
          theme(
            plot.background = element_blank(),
            plot.title = black_bold_16(),
            panel.background = element_rect(
              fill = "transparent",
              color = "black",
              size = 1
            ),
            axis.title.x = element_blank(),
            axis.title.y = black_bold_16(),
            axis.text.x = black_bold_16_45(),
            axis.text.y = black_bold_16()
          )
      }
      lapply(c(1, 2), contribplot)
    }
    plot.pca.factomineR()
  }
  plot.primary.PCA()
  
  plot.metastatic.PCA <- function() {
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c( # "Primary Tumor (TCGA)",
        "Metastatic Tumor (TCGA)"
      ),
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.list, "sample.type", "primary.site")]
    EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
    ## remove the last two columns
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset) - 2)]
    rownames(df1) <- NULL
    plot.pca.factomineR <- function() {
      res.pca <- PCA(df1,
                     scale.unit = TRUE,
                     ncp = 10,
                     graph = FALSE
      )
      
      biplot <- fviz_pca_biplot(res.pca,
                                axes = c(1, 2),
                                labelsize = 5,
                                col.ind = EIF.TCGA.RNAseq.anno.subset$primary.site,
                                title = "PCA - Biplot (Metastatic tumors)",
                                palette = col_vector,
                                pointshape = 20,
                                pointsize = 0.75,
                                # addEllipses = TRUE,
                                label = "var",
                                col.var = "black",
                                repel = TRUE
      ) +
        xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16(),
          legend.title = element_blank(),
          legend.position = c(0, .625),
          legend.background = element_blank(),
          legend.text = black_bold_16()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = "EIFPCAmetastatic.pdf",
        plot = biplot,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      
      eig <- fviz_eig(res.pca,
                      labelsize = 6,
                      geom = "bar",
                      width = 0.7,
                      addlabels = TRUE
      ) +
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16()
        )
      print(eig)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = "EIFeigmetastatic.pdf",
        plot = eig,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      
      var <- get_pca_var(res.pca)
      pdf(file.path(
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = "EIFcormetastatic.pdf"
      ),
      width = 9,
      height = 9,
      useDingbats = FALSE
      )
      corrplot(var$cos2, # cos2 is better than contribute
               title = "PCA (Metastatic tumors)",
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      dev.off()
      corrplot(var$cos2, # cos2 is better than contribute
               title = "PCA (Metastatic tumors)",
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      
      contribplot <- function(x) {
        fviz_contrib(res.pca,
                     choice = "var",
                     axes = x,
                     top = 10,
                     fill = "lightblue",
                     color = "black"
        ) +
          theme_minimal() +
          theme(
            plot.background = element_blank(),
            plot.title = black_bold_16(),
            panel.background = element_rect(
              fill = "transparent",
              color = "black",
              size = 1
            ),
            axis.title.x = element_blank(),
            axis.title.y = black_bold_16(),
            axis.text.x = black_bold_16_45(),
            axis.text.y = black_bold_16()
          )
      }
      lapply(c(1, 2), contribplot)
    }
    plot.pca.factomineR()
  }
  plot.metastatic.PCA()
}
plot.EIF.TCGA.PCA.all.tumor(c(
  "EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))

plot.EIF.TCGA.PCA.all.tumor(c(
  "EIF4G1","EIF4G2", "EIF4G3",
  "EIF4A1","EIF4A2", 
  "EIF4E", "EIF4E2", "EIF4E3", 
  "EIF3D", "EIF4EBP1", 
  "EIF4H", "EIF4B", "MYC", "JUN"
))


plot.EIF.GTEX.PCA.all.tissue <- function(EIF.list) {
  tissue.GTEX.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX.anno <- read_tsv(file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt"))
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno)
    subset <- TCGA.GTEX.anno[, c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    TCGA.GTEX <- as.data.frame(TCGA.GTEX)
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE
    )
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  
  get.EIF.TCGA.GTEX <- function(EIF.list) {
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[, c(
      EIF.list,
      "sample.type",
      "primary.site"
    ),
    drop = FALSE
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0,
    ]
    # EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site %in% c("Brain","Blood","Pancreas"), ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        "Normal Tissue",
        "Solid Tissue Normal"
      ),
    ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c(
        "Normal Tissue",
        "Solid Tissue Normal",
        "Primary Tumor",
        "Metastatic"
      ),
      labels = c(
        "Healthy Tissue (GTEx)",
        "Adjacent Normal Tissue (TCGA)",
        "Primary Tumor (TCGA)",
        "Metastatic Tumor (TCGA)"
      )
    )
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)
  
  plot.PCA <- function() {
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type == "Healthy Tissue (GTEx)",
    ]
    EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
    ## remove the last two columns
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset) - 2)]
    rownames(df1) <- NULL
    
    plot.pca.factomineR <- function() {
      res.pca <- PCA(df1,
                     scale.unit = TRUE,
                     ncp = 10,
                     graph = FALSE
      )
      
      biplot <- fviz_pca_biplot(res.pca,
                                axes = c(1, 2),
                                labelsize = 5,
                                col.ind = EIF.TCGA.RNAseq.anno.subset$primary.site,
                                title = "PCA - Biplot (Healthy tissues)",
                                palette = col_vector, ,
                                pointshape = 20,
                                pointsize = 0.75,
                                # addEllipses = TRUE,
                                label = "var",
                                col.var = "black",
                                repel = TRUE
      ) +
        xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16(),
          legend.title = element_blank(),
          legend.position = c(0, .625),
          legend.background = element_blank(),
          legend.text = black_bold_16()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "GTEX"),
        filename = "EIFPCAGTEX.pdf",
        plot = biplot,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      
      eig <- fviz_eig(res.pca,
                      labelsize = 6,
                      geom = "bar",
                      width = 0.7,
                      addlabels = TRUE
      ) +
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16()
        )
      print(eig)
      
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "GTEX"),
        filename = "EIFeigGTEX.pdf",
        plot = eig,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      var <- get_pca_var(res.pca)
      # fviz_pca_var(res.pca, col.var="contrib")
      
      pdf(file.path(
        path = file.path(output.directory, "PCA", "GTEX"),
        filename = "EIFcorGTEX.pdf"
      ),
      width = 9,
      height = 9,
      useDingbats = FALSE
      )
      corrplot(var$cos2, # cos2 is better than contribute
               title = "PCA (Healthy tissues)",
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      dev.off()
      corrplot(var$cos2, # cos2 is better than contribute
               title = "PCA (Healthy tissues)",
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      
      contribplot <- function(x) {
        fviz_contrib(res.pca,
                     choice = "var",
                     axes = x,
                     top = 10,
                     fill = "lightblue",
                     color = "black"
        ) +
          theme_minimal() +
          theme(
            plot.background = element_blank(),
            plot.title = black_bold_16(),
            panel.background = element_rect(
              fill = "transparent",
              color = "black",
              size = 1
            ),
            axis.title.x = element_blank(),
            axis.title.y = black_bold_16(),
            axis.text.x = black_bold_16_45(),
            axis.text.y = black_bold_16()
          )
      }
      lapply(c(1, 2), contribplot)
    }
    plot.pca.factomineR()
  }
  plot.PCA()
}
plot.EIF.GTEX.PCA.all.tissue(c(
  "EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))


plot.EIF.TCGA.GTEX.PCA.all.tumor.tissue <- function(EIF.list) {
  tissue.GTEX.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX.anno <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names
    subset <- TCGA.GTEX.anno[, c("_sample_type", "primary disease or tissue", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.disease", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    TCGA.GTEX <- as.data.frame(TCGA.GTEX)
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE
    )
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)]
    # TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()

  get.EIF.TCGA.GTEX <- function(EIF.list) {
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[, c(
      EIF.list,
      "sample.type",
      "primary.disease",
      "primary.site"
    ),
    drop = FALSE
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0,
    ]
    # EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site == "Brain", ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        "Normal Tissue"
      ),
    ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c(
        "Normal Tissue",
        # "Solid Tissue Normal",
        "Primary Tumor",
        "Metastatic"
      ),
      labels = c(
        "Healthy Tissue (GTEx)",
        # "Adjacent Normal Tissue (TCGA)",
        "Primary Tumor (TCGA)",
        "Metastatic Tumor (TCGA)"
      )
    )
    # EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)

  plot.PCA.prcomp <- function() {
    # the variables should be scaled to have unit variance
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset) - 3)]
    PCA <- prcomp(df1, center = TRUE, scale = TRUE)
    # Extract PC axes for plotting
    PCAvalues <- data.frame(
      Sample.type = EIF.TCGA.RNAseq.anno.subset$sample.type,
      PCA$x
    )
    # Extract loadings of the variables
    PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
    # Plot
    gr <- factor(EIF.TCGA.RNAseq.anno.subset$sample.type)
    pca3d(PCA, group = gr)

    plot.PCA.3D <- function() {
      get_colors <- function(groups, group.col = palette()) {
        groups <- as.factor(groups)
        ngrps <- length(levels(groups))
        if (ngrps > length(group.col)) {
          group.col <- rep(group.col, ngrps)
        }
        color <- group.col[as.numeric(groups)]
        names(color) <- as.vector(groups)
        return(color)
      }
      cols <- get_colors(PCAvalues$Sample.type, brewer.pal(n = 4, name = "Dark2"))
      plot3d(PCAvalues[, 2:4], col = cols)
      legend3d("topright",
        legend = levels(PCAvalues$Sample.type),
        pch = 16, col = brewer.pal(n = 4, name = "Dark2"),
        cex = 1
      )
    }
    plot.PCA.3D()
    p <- ggplot(
      PCAvalues,
      aes(
        x = PC1,
        y = PC2,
        colour = Sample.type
      )
    ) +
      # xlim(-5, 8) +
      geom_point(size = 0.2) +
      geom_segment(
        data = PCAloadings,
        aes(
          x = 0,
          y = 0,
          xend = (PC1 * 5),
          yend = (PC2 * 5)
        ),
        arrow = arrow(length = unit(1 / 3, "picas")),
        color = "black"
      ) +
      annotate("text",
        size = 6,
        fontface = "bold",
        x = (PCAloadings$PC1 * 5),
        y = (PCAloadings$PC2 * 5),
        label = PCAloadings$Variables
      ) +
      # stat_n_text(geom = "label") +
      ggtitle("Principal Component Analysis") +
      theme(
        plot.background = element_blank(),
        plot.title = black_bold_16(),
        panel.background = element_rect(
          fill = "transparent",
          color = "black",
          size = 1
        ),
        axis.title = black_bold_16(),
        axis.text = black_bold_16(),
        legend.title = element_blank(),
        legend.position = c(0.3, 0.93),
        legend.background = element_blank(),
        legend.text = black_bold_16(),
        legend.key = element_blank()
      )
    print(p)
  }
  plot.PCA.prcomp()

  plot.pca.factomineR <- function() {
    df1 <- EIF.TCGA.RNAseq.anno.subset
    res.pca <- PCA(df1[1:(length(df1) - 3)],
      scale.unit = TRUE,
      ncp = 14,
      graph = FALSE
    )
    biplot <- fviz_pca_biplot(res.pca,
      axes = c(1, 2),
      labelsize = 5,
      col.ind = df1$sample.type,
      palette = c("#D55E00", "#009E73", "#CC79A7", "#0072B2"),
      # palette    = c("#CC79A7","#0072B2"),
      pointshape = 20,
      pointsize = 0.75,
      # addEllipses = TRUE,
      title = "PCA - Biplot (Healthy Tissues + Tumors)",
      label = "var",
      col.var = "black",
      repel = TRUE
    ) +
      xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
      # xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
      theme_classic() +
      theme(
        plot.background = element_blank(),
        plot.title = black_bold_16(),
        panel.background = element_rect(
          fill = "transparent",
          color = "black",
          size = 1
        ),
        axis.title.x = black_bold_16(),
        axis.title.y = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        legend.title = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.text = black_bold_16()
      )
    print(biplot)
    ggplot2::ggsave(
      path = file.path(output.directory, "PCA", "All"),
      filename = "EIFPCAall.pdf",
      plot = biplot,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )

    plot.selected.PCA <- function(sample, color) {
      test <- df1[df1$sample.type == sample, ]
      sample.type <- levels(df1$sample.type)
      biplot <- fviz_pca_biplot(res.pca,
        axes = c(1, 2),
        labelsize = 5,
        col.ind = df1$sample.type,
        # palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"),
        palette = color,
        select.ind = list(name = row.names(test)),
        pointshape = 20,
        pointsize = 0.75,
        # addEllipses = TRUE,
        title = "PCA - Biplot (Healthy Tissues + Tumors)",
        label = "var",
        col.var = "black",
        repel = TRUE
      ) +
        xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
        # xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
        theme_classic() +
        # scale_alpha_manual(values=c(0.1, 0.1, 0.1, 0.1),guide=F)+
        # scale_x_continuous(breaks = seq(-6, 8, 2), limits=c(-5, 8)) +
        # scale_y_continuous(breaks = seq(-4, 6, 2), limits=c(-4, 7)) +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16(),
          legend.title = element_blank(),
          legend.position = c(0, 0),
          legend.justification = c(0, 0),
          legend.background = element_blank(),
          legend.text = black_bold_16()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "All"),
        filename = paste0("EIFPCAall", sample, ".pdf"),
        plot = biplot,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    plot.selected.PCA("Healthy Tissue (GTEx)", "#D55E00")
    plot.selected.PCA("Primary Tumor (TCGA)", "#009E73")
    plot.selected.PCA("Metastatic Tumor (TCGA)", "#CC79A7")

    plot.selected.healthy.color.PCA <- function(sample, color) {
      test <- df1[df1$sample.type == sample, ]
      sample.type <- levels(df1$sample.type)
      biplot <- fviz_pca_biplot(res.pca,
        axes = c(1, 2),
        labelsize = 5,
        col.ind = df1$primary.site,
        # palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"),
        palette = color,
        select.ind = list(name = row.names(test)),
        pointshape = 20,
        pointsize = 0.75,
        # addEllipses = TRUE,
        title = "PCA - Biplot (Healthy Tissues + Tumors)",
        label = "var",
        col.var = "black",
        repel = TRUE
      ) +
        xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
        # xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
        theme_classic() +
        # scale_alpha_manual(values=c(0.1, 0.1, 0.1, 0.1),guide=F)+
        # scale_x_continuous(breaks = seq(-6, 8, 2), limits=c(-5, 8)) +
        # scale_y_continuous(breaks = seq(-4, 6, 2), limits=c(-4, 7)) +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16(),
          legend.title = element_blank(),
          legend.position = c(0, .625),
          legend.justification = c(0, 0),
          legend.background = element_blank(),
          legend.text = black_bold_12()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "All"),
        filename = paste0("EIFPCAall", sample, "color.pdf"),
        plot = biplot,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    plot.selected.healthy.color.PCA("Healthy Tissue (GTEx)", col_vector)

    plot.selected.tumor.color.PCA <- function(sample, color) {
      test <- df1[df1$sample.type == sample, ]
      # sample.type <- levels(df1$sample.type)
      test$primary.disease <- as.factor(test$primary.disease)
      biplot <- fviz_pca_biplot(res.pca,
        axes = c(1, 2),
        labelsize = 5,
        col.ind = df1$primary.disease,
        # palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"),
        palette = color,
        select.ind = list(name = row.names(test)),
        pointshape = 20,
        pointsize = 0.75,
        # addEllipses = TRUE,
        title = "PCA - Biplot (Healthy Tissues + Tumors)",
        label = "var",
        col.var = "black",
        repel = TRUE
      ) +
        xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
        # xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
        theme_classic() +
        # scale_alpha_manual(values=c(0.1, 0.1, 0.1, 0.1),guide=F)+
        # scale_x_continuous(breaks = seq(-6, 8, 2), limits=c(-5, 8)) +
        # scale_y_continuous(breaks = seq(-4, 6, 2), limits=c(-4, 7)) +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16(),
          legend.title = element_blank(),
          legend.position = c(0, .625),
          legend.justification = c(0, 0),
          legend.background = element_blank(),
          legend.text = black_bold_12()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "All"),
        filename = paste0("EIFPCAall", sample, "color.pdf"),
        plot = biplot,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    plot.selected.tumor.color.PCA("Primary Tumor (TCGA)", col_vector)
    plot.selected.tumor.color.PCA("Metastatic Tumor (TCGA)", col_vector)

    eig <- fviz_eig(res.pca,
      labelsize = 6,
      geom = "bar",
      width = 0.7,
      addlabels = TRUE
    ) +
      # geom_text(aes(label = res.pca$eig, size = 18)) +
      theme_classic() +
      theme(
        plot.background = element_blank(),
        plot.title = black_bold_16(),
        panel.background = element_rect(
          fill = "transparent",
          color = "black",
          size = 1
        ),
        axis.title.x = black_bold_16(),
        axis.title.y = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16()
      )
    print(eig)
    ggplot2::ggsave(
      path = file.path(output.directory, "PCA", "All"),
      filename = "EIFPCAeig.pdf",
      plot = eig,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )

    var <- get_pca_var(res.pca)
    corrplot(var$contrib, is.corr = FALSE)

    pdf(file.path(
      path = file.path(output.directory, "PCA", "All"),
      filename = "EIFPCAcor.pdf"
    ),
    width = 9,
    height = 9,
    useDingbats = FALSE
    )
    corrplot(var$cos2, # cos2 is better than contribute
      is.corr = FALSE,
      tl.cex = 1.5,
      number.cex = 1.5,
      method = "color",
      addgrid.col = "gray",
      addCoef.col = "black",
      tl.col = "black"
    )
    dev.off()
    corrplot(var$cos2, # cos2 is better than contribute
      is.corr = FALSE,
      tl.cex = 1.5,
      number.cex = 1.5,
      method = "color",
      addgrid.col = "gray",
      addCoef.col = "black",
      tl.col = "black"
    )

    contribplot <- function(x) {
      fviz_contrib(res.pca,
        choice = "var",
        axes = x,
        top = 10,
        fill = "lightblue",
        color = "black"
      ) +
        theme_minimal() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = element_blank(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16_45(),
          axis.text.y = black_bold_16()
        )
    }
    lapply(c(1, 2), contribplot)
  }
  plot.pca.factomineR()
}
plot.EIF.TCGA.GTEX.PCA.all.tumor.tissue(c(
  "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))


plot.EIF.TCGA.GTEX.PCA.all.tumor.tissue(c(
  "EIF4G1","EIF4G2", "EIF4G3",
  "EIF4A1","EIF4A2", 
  "EIF4E", "EIF4E2", "EIF4E3", 
  "EIF3D", "EIF4EBP1", 
  "EIF4H", "EIF4B", "MYC", "JUN"
))


plot.EIF.TCGA.GTEX.PCA.all.tumor.tissue2 <- function(EIF.list) {
  tissue.GTEX.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX.anno <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names
    subset <- TCGA.GTEX.anno[, c("_sample_type", "primary disease or tissue", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.disease", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    TCGA.GTEX <- as.data.frame(TCGA.GTEX)
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
      ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE
    )
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)]
    # TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  
  get.EIF.TCGA.GTEX <- function(EIF.list) {
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[, c(
      EIF.list,
      "sample.type",
      "primary.disease",
      "primary.site"
    ),
    drop = FALSE
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0,
      ]
    # EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site == "Brain", ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        "Normal Tissue"
      ),
      ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c(
        "Normal Tissue",
        # "Solid Tissue Normal",
        "Primary Tumor",
        "Metastatic"
      ),
      labels = c(
        "Healthy Tissue (GTEx)",
        # "Adjacent Normal Tissue (TCGA)",
        "Primary Tumor (TCGA)",
        "Metastatic Tumor (TCGA)"
      )
    )
    # EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)
  
  plot.pca.factomineR <- function() {
    df1 <- EIF.TCGA.RNAseq.anno.subset
    res.pca <- PCA(df1[1:(length(df1) - 3)],
                   scale.unit = TRUE,
                   ncp = 14,
                   graph = FALSE
    )
    var <- fviz_pca_var(res.pca)
    print(var)
    
    biplo <- fviz_pca_var(res.pca)
    print(biplo)
    
    biplot <- fviz_pca_biplot(res.pca,
                              axes = c(1, 2),
                              labelsize = 5,
                              col.ind = df1$sample.type,
                              palette = c("#D55E00", "#009E73", "#CC79A7", "#0072B2"),
                              # palette    = c("#CC79A7","#0072B2"),
                              pointshape = 20,
                              pointsize = 0.75,
                              # addEllipses = TRUE,
                              title = "PCA - Biplot (Healthy Tissues + Tumors)",
                              label = "var",
                              col.var = "black",
                              repel = TRUE
    ) +
      #xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
      # xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
      theme_classic() +
      theme(
        plot.background = element_blank(),
        plot.title = black_bold_16(),
        panel.background = element_rect(
          fill = "transparent",
          color = "black",
          size = 1
        ),
        axis.title.x = black_bold_16(),
        axis.title.y = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        legend.title = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.text = black_bold_16()
      )
    print(biplot)
    ggplot2::ggsave(
      path = file.path(output.directory, "PCA", "All"),
      filename = "EIFPCAall.pdf",
      plot = biplot,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )
    var <- get_pca_var(res.pca)
    corrplot(var$contrib, is.corr = FALSE)
    
    pdf(file.path(
      path = file.path(output.directory, "PCA", "All"),
      filename = "EIFPCAcor.pdf"
    ),
    width = 9,
    height = 9,
    useDingbats = FALSE
    )
    corrplot(var$cos2, # cos2 is better than contribute
             is.corr = FALSE,
             tl.cex = 1.5,
             number.cex = 1.5,
             method = "color",
             addgrid.col = "gray",
             addCoef.col = "black",
             tl.col = "black"
    )
    dev.off()
    corrplot(var$cos2, # cos2 is better than contribute
             is.corr = FALSE,
             tl.cex = 1.5,
             number.cex = 1.5,
             method = "color",
             addgrid.col = "gray",
             addCoef.col = "black",
             tl.col = "black"
    )
  }
  plot.pca.factomineR()
}
plot.EIF.TCGA.GTEX.PCA.all.tumor.tissue2(c(
  "EIF4G1",#"EIF4G2", "EIF4G2",
  "EIF4A1",#"EIF4A2", 
  "EIF4E", #"EIF4E2", #"EIF4E3", 
  #"EIF3D", 
  "EIF4EBP1", "PABPC1", "MKNK1", "MKNK2",
  "EIF4B", "EIF4H", 
  "MYC", "JUN"
))

plot.EIF.TCGA.GTEX.PCA.each.tumor <- function(EIF.list, tissue) {
  tissue.GTEX.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX.anno <- read_tsv(
      file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt")
    )
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    # TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names
    subset <- TCGA.GTEX.anno[, c("_sample_type", "_primary_site", "primary disease or tissue"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site", "primary.disease")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    TCGA.GTEX <- as.data.frame(TCGA.GTEX)
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE
    )
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)]
    # TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  get.EIF.TCGA.GTEX <- function(EIF.list) {
    # EIF.list <- c("EIF4E", "EIF4G1", "EIF4G2", "EIF4A1","EIF4EBP1", "PABPC1",
    #              "MKNK1","MKNK2", "MTOR", "RPTOR", "RPS6KB1","MYC")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[, c(
      EIF.list,
      "sample.type",
      "primary.site", "primary.disease"
    ),
    drop = FALSE
    ]

    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[c(EIF.list, "sample.type", "primary.site", "primary.disease")]

    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0,
    ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        "Normal Tissue",
        "Solid Tissue Normal"
      ),
    ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c(
        "Normal Tissue",
        "Solid Tissue Normal",
        "Primary Tumor",
        "Metastatic"
      ),
      labels = c(
        "Healthy Tissue (GTEx)",
        "Adjacent Normal Tissue (TCGA)",
        "Primary Tumor (TCGA)",
        "Metastatic Tumor (TCGA)"
      )
    )
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset$primary.disease <- as.factor(EIF.TCGA.RNAseq.anno.subset$primary.disease)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno <- get.EIF.TCGA.GTEX(EIF.list)

  EIF.PCA.tissue <- function(x) {
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno[
      EIF.TCGA.RNAseq.anno$primary.site == x,
    ]
    print(summary(EIF.TCGA.RNAseq.anno.subset))
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset) - 3)]
    rownames(df1) <- NULL
    plot.pca.factomineR <- function(x) {
      res.pca <- PCA(df1,
        scale.unit = TRUE,
        ncp = 10,
        graph = FALSE
      )

      biplot <- fviz_pca_biplot(res.pca,
        axes = c(1, 2),
        labelsize = 5,
        col.ind = EIF.TCGA.RNAseq.anno.subset$sample.type,
        palette = c("#D55E00", "#CC79A7", "#009E73", "#0072B2"),
        pointshape = 20,
        pointsize = 0.75,
        # addEllipses = TRUE, ellipse.level = 0.9,
        label = "var",
        col.var = "black",
        repel = TRUE,
        title = paste0("PCA - Biplot (", x, ")")
      ) +
        # scale_x_continuous(breaks = seq(-6, 12, 2), limits=c(-5, 12)) +
        # scale_y_continuous(breaks = seq(-4, 10, 2), limits=c(-4, 10)) +
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16(),
          legend.title = element_blank(),
          legend.position = c(0, 0),
          legend.justification = c(0, 0),
          legend.background = element_blank(),
          legend.text = black_bold_16()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "Lung"),
        filename = paste0(x, "EIFPCA.pdf"),
        plot = biplot,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )

      eig <- fviz_eig(res.pca,
        labelsize = 6,
        geom = "bar",
        width = 0.7,
        addlabels = TRUE
      ) +
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = black_bold_16(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16(),
          axis.text.y = black_bold_16()
        )
      print(eig)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "Lung"),
        filename = paste0(x, "EIFeig.pdf"),
        plot = eig,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )

      contribplot <- function(x) {
        p <- fviz_contrib(res.pca,
          choice = "var",
          axes = x,
          top = 10,
          fill = "lightblue",
          color = "black"
        ) +
          theme_classic() +
          theme(
            plot.background = element_blank(),
            plot.title = black_bold_16(),
            panel.background = element_rect(
              fill = "transparent",
              color = "black",
              size = 1
            ),
            axis.title.x = element_blank(),
            axis.title.y = black_bold_16(),
            axis.text.x = black_bold_16_45(),
            axis.text.y = black_bold_16()
          )
        print(p)
        ggplot2::ggsave(
          path = file.path(output.directory, "PCA", "Lung"),
          filename = paste0("EIFcontri", x, ".pdf"),
          plot = p,
          width = 8,
          height = 4,
          useDingbats = FALSE
        )
      }
      lapply(c(1, 2), contribplot)

      var <- get_pca_var(res.pca)
      # corrplot(var$contrib, is.corr=FALSE)
      # fviz_pca_var(res.pca, col.var="contrib")
      pdf(file.path(
        path = file.path(output.directory, "PCA", "Lung"),
        filename = paste0(x, "EIFPCAcor.pdf")
      ),
      width = 9,
      height = 9,
      useDingbats = FALSE
      )
      corrplot(var$cos2, # cos2 is better than contribute
        is.corr = FALSE,
        tl.cex = 1.5,
        number.cex = 1.5,
        method = "color",
        addgrid.col = "gray",
        addCoef.col = "black",
        tl.col = "black"
      )
      dev.off()
    }
    plot.pca.factomineR(x)
  }
  # lapply(disease.list, EIF.PCA.tissue)
  EIF.PCA.tissue(tissue)
}
plot.EIF.TCGA.GTEX.PCA.each.tumor(
  c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
    "PABPC1", "MKNK1", "MKNK2"#, "MYC"
    ),
  "Skin"
)

plot.EIF.CPTAC.PCA.LUAD <- function(EIF.list) {
  CPTAC.LUAD.Sample <- read_excel(
    file.path(data.file.directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx")
  )
  CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample[, c("Aliquot (Specimen Label)", "Type")]
  CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample.ID[
    !duplicated(CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`),
  ]
  row.names(CPTAC.LUAD.Sample.ID) <- CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`
  CPTAC.LUAD.Sample.ID <- as.data.frame(CPTAC.LUAD.Sample.ID)
  CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)` <- NULL

  CPTAC.LUAD.Proteomics <- fread(
    file.path(data.file.directory, "CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv"),
    data.table = FALSE
  )
  EIF.CPTAC.LUAD.Proteomics <- CPTAC.LUAD.Proteomics[CPTAC.LUAD.Proteomics$Gene %in% EIF.list, ]
  row.names(EIF.CPTAC.LUAD.Proteomics) <- EIF.CPTAC.LUAD.Proteomics$Gene
  EIF.CPTAC.LUAD.Proteomics <- select(EIF.CPTAC.LUAD.Proteomics, -contains("Unshared"))

  EIF.CPTAC.LUAD.Proteomics$Gene <- NULL
  EIF.CPTAC.LUAD.Proteomics <- EIF.CPTAC.LUAD.Proteomics[1:(length(EIF.CPTAC.LUAD.Proteomics) - 6)]
  EIF.CPTAC.LUAD.Proteomics.t <- data.table::transpose(EIF.CPTAC.LUAD.Proteomics)
  rownames(EIF.CPTAC.LUAD.Proteomics.t) <- colnames(EIF.CPTAC.LUAD.Proteomics)
  colnames(EIF.CPTAC.LUAD.Proteomics.t) <- rownames(EIF.CPTAC.LUAD.Proteomics)
  rownames(EIF.CPTAC.LUAD.Proteomics.t) <- sub(" Log Ratio", "", rownames(EIF.CPTAC.LUAD.Proteomics.t))
  EIF.CPTAC.LUAD.Proteomics.Sampletype <- merge(EIF.CPTAC.LUAD.Proteomics.t,
    CPTAC.LUAD.Sample.ID,
    by    = "row.names",
    all.x = TRUE
  )
  rownames(EIF.CPTAC.LUAD.Proteomics.Sampletype) <- EIF.CPTAC.LUAD.Proteomics.Sampletype$Row.names
  EIF.CPTAC.LUAD.Proteomics.Sampletype$Row.names <- NULL
  EIF.CPTAC.LUAD.Proteomics.Sampletype$Type <- factor(EIF.CPTAC.LUAD.Proteomics.Sampletype$Type,
    levels = c("Normal", "Tumor"),
    labels = c("Adjacent Normal Tissue (CPTAC)", "Primary Tumor (CPTAC)")
  )
  EIF.CPTAC.LUAD.Proteomics.Sampletype <- EIF.CPTAC.LUAD.Proteomics.Sampletype[!is.na(EIF.CPTAC.LUAD.Proteomics.Sampletype$Type), ]
  EIF.CPTAC.LUAD.Proteomics.Sampletype <- EIF.CPTAC.LUAD.Proteomics.Sampletype[,c(EIF.list, "Type")]


  df1 <- EIF.CPTAC.LUAD.Proteomics.Sampletype[1:(length(EIF.CPTAC.LUAD.Proteomics.Sampletype) - 1)]
  rownames(df1) <- NULL

  # TODO: 'nb' is defined here but not used
  nb <- missMDA::estim_ncpPCA(df1, ncp.max = 5)
  res.comp <- missMDA::imputePCA(df1, ncp = 2)

  # nb <- missMDA::estim_ncpMCA(df1, ncp.max = 5)
  # res.comp <- missMDA::MIPCA(df1,ncp = nb$ncp, nboot = 1000)

  res.pca <- PCA(res.comp$completeObs,
    scale.unit = TRUE,
    ncp = 10,
    graph = FALSE
  )


  biplot <- fviz_pca_biplot(res.pca,
    axes = c(1, 2),
    labelsize = 5,
    col.ind = EIF.CPTAC.LUAD.Proteomics.Sampletype$Type,
    palette = c("#D55E00", "#009E73"),
    pointshape = 20,
    pointsize = 0.75,
    title = "PCA - Biplot (LUAD)",
    label = "var",
    col.var = "black",
    repel = TRUE
  ) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      plot.title = black_bold_16(),
      panel.background = element_rect(
        fill = "transparent",
        color = "black",
        size = 1
      ),
      axis.title.x = black_bold_16(),
      axis.title.y = black_bold_16(),
      axis.text.x = black_bold_16(),
      axis.text.y = black_bold_16(),
      legend.title = element_blank(),
      legend.position = c(0, 0),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.text = black_bold_16()
    )
  print(biplot)
  ggplot2::ggsave(
    path = file.path(output.directory, "PCA", "Lung"),
    filename = "EIFLUADPCA.pdf",
    plot = biplot,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
  eig <- fviz_eig(res.pca,
    labelsize = 6,
    geom = "bar",
    width = 0.7,
    addlabels = TRUE
  ) +
    # geom_text(aes(label = res.pca$eig, size = 18)) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      plot.title = black_bold_16(),
      panel.background = element_rect(
        fill = "transparent",
        color = "black",
        size = 1
      ),
      axis.title.x = black_bold_16(),
      axis.title.y = black_bold_16(),
      axis.text.x = black_bold_16(),
      axis.text.y = black_bold_16()
    )
  print(eig)
  ggplot2::ggsave(
    path = file.path(output.directory, "PCA", "Lung"),
    filename = "EIFLUADEig.pdf",
    plot = eig,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
  var <- get_pca_var(res.pca)
  # fviz_pca_var(res.pca, col.var="contrib")
  pdf(file.path(
    path = file.path(output.directory, "PCA", "Lung"),
    filename = "EIFLUADcor.pdf"
  ),
  width = 9,
  height = 9,
  useDingbats = FALSE
  )
  corrplot(var$cos2, # cos2 is better than contribute
    is.corr = FALSE,
    tl.cex = 1.5,
    number.cex = 1.5,
    method = "color",
    addgrid.col = "gray",
    addCoef.col = "black",
    tl.col = "black"
  )
  dev.off()
  corrplot(var$cos2, # cos2 is better than contribute
    is.corr = FALSE,
    tl.cex = 1.5,
    number.cex = 1.5,
    method = "color",
    addgrid.col = "gray",
    addCoef.col = "black",
    tl.col = "black"
  )
}
plot.EIF.CPTAC.PCA.LUAD(c("EIF4E", "EIF4G1", "EIF4A1", "PABPC1", 
                          "MKNK1", "MKNK2", "EIF4EBP1"))


## Figure 5 ##
######################################
#### RNA protein correlation CCLE ####
######################################

plot.EIF.cor.CCLE <- function() {
  get_EIF_CCLE_RNA <- function(){
    CCLE_RNA <- fread(
      file.path(data.file.directory, "CCLE_expression_full.csv"),
      data.table = FALSE
    )
    EIF_CCLE_RNA <- CCLE_RNA[c('EIF4G1 (ENSG00000114867)', 
                               'EIF4A1 (ENSG00000161960)',
                               'EIF4E (ENSG00000151247)',
                               'EIF4EBP1 (ENSG00000187840)', "V1")]
    
    CCLE_Anno <- fread(
      file.path(data.file.directory, "sample_info.csv"),
      data.table = FALSE
    )
    CCLE_Anno <- CCLE_Anno[, 1:2]
    CCLE_RNA_Anno <- merge(EIF_CCLE_RNA, CCLE_Anno, by.x = "V1", by.y = "DepMap_ID", all=T)
    CCLE_RNA_Anno <- na.omit(CCLE_RNA_Anno)
    CCLE_RNA_Anno$V1 <- NULL
    colnames(CCLE_RNA_Anno) <- c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "Celline")
    return(CCLE_RNA_Anno)
    
  }
  EIF_CCLE_RNA <- get_EIF_CCLE_RNA()
  
  
  get_EIF_CCLE_Pro <- function(){
    CCLE_PRO <- fread(
      file.path(data.file.directory, "protein_quant_current_normalized.csv"),
      data.table = FALSE
    )
    EIF_CCLE_PRO <- CCLE_PRO[CCLE_PRO$Gene_Symbol %in% c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1"), ]
    EIF_CCLE_PRO <- EIF_CCLE_PRO[, -grep("_Peptides", colnames(EIF_CCLE_PRO))]
    EIF_CCLE_PRO <-  EIF_CCLE_PRO[-2, , drop = FALSE]
    row.names(EIF_CCLE_PRO) <- EIF_CCLE_PRO$Gene_Symbol
    EIF_CCLE_PRO$Gene_Symbol <- NULL
    EIF_CCLE_PRO <-  EIF_CCLE_PRO[ , -(1:5) , drop = FALSE]
    EIF_CCLE_PRO.T <- t(EIF_CCLE_PRO)
    EIF_CCLE_PRO.T <- as.data.frame(EIF_CCLE_PRO.T)
    EIF_CCLE_PRO.T$Celline <- sub("\\_.*", "", row.names(EIF_CCLE_PRO.T))
    EIF_CCLE_PRO.T$Type <- sub(".*_ *(.*?) *_.*", "\\1", row.names(EIF_CCLE_PRO.T))
    EIF_CCLE_PRO.T$Type <- as.factor(EIF_CCLE_PRO.T$Type)
    return(EIF_CCLE_PRO.T)
  }
  EIF_CCLE_Pro <- get_EIF_CCLE_Pro()
  
  CCLE_RNA_Pro <- merge(EIF_CCLE_RNA, EIF_CCLE_Pro, by = "Celline",  suffixes = c("RNA","Pro"),all=T)
  CCLE_RNA_Pro <- na.omit(CCLE_RNA_Pro)
  
  Scatter.plot <- function(x){
    p1 <- ggscatter(CCLE_RNA_Pro, 
                    x = paste0(x, "Pro"), 
                    y = paste0(x, "RNA"), #color = "Type", 
                    add = "reg.line", #conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson", title = x,
                    xlab = "Protein expresion", ylab = "RNA expression")+
      stat_cor(aes(color = Type),     
               label.x.npc = 0.5, 
               label.y.npc = 0.5, hjust = 0
      )+           # Add correlation coefficient
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) 
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "LUAD"),
      filename = paste(x, "corCCLE.pdf"),
      plot = p1,
      width = 3,
      height = 3,
      useDingbats = FALSE
    )
    
    p2 <- ggscatter(CCLE_RNA_Pro, 
                    x = paste0(x, "Pro"), 
                    y = paste0(x, "RNA"), color = "Type", 
                    add = "reg.line", #conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson", #title = x,
                    #xlab = "Protein expresion", 
                    #ylab = "RNA expression"
    )+
      stat_cor(aes(color = Type),     
               label.x.npc = 0.5, 
               label.y.npc = 0.5, hjust = 0
      )+           # Add correlation coefficient
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        #legend.position = c(0.8, 0.2),
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) 
    print(p2)
    ggplot2::ggsave(
      path = file.path(output.directory, "LUAD"),
      filename = paste(x, "pro corCCLE.pdf"),
      plot = p2,
      width = 6,
      height = 3,
      useDingbats = FALSE
    )
  }
  lapply(c("EIF4G1", "EIF4A1", "EIF4E"), Scatter.plot)
}

plot.EIF.cor.CCLE()
  
#####################################
## Heatmap of correlation analysis ##
#####################################

plot.bargraph.CORs <- function(
  EIF.cor.tumor,
  EIF.cor.normal,
  tumor.label,
  normal.label,
  gene,
  posCORs,
  negCORs,
  output.poscor.filename,
  output.negcor.filename,
  coord_flip.ylim) {
  EIF.cor.tumor$label <- tumor.label
  EIF.cor.tumor$gene <- row.names(EIF.cor.tumor)
  EIF.cor.normal$label <- normal.label
  EIF.cor.normal$gene <- row.names(EIF.cor.normal)
  EIF.cor <- rbind(EIF.cor.tumor, EIF.cor.normal, make.row.names = F)
  
  EIF.cor$label <- factor(EIF.cor$label, levels = c(tumor.label, normal.label))
  EIF.cor$gene <- factor(EIF.cor$gene, levels = c("EIF4EBP1", "EIF4A1", "EIF4G1", "EIF4E"))
  levels(EIF.cor$gene)
  
  p1 <- ggplot(
    data = EIF.cor,
    aes(
      x = gene,
      y = posCORs,
      fill = label
    ),
    color = label
  ) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = posCORs),
              position = position_dodge(width = 0.9),
              size = 3.5
    ) +
    scale_fill_manual(values = c("#CC79A7", "#0072B2", "#E69F00", "#009E73", "#D55E00")) + # for color-blind palettes
    labs(y = paste("number of positively correlating genes")) +
    coord_flip(ylim = c(0, coord_flip.ylim)) +
    guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
    theme_bw() +
    theme(
      plot.title = black_bold_18(),
      axis.title.x = black_bold_18(),
      axis.title.y = element_blank(),
      axis.text.x = black_bold_18(),
      axis.text.y = black_bold_18(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_18(),
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_18()
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output.directory, "Heatmap"),
    filename = paste(output.poscor.filename),
    plot = p1,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
  
  p2 <- ggplot(
    data = EIF.cor,
    aes(
      x = gene,
      y = negCORs,
      fill = label
    ),
    color = label
  ) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = negCORs),
              position = position_dodge(width = 0.9),
              size = 3.5
    ) +
    scale_fill_manual(values = c("#CC79A7", "#0072B2", "#E69F00", "#009E73", "#D55E00")) + # for color-blind palettes
    labs(y = paste("number of negatively correlating genes")) +
    coord_flip(ylim = c(0, coord_flip.ylim)) +
    guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
    theme_bw() +
    theme(
      plot.title = black_bold_18(),
      axis.title.x = black_bold_18(),
      axis.title.y = element_blank(),
      axis.text.x = black_bold_18(),
      axis.text.y = black_bold_18(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_18(),
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_18()
    )
  print(p2)
  ggplot2::ggsave(
    path = file.path(output.directory, "Heatmap"),
    filename = output.negcor.filename,
    plot = p2,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
}


### find posCOR and negCOR in the overlapping CORs from all cancer cases
plot.Venn.all <- function() {
  Data <- read_tsv(file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt"))
  Data <- Data[!duplicated(Data$sample), ]
  Data <- na.omit(Data)
  row.names(Data) <- Data$sample
  Data <- as.data.frame(Data)
  Data$sample <- NULL
  Sample.ID <- row.names(Data)
  subset <- as.data.frame(Data$`_sample_type`)
  row.names(subset) <- row.names(Data)
  colnames(subset) <- "sample.type"
  tissue.GTEX.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    TCGA.GTEX <- as.data.frame(TCGA.GTEX)
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX <- as.data.frame(TCGA.GTEX)
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE
    )
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)]
    TCGA.GTEX.Lung.sampletype <- na.omit(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  gene.name <- names(TCGA.GTEX.sampletype)
  gene.name <- gene.name [!gene.name %in% c("Row.names", "sample.type")]

  ### TO BE CONTINUED
  EIF.correlation <- function(y, z) {
    TCGA.GTEX.tumor.lung <- TCGA.GTEX.sampletype[
      TCGA.GTEX.sampletype$sample.type %in% y,
    ]
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.tumor.lung[[x]],
        TCGA.GTEX.tumor.lung[[y]],
        method = "pearson"
      )
      res <- data.frame(x,
        y,
        result[c(
          "estimate",
          "p.value",
          "statistic",
          "method"
        )],
        stringsAsFactors = FALSE
      )
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(
        rbind.data.frame,
        lapply(gene.name,
          correlation.coefficient,
          y = x
        )
      )
      rownames(cor.data) <- cor.data[, 1]
      # cor.data1 <- cor.data[cor.data[, "p.value"] <= 0.05,]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    plot.pos.Venn <- function() {
      c4 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05
      )
      summary(c4)
      b <- vennCounts(c4)
      colnames(b) <- c(
        "EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts"
      )
      vennDiagram(b)
      pos.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], # EIF4E
        EIF4G1      = b[5, "Counts"], # EIF4G1
        EIF4A1      = b[3, "Counts"], # EIF4A1
        EIF4EBP1    = b[2, "Counts"], # EIF4EBP1
        "EIF4E&EIF4G1"   = b[13, "Counts"],
        "EIF4E&EIF4A1"   = b[11, "Counts"],
        "EIF4E&EIF4EBP1" = b[10, "Counts"],
        "EIF4G1&EIF4A1"   = b[7, "Counts"],
        "EIF4G1&EIF4EBP1"   = b[6, "Counts"],
        "EIF4A1&EIF4EBP1"   = b[4, "Counts"],
        "EIF4E&EIF4G1&EIF4A1" = b[15, "Counts"],
        "EIF4E&EIF4G1&EIF4EBP1" = b[14, "Counts"],
        "EIF4E&EIF4A1&EIF4EBP1" = b[12, "Counts"],
        "EIF4G1&EIF4A1&EIF4EBP1" = b[8, "Counts"],
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]
      ))
      p2 <- plot(pos.Venn2,
        # key = TRUE,
        main = paste(z, "posCOR"),
        lwd = 0,
        fill = c("#999999", "#009E73", "#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels = list(
          labels = c(
            "EIF4E",
            "EIF4G1",
            "EIF4A1",
            "EIF4EBP1"
          ),
          cex = 1.25
        )
      )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "Heatmap"),
        filename = paste("all", z, "pos4Venn.pdf"),
        plot = p2,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    plot.pos.Venn()
    plot.neg.Venn <- function() {
      c4 <- cbind(
        EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05
      )
      summary(c4)
      b <- vennCounts(c4)
      colnames(b) <- c(
        "EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts"
      )
      vennDiagram(b)
      neg.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], # EIF4E
        EIF4G1      = b[5, "Counts"], # EIF4G1
        EIF4A1      = b[3, "Counts"], # EIF4A1
        EIF4EBP1    = b[2, "Counts"], # EIF4EBP1
        "EIF4E&EIF4G1"   = b[13, "Counts"],
        "EIF4E&EIF4A1"   = b[11, "Counts"],
        "EIF4E&EIF4EBP1" = b[10, "Counts"],
        "EIF4G1&EIF4A1"   = b[7, "Counts"],
        "EIF4G1&EIF4EBP1"   = b[6, "Counts"],
        "EIF4A1&EIF4EBP1"   = b[4, "Counts"],
        "EIF4E&EIF4G1&EIF4A1" = b[15, "Counts"],
        "EIF4E&EIF4G1&EIF4EBP1" = b[14, "Counts"],
        "EIF4E&EIF4A1&EIF4EBP1" = b[12, "Counts"],
        "EIF4G1&EIF4A1&EIF4EBP1" = b[8, "Counts"],
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]
      ))
      p2 <- plot(neg.Venn2,
        # key = TRUE,
        main = paste(z, "negCOR"),
        lwd = 0,
        fill = c("#999999", "#009E73", "#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels = list(
          labels = c(
            "EIF4E",
            "EIF4G1",
            "EIF4A1",
            "EIF4EBP1"
          ),
          cex = 1.25
        )
      )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "Heatmap"),
        filename = paste("all", z, "neg4Venn.pdf"),
        plot = p2,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    plot.neg.Venn()

    c4 <- cbind(
      EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
      EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
      EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
      EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05
    )
    colnames(c4) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    df <- as.data.frame(summary(c4))
    df1 <- df[df$Freq %like% "TRUE", ]
    df1$Var1 <- NULL
    df1$Var2 <- gsub(" ", "", df1$Var2)
    row.names(df1) <- df1$Var2
    df1$Var2 <- NULL
    df1$Freq <- gsub("TRUE :", "", df1$Freq)
    df1$Freq <- as.numeric(df1$Freq)
    colnames(df1) <- "posCORs"

    c5 <- cbind(
      EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
      EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
      EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
      EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05
    )
    colnames(c5) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    dt <- as.data.frame(summary(c5))
    dt1 <- dt[dt$Freq %like% "TRUE", ]
    dt1$Var1 <- NULL
    dt1$Var2 <- gsub(" ", "", dt1$Var2)
    row.names(dt1) <- dt1$Var2
    dt1$Var2 <- NULL
    dt1$Freq <- gsub("TRUE :", "", dt1$Freq)
    dt1$Freq <- as.numeric(dt1$Freq)
    colnames(dt1) <- "negCORs"
    df2 <- cbind(df1, dt1)
    return(df2)
  }
  all.sample.type <- levels(subset$sample.type)
  all.tumor.type <- all.sample.type [!all.sample.type %in% c(
    "Cell Line",
    "Normal Tissue",
    "Solid Tissue Normal"
  )]

  EIF.cor.tumor <- EIF.correlation(y = all.tumor.type, z = "tumor")
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"), z = "normal")

  plot.bargraph.CORs(
    EIF.cor.tumor = EIF.cor.tumor,
    EIF.cor.normal = EIF.cor.normal,
    tumor.label = "tumor",
    normal.label = "normal",
    gene = gene,
    posCORs = posCORs,
    negCORs = negCORs,
    output.poscor.filename = paste("all posCORs.pdf"),
    output.negcor.filename = paste("all negCORs.pdf"),
    coord_flip.ylim = 14000)
}

### find posCOR and negCOR in the overlapping CORs from lung cancer cases
plot.Venn.lung <- function(x) {
  # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
  Sampletype <- read_tsv(file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt"))
  tissue.GTEX.TCGA.gene <- function(x) {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    )
    Lung <- Sampletype[Sampletype$`_primary_site` == x, ]
    Lung.ID <- as.vector(Lung$sample)
    Lung.ID <- na.omit(Lung.ID) # NA in the vector
    TCGA.GTEX.Lung <- TCGA.GTEX %>% select("sample", all_of(Lung.ID))
    TCGA.GTEX.Lung <- TCGA.GTEX.Lung[
      !duplicated(TCGA.GTEX.Lung$sample),
      !duplicated(colnames(TCGA.GTEX.Lung))
    ]
    row.names(TCGA.GTEX.Lung) <- TCGA.GTEX.Lung$sample
    TCGA.GTEX.Lung <- as.data.frame(TCGA.GTEX.Lung)
    TCGA.GTEX.Lung$sample <- NULL
    TCGA.GTEX.Lung.t <- data.table::transpose(TCGA.GTEX.Lung)
    rownames(TCGA.GTEX.Lung.t) <- colnames(TCGA.GTEX.Lung)
    colnames(TCGA.GTEX.Lung.t) <- rownames(TCGA.GTEX.Lung)

    Lung <- Lung[!duplicated(Lung$sample), ]
    Lung <- na.omit(Lung)
    row.names(Lung) <- Lung$sample
    Lung <- as.data.frame(Lung)
    Lung$sample <- NULL
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.Lung.t,
      Lung,
      by    = "row.names",
      all.x = TRUE
    )
    TCGA.GTEX.Lung.sampletype <- as.data.frame(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype.lung <- tissue.GTEX.TCGA.gene(x)
  row.names(TCGA.GTEX.sampletype.lung) <- TCGA.GTEX.sampletype.lung$Row.names
  TCGA.GTEX.sampletype.lung$Row.names <- NULL
  TCGA.GTEX.sampletype.lung$`_sample_type` <- as.factor(
    TCGA.GTEX.sampletype.lung$`_sample_type`
  )
  geneID <- colnames(Sampletype[Sampletype$`_primary_site` == x, ])
  TCGA.GTEX.lung <- TCGA.GTEX.sampletype.lung[
    ,
    !names(TCGA.GTEX.sampletype.lung) %in% geneID
  ]
  gene.name <- names(TCGA.GTEX.lung)

  EIF.correlation <- function(y, z) {
    TCGA.GTEX.subset.lung <- TCGA.GTEX.sampletype.lung[
      TCGA.GTEX.sampletype.lung$`_sample_type` %in% y,
    ]
    TCGA.GTEX.subset.lung$`_sample_type` <- droplevels(
      TCGA.GTEX.subset.lung$`_sample_type`
    )
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.subset.lung[[x]],
        TCGA.GTEX.subset.lung[[y]],
        method = "pearson"
      )
      res <- data.frame(x,
        y,
        result[c(
          "estimate",
          "p.value",
          "statistic",
          "method"
        )],
        stringsAsFactors = FALSE
      )
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(
        rbind.data.frame,
        lapply(gene.name,
          correlation.coefficient,
          y = x
        )
      )
      rownames(cor.data) <- cor.data[, 1]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")

    plot.pos.Venn <- function(y) {
      c4 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05
      )
      b <- vennCounts(c4)
      colnames(b) <- c(
        "EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts"
      )
      vennDiagram(b)
      pos.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], # EIF4E
        EIF4G1      = b[5, "Counts"], # EIF4G1
        EIF4A1      = b[3, "Counts"], # EIF4A1
        EIF4EBP1    = b[2, "Counts"], # EIF4EBP1
        "EIF4E&EIF4G1"   = b[13, "Counts"],
        "EIF4E&EIF4A1"   = b[11, "Counts"],
        "EIF4E&EIF4EBP1" = b[10, "Counts"],
        "EIF4G1&EIF4A1"   = b[7, "Counts"],
        "EIF4G1&EIF4EBP1"   = b[6, "Counts"],
        "EIF4A1&EIF4EBP1"   = b[4, "Counts"],
        "EIF4E&EIF4G1&EIF4A1" = b[15, "Counts"],
        "EIF4E&EIF4G1&EIF4EBP1" = b[14, "Counts"],
        "EIF4E&EIF4A1&EIF4EBP1" = b[12, "Counts"],
        "EIF4G1&EIF4A1&EIF4EBP1" = b[8, "Counts"],
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]
      ))
      p2 <- plot(pos.Venn2,
        # key = TRUE,
        lwd = 0,
        fill = c("#999999", "#009E73", "#56B4E9", "#E69F00"),
        main = paste(x, z, "posCOR"),
        quantities = list(cex = 1.25),
        labels = list(
          labels = c(
            "EIF4E",
            "EIF4G1",
            "EIF4A1",
            "EIF4EBP1"
          ),
          cex = 1.25
        )
      )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "Heatmap"),
        filename = paste(x, z, "pos4Venn.pdf"),
        plot = p2,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    plot.pos.Venn(y)

    plot.neg.Venn <- function(y) {
      c4 <- cbind(
        EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05
      )
      b <- vennCounts(c4)
      colnames(b) <- c(
        "EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts"
      )
      vennDiagram(b)
      neg.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], # EIF4E
        EIF4G1      = b[5, "Counts"], # EIF4G1
        EIF4A1      = b[3, "Counts"], # EIF4A1
        EIF4EBP1    = b[2, "Counts"], # EIF4EBP1
        "EIF4E&EIF4G1"   = b[13, "Counts"],
        "EIF4E&EIF4A1"   = b[11, "Counts"],
        "EIF4E&EIF4EBP1" = b[10, "Counts"],
        "EIF4G1&EIF4A1"   = b[7, "Counts"],
        "EIF4G1&EIF4EBP1"   = b[6, "Counts"],
        "EIF4A1&EIF4EBP1"   = b[4, "Counts"],
        "EIF4E&EIF4G1&EIF4A1" = b[15, "Counts"],
        "EIF4E&EIF4G1&EIF4EBP1" = b[14, "Counts"],
        "EIF4E&EIF4A1&EIF4EBP1" = b[12, "Counts"],
        "EIF4G1&EIF4A1&EIF4EBP1" = b[8, "Counts"],
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]
      ))
      p2 <- plot(neg.Venn2,
        # key = TRUE,
        lwd = 0,
        fill = c("#999999", "#009E73", "#56B4E9", "#E69F00"),
        main = paste(x, z, "negCOR"),
        quantities = list(cex = 1.25),
        labels = list(
          labels = c(
            "EIF4E",
            "EIF4G1",
            "EIF4A1",
            "EIF4EBP1"
          ),
          cex = 1.25
        )
      )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "Heatmap"),
        filename = paste(x, z, "neg4Venn.pdf"),
        plot = p2,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    plot.neg.Venn(y)

    c4 <- cbind(
      EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
      EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
      EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
      EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05
    )
    colnames(c4) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    df <- as.data.frame(summary(c4))
    df1 <- df[df$Freq %like% "TRUE", ]
    df1$Var1 <- NULL
    df1$Var2 <- gsub(" ", "", df1$Var2)
    row.names(df1) <- df1$Var2
    df1$Var2 <- NULL
    df1$Freq <- gsub("TRUE :", "", df1$Freq)
    df1$Freq <- as.numeric(df1$Freq)
    colnames(df1) <- "posCORs"

    c5 <- cbind(
      EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
      EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
      EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
      EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05
    )
    colnames(c5) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    dt <- as.data.frame(summary(c5))
    dt1 <- dt[dt$Freq %like% "TRUE", ]
    dt1$Var1 <- NULL
    dt1$Var2 <- gsub(" ", "", dt1$Var2)
    row.names(dt1) <- dt1$Var2
    dt1$Var2 <- NULL
    dt1$Freq <- gsub("TRUE :", "", dt1$Freq)
    dt1$Freq <- as.numeric(dt1$Freq)
    colnames(dt1) <- "negCORs"
    df2 <- cbind(df1, dt1)
    return(df2)
  }
  EIF.cor.tumor <- EIF.correlation(
    y = c(
      "Primary Tumor",
      "Metastatic",
      "Recurrent Tumor"
    ),
    z = "tumor"
  )
  EIF.cor.normal <- EIF.correlation(
    y = c("Normal Tissue"),
    z = "normal"
  )
  plot.bargraph.CORs(
    EIF.cor.tumor = EIF.cor.tumor,
    EIF.cor.normal = EIF.cor.normal,
    tumor.label = paste(x, "tumor"),
    normal.label = paste("Normal", x),
    gene = gene,
    posCORs = posCORs,
    negCORs = negCORs,
    output.poscor.filename = paste(x, "posCORs barplot.pdf"),
    output.negcor.filename = paste(x, "negCORs barplot.pdf"),
    coord_flip.ylim = 15000)
}

### plot heatmapand pathway analysis on clusters from all cancer cases ###
plot.heatmap.total <- function() {
  Data <- read_tsv(file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt"))
  Data <- Data[!duplicated(Data$sample), ]
  Data <- na.omit(Data)
  row.names(Data) <- Data$sample
  Data <- as.data.frame(Data)
  Data$sample <- NULL
  Sample.ID <- row.names(Data)
  subset <- as.data.frame(Data$`_sample_type`)
  row.names(subset) <- row.names(Data)
  colnames(subset) <- "sample.type"
  tissue.GTEX.TCGA.gene <- function() {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    ) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[
      !duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))
    ]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX <- as.data.frame(TCGA.GTEX)
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[, colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE
    )
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)]
    TCGA.GTEX.Lung.sampletype <- na.omit(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  gene.name <- names(TCGA.GTEX.sampletype)
  gene.name <- gene.name [!gene.name %in% c("Row.names", "sample.type")]

  ### TO BE CONTINUED
  EIF.correlation <- function(y, z) {
    TCGA.GTEX.tumor.lung <- TCGA.GTEX.sampletype[
      TCGA.GTEX.sampletype$sample.type %in% y,
    ]
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.tumor.lung[[x]],
        TCGA.GTEX.tumor.lung[[y]],
        method = "pearson"
      )
      res <- data.frame(x,
        y,
        result[c(
          "estimate",
          "p.value",
          "statistic",
          "method"
        )],
        stringsAsFactors = FALSE
      )
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(
        rbind.data.frame,
        lapply(gene.name,
          correlation.coefficient,
          y = x
        )
      )
      rownames(cor.data) <- cor.data[, 1]
      # cor.data1 <- cor.data[cor.data[, "p.value"] <= 0.05,]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    cor.data <- cbind(
      setNames(data.frame(EIF4E.cor[, c(3, 4)]), c("EIF4E", "EIF4E.p")),
      setNames(data.frame(EIF4G1.cor[, c(3, 4)]), c("EIF4G1", "EIF4G1.p")),
      setNames(data.frame(EIF4A1.cor[, c(3, 4)]), c("EIF4A1", "EIF4A1.p")),
      setNames(data.frame(EIF4EBP1.cor[, c(3, 4)]), c("EIF4EBP1", "EIF4EBP1.p"))
    )
    return(cor.data)
  }
  all.sample.type <- levels(subset$sample.type)
  all.tumor.type <- all.sample.type [!all.sample.type %in% c(
    "Cell Line",
    "Normal Tissue",
    "Solid Tissue Normal"
  )]
  EIF.cor.tumor <- EIF.correlation(y = all.tumor.type, z = "tumor")
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"), z = "normal")
  cor.data <- cbind(
    setNames(
      data.frame(EIF.cor.tumor[1:8]),
      c(
        "EIF4E.tumor", "EIF4E.p.tumor",
        "EIF4G1.tumor", "EIF4G1.p.tumor",
        "EIF4A1.tumor", "EIF4A1.p.tumor",
        "EIF4EBP1.tumor", "EIF4EBP1.p.tumor"
      )
    ),
    setNames(
      data.frame(EIF.cor.normal[1:8]),
      c(
        "EIF4E.normal", "EIF4E.p.normal",
        "EIF4G1.normal", "EIF4G1.p.normal",
        "EIF4A1.normal", "EIF4A1.p.normal",
        "EIF4EBP1.normal", "EIF4EBP1.p.normal"
      )
    )
  )
  DF <- as.matrix(na.omit(cor.data[
    cor.data$EIF4E.tumor > 0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4E.tumor < -0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor > 0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor < -0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor > 0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor < -0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor > 0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor < -0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4E.normal > 0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4E.normal < -0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4G1.normal > 0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4G1.normal < -0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal > 0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal < -0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal > 0.3 & cor.data$EIF4EBP1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal < -0.3 & cor.data$EIF4EBP1.p.normal <= 0.05,
  ]))
  DF <- DF[, c(1, 3, 5, 7, 9, 11, 13, 15)]
  pheatmap::pheatmap(DF,
    # main = "Correlation Coefficient Heatmap",
    # annotation_row = my_sample_col[,"_sample_type",drop=FALSE],
    angle_col = c("0"),
    fontsize = 12,
    fontface = "bold",
    color = colorRampPalette(rev(brewer.pal(
      n = 7,
      name = "RdYlBu"
    )))(100),
    show_rownames = FALSE,
    show_colnames = FALSE
  )
  # DF_scaled = t(scale(t(DF)))
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options
  ht1 <- Heatmap(DF,
    name = "Correlation Coefficient Heatmap (All)",
    heatmap_legend_param = list(
      labels_gp = gpar(font = 15),
      legend_width = unit(6, "cm"),
      direction = "horizontal"
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    bottom_annotation = HeatmapAnnotation(
      annotation_legend_param = list(direction = "horizontal"),
      type = c(
        "tumor", "tumor", "tumor", "tumor",
        "normal", "normal", "normal", "normal"
      ),
      col = list(type = c(
        "normal" = "royalblue",
        "tumor" = "pink"
      )),
      cn = anno_text(gsub("\\..*", "", colnames(DF)),
        location = 0,
        rot = 0,
        just = "center",
        gp = gpar(
          fontsize = 15,
          fontface = "bold"
        )
      )
    ),
    # cluster_rows  = as.dendrogram(hclust(dist(DF))),
    # row_split     = 3,
    row_km = 3,
    # row_km_repeats = 100,
    row_title = "cluster_%s",
    row_title_gp = gpar(
      fontsize = 15,
      fontface = "bold"
    ),
    border = TRUE,
    col = circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "#EEEEEE", "red")
    )
  )
  ht <- draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )

  pdf(file.path(
    path = file.path(output.directory, "Heatmap"),
    filename = "all tumors heatmap.pdf"
  ),
  width = 8,
  height = 8,
  useDingbats = FALSE
  )
  ht <- draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )
  dev.off()

  ## try to extract clusters from heatmap
  # Saving row names of cluster one
  plot.cluster.pathway <- function() {
    cluster.geneID.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]], ])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez <- mapIds(org.Hs.eg.db,
        keys = c1$V1,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )
      # c1 <- c1[!is.na(c1)]
      c1 <- na.omit(c1)
      return(c1$entrez)
    }
    cluster.num <- as.character(c(1:3))
    names(cluster.num) <- paste("cluster", 1:3)
    cluster.data <- lapply(cluster.num, cluster.geneID.list)
    ck.GO <- compareCluster(
      geneCluster = cluster.data,
      fun = "enrichGO",
      OrgDb = "org.Hs.eg.db"
    )
    ck.KEGG <- compareCluster(
      geneCluster = cluster.data,
      fun = "enrichKEGG"
    )
    ck.REACTOME <- compareCluster(
      geneCluster = cluster.data,
      fun = "enrichPathway"
    )
    head(as.data.frame(ck.REACTOME))
    p1 <- dotplot(ck.GO,
      title = "The Most Enriched GO Pathways",
      showCategory = 8,
      font.size = 18,
      includeAll = FALSE
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = black_bold_16(),
        legend.text = black_bold_16(),
        strip.text = black_bold_16()
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste("all tumors GO.pdf"),
      plot = p1,
      width = 10,
      height = 8,
      useDingbats = FALSE
    )

    p2 <- dotplot(ck.KEGG,
      title = "The Most Enriched KEGG Pathways",
      showCategory = 8,
      font.size = 18,
      includeAll = FALSE
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = black_bold_16(),
        legend.text = black_bold_16(),
        strip.text = black_bold_16()
      )
    print(p2)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste("all tumors KEGG.pdf"),
      plot = p2,
      width = 12,
      height = 8,
      useDingbats = FALSE
    )

    p3 <- dotplot(ck.REACTOME,
      title = "The Most Enriched REACTOME Pathways",
      showCategory = 8,
      font.size = 16,
      includeAll = FALSE
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = black_bold_16(),
        legend.text = black_bold_16(),
        strip.text = black_bold_16()
      )
    print(p3)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste("all tumors REACTOME.pdf"),
      plot = p3,
      width = 12,
      height = 8,
      useDingbats = FALSE
    )
  }
  plot.cluster.pathway()
}

### plot heatmapand pathway analysis on clusters from lung cancer cases ###
plot.heatmap.lung <- function(x) {
  # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
  Sampletype <- read_tsv(file.path(data.file.directory, "TcgaTargetGTEX_phenotype.txt"))
  tissue.GTEX.TCGA.gene <- function(x) {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      file.path(data.file.directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"),
      data.table = FALSE
    )
    Lung <- Sampletype[Sampletype$`_primary_site` == x, ]
    Lung.ID <- as.vector(Lung$sample)
    Lung.ID <- na.omit(Lung.ID) # NA in the vector
    TCGA.GTEX.Lung <- TCGA.GTEX %>% select("sample", all_of(Lung.ID))
    TCGA.GTEX.Lung <- TCGA.GTEX.Lung[
      !duplicated(TCGA.GTEX.Lung$sample),
      !duplicated(colnames(TCGA.GTEX.Lung))
    ]
    row.names(TCGA.GTEX.Lung) <- TCGA.GTEX.Lung$sample
    TCGA.GTEX.Lung <- as.data.frame(TCGA.GTEX.Lung)
    TCGA.GTEX.Lung$sample <- NULL
    TCGA.GTEX.Lung.t <- data.table::transpose(TCGA.GTEX.Lung)
    rownames(TCGA.GTEX.Lung.t) <- colnames(TCGA.GTEX.Lung)
    colnames(TCGA.GTEX.Lung.t) <- rownames(TCGA.GTEX.Lung)

    Lung <- Lung[!duplicated(Lung$sample), ]
    Lung <- na.omit(Lung)
    row.names(Lung) <- Lung$sample
    Lung <- as.data.frame(Lung)
    Lung$sample <- NULL
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.Lung.t,
      Lung,
      by    = "row.names",
      all.x = TRUE
    )
    TCGA.GTEX.Lung.sampletype <- as.data.frame(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype.lung <- tissue.GTEX.TCGA.gene(x)
  row.names(TCGA.GTEX.sampletype.lung) <- TCGA.GTEX.sampletype.lung$Row.names
  TCGA.GTEX.sampletype.lung <- as.data.frame(TCGA.GTEX.sampletype.lung)
  TCGA.GTEX.sampletype.lung$Row.names <- NULL
  TCGA.GTEX.sampletype.lung$`_sample_type` <- as.factor(
    TCGA.GTEX.sampletype.lung$`_sample_type`
  )
  # remove Solid Tissue Normal data
  # TCGA.GTEX.sampletype.lung <- TCGA.GTEX.sampletype.lung[
  #  !(TCGA.GTEX.sampletype.lung$`_sample_type` %in% "Solid Tissue Normal"), ]

  geneID <- colnames(Sampletype[Sampletype$`_primary_site` == x, ])
  TCGA.GTEX.lung <- TCGA.GTEX.sampletype.lung[
    ,
    !names(TCGA.GTEX.sampletype.lung) %in% geneID
  ]
  gene.name <- names(TCGA.GTEX.lung)
  EIF.correlation <- function(y) {
    TCGA.GTEX.subset.lung <- TCGA.GTEX.sampletype.lung[
      TCGA.GTEX.sampletype.lung$`_sample_type` %in% y,
    ]
    TCGA.GTEX.subset.lung$`_sample_type` <- droplevels(
      TCGA.GTEX.subset.lung$`_sample_type`
    )
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.subset.lung[[x]],
        TCGA.GTEX.subset.lung[[y]],
        method = "pearson"
      )
      res <- data.frame(x,
        y,
        result[c(
          "estimate",
          "p.value",
          "statistic",
          "method"
        )],
        stringsAsFactors = FALSE
      )
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(
        rbind.data.frame,
        lapply(gene.name,
          correlation.coefficient,
          y = x
        )
      )
      rownames(cor.data) <- cor.data[, 1]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")

    plot.pos.Venn <- function(y) {
      c3 <- cbind(
        EIF4E.cor$estimate > 0.3,
        EIF4G1.cor$estimate > 0.3,
        EIF4A1.cor$estimate > 0.3
      )
      a <- vennCounts(c3)
      colnames(a) <- c(
        "EIF4E",
        "EIF4G1",
        "EIF4A1",
        "Counts"
      )
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(
        A = a[5, "Counts"],
        B = a[3, "Counts"],
        C = a[2, "Counts"],
        "A&B" = a[7, "Counts"],
        "A&C" = a[6, "Counts"],
        "B&C" = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]
      ))
      p1 <- plot(pos.Venn,
        # key = TRUE,
        lwd = 0,
        border = "black",
        fill = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.5),
        # main       = y,
        labels = list(
          labels = c(
            "EIF4E posCOR",
            "EIF4G1 posCOR",
            "EIF4A1 posCOR"
          ),
          cex = 1.5
        )
      )
      print(p1)
      ggplot2::ggsave(
        path = file.path(output.directory, "Heatmap"),
        filename = paste(x, "Venn.pdf"),
        plot = p1,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
    }
    #plot.pos.Venn(y)
    cor.data <- cbind(
      setNames(data.frame(EIF4E.cor[, c(3, 4)]), c("EIF4E", "EIF4E.p")),
      setNames(data.frame(EIF4G1.cor[, c(3, 4)]), c("EIF4G1", "EIF4G1.p")),
      setNames(data.frame(EIF4A1.cor[, c(3, 4)]), c("EIF4A1", "EIF4A1.p")),
      setNames(data.frame(EIF4EBP1.cor[, c(3, 4)]), c("EIF4EBP1", "EIF4EBP1.p"))
    )
    return(cor.data)
  }
  EIF.cor.tumor <- EIF.correlation(y = c(
    "Primary Tumor",
    "Metastatic",
    "Recurrent Tumor"
  ))
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"))
  cor.data <- cbind(
    setNames(
      data.frame(EIF.cor.tumor[1:8]),
      c(
        "EIF4E.tumor", "EIF4E.p.tumor",
        "EIF4G1.tumor", "EIF4G1.p.tumor",
        "EIF4A1.tumor", "EIF4A1.p.tumor",
        "EIF4EBP1.tumor", "EIF4EBP1.p.tumor"
      )
    ),
    setNames(
      data.frame(EIF.cor.normal[1:8]),
      c(
        "EIF4E.normal", "EIF4E.p.normal",
        "EIF4G1.normal", "EIF4G1.p.normal",
        "EIF4A1.normal", "EIF4A1.p.normal",
        "EIF4EBP1.normal", "EIF4EBP1.p.normal"
      )
    )
  )
  DF <- as.matrix(na.omit(cor.data[
    cor.data$EIF4E.tumor > 0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4E.tumor < -0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor > 0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor < -0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor > 0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor < -0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor > 0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor < -0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4E.normal > 0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4E.normal < -0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4G1.normal > 0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4G1.normal < -0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal > 0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal < -0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal > 0.3 & cor.data$EIF4EBP1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal < -0.3 & cor.data$EIF4EBP1.p.normal <= 0.05,
  ]))
  DF <- DF[, c(1, 3, 5, 7, 9, 11, 13, 15)]
  # DF <- as.matrix(na.omit(cor.data))
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options)
  pheatmap::pheatmap(DF,
    # main = "Correlation Coefficient Heatmap",
    angle_col = c("0"),
    fontsize = 12,
    fontface = "bold",
    color = colorRampPalette(rev(brewer.pal(
      n = 7,
      name = "RdYlBu"
    )))(100),
    show_rownames = FALSE,
    show_colnames = TRUE
  )
  ht1 <- Heatmap(DF, # km=3,
    name = "Correlation Coefficient Heatmap(Lung)",
    heatmap_legend_param = list(
      direction = "horizontal",
      legend_width = unit(6, "cm")
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    bottom_annotation = HeatmapAnnotation(
      annotation_legend_param = list(direction = "horizontal"),
      type = c(
        "tumor", "tumor", "tumor", "tumor",
        "normal", "normal", "normal", "normal"
      ),
      col = list(type = c(
        "normal" = "royalblue",
        "tumor" = "pink"
      )),
      cn = anno_text(gsub("\\..*", "", colnames(DF)),
        location = 0,
        rot = 0,
        just = "center",
        gp = gpar(
          fontsize = 15,
          fontface = "bold"
        )
      )
    ),
    row_km = 3,
    # row_km_repeats = 100,
    row_title = "cluster_%s",
    row_title_gp = gpar(
      fontsize = 15,
      fontface = "bold"
    ),
    border = TRUE,
    col = circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "#EEEEEE", "red")
    )
  )
  ht <- draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )

  pdf(file.path(
    path = file.path(output.directory, "Heatmap"),
    filename = paste(x, "tumors heatmap.pdf")
  ),
  width = 8,
  height = 8,
  useDingbats = FALSE
  )
  ht <- draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )
  dev.off()
  ## try to extract clusters from heatmap
  # Saving row names of cluster one
  plot.cluster.pathway <- function() {
    cluster.geneID.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]], ])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez <- mapIds(org.Hs.eg.db,
        keys = c1$V1,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )
      c1 <- na.omit(c1)
      return(c1$entrez)
    }
    cluster.num <- as.character(c(1:3))
    names(cluster.num) <- paste("cluster", 1:3)
    cluster.data <- lapply(cluster.num, cluster.geneID.list)
    ck.GO <- compareCluster(
      geneCluster = cluster.data,
      fun = "enrichGO",
      OrgDb = "org.Hs.eg.db"
    )
    ck.KEGG <- compareCluster(
      geneCluster = cluster.data,
      fun = "enrichKEGG"
    )
    ck.REACTOME <- compareCluster(
      geneCluster = cluster.data,
      fun = "enrichPathway"
    )
    p1 <- dotplot(ck.GO,
      title = "The Most Enriched GO Pathways",
      showCategory = 8,
      font.size = 18,
      includeAll = FALSE
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = black_bold_16(),
        legend.text = black_bold_16(),
        strip.text = black_bold_16()
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste(x, "Go.pdf"),
      plot = p1,
      width = 10,
      height = 8,
      useDingbats = FALSE
    )
    p2 <- dotplot(ck.KEGG,
      title = "The Most Enriched KEGG Pathways",
      showCategory = 8,
      font.size = 18,
      includeAll = FALSE
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = black_bold_16(),
        legend.text = black_bold_16(),
        strip.text = black_bold_16()
      )
    print(p2)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste(x, "KEGG.pdf"),
      plot = p2,
      width = 12,
      height = 8,
      useDingbats = FALSE
    )
    p3 <- dotplot(ck.REACTOME,
      title = "The Most Enriched REACTOME Pathways",
      showCategory = 8,
      font.size = 16,
      includeAll = FALSE
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16(),
        axis.title = black_bold_16(),
        axis.text.x = black_bold_16(),
        axis.text.y = black_bold_16(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = black_bold_16(),
        legend.text = black_bold_16(),
        strip.text = black_bold_16()
      )
    print(p3)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste(x, "REACTOME.pdf"),
      plot = p3,
      width = 14,
      height = 8,
      useDingbats = FALSE
    )
  }
  plot.cluster.pathway()

  ### select posCOR for gene expression data ###
  DF.tumor <- as.matrix(na.omit(cor.data[cor.data$EIF4E.tumor > 0.3 |
    cor.data$EIF4G1.tumor > 0.3 |
    cor.data$EIF4A1.tumor > 0.3, ]))
  plot.cluster.heatmap <- function() {
    c1 <- t(t(row.names(DF[row_order(ht1)[["1"]], ])))
    c1 <- as.data.frame(c1)
    c1$V1 <- as.character(c1$V1)
    gene.list <- c1$V1 # row.names(DF.tumor)
    TCGA.GTEX.lung.genelist <- TCGA.GTEX.lung[, colnames(TCGA.GTEX.lung) %in% gene.list]
    DF2 <- as.matrix(na.omit(TCGA.GTEX.lung.genelist))
    DF2.t <- t(DF2)

    sample.ID <- row.names(DF2)
    sample.ID.type <- Sampletype[Sampletype$sample %in% sample.ID, ]

    my_sample_col <- sample.ID.type[
      ,
      colnames(sample.ID.type) %in% c("sample", "_sample_type")
    ]
    row.names(my_sample_col) <- my_sample_col$sample
    my_sample_col$sample <- NULL
    my_sample_col$`_sample_type` <- as.factor(my_sample_col$`_sample_type`)
    my_sample_col <- as.data.frame(my_sample_col)

    breaksList <- seq(-3, 3, by = 1)
    p4 <- pheatmap::pheatmap(
      DF2.t,
      main = "Gene expression Heatmap",
      annotation_col = my_sample_col[, "_sample_type", drop = FALSE],
      scale = "row",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete", # has to be complete
      treeheight_row = 0, treeheight_col = 0,
      cutree_cols = 2,
      angle_col = c("0"),
      fontsize = 12,
      fontface = "bold",
      color = colorRampPalette(rev(brewer.pal(
        n = 7,
        name = "RdYlBu"
      )))(length(breaksList)),
      breaks = breaksList,
      show_rownames = FALSE,
      show_colnames = FALSE
    )
    print(p4)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste(x, "cluster3 heatmap.pdf"),
      plot = p4,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )
  }
}

###
plot.Venn.all()
plot.Venn.lung(x = "Lung")
plot.heatmap.total()
plot.heatmap.lung(x = "Lung")


## Figure 6 ##
######################################
#### RNA protein correlation LUAD ####
######################################
EIF.pro.correlation <- function(y) {
  LUAD.Pro <- read_excel(
    file.path(data.file.directory, "Protein.xlsx"), col_names = FALSE
  )
  LUAD.Pro <- as.data.frame(LUAD.Pro)
  LUAD.Pro$...1 <- make.unique(LUAD.Pro$...1)
  row.names(LUAD.Pro) <- LUAD.Pro$...1
  LUAD.Pro <- LUAD.Pro[,-(1:1), drop = FALSE]
  LUAD.Pro.T <- t(LUAD.Pro)
  
  Pro_List <- colnames(LUAD.Pro.T) 
  Pro_List <- Pro_List [! Pro_List %in% c("Type", "Sample")]
  LUAD.Pro.T <- as.data.frame(LUAD.Pro.T, stringsAsFactors = FALSE)
  
  LUAD.Pro.T[Pro_List] <- sapply(LUAD.Pro.T[Pro_List], as.numeric)
  
  LUAD.Pro <- LUAD.Pro.T
  
  Scatter.plot <- function(x, y, z){
    #colnames(LUAD.Pro)[1] <- "Type"
    LUAD.Pro <- LUAD.Pro[
      LUAD.Pro$Type %in% "Tumor",
    ]
    p1 <- ggscatter(LUAD.Pro, 
                    x = x, 
                    y = y, #color = "black",
                    add = "reg.line", #conf.int = TRUE, 
                    add.params = list(color = "black", fill = "lightgray"), 
                    cor.coef = TRUE, cor.method = "pearson",
                    color = z,
                    xlab = paste("log2(", x, "protein ratio)"), 
                    ylab = paste("log2(", y, "protein ratio)")) + 
      # scale_y_continuous(breaks= scales::pretty_breaks())+
      scale_y_continuous(
        breaks = get_breaks(by = 1, from = -1),
        limits = c(-1, 2)) + # for 3G
      # Add correlation coefficient
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) #+ 
    #stat_cor(aes(color = Type),
    #         label.x.npc = 0, 
    #         label.y.npc = 0.9, 
    #         hjust = 0) 
    print(p1)
    
    ggplot2::ggsave(
      path = file.path(output.directory, "LUAD"),
      filename = paste(x, y, "cor.pdf"),
      plot = p1,
      width = 3,
      height = 3,
      useDingbats = FALSE
    )
  }
  
  Scatter.plot(x = "EIF4E",
               y = "EIF4G1", z = "dark red")
  Scatter.plot(x = "EIF4G1",
               y = "EIF4A1", z = "dark green")
  Scatter.plot(x = "EIF4A1",
               y = "EIF4E", z = "dark blue")
  
  ## cell division
  Scatter.plot(x = "EIF4G1",
               y = "CKAP2", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "CKAP2", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "CKAP2", z = "dark blue")
  
  Scatter.plot(x = "EIF4G1",
               y = "CCNA2", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "CCNA2", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "CCNA2", z = "dark blue")
  
  Scatter.plot(x = "EIF4G1",
               y = "ERCC6L", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "ERCC6L", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "ERCC6L", z = "dark blue")
  
  
  Scatter.plot(x = "EIF4G1",
               y = "MCM7", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "MCM7", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "MCM7", z = "dark blue")
  
  ## translation
  Scatter.plot(x = "EIF4G1",
               y = "RPS2", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "RPS2", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "RPS2", z = "dark blue")
  
  Scatter.plot(x = "EIF4G1",
               y = "EIF3B", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "EIF3B", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "EIF3B", z = "dark blue")
  
  Scatter.plot(x = "EIF4G1",
               y = "EIF3G", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "EIF3G", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "EIF3G", z = "dark blue")
  
  Scatter.plot(x = "EIF4G1",
               y = "EIF2S3", z = "dark green")
  Scatter.plot(x = "EIF4E",
               y = "EIF2S3", z = "dark red")
  Scatter.plot(x = "EIF4A1",
               y = "EIF2S3", z = "dark blue")
}    
EIF.pro.correlation()

######################################
## Boxplots for phosphor-proteomics ##
######################################
plot.EIF4.CPTAC.pro.LUAD <- function() {
  get.clinic.data <- function() {
    CPTAC.LUAD.Sample <- read_excel(
      file.path(data.file.directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx")
    )
    CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample[, c("Aliquot (Specimen Label)", "Type")]
    CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample.ID[
      !duplicated(CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`),
    ]
    row.names(CPTAC.LUAD.Sample.ID) <- CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`
    CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)` <- NULL
    CPTAC.LUAD.clinic <- read_excel(
      file.path(data.file.directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx"),
      sheet = 2
    )
    CPTAC.LUAD.clinic.Sampletype <- merge(CPTAC.LUAD.clinic,
      CPTAC.LUAD.Sample,
      by.x = "case_id",
      by.y = "Participant ID (case_id)"
    )
    # colnames(EIF.CPTAC.LUAD.clinic.Sampletype)
    CPTAC.LUAD.clinic.Sampletype <- CPTAC.LUAD.clinic.Sampletype[
      ,
      c("tumor_stage_pathological", "Aliquot (Specimen Label)", "Type")
    ]
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$Type == "Normal"] <- "Normal"
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological %in% c("Stage III", "Stage IIIA", "Stage IIIB")] <- "Stage III"
    return(CPTAC.LUAD.clinic.Sampletype)
  }
  CPTAC.LUAD.clinic.Sampletype <- get.clinic.data()

  get.proteomics.data <- function(x) {
    CPTAC.LUAD.Proteomics <- fread(
      file.path(data.file.directory, "CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv"),
      data.table = FALSE
    )
    EIF.CPTAC.LUAD.Proteomics <- CPTAC.LUAD.Proteomics[
      CPTAC.LUAD.Proteomics$Gene %in% x,
    ]
    row.names(EIF.CPTAC.LUAD.Proteomics) <- EIF.CPTAC.LUAD.Proteomics$Gene
    EIF.CPTAC.LUAD.Proteomics <- select(EIF.CPTAC.LUAD.Proteomics, -contains("Unshared"))

    EIF.CPTAC.LUAD.Proteomics$Gene <- NULL
    EIF.CPTAC.LUAD.Proteomics <- EIF.CPTAC.LUAD.Proteomics[1:(length(EIF.CPTAC.LUAD.Proteomics) - 6)]
    EIF.CPTAC.LUAD.Proteomics.t <- data.table::transpose(EIF.CPTAC.LUAD.Proteomics)
    rownames(EIF.CPTAC.LUAD.Proteomics.t) <- colnames(EIF.CPTAC.LUAD.Proteomics)
    colnames(EIF.CPTAC.LUAD.Proteomics.t) <- rownames(EIF.CPTAC.LUAD.Proteomics)
    rownames(EIF.CPTAC.LUAD.Proteomics.t) <- sub(" Log Ratio", "", rownames(EIF.CPTAC.LUAD.Proteomics.t))

    EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype <- merge(EIF.CPTAC.LUAD.Proteomics.t,
      CPTAC.LUAD.clinic.Sampletype,
      by.x  = "row.names",
      by.y  = "Aliquot (Specimen Label)"
    )
    row.names(EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype) <- EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype$Row.names
    EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype$Row.names <- NULL
    EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype$tumor_stage_pathological <- as.factor(
      EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype$tumor_stage_pathological
    )
    return(EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype)
  }
  EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype <- get.proteomics.data(
    x = c("EIF4G1", "EIF4A1", "EIF4EBP1", "EIF4E", "EIF4E2")
  )

  
  
  
  df <- melt(EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype,
    variable.name = "Gene"
  )
  # df$tumor_stage_pathological <- as.factor(df$tumor_stage_pathological)
  # levels(df$tumor_stage_pathological)

  pro.plot <- function() {
    p1 <- ggplot(
      data = df,
      aes(
        x = tumor_stage_pathological,
        y = 2**value
      )
    ) + # ylim(0, 2)+
      geom_boxplot(
        data = df,
        aes(fill = Gene),
        # alpha         = 0,
        # size     = .75,
        # width    = 1,
        outlier.color = NA,
        position = position_dodge(width = .9)
      ) +
      stat_n_text(
        size = 4,
        fontface = "bold",
        angle = 90,
        hjust = 0
      ) +
      labs(x = NULL, y = NULL) +
      # scale_color_manual(values = col_vector)+
      # scale_fill_manual(values = col_vector) +
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_discrete(drop = F) +
      facet_wrap(~Gene,
        scales = "free_y",
        nrow = 4, ncol = 7
      ) +
      facet_rep_wrap(~Gene,
        scales = "free_y",
        nrow = 4,
        ncol = 7,
        repeat.tick.labels = c("x", "y")
      ) +
      # ggplot2::facet_grid_sc(cols = vars(Gene),
      #              scales = list(y = scales_y))+
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12_45(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) +
      stat_compare_means(
        comparisons = list(
          c("Normal", "Stage I"),
          c("Normal", "Stage II"),
          c("Normal", "Stage III")
        ),
        method = "t.test",
        label = "p.signif", # label.y = 1,
        size = 4
      )
    print(p1)

    p2 <- ggplot(
      data = df,
      aes(
        x = Type, # tumor_stage_pathological,
        y = 2**value
      )
    ) + # ylim(0, 2)+
      geom_boxplot(
        data = df,
        aes(fill = Gene),
        # alpha         = 0,
        # size     = .75,
        # width    = 1,
        outlier.color = NA,
        position = position_dodge(width = .9)
      ) +
      stat_n_text(
        size = 4,
        fontface = "bold",
        angle = 90,
        hjust = 0
      ) +
      labs(x = NULL, y = NULL) +
      # scale_color_manual(values = col_vector)+
      # scale_fill_manual(values = col_vector) +
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_discrete(drop = F) +
      facet_wrap(~Gene,
        scales = "free_y",
        nrow = 4, ncol = 7
      ) +
      facet_rep_wrap(~Gene,
        scales = "free_y",
        nrow = 4,
        ncol = 7,
        repeat.tick.labels = c("x", "y")
      ) +
      # ggplot2::facet_grid_sc(cols = vars(Gene),
      #              scales = list(y = scales_y))+
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12_45(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) +
      stat_compare_means(
        comparisons = list(c("Normal", "Tumor")),
        method = "t.test",
        label = "p.signif", # label.y = 1,
        size = 4
      )
    print(p2)
  }
  pro.plot()

  nor.pro.plot <- function() {
    gene <- levels(df$Gene)
    for (x in gene) {
      df1.5 <- df[df$Gene == x, ]
      means <- aggregate(value ~ tumor_stage_pathological, df1.5, median)
      df1.5 <- df[df$Gene == x, ]

      stage <- levels(df1.5$tumor_stage_pathological)
      for (val in stage) {
        y <- means[means$tumor_stage_pathological == "Normal", ]$value
        # df1.5$nor <-  NULL
        df1.5$nor[df1.5$tumor_stage_pathological == val] <- df1.5$value[df1.5$tumor_stage_pathological == val] - y
      }
      df$nor[df$Gene == x] <- df1.5$nor[df1.5$Gene == x]
    }

    df <- na.omit(df)
    dataMedian <- summarise(group_by(df, Gene, Type), MD = 2**median(nor))

    p1 <- ggplot(
      data = df,
      aes(
        x = tumor_stage_pathological,
        y = 2**nor
      )
    ) + # ylim(0, 2)+
      geom_boxplot(
        data = df,
        aes(fill = Gene),
        # alpha         = 0,
        # size     = .75,
        # width    = 1,
        outlier.color = NA,
        position = position_dodge(width = .9)
      ) +
      stat_n_text(
        size = 4,
        fontface = "bold",
        angle = 90,
        hjust = 0
      ) +
      labs(x = NULL, y = NULL) +
      # scale_color_manual(values = col_vector)+
      # scale_fill_manual(values = col_vector) +
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_discrete(drop = F) +
      facet_wrap(~Gene,
        # scales = "free_y",
        nrow = 4, ncol = 7
      ) +
      facet_rep_wrap(~Gene,
        # scales = 'free_y',
        nrow = 4,
        ncol = 7,
        repeat.tick.labels = c("x", "y")
      ) +
      # ggplot2::facet_grid_sc(cols = vars(Gene),
      #              scales = list(y = scales_y))+
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12_45(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) +
      stat_compare_means(
        comparisons = list(
          c("Normal", "Stage I"),
          c("Normal", "Stage II"),
          c("Normal", "Stage III")
        ),
        method = "t.test",
        label = "p.signif", # label.y = 1,
        size = 4
      )
    print(p1)

    p2 <- ggplot(
      data = df,
      aes(
        x = Type, # tumor_stage_pathological,
        y = 2**nor
      )
    ) + # ylim(0, 2)+
      geom_boxplot(
        data = df,
        aes(fill = Gene),
        # alpha         = 0,
        # size     = .75,
        # width    = 1,
        outlier.color = NA,
        position = position_dodge(width = .9)
      ) +
      geom_text(
        data = dataMedian, aes(
          x = Type,
          y = MD,
          label = MD
        ),
        position = position_dodge(width = 0.8), size = 3, vjust = -0.5
      ) +
      stat_n_text(
        size = 4,
        fontface = "bold",
        angle = 90,
        hjust = 0
      ) +
      labs(x = NULL, y = NULL) +
      # scale_color_manual(values = col_vector)+
      # scale_fill_manual(values = col_vector) +
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_discrete(drop = F) +
      facet_wrap(~Gene,
        # scales = "free_y",
        nrow = 4, ncol = 7
      ) +
      facet_rep_wrap(~Gene,
        # scales = 'free_y',
        nrow = 4,
        ncol = 7,
        repeat.tick.labels = c("x", "y")
      ) +
      # ggplot2::facet_grid_sc(cols = vars(Gene),
      #              scales = list(y = scales_y))+
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12_45(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) +
      stat_compare_means(
        comparisons = list(c("Normal", "Tumor")),
        method = "t.test",
        label = "p.signif", # label.y = 1,
        size = 4
      )
    print(p2)
  }
  nor.pro.plot()

  get.proteomics.phosproteomics.data <- function() {
    CPTAC.LUAD.Phosproteomics <- fread(
      file.path(data.file.directory, "CPTAC3_Lung_Adeno_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv"),
      data.table = FALSE
    )

    get.phos.data <- function(x) {
      EIF.CPTAC.LUAD.Phosproteomics <- CPTAC.LUAD.Phosproteomics[
        CPTAC.LUAD.Phosproteomics$Gene %in% x,
      ]
      colnames(EIF.CPTAC.LUAD.Phosproteomics) <- sub(
        " Log Ratio",
        "",
        colnames(EIF.CPTAC.LUAD.Phosproteomics)
      )
      row.names(EIF.CPTAC.LUAD.Phosproteomics) <- EIF.CPTAC.LUAD.Phosproteomics$Phosphosite
      EIF.CPTAC.LUAD.Phosproteomics$Phosphosite <- NULL
      EIF.CPTAC.LUAD.Phosproteomics <- EIF.CPTAC.LUAD.Phosproteomics[1:(length(EIF.CPTAC.LUAD.Phosproteomics) - 3)]
      EIF.CPTAC.LUAD.Phosproteomics.t <- t(EIF.CPTAC.LUAD.Phosproteomics)
      colnames(EIF.CPTAC.LUAD.Phosproteomics.t) <- gsub(
        ".*:", paste0(x, ":"),
        colnames(EIF.CPTAC.LUAD.Phosproteomics.t)
      )
      return(EIF.CPTAC.LUAD.Phosproteomics.t)
    }
    phos.data <- lapply(
      c("EIF4G1", "EIF4EBP1", "EIF4A1", "EIF4E2"),
      get.phos.data
    )
    phos.data.all <- do.call("cbind", phos.data)

    combine.pho.phos.data <- function() {
      EIF.CPTAC.LUAD.Phosproteomics.Sampletype <- merge(
        phos.data.all,
        CPTAC.LUAD.clinic.Sampletype,
        by.x = "row.names",
        by.y = "Aliquot (Specimen Label)"
      )
      rownames(EIF.CPTAC.LUAD.Phosproteomics.Sampletype) <- EIF.CPTAC.LUAD.Phosproteomics.Sampletype$Row.names
      EIF.CPTAC.LUAD.Phosproteomics.Sampletype$Row.names <- NULL

      EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype$tumor_stage_pathological <- NULL
      EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype$Type <- NULL

      EIF.CPTAC.LUAD.phos.proteomics.Sampletype <- merge(
        EIF.CPTAC.LUAD.Phosproteomics.Sampletype,
        EIF.CPTAC.LUAD.Proteomics.clinic.Sampletype,
        by  = "row.names"
      )
      row.names(EIF.CPTAC.LUAD.phos.proteomics.Sampletype) <- EIF.CPTAC.LUAD.phos.proteomics.Sampletype$Row.names
      EIF.CPTAC.LUAD.phos.proteomics.Sampletype$Row.names <- NULL
      return(EIF.CPTAC.LUAD.phos.proteomics.Sampletype)
    }
    EIF.CPTAC.LUAD.phos.proteomics.Sampletype <- combine.pho.phos.data()

    df2 <- reshape2::melt(EIF.CPTAC.LUAD.phos.proteomics.Sampletype,
      variable.name = "Gene"
    )
    # df2$Gene <- gsub(".*:",paste0("EIF4G1",":"),df2$Gene )
    df2$tumor_stage_pathological <- as.factor(df2$tumor_stage_pathological)
    df2$Gene <- as.factor(df2$Gene)
    levels(df2$Gene)
    df2$Gene <- ordered(df2$Gene,
      levels = c(
        "EIF4G1", "EIF4G1:t212", "EIF4G1:t212t214", "EIF4G1:t212t218",
        "EIF4G1:t214", "EIF4G1:t214t218",
        "EIF4G1:t218", "EIF4G1:s503", "EIF4G1:t654",
        "EIF4G1:s711", "EIF4G1:s1035", "EIF4G1:s1068",
        "EIF4G1:t1080s1087", "EIF4G1:t1080s1099",
        "EIF4G1:s1084", "EIF4G1:s1084s1087", "EIF4G1:s1084s1099",
        "EIF4G1:s1087",
        "EIF4G1:s1099", "EIF4G1:s1151s1154", "EIF4G1:s1152s1154",
        "EIF4G1:s1201", "EIF4G1:s1206",
        "EIF4G1:s1216", "EIF4G1:t1218",
        "EIF4G1:s1238", "EIF4G1:s1245",
        "EIF4EBP1", "EIF4EBP1:s35t46",
        "EIF4EBP1:t36t46", "EIF4EBP1:t37t46",
        "EIF4EBP1:s65", "EIF4EBP1:s65t68", "EIF4EBP1:t70",
        "EIF4EBP1:s83", "EIF4EBP1:s85", "EIF4EBP1:s86",
        "EIF4EBP1:s94", "EIF4EBP1:s96", "EIF4EBP1:s101",
        "EIF4A1", "EIF4A1:s78", "EIF4A1:y70", "EIF4A1:s205", "EIF4A1:t207",
        "EIF4E",
        "EIF4E2", "EIF4E2:s13"
      )
    )
    df2$Gene <- droplevels(df2$Gene)
    df2 <- df2[!is.na(df2$Gene), ]
    df2 <- df2[!is.na(df2$value), ]
    return(df2)
  }
  df2 <- get.proteomics.phosproteomics.data()

  phos.plot <- function() {
    color <- function(x) {
      n <- x
      qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
      col_vector <- unlist(mapply(
        brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      col <- sample(col_vector, n)
      return(col)
    }
    col_vector <- color(48)

    p2 <- ggplot(
      data = df2,
      aes(
        x = tumor_stage_pathological,
        y = 2**value
      )
    ) + # ylim(0, 2)+
      geom_boxplot(
        data = df2,
        aes(fill = Gene),
        # alpha         = 0,
        # size     = .75,
        # width    = 1,
        outlier.color = "black",
        position = position_dodge(width = .9)
      ) +
      stat_n_text(
        size = 4,
        fontface = "bold",
        angle = 90,
        hjust = 0
      ) +
      labs(x = NULL, y = NULL) +
      # scale_color_manual(values = col_vector)+
      # scale_fill_manual(values = col_vector) +
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_discrete(drop = F) +
      facet_wrap(~Gene, scales = "free_y", ncol = 7) +
      facet_rep_wrap(~Gene,
        scales = "free_y",
        ncol = 7,
        repeat.tick.labels = c("x", "y")
      ) +
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12_45(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) +
      stat_compare_means(
        comparisons = list(
          c("Normal", "Stage I"),
          c("Normal", "Stage II"),
          c("Normal", "Stage III")
        ),
        method = "t.test",
        label = "p.signif", # label.y = 1,
        size = 4
      )
    print(p2)

    p3 <- function(x) {
      mydf <- df2[df2$Gene %in% x, ]
      sts <- boxplot.stats(mydf$value)$stats # Compute lower and upper whisker limits
      p2 <- ggplot(
        data = mydf,
        aes(
          x = tumor_stage_pathological,
          y = 2**value
        )
      ) + # ylim(0, 2)+
        geom_boxplot(
          data = mydf,
          aes(fill = Gene),
          # alpha         = 0,
          # size     = .75,
          # width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        # coord_cartesian(ylim = c(sts[2]/2, max(sts)*2))+
        stat_n_text(
          size = 4,
          fontface = "bold",
          angle = 0,
          hjust = 0.5
        ) +
        labs(x = NULL, y = NULL) +

        # coord_cartesian(ylim = ylim1*1.2)+
        scale_x_discrete(limits = c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV")) +
        # scale_fill_manual(values = col_vector) +
        # scale_fill_brewer(palette="Dark2")+
        scale_fill_discrete(drop = F) +
        ggplot2::facet_grid(. ~ Gene) +
        facet_wrap(~Gene, scales = "free_y") +
        # facet_rep_wrap(~ Gene,
        #               scales = 'free_y',
        #               ncol   = 7,
        #               repeat.tick.labels = c('x','y')) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title.x = black_bold_12(),
          axis.title.y = black_bold_12(),
          axis.text.x = black_bold_12_45(),
          axis.text.y = black_bold_12(),
          # axis.line.x      = element_line(color = "black"),
          # axis.line.y      = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.text = black_bold_12(),
          strip.background = element_rect(fill = "white")
        ) +
        stat_compare_means(
          comparisons = list(
            c("Normal", "Stage I"),
            c("Normal", "Stage II"),
            c("Normal", "Stage III")
          ),
          method = "t.test",
          label = "p.signif", # label.y = 1,
          size = 4
        )
      print(p2)
    }
    p3(x = "EIF4EBP1:s65")
    lapply(
      c(
        "EIF4G1", "EIF4G1:t212", "EIF4G1:t212t214", "EIF4G1:t212t218",
        "EIF4G1:t214", "EIF4G1:t214t218",
        "EIF4G1:t218", "EIF4G1:s503", "EIF4G1:t654",
        "EIF4G1:s711", "EIF4G1:s1035", "EIF4G1:s1068",
        "EIF4G1:t1080s1087", "EIF4G1:t1080s1099",
        "EIF4G1:s1084", "EIF4G1:s1084s1087", "EIF4G1:s1084s1099",
        "EIF4G1:s1087",
        "EIF4G1:s1099", "EIF4G1:s1151s1154", "EIF4G1:s1152s1154",
        "EIF4G1:s1201", "EIF4G1:s1206",
        "EIF4G1:s1216", "EIF4G1:t1218",
        "EIF4G1:s1238", "EIF4G1:s1245",
        "EIF4EBP1", "EIF4EBP1:s35t46",
        "EIF4EBP1:t36t46", "EIF4EBP1:t37t46",
        "EIF4EBP1:s65", "EIF4EBP1:s65t68", "EIF4EBP1:t70",
        "EIF4EBP1:s83", "EIF4EBP1:s85", "EIF4EBP1:s86",
        "EIF4EBP1:s94", "EIF4EBP1:s96", "EIF4EBP1:s101",
        "EIF4A1", "EIF4A1:s78", "EIF4A1:y70", "EIF4A1:s205", "EIF4A1:t207",
        "EIF4E",
        "EIF4E2", "EIF4E2:s13"
      ),
      p3
    )
  }
  phos.plot()

  nor.phos.plot <- function() {
    gene <- levels(df2$Gene)
    for (x in gene) {
      df1.5 <- df2[df2$Gene == x, ]
      avg <- aggregate(value ~ tumor_stage_pathological, df1.5, median)
      df1.5 <- df2[df2$Gene == x, ]

      stage <- levels(df1.5$tumor_stage_pathological)
      for (val in stage) {
        y <- avg[avg$tumor_stage_pathological == "Normal", ]$value
        # df1.5$nor <-  NULL
        df1.5$nor[df1.5$tumor_stage_pathological == val] <- df1.5$value[df1.5$tumor_stage_pathological == val] - y
      }
      df2$nor[df2$Gene == x] <- df1.5$nor[df1.5$Gene == x]
    }
    df2 <- na.omit(df2)
    dataMedian <- summarise(group_by(df2, Gene, Type), MD = 2**median(nor))

    color <- function(x) {
      n <- x
      qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
      col_vector <- unlist(mapply(
        brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      col <- sample(col_vector, n)
      return(col)
    }
    col_vector <- color(48)
    p1 <- ggplot(
      data = df2,
      aes(
        x = tumor_stage_pathological,
        y = 2**nor
      )
    ) + # ylim(0, 5)+
      geom_boxplot(
        data = df2,
        aes(fill = Gene),
        # alpha         = 0,
        # size     = .75,
        # width    = 1,
        outlier.shape = NA,
        position = position_dodge(width = .9)
      ) +
      stat_n_text(
        size = 4,
        fontface = "bold",
        angle = 90,
        vjust = 0.5,
        hjust = 0
      ) +
      labs(x = NULL, y = NULL) +
      # scale_color_manual(values = col_vector)+
      # scale_fill_manual(values = col_vector) +
      # scale_fill_brewer(palette="Dark2")+
      scale_fill_discrete(drop = F) +
      coord_cartesian(ylim = c(-5, 12)) +
      facet_wrap(~Gene,
        # scales = "free_y",
        ncol = 8
      ) +
      facet_rep_wrap(~Gene,
        # scales = 'free_y',
        ncol = 8,
        repeat.tick.labels = c("x", "y")
      ) +
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12_45(),
        axis.text.y = black_bold_12(),
        # axis.line.x      = element_line(color = "black"),
        # axis.line.y      = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.text = black_bold_12(),
        strip.background = element_rect(fill = "white")
      ) +
      stat_compare_means(
        comparisons = list(
          c("Normal", "Stage I"),
          c("Normal", "Stage II"),
          c("Normal", "Stage III")
        ),
        method = "t.test",
        label = "p.signif", label.y = c(7, 9, 11),
        size = 4
      )
    print(p1)

    p3 <- function(x) {
      mydf <- df2[df2$Gene %in% x, ]
      sts <- boxplot.stats(mydf$value)$stats # Compute lower and upper whisker limits
      hline <- dataMedian$MD[dataMedian$Gene == x & dataMedian$Type == "Tumor"]
      hline <- round(hline, digits = 3)
      p2 <- ggplot(
        data = mydf,
        aes(
          x = tumor_stage_pathological,
          y = 2**nor
        )
      ) + # ylim(0, 2)+
        geom_boxplot(
          data = mydf,
          aes(fill = Gene),
          # alpha         = 0,
          # size     = .75,
          # width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        geom_hline(yintercept = hline, colour = "dark red", linetype = "dashed") +
        annotate("text",
          label = hline,
          x = "Stage IV",
          y = hline,
          vjust = -1,
          size = 4,
          colour = "dark red"
        ) +
        stat_n_text(
          size = 4,
          fontface = "bold", y.pos = 0.4,
          angle = 0,
          hjust = 0.5
        ) +
        labs(x = NULL, y = NULL) +
        coord_cartesian(ylim = c(0.4, 2.2)) +
        scale_y_continuous(breaks = seq(0.4, 2.2, by = 0.3)) +
        scale_x_discrete(limits = c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV")) +
        # scale_fill_manual(values = col_vector) +
        # scale_fill_brewer(palette="Dark2")+
        scale_fill_discrete(drop = F) +
        ggplot2::facet_grid(. ~ Gene) +
        facet_wrap(~Gene) +
        # facet_rep_wrap(~ Gene,
        #               scales = 'free_y',
        #               ncol   = 7,
        #               repeat.tick.labels = c('x','y')) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title.x = black_bold_12(),
          axis.title.y = black_bold_12(),
          axis.text.x = black_bold_12_45(),
          axis.text.y = black_bold_12(),
          # axis.line.x      = element_line(color = "black"),
          # axis.line.y      = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.text = black_bold_12(),
          strip.background = element_rect(fill = "white")
        ) +
        stat_compare_means(
          comparisons = list(
            c("Normal", "Stage I"),
            c("Normal", "Stage II"),
            c("Normal", "Stage III")
          ),
          method = "t.test",
          label = "p.signif", label.y = c(1.7, 1.9, 2.1),
          size = 4
        )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "CPTAC"),
        filename = paste0(str_remove(x, ":"), "pro.pdf"),
        plot = p2,
        width = 2.5,
        height = 3,
        useDingbats = FALSE
      )
    }
    p3(x = "EIF4A1:s78")
    lapply(c("EIF4G1", "EIF4EBP1", "EIF4A1", "EIF4E", "EIF4E2"), p3)

    p4 <- function(x) {
      mydf <- df2[df2$Gene %in% x, ]
      hline <- dataMedian$MD[dataMedian$Gene == x & dataMedian$Type == "Tumor"]
      hline <- round(hline, digits = 3)
      p2 <- ggplot(
        data = mydf,
        aes(
          x = tumor_stage_pathological,
          y = 2**nor
        )
      ) + # ylim(0, 2)+
        geom_boxplot(
          data = mydf,
          aes(fill = Gene),
          # alpha         = 0,
          # size     = .75,
          # width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        geom_hline(yintercept = hline, colour = "dark red", linetype = "dashed") +
        annotate("text",
          label = hline,
          x = "Stage IV",
          y = hline,
          vjust = -1,
          size = 4,
          colour = "dark red"
        ) +
        stat_n_text(
          size = 4,
          fontface = "bold", y.pos = -1,
          angle = 0,
          hjust = 0.5
        ) +
        labs(x = NULL, y = NULL) +
        coord_cartesian(ylim = c(-1, 11)) +
        scale_y_continuous(breaks = seq(-1, 11, by = 2)) +
        scale_x_discrete(limits = c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV")) +
        # scale_fill_manual(values = col_vector) +
        # scale_fill_brewer(palette="Dark2")+
        scale_fill_discrete(drop = F) +
        ggplot2::facet_grid(. ~ Gene) +
        facet_wrap(~Gene) +
        # facet_rep_wrap(~ Gene,
        #               scales = 'free_y',
        #               ncol   = 7,
        #               repeat.tick.labels = c('x','y')) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title.x = black_bold_12(),
          axis.title.y = black_bold_12(),
          axis.text.x = black_bold_12_45(),
          axis.text.y = black_bold_12(),
          # axis.line.x      = element_line(color = "black"),
          # axis.line.y      = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.text = black_bold_12(),
          strip.background = element_rect(fill = "white")
        ) +
        stat_compare_means(
          comparisons = list(
            c("Normal", "Stage I"),
            c("Normal", "Stage II"),
            c("Normal", "Stage III")
          ),
          method = "t.test",
          label = "p.signif", label.y = c(7.2, 8.4, 9.6),
          size = 4
        )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "CPTAC"),
        filename = paste0(str_remove(x, ":"), "phospro.pdf"),
        plot = p2,
        width = 2.5,
        height = 3,
        useDingbats = FALSE
      )
    }
    p4(x = "EIF4EBP1:s65")
    lapply(
      c( # "EIF4G1",
        "EIF4G1:t212", "EIF4G1:t212t214", "EIF4G1:t212t218",
        "EIF4G1:t214", "EIF4G1:t214t218",
        "EIF4G1:t218", "EIF4G1:s503", "EIF4G1:t654",
        "EIF4G1:s711", "EIF4G1:s1035", "EIF4G1:s1068",
        "EIF4G1:t1080s1087", "EIF4G1:t1080s1099",
        "EIF4G1:s1084", "EIF4G1:s1084s1087", "EIF4G1:s1084s1099",
        "EIF4G1:s1087",
        "EIF4G1:s1099", "EIF4G1:s1151s1154", "EIF4G1:s1152s1154",
        "EIF4G1:s1201", "EIF4G1:s1206",
        "EIF4G1:s1216", "EIF4G1:t1218",
        "EIF4G1:s1238", "EIF4G1:s1245",
        # "EIF4EBP1",
        "EIF4EBP1:s35t46", "EIF4EBP1:t36t46", "EIF4EBP1:t37t46",
        "EIF4EBP1:s65", "EIF4EBP1:s65t68", "EIF4EBP1:t70",
        "EIF4EBP1:s83", "EIF4EBP1:s85", "EIF4EBP1:s86",
        "EIF4EBP1:s94", "EIF4EBP1:s96", "EIF4EBP1:s101",
        # "EIF4A1","EIF4A1:s78",
        "EIF4A1:y70", "EIF4A1:s205", "EIF4A1:t207",
        # "EIF4E",
        "EIF4E2", "EIF4E2:s13"
      ),
      p4
    )
  }
  nor.phos.plot()
}
#plot.EIF4.CPTAC.pro.LUAD()

plot.EIF4.CPTAC.pro.LUAD2 <- function() {
  get.clinic.data <- function() {
    CPTAC.LUAD.Sample <- read_excel(
      file.path(data.file.directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx")
    )
    CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample[, c("Aliquot (Specimen Label)", "Type")]
    CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample.ID[
      !duplicated(CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`),
      ]
    row.names(CPTAC.LUAD.Sample.ID) <- CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`
    CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)` <- NULL
    CPTAC.LUAD.clinic <- read_excel(
      file.path(data.file.directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx"),
      sheet = 2
    )
    CPTAC.LUAD.clinic.Sampletype <- merge(CPTAC.LUAD.clinic,
                                          CPTAC.LUAD.Sample,
                                          by.x = "case_id",
                                          by.y = "Participant ID (case_id)"
    )
    # colnames(EIF.CPTAC.LUAD.clinic.Sampletype)
    CPTAC.LUAD.clinic.Sampletype <- CPTAC.LUAD.clinic.Sampletype[
      ,
      c("tumor_stage_pathological", "Aliquot (Specimen Label)", "Type")
      ]
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$Type == "Normal"] <- "Normal"
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
    CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological[CPTAC.LUAD.clinic.Sampletype$tumor_stage_pathological %in% c("Stage III", "Stage IIIA", "Stage IIIB")] <- "Stage III"
    return(CPTAC.LUAD.clinic.Sampletype)
  }
  CPTAC.LUAD.clinic.Sampletype <- get.clinic.data()
  
  get.proteomics.data <- function(EIF_list) {
    LUAD.Proteomics <- read_excel(
      file.path(data.file.directory, "Protein.xlsx"), col_names = FALSE
    )
    #EIF_list <- c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1","Type","Sample")
    EIF.LUAD.Proteomics <- LUAD.Proteomics[LUAD.Proteomics$...1 %in% c(EIF_list, "Type","Sample"), ]
    EIF.LUAD.Proteomics.T <- t(EIF.LUAD.Proteomics)
    colnames(EIF.LUAD.Proteomics.T) <- EIF.LUAD.Proteomics.T[1, ]
    EIF.LUAD.Proteomics.T <- EIF.LUAD.Proteomics.T[-1, ] 

    EIF.LUAD.Proteomics.T <- as.data.frame(EIF.LUAD.Proteomics.T, stringsAsFactors = FALSE)

    EIF.LUAD.Proteomics.T[EIF_list] <- sapply(EIF.LUAD.Proteomics.T[EIF_list], as.numeric)
    return(EIF.LUAD.Proteomics.T)
  }
  EIF.LUAD.Proteomics <- get.proteomics.data(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1","AKT1","RPS6KB1","MKNK1","MKNK2","EIF4B", "EIF4H", "MAPK8","MAPK14","RPTOR"))
  
  
  get.rna.data <- function(EIF_list) {
    LUAD.RNA <- read_excel(
      file.path(data.file.directory, "RNA.xlsx"), col_names = FALSE
    )
    #EIF_list <- c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1","Type","Sample")
    EIF.LUAD.RNA <- LUAD.RNA[LUAD.RNA$...1 %in% c(EIF_list, "Type","Sample"), ]
    EIF.LUAD.RNA.T <- t(EIF.LUAD.RNA)
    colnames(EIF.LUAD.RNA.T) <- EIF.LUAD.RNA.T[1, ]
    EIF.LUAD.RNA.T <- EIF.LUAD.RNA.T[-1, ] 
    EIF.LUAD.RNA.T <- as.data.frame(EIF.LUAD.RNA.T, stringsAsFactors = FALSE)
    EIF.LUAD.RNA.T[EIF_list] <- sapply(EIF.LUAD.RNA.T[EIF_list], as.numeric)
    
    return(EIF.LUAD.RNA.T)
  }
  EIF.LUAD.RNA <- get.rna.data(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1"))
  
  get.phos.data <- function(EIF_list) {
    #EIF_list <- c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1","AKT1","RPS6KB1")
    LUAD.Phos <- read_excel(
      file.path(data.file.directory, "Phos.xlsx"), col_names = FALSE
    )
    #EIF_list <- c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1","Type","Sample")
    EIF.LUAD.Phos <- LUAD.Phos[LUAD.Phos$...1 %in% c(EIF_list, "Type","Sample"), ]
    EIF.LUAD.Phos <- as.data.frame(EIF.LUAD.Phos)
    row.names(EIF.LUAD.Phos) <- paste(EIF.LUAD.Phos$...1,EIF.LUAD.Phos$...2)
    EIF.LUAD.Phos <- EIF.LUAD.Phos[,-(1:2), drop = FALSE]

    EIF.LUAD.Phos.T <- t(EIF.LUAD.Phos)
    Phos_List <- colnames(EIF.LUAD.Phos.T) 
    Phos_List <- Phos_List [! Phos_List %in% c("Type na", "Sample na")]
    EIF.LUAD.Phos.T <- as.data.frame(EIF.LUAD.Phos.T, stringsAsFactors = FALSE)
    
    EIF.LUAD.Phos.T[Phos_List] <- sapply(EIF.LUAD.Phos.T[Phos_List], as.numeric)
    return(EIF.LUAD.Phos.T)
  }
  EIF.LUAD.Phos <- get.phos.data(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1","AKT1","MTOR","RPS6KB1", "MKNK1","MKNK2","EIF4B", "EIF4H", "MAPK8","MAPK14","RPTOR"))
  
  
  plot.EIF.cor.LUAD <- function(){
    EIF.LUAD.Proteomics <- get.proteomics.data(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1"))
    EIF.LUAD.Proteomics.Tumor <- EIF.LUAD.Proteomics[EIF.LUAD.Proteomics$Type == "Tumor", ] 
    EIF.LUAD.Proteomics.NAT <- EIF.LUAD.Proteomics[EIF.LUAD.Proteomics$Type == "NAT", ] 
    
    EIF.LUAD.RNA.Tumor <- EIF.LUAD.RNA[EIF.LUAD.RNA$Type == "Tumor", ] 
    EIF.LUAD.RNA.NAT <- EIF.LUAD.RNA[EIF.LUAD.RNA$Type == "NAT", ] 
    
    tumor.total <- merge(EIF.LUAD.Proteomics.Tumor, EIF.LUAD.RNA.Tumor, by=c("Sample","Type"), suffixes = c(".pro",".rna"))
    NAT.total <- merge(EIF.LUAD.Proteomics.NAT, EIF.LUAD.RNA.NAT, by=c("Sample","Type"), suffixes = c(".pro",".rna"))
    
    Scatter.plot <- function(x){
      p1 <- ggscatter(tumor.total, 
                      x = paste0(x, ".pro"), 
                      y = paste0(x, ".rna"), #color = "Type",
                      add = "reg.line", #conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson", title = x,
                      xlab = "Protein expresion", ylab = "RNA expression")+
        theme_bw() +
        theme(
          plot.title = black_bold_12(),
          axis.title.x = black_bold_12(),
          axis.title.y = black_bold_12(),
          axis.text.x = black_bold_12(),
          axis.text.y = black_bold_12(),
          # axis.line.x      = element_line(color = "black"),
          # axis.line.y      = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.text = black_bold_12(),
          strip.background = element_rect(fill = "white")
        ) 
      print(p1)
      ggplot2::ggsave(
                      path = file.path(output.directory, "LUAD"),
                      filename = paste(x, "corLUAD.pdf"),
                      plot = p1,
                      width = 3,
                      height = 3,
                      useDingbats = FALSE
                    )
    }
    Scatter.plot("EIF4G1")
    lapply(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1"), Scatter.plot)
  }
  plot.EIF.cor.LUAD()
  
  get.proteomics.phosproteomics.data <- function() {
    combine.pho.phos.data <- function() {
      EIF.LUAD.Phos.Sampletype <- merge(
        EIF.LUAD.Phos,
        CPTAC.LUAD.clinic.Sampletype,
        by.x = "Sample na",
        by.y = "Aliquot (Specimen Label)"
      )
      
      EIF.LUAD.Phos.Sampletype$`Type na` <- NULL
      EIF.LUAD.Phos.Sampletype$Type <- NULL
      
      #EIF.LUAD.Proteomics$Type <- NULL
      
      EIF.LUAD.phos.pro.Sampletype <- merge(
        EIF.LUAD.Phos.Sampletype,
        EIF.LUAD.Proteomics,
        by.x = "Sample na",
        by.y = "Sample"
      )
      row.names(EIF.LUAD.phos.pro.Sampletype) <- EIF.LUAD.phos.pro.Sampletype$Row.names
      EIF.LUAD.phos.pro.Sampletype$Row.names <- NULL
      return(EIF.LUAD.phos.pro.Sampletype)
    }
    EIF.LUAD.phos.pro.Sampletype <- combine.pho.phos.data()
    
    df2 <- reshape2::melt(EIF.LUAD.phos.pro.Sampletype,
                          variable.name = "Gene"
    )
    # df2$Gene <- gsub(".*:",paste0("EIF4G1",":"),df2$Gene )
    df2$tumor_stage_pathological <- as.factor(df2$tumor_stage_pathological)
    df2$Gene <- as.factor(df2$Gene)
    levels(df2$Gene)
    df2$Gene <- droplevels(df2$Gene)
    df2 <- df2[!is.na(df2$Gene), ]
    df2 <- df2[!is.na(df2$value), ]
    return(df2)
  }
  df2 <- get.proteomics.phosproteomics.data()
  df2$Gene <- factor(as.character(df2$Gene))
  levels(df2$Gene)
  
  
  nor.phos.plot <- function() {
    gene <- levels(df2$Gene)
    for (x in gene) {
      df1.5 <- df2[df2$Gene == x, ]
      avg <- aggregate(value ~ tumor_stage_pathological, df1.5, median)
      df1.5 <- df2[df2$Gene == x, ]
      
      stage <- levels(df1.5$tumor_stage_pathological)
      for (val in stage) {
        y <- avg[avg$tumor_stage_pathological == "Normal", ]$value
        # df1.5$nor <-  NULL
        df1.5$nor[df1.5$tumor_stage_pathological == val] <- df1.5$value[df1.5$tumor_stage_pathological == val] - y
      }
      df2$nor[df2$Gene == x] <- df1.5$nor[df1.5$Gene == x]
    }
    df2 <- na.omit(df2)
    dataMedian <- summarise(group_by(df2, Gene, Type), MD = 2**median(nor))
    
    color <- function(x) {
      n <- x
      qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
      col_vector <- unlist(mapply(
        brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      col <- sample(col_vector, n)
      return(col)
    }
    col_vector <- color(48)
    
    
    p3 <- function(x) {
      mydf <- df2[df2$Gene %in% x, ]
      sts <- boxplot.stats(mydf$value)$stats # Compute lower and upper whisker limits
      hline <- dataMedian$MD[dataMedian$Gene == x & dataMedian$Type == "Tumor"]
      hline <- round(hline, digits = 3)
      p2 <- ggplot(
        data = mydf,
        aes(
          x = tumor_stage_pathological,
          y = 2**nor
        )
      ) + #ylab("Normalized peptide counts") +
        geom_boxplot(
          data = mydf,
          aes(fill = Gene),
          # alpha         = 0,
          # size     = .75,
          # width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        geom_hline(yintercept = hline, colour = "dark red", linetype = "dashed") +
        annotate("text",
                 label = hline,
                 x = "Stage IV",
                 y = hline,
                 vjust = -1,
                 size = 4,
                 colour = "dark red"
        ) +
        stat_n_text(
          size = 4,
          fontface = "bold", y.pos = 0,
          angle = 0,
          hjust = 0.5
        ) +
        labs(x = NULL, y = "Normalized peptide ratio") +
        #coord_cartesian(ylim = c(-0.5, 7.5)) +
        #scale_y_continuous(breaks = seq(0, 7.5, by = 1)) +
        coord_cartesian(ylim = c(0, 3)) +
        scale_y_continuous(breaks = seq(0, 3, by = 1)) +
        scale_x_discrete(limits = c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV")) +
        # scale_fill_manual(values = col_vector) +
        # scale_fill_brewer(palette="Dark2")+
        scale_fill_discrete(drop = F) +
        ggplot2::facet_grid(. ~ Gene) +
        facet_wrap(~Gene) +
        # facet_rep_wrap(~ Gene,
        #               scales = 'free_y',
        #               ncol   = 7,
        #               repeat.tick.labels = c('x','y')) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title.x = black_bold_12(),
          axis.title.y = black_bold_12(),
          axis.text.x = black_bold_12_45(),
          axis.text.y = black_bold_12(),
          # axis.line.x      = element_line(color = "black"),
          # axis.line.y      = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.text = black_bold_12(),
          strip.background = element_rect(fill = "white")
        ) +
        stat_compare_means(
          comparisons = list(
            c("Normal", "Stage I"),
            c("Normal", "Stage II"),
            c("Normal", "Stage III")
          ),
          method = "t.test",
          label = "p.signif", 
          label.y = c(5, 6, 7),
          #label.y = c(2.4, 2.7, 3),
          size = 4
        )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "CPTAC"),
        filename = paste0(str_remove(x, ":"), "pro.pdf"),
        plot = p2,
        width = 3,
        height = 3,
        useDingbats = FALSE
      )
    }
    List <- levels(df2$Gene)
    #toMatch <- c("EIF4G1","EIF4E", "EIF4A1", "MKNK1", "MKNK2")
    #toMatch <- c("AKT1","MTOR")
    toMatch <- c("EIF4G1","EIF4E", "EIF4A1","MTOR", "AKT1", "EIF4EBP1","MKNK1","MKNK2")
    List2 <- List[grepl(paste(toMatch,collapse="|"), List)]
    lapply(List2, p3)
    #lapply(c("EIF4B","EIF4H", "AKT1"), p3)
    
    p4 <- function(x) {
      mydf <- df2[df2$Gene %in% x, ]
      hline <- dataMedian$MD[dataMedian$Gene == x & dataMedian$Type == "Tumor"]
      hline <- round(hline, digits = 3)
      p2 <- ggplot(
        data = mydf,
        aes(
          x = tumor_stage_pathological,
          y = 2**nor
        )
      ) + # ylim(0, 2)+
        geom_boxplot(
          data = mydf,
          aes(fill = Gene),
          # alpha         = 0,
          # size     = .75,
          # width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        geom_hline(yintercept = hline, colour = "dark red", linetype = "dashed") +
        annotate("text",
                 label = hline,
                 x = "Stage IV",
                 y = hline,
                 vjust = -1,
                 size = 4,
                 colour = "dark red"
        ) +
        stat_n_text(
          size = 4,
          fontface = "bold", y.pos = -1,
          angle = 0,
          hjust = 0.5
        ) +
        labs(x = NULL, y = "Normalized peptide ratio") +
        coord_cartesian(ylim = c(-1, 15)) +
        scale_y_continuous(breaks = seq(-1, 15, by = 2)) +
        scale_x_discrete(limits = c("Normal", "Stage I", "Stage II", "Stage III", "Stage IV")) +
        # scale_fill_manual(values = col_vector) +
        # scale_fill_brewer(palette="Dark2")+
        scale_fill_discrete(drop = F) +
        ggplot2::facet_grid(. ~ Gene) +
        facet_wrap(~Gene) +
        # facet_rep_wrap(~ Gene,
        #               scales = 'free_y',
        #               ncol   = 7,
        #               repeat.tick.labels = c('x','y')) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title.x = black_bold_12(),
          axis.title.y = black_bold_12(),
          axis.text.x = black_bold_12_45(),
          axis.text.y = black_bold_12(),
          # axis.line.x      = element_line(color = "black"),
          # axis.line.y      = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.text = black_bold_12(),
          strip.background = element_rect(fill = "white")
        ) +
        stat_compare_means(
          comparisons = list(
            c("Normal", "Stage I"),
            c("Normal", "Stage II"),
            c("Normal", "Stage III")
          ),
          method = "t.test",
          label = "p.signif", label.y = c(12.4, 13.6, 14.8),
          size = 4
        )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "CPTAC"),
        filename = paste0(str_remove(x, ":"), "phospro.pdf"),
        plot = p2,
        width = 3,
        height = 3,
        useDingbats = FALSE
      )
    }
    List <- levels(df2$Gene)
    #toMatch <- c("EIF4G1", "EIF4E", "EIF4A1", "MKNK1","MKNK2")
    toMatch <- c("EIF4G1", "EIF4B","EIF4H")
    List2 <- List[grepl(paste(toMatch,collapse="|"), List)]
    lapply(List2, p4)
  }
  nor.phos.plot()
}
plot.EIF4.CPTAC.pro.LUAD2()



  
  
  
# Figure 1
plot.bargraph.EIF.CNV.sum(c("TP53", "EIF4A1","EIF4A2", "EIF4E", "EIF4E2", "EIF4E3", 
                            "MYC", "EIF3D", "EIF4EBP1", "EIF4G1","EIF4G2", "EIF4G3",
                            "EIF4H", "EIF4B"))


plot.matrix.EIF.CNV.corr(c("TP53", "EIF4A1","EIF4A2", "EIF4E", "EIF4E2", "EIF4E3", 
                           "MYC", "EIF3D", "EIF4EBP1", "EIF4G1", "EIF4G2", "EIF4G3",
                           "EIF4H", "EIF4B"))

lapply(c("TP53", "EIF4A1","EIF4A2", "EIF4E", "EIF4E2", "EIF4E3", 
         "MYC", "EIF3D", "EIF4EBP1", "EIF4G1","EIF4G2", "EIF4G3",
         "EIF4H","EIF4B"), plot.bargraph.EIF.CNV.TCGA)


plot.boxgraph.EIF.RNAseq.TCGA (c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
                                 "EIF4G1","EIF4G2","EIF4G3",
                                 "EIF4B","EIF4H",
                                 "EIF4A1","EIF4A2","EIF3D","TP53","MYC"))

lapply(c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
         "EIF4G1","EIF4G2","EIF4G3",
         "EIF4B","EIF4H",
         "EIF4A1","EIF4A2","EIF3D","TP53","MYC"),
      plot.boxgraph.EIF.CNVratio.TCGA)

plot.boxgraph.EIF.CNVratio.TCGA("EIF4G1")


# Figure 2
plot.boxgraph.EIF.RNAseq.TCGA.GTEX(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1"))
plot.boxgraph.EIF.ratio.TCGA.GTEX(c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1"))
plot.violingraph.EIF.RNAseq.TCGA(c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1"))
plot.violingraph.EIF.RNAseq.TCGA(c("EIF4A1", "EIF4A2", 
                                   "EIF4E", "EIF4E2", "EIF4E3", 
                                   "EIF3D", "EIF4EBP1", 
                                   "EIF4G1", "EIF4G2", "EIF4G3",
                                   "EIF4H", "EIF4B", "MYC", "TP53"))


plot.violingraph.EIF.ratio.TCGA(c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1"))


# Figure 3
plot.EIF.TCGA.GTEX.PCA.all.tumor.tissue(c(
  "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))


plot.EIF.TCGA.GTEX.PCA.all.tumor.tissue(c(
  "EIF4G1","EIF4G2", "EIF4G3",
  "EIF4A1","EIF4A2", 
  "EIF4E", "EIF4E2", "EIF4E3", 
  "EIF3D", "EIF4EBP1", 
  "EIF4H", "EIF4B", "MYC", "JUN"
))


plot.EIF.TCGA.GTEX.PCA.each.tumor(
  c(
    "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
    "PABPC1", "MKNK1", "MKNK2"#, "MYC"
  ),
  "Skin"
)

plot.EIF.TCGA.PCA.all.tumor(c(
  "EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))

plot.EIF.TCGA.PCA.all.tumor(c(
  "EIF4G1","EIF4G2", "EIF4G3",
  "EIF4A1","EIF4A2", 
  "EIF4E", "EIF4E2", "EIF4E3", 
  "EIF3D", "EIF4EBP1", 
  "EIF4H", "EIF4B", "MYC", "JUN"
))

plot.EIF.GTEX.PCA.all.tissue(c(
  "EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))

plot.EIF.CPTAC.PCA.LUAD()


# Figure 4
plot.km.EIF.all.tumors("HNRNPL")
lapply(
  c(
    "EIF4G1","EIF4G2", "EIF4G3",
    "EIF4A1","EIF4A2", 
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF3D", "EIF3E",
    "EIF4EBP1", "EIF4EBP2", 
    "EIF4H", "EIF4B", "MYC",
    "PABPC1", "MKNK1", "MKNK2"
  ), plot.km.EIF.all.tumors
)

plot.km.EIF.each.tumor(
  EIF = "EIF4G1",
  tumor = c("lung squamous cell carcinoma")
)

lapply(
  c(
    "EIF4G1","EIF4G2", "EIF4G3",
    "EIF4A1","EIF4A2", 
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF3D", "EIF3E","EIF4EBP1", "EIF4EBP2", 
    "EIF4H", "EIF4B", "MYC",
    "PABPC1", "MKNK1", "MKNK2"
  ),
  plot.km.EIF.each.tumor,
  tumor = "lung adenocarcinoma"
)
plot.coxph.EIF.all.tumors()
plot.coxph.EIF.each.tumor(c("lung adenocarcinoma"))


# Figure 5
plot.Venn.all()
plot.Venn.lung(x = "Lung")
plot.heatmap.total()
plot.heatmap.lung(x = "Lung")


# Figure 6
plot.EIF4.CPTAC.pro.LUAD()
