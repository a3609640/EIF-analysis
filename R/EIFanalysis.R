library(AnnotationDbi)
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
library(glmnet)
library(gplots)
library(gridExtra)
library(igraph)
library(KEGG.db)
library(limma)
library(missMDA)
library(org.Hs.eg.db)
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

##############################################
## boxplot for EIF expression across tumors ##
##############################################
plot.boxgraph.EIF.RNAseq.TCGA <- function (EIF.gene) {
  pan.TCGA.gene <- function(){
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.pancancer <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
      TCGA.sampletype,
      by    = "row.names",
      all.x = TRUE)
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    sample.type.list <- levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    cancer.type.list <- levels(TCGA.RNAseq.anno$primary.disease)
    return(TCGA.RNAseq.sampletype)
  }
  TCGA.RNAseq.anno <- pan.TCGA.gene()
  
  pancancer.TCGA.EIF <- function(){
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno [
      !TCGA.RNAseq.anno$sample.type %in% "Solid Tissue Normal", ]
    row.names(TCGA.RNAseq.anno.subset) <- TCGA.RNAseq.anno.subset$Row.names
    TCGA.RNAseq.anno.subset$Row.names <- NULL
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno.subset[ ,
      colnames(TCGA.RNAseq.anno.subset) %in% c(EIF.gene, 
                                               "sample.type",
                                               "primary.disease")]
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset [
    #  EIF.TCGA.RNAseq.anno.subset$primary.disease %in% c("lung squamous cell carcinoma","lung adenocarcinoma"), ]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) + 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) -1)
    EIF.TCGA.RNAseq.anno.subset$`EIF4G1-(EIF4E+EIF4EBP1)` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4G1) - 2**(EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1`) +1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.gene, 
                                                                 "EIF4E+EIF4EBP1",
                                                                 #"EIF4G1-(EIF4E+EIF4EBP1)",
                                                                 "sample.type",
                                                                 "primary.disease")]
    #EIF.TCGA.RNAseq.anno.subset$delta <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4G1 - 2**EIF.TCGA.RNAseq.anno.subset$EIF4E - 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -1)
    EIF.TCGA.RNAseq.mean <- aggregate(EIF.TCGA.RNAseq.anno.subset[, 1:5], list(EIF.TCGA.RNAseq.anno.subset$primary.disease), median)
    cor_5 <- rcorr(as.matrix(EIF.TCGA.RNAseq.mean[ ,2:6]))
    M <- cor_5$r
    p_mat <- cor_5$P
    pdf(file.path(
      path        = "~/Documents/EIF_output/Expression", 
      filename    = "EIFmediancor.pdf"), 
      width       = 6, 
      height      = 6, 
      useDingbats = FALSE)
    corrplot(
      M, 
      method      = "color", 
      tl.cex      = 1, 
      tl.srt      = 45,
      number.cex  = 1, 
      addgrid.col = "gray",
      addCoef.col = "black", 
      tl.col      = "black",
      #type        = "upper", 
      order       = "FPC", 
      p.mat       = p_mat, 
      sig.level   = 0.05, #insig = "blank" 
    )
    dev.off()
    
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    # correlation plot
    my_data <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.gene, "EIF4E+EIF4EBP1")]
    cor_5 <- rcorr(as.matrix(my_data))
    M <- cor_5$r
    p_mat <- cor_5$P
    
    pdf(file.path(
      path        = "~/Documents/EIF_output/Expression", 
      filename    = "EIFsumPCAcor.pdf"), 
      width       = 6, 
      height      = 6, 
      useDingbats = FALSE)
    corrplot(
      M, 
      method      = "color", 
      tl.cex      = 1, tl.srt = 45,
      number.cex  = 1, 
      addgrid.col = "gray",
      addCoef.col = "black", 
      tl.col      = "black",
      #type        = "upper", 
      order       = "FPC", 
      p.mat       = p_mat, 
      sig.level   = 0.05, #insig = "blank" 
    )
    dev.off()

    
    EIF.TCGA.RNAseq.anno.subset.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long$sample.type <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$sample.type)
    EIF.TCGA.RNAseq.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.disease)
    return(EIF.TCGA.RNAseq.anno.subset.long)}
  pancancer.TCGA.EIF.long <- pancancer.TCGA.EIF()
  
  # reorder bars by explicitly ordering factor levels
  make.plot <- function () {
    mean <- within(pancancer.TCGA.EIF.long[
      pancancer.TCGA.EIF.long$variable == "EIF4E", ], # TCGAstudy is one column in df2
    primary.disease <- reorder(primary.disease, value, median))
    mean$primary.disease <- as.factor(mean$primary.disease)
    neworder <- levels(mean$primary.disease)
    x.ordered <- factor(pancancer.TCGA.EIF.long$primary.disease, 
                        levels = neworder)
    pancancer.TCGA.EIF.long1 <- pancancer.TCGA.EIF.long[
      pancancer.TCGA.EIF.long$variable %in% EIF.gene, ]
    f1 <- factor(pancancer.TCGA.EIF.long1$primary.disease)
    f.ordered1 <- fct_rev(f1)
    p1 <- ggplot(data = pancancer.TCGA.EIF.long1,
                 aes(x     = f.ordered1,  
                    #x     = x.ordered, # order primary disease
                     y     = 2**value,
                     color = variable)) +
      scale_y_continuous(trans = log2_trans(), 
                         labels = label_comma()) +
      stat_n_text(size     = 5, 
                  fontface = "bold", 
                  hjust    = 0) +
      geom_boxplot(
        alpha    = .01,
        #size     = .75,
        #width    = 1,
        position = position_dodge(width = .9)
      ) +
      scale_color_manual(
        values = c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00"),
        breaks = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1"),
        labels = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1")) + #for color-blind palettes
      labs(x = "primary disease",
           y = paste("normalized RNA counts")) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_12,
        axis.title.x         = black_bold_tahoma_12,
        axis.title.y         = element_blank(),
        axis.text.x          = black_bold_tahoma_12,
        axis.text.y          = black_bold_tahoma_12,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        legend.title         = element_blank(),
        legend.text          = black_bold_tahoma_12,
        legend.position      = "top",
        legend.justification = "left",
        legend.box           = "horizontal", 
        strip.text           = black_bold_tahoma_12
        )
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/Expression/TCGA", 
      filename    = "EIFexpressionTCGA.pdf", 
      plot        = p1,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    
    pancancer.TCGA.EIF.long2 <- pancancer.TCGA.EIF.long[
      pancancer.TCGA.EIF.long$variable %in% c("EIF4G1","EIF4E+EIF4EBP1"), ]
    f2 <- factor(pancancer.TCGA.EIF.long2$primary.disease)
    f.ordered2 <- fct_rev(f2)
    p2 <- ggplot(data = pancancer.TCGA.EIF.long2,
      aes(
        x     = f.ordered2,  
        #x     = x.ordered, # order primary disease
        y     = 2**value,
        color = variable)) +
      stat_n_text(size = 5, fontface = "bold", hjust = 0) +
      stat_summary(
        fun.y = "median",
        geom  = 'line', # alpha = 0.5,
        aes(group = variable, colour = variable),
            position = position_dodge(width = 0.5)) +
      scale_y_continuous(trans = log2_trans(), labels = label_comma()) +
      geom_boxplot(aes(colour = variable),
        alpha    = .01,
        #size     = .75,
        #width    = 1,
        position = position_dodge(width = 0.5)
      ) +
    scale_color_manual(
      values = c("#009E73","#CC79A7"),#for color-blind palettes
      breaks = c("EIF4E+EIF4EBP1", "EIF4G1"),
      labels = c("EIF4E+EIF4EBP1", "EIF4G1")) + 
    labs(x = "primary disease",
         y = paste("normalized RNA counts")) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title           = black_bold_tahoma_12,
      axis.title.x         = black_bold_tahoma_12,
      axis.title.y         = element_blank(),
      axis.text.x          = black_bold_tahoma_12,
      axis.text.y          = black_bold_tahoma_12,
      axis.line.x          = element_line(color = "black"),
      axis.line.y          = element_line(color = "black"),
      panel.grid           = element_blank(),
      legend.title         = element_blank(),
      legend.text          = black_bold_tahoma_12,
      legend.position      = "top",
      legend.justification = "left",
      legend.box           = "horizontal", 
      strip.text           = black_bold_tahoma_12
    ) 
    print(p2)
    ggsave(
    path        = "~/Documents/EIF_output/Expression/TCGA", 
    filename    = "EIFsumexpressionTCGA.pdf", 
    plot        = p2,
    width       = 7, 
    height      = 8, 
    useDingbats = FALSE)
    }
  make.plot()
  
}
plot.boxgraph.EIF.RNAseq.TCGA(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1"))

plot.boxgraph.EIF.RNAseq.GTEX <- function (EIF.gene) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE)
    # check the name of the last column
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  GTEX.RNAseq.anno <- tissue.GTEX.TCGA.gene()
  
  GTEX.EIF <- function(){
    TCGA.RNAseq.anno.subset <- GTEX.RNAseq.anno [
      GTEX.RNAseq.anno$sample.type %in% "Normal Tissue", ]
    row.names(TCGA.RNAseq.anno.subset) <- TCGA.RNAseq.anno.subset$Row.names
    TCGA.RNAseq.anno.subset$Row.names <- NULL
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno.subset[ ,
      colnames(TCGA.RNAseq.anno.subset) %in% c(EIF.gene, 
                                               "sample.type",
                                               "primary.site")]
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset [
    #  EIF.TCGA.RNAseq.anno.subset$primary.disease %in% c("lung squamous cell carcinoma","lung adenocarcinoma"), ]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) + 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) -1)
    EIF.TCGA.RNAseq.anno.subset$`EIF4G1-(EIF4E+EIF4EBP1)` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4G1) - 2**(EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1`) +1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.gene, 
      "EIF4E+EIF4EBP1",
      #"EIF4G1-(EIF4E+EIF4EBP1)",
      "sample.type",
      "primary.site")]
    #EIF.TCGA.RNAseq.anno.subset$delta <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4G1 - 2**EIF.TCGA.RNAseq.anno.subset$EIF4E - 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -1
    
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long$sample.type <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$sample.type)
    EIF.TCGA.RNAseq.anno.subset.long$primary.site <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.site)
    EIF.TCGA.RNAseq.anno.subset.long <- EIF.TCGA.RNAseq.anno.subset.long[
      !EIF.TCGA.RNAseq.anno.subset.long$value == 0, ]
    return(EIF.TCGA.RNAseq.anno.subset.long)}
  GTEX.EIF.long <- GTEX.EIF()
  # reorder bars by explicitly ordering factor levels
  make.plot <- function () {
    make.plot1 <- function(){
      GTEX.EIF.long1 <- GTEX.EIF.long[
        GTEX.EIF.long$variable %in% EIF.gene, ]
      f1 <- factor(GTEX.EIF.long1$primary.site)
      f.ordered1 <- fct_rev(f1)
    p1 <- ggplot(data = GTEX.EIF.long1,
      aes(x     = f.ordered1,  
        #x     = x.ordered, # order primary disease
        y     = 2**value,
        color = variable)) +
      scale_y_continuous(trans = log2_trans(), 
        labels = label_comma()) +
      stat_n_text(size     = 5, 
        fontface = "bold", 
        hjust    = 0) +
      geom_boxplot(
        alpha    = .01,
        #size     = .75,
        #width    = 1,
        position = position_dodge(width = .9)
      ) +
      scale_color_manual(
        values = c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00"),
        breaks = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1"),
        labels = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1")) + #for color-blind palettes
      labs(x = "primary disease",
        y = paste("normalized RNA counts")) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_12,
        axis.title.x         = black_bold_tahoma_12,
        axis.title.y         = element_blank(),
        axis.text.x          = black_bold_tahoma_12,
        axis.text.y          = black_bold_tahoma_12,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        legend.title         = element_blank(),
        legend.text          = black_bold_tahoma_12,
        legend.position      = "top",
        legend.justification = "left",
        legend.box           = "horizontal", 
        strip.text           = black_bold_tahoma_12
      )
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/Expression/GTEX", 
      filename    = "EIFexpressionGTEX.pdf", 
      plot        = p1,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
  }
    make.plot1()
    
    make.plot2 <- function(){
      GTEX.EIF.long2 <- GTEX.EIF.long[
      GTEX.EIF.long$variable %in% c("EIF4G1","EIF4E+EIF4EBP1"), ]
    f2 <- factor(GTEX.EIF.long2$primary.site)
    f.ordered2 <- fct_rev(f2)
    p2 <- ggplot(data = GTEX.EIF.long2,
      aes(
        x     = f.ordered2,  
        #x     = x.ordered, # order primary disease
        y     = 2**value,
        color = variable)) +
      stat_n_text(size = 5, fontface = "bold", hjust = 0) +
      stat_summary(
        fun.y = "median",
        geom  = 'line', # alpha = 0.5,
        aes(group = variable, colour = variable),
        position = position_dodge(width = 0.5)) +
      scale_y_continuous(trans = log2_trans(), labels = label_comma()) +
      geom_boxplot(aes(colour = variable),
        alpha    = .01,
        #size     = .75,
        #width    = 1,
        position = position_dodge(width = 0.5)
      ) +
      scale_color_manual(
        values = c("#009E73","#CC79A7"),#for color-blind palettes
        breaks = c("EIF4E+EIF4EBP1", "EIF4G1"),
        labels = c("EIF4E+EIF4EBP1", "EIF4G1")) + 
      labs(x = "primary disease",
        y = paste("normalized RNA counts")) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_12,
        axis.title.x         = black_bold_tahoma_12,
        axis.title.y         = element_blank(),
        axis.text.x          = black_bold_tahoma_12,
        axis.text.y          = black_bold_tahoma_12,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        legend.title         = element_blank(),
        legend.text          = black_bold_tahoma_12,
        legend.position      = "top",
        legend.justification = "left",
        legend.box           = "horizontal", 
        strip.text           = black_bold_tahoma_12
      ) 
    print(p2)
    ggsave(
      path        = "~/Documents/EIF_output/Expression/GTEX", 
      filename    = "EIFsumexpressionGTEX.pdf", 
      plot        = p2,
      width       = 7, 
      height      = 8, 
      useDingbats = FALSE)
    }
    make.plot2()
  }
  make.plot()
  
}
plot.boxgraph.EIF.RNAseq.GTEX(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1"))

plot.boxgraph.EIF.RNAseq.TCGA.GTEX <- function (EIF.gene) {
  pan.TCGA.gene <- function(){
    TCGA.pancancer <- fread(
      #"~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count",
      data.table = FALSE)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
  

    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)

    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    subset <- TCGA.sampletype[ ,c("_sample_type", 
                                  "primary disease or tissue",
                                  "_primary_site",
                                  "_study"), 
                                  drop = FALSE]
    row.names(subset) <- row.names(TCGA.sampletype)
    colnames(subset) <- c("sample.type", "primary.disease","primary.site","study")
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    subset,
                                    by    = "row.names",
                                    all.x = TRUE)
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)

    
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[ ,
      colnames(TCGA.RNAseq.anno) %in% c(EIF.gene, 
        "sample.type","primary.disease","primary.site","study")]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) + 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) -1)
    EIF.TCGA.RNAseq.anno.subset$`EIF4G1+EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4G1) + 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) -1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.gene, 
      "EIF4E+EIF4EBP1","EIF4G1+EIF4EBP1", "sample.type","primary.disease","primary.site","study")]
    #EIF.TCGA.RNAseq.anno.subset$delta <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4G1 - 2**EIF.TCGA.RNAseq.anno.subset$EIF4E - 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -1)
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long <- EIF.TCGA.RNAseq.anno.subset.long[
      !EIF.TCGA.RNAseq.anno.subset.long$value == 0, ]
    EIF.TCGA.RNAseq.anno.subset.long$sample.type <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$sample.type)
    EIF.TCGA.RNAseq.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.disease)
    return(EIF.TCGA.RNAseq.anno.subset.long)
  }
  TCGA.RNAseq.anno <- pan.TCGA.gene()
  
  pancancer.TCGA.EIF <- function(){
    TCGA.RNAseq.anno <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "TCGA", ]
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno [
      !TCGA.RNAseq.anno$sample.type %in% c("Solid Tissue Normal"), ]
            return(TCGA.RNAseq.anno.subset)}
  pancancer.TCGA.EIF.long <- pancancer.TCGA.EIF()
  
  GTEX.EIF <- function(){
    TCGA.RNAseq.anno <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "GTEX", ]
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[
      TCGA.RNAseq.anno$sample.type %in% "Normal Tissue", ]
    return(TCGA.RNAseq.anno.subset)}
  GTEX.EIF.long <- GTEX.EIF()
  
  # reorder bars by explicitly ordering factor levels
  make.plot <- function () {
    make.plot1 <- function(){
      pancancer.TCGA.EIF.long1 <- pancancer.TCGA.EIF.long[
        pancancer.TCGA.EIF.long$variable %in% EIF.gene, ]
      f1 <- factor(pancancer.TCGA.EIF.long1$primary.disease)
      f.ordered1 <- fct_rev(f1)
      p1 <- ggplot(data = pancancer.TCGA.EIF.long1,
        aes(x     = f.ordered1,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) + 
        scale_y_continuous(trans = log2_trans(), limits = c(2**4,2**17),
          labels = label_comma()) +
        #stat_n_text(size     = 5, 
        #  fontface = "bold", 
        #  hjust    = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        scale_color_manual(
          values = c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00"),
          breaks = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1"),
          labels = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        )
      print(p1)
      
      GTEX.EIF.long1 <- GTEX.EIF.long[
        GTEX.EIF.long$variable %in% EIF.gene, ]
      f2 <- factor(GTEX.EIF.long1$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = GTEX.EIF.long1,
        aes(x     = f.ordered2,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) +
        scale_y_continuous(trans = log2_trans(), limits = c(2**4,2**17),
          labels = label_comma()) +
        #stat_n_text(size     = 5, 
        #  fontface = "bold", 
        #  hjust    = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        scale_color_manual(
          values = c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00"),
          breaks = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1"),
          labels = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        )
      print(p2)
      
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFexpressionTCGAGTEX.pdf", 
        plot        = p,
        width       = 8, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot1()
    
    make.plot2 <- function(){
      pancancer.TCGA.EIF.long2 <- pancancer.TCGA.EIF.long[
        pancancer.TCGA.EIF.long$variable %in% c("EIF4E+EIF4EBP1","EIF4G1"), ]
      f1 <- factor(pancancer.TCGA.EIF.long2$primary.disease)
      f.ordered1 <- fct_rev(f1)
      p1 <- ggplot(data = pancancer.TCGA.EIF.long2,
        aes(
          x     = f.ordered1,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        stat_summary(
          fun.y = "median",
          geom  = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.5)) +
        scale_y_continuous(trans = log2_trans(), limits = c(2**9, 2**16),
          labels = label_comma()) +
        geom_boxplot(aes(colour = variable),
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.5)
        ) +
        scale_color_manual(
          values = c("#009E73","#CC79A7"),#for color-blind palettes
          breaks = c("EIF4E+EIF4EBP1", "EIF4G1"),
          labels = c("EIF4E+EIF4EBP1", "EIF4G1")) + 
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        ) 
      print(p1)
      
      GTEX.EIF.long2 <- GTEX.EIF.long[
        GTEX.EIF.long$variable %in% c("EIF4G1","EIF4E+EIF4EBP1"), ]
      f2 <- factor(GTEX.EIF.long2$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = GTEX.EIF.long2,
        aes(
          x     = f.ordered2,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        stat_summary(
          fun.y = "median",
          geom  = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.5)) +
        scale_y_continuous(trans = log2_trans(), limits = c(2**9, 2**16),
          labels = label_comma()) +
        geom_boxplot(aes(colour = variable),
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.5)
        ) +
        scale_color_manual(
          values = c("#009E73","#CC79A7"),#for color-blind palettes
          breaks = c("EIF4E+EIF4EBP1", "EIF4G1"),
          labels = c("EIF4E+EIF4EBP1", "EIF4G1")) + 
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        ) 
      print(p2)
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFsumexpressionTCGAGTEX.pdf", 
        plot        = p,
        width       = 7, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot2()
    
    make.plot3 <- function(){
      pancancer.TCGA.EIF.long1 <- pancancer.TCGA.EIF.long[
        pancancer.TCGA.EIF.long$variable %in% EIF.gene, ]
      f1 <- factor(pancancer.TCGA.EIF.long1$primary.disease)
      f.ordered1 <- fct_rev(f1)
      p1 <- ggplot(data = pancancer.TCGA.EIF.long1,
        aes(x     = f.ordered1,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) + 
        scale_y_continuous(trans = log2_trans(), limits = c(2**4,2**17),
          labels = label_comma()) +
        #stat_n_text(size     = 5, 
        #  fontface = "bold", 
        #  hjust    = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        scale_color_manual(
          values = c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00"),
          breaks = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1"),
          labels = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        )
      print(p1)
      
      pancancer.TCGA.EIF.long2 <- pancancer.TCGA.EIF.long[
        pancancer.TCGA.EIF.long$variable %in% c("EIF4E+EIF4EBP1","EIF4G1"), ]
      f2 <- factor(pancancer.TCGA.EIF.long2$primary.disease)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = pancancer.TCGA.EIF.long2,
        aes(
          x     = f.ordered2,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        stat_summary(
          fun.y = "median",
          geom  = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.5)) +
        scale_y_continuous(trans = log2_trans(), limits = c(2**4, 2**17),
          labels = label_comma()) +
        geom_boxplot(aes(colour = variable),
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.5)
        ) +
        scale_color_manual(
          values = c("#009E73","#CC79A7"),#for color-blind palettes
          breaks = c("EIF4E+EIF4EBP1", "EIF4G1"),
          labels = c("EIF4E+EIF4EBP1", "EIF4G1")) + 
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        ) 
      print(p2)
      
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFexprTCGA.pdf", 
        plot        = p,
        width       = 8, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot3()
    
    make.plot4 <- function(){
      GTEX.EIF.long1 <- GTEX.EIF.long[
        GTEX.EIF.long$variable %in% EIF.gene, ]
      f1 <- factor(GTEX.EIF.long1$primary.site)
      f.ordered1 <- fct_rev(f1)
      p1 <- ggplot(data = GTEX.EIF.long1,
        aes(x     = f.ordered1,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) +
        scale_y_continuous(trans = log2_trans(), limits = c(2**4,2**17),
          labels = label_comma()) +
        #stat_n_text(size     = 5, 
        #  fontface = "bold", 
        #  hjust    = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = .9)
        ) +
        scale_color_manual(
          values = c("#0072B2","#009E73","#D55E00","#CC79A7","#E69F00"),
          breaks = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1"),
          labels = c("EIF4E","EIF4EBP1","EIF4G1","EIF4A1")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        )
      print(p1)
      
      GTEX.EIF.long2 <- GTEX.EIF.long[
        GTEX.EIF.long$variable %in% c("EIF4G1","EIF4E+EIF4EBP1"), ]
      f2 <- factor(GTEX.EIF.long2$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = GTEX.EIF.long2,
        aes(
          x     = f.ordered2,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        stat_summary(
          fun.y = "median",
          geom  = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.5)) +
        scale_y_continuous(trans = log2_trans(), limits = c(2**4, 2**17),
          labels = label_comma()) +
        geom_boxplot(aes(colour = variable),
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.5)
        ) +
        scale_color_manual(
          values = c("#009E73","#CC79A7"),#for color-blind palettes
          breaks = c("EIF4E+EIF4EBP1", "EIF4G1"),
          labels = c("EIF4E+EIF4EBP1", "EIF4G1")) + 
        labs(x = "primary disease",
          y = paste("normalized RNA counts")) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.text          = black_bold_tahoma_12,
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          strip.text           = black_bold_tahoma_12
        ) 
      print(p2)
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFexprGTEX.pdf", 
        plot        = p,
        width       = 7, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot4()
  }
  make.plot()
  
}
plot.boxgraph.EIF.RNAseq.TCGA.GTEX(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1"))

##########################################
## boxplot for EIF ratios across tumors ##
##########################################
plot.boxgraph.EIF.ratio.TCGA <- function (EIF.gene) {
  pan.TCGA.gene <- function(){
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.pancancer <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
      TCGA.sampletype,
      by    = "row.names",
      all.x = TRUE)
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    sample.type.list <- levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    cancer.type.list <- levels(TCGA.RNAseq.anno$primary.disease)
    return(TCGA.RNAseq.anno)
  }
  TCGA.RNAseq.anno <- pan.TCGA.gene()
  
  pancancer.TCGA.EIF.ratio <- function () {
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[
      !TCGA.RNAseq.anno$sample.type %in% "Solid Tissue Normal", ]
    row.names(TCGA.RNAseq.anno.subset) <- TCGA.RNAseq.anno.subset$Row.names
    TCGA.RNAseq.anno.subset$Row.names <- NULL
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno.subset[ ,
      colnames(TCGA.RNAseq.anno.subset) %in% c(EIF.gene, 
                                               "sample.type",
                                               "primary.disease")]
    EIF.TCGA.RNAseq.anno.subset$sum <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - 2)
    EIF.TCGA.GTEX.score <- EIF.TCGA.RNAseq.anno.subset
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4G1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4EBP1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    #EIF.TCGA.GTEX.score$`PABPC1:\nEIF4E` <- 
    #  (EIF.TCGA.RNAseq.anno.subset$PABPC1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4E:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4E - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4EBP1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - EIF.TCGA.RNAseq.anno.subset$sum)
    ratio <- c("EIF4A1:\nEIF4E", "EIF4A1:\nEIF4EBP1", "EIF4G1:\nEIF4E", 
               "EIF4G1:\nEIF4EBP1", "EIF4EBP1:\nEIF4E", "EIF4A1:\nEIF4G1",
               "EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1",
               "EIF4E:\nEIF4E+EIF4EBP1", "EIF4EBP1:\nEIF4E+EIF4EBP1")
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.score[ ,c(ratio,
                                                   "sample.type",
                                                   "primary.disease"), 
                                              drop = FALSE]
    counts <- as.data.frame(table(EIF.TCGA.GTEX.score$primary.disease))
    colnames(counts) <- c("primary.disease","Freq")
    EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score)
    # reorder bars by explicitly ordering factor levels
    EIF.TCGA.GTEX.score.long$primary.disease <- as.factor(
      EIF.TCGA.GTEX.score.long$primary.disease)
    return(EIF.TCGA.GTEX.score.long)
    }
  pancancer.TCGA.EIF.ratio.long <- pancancer.TCGA.EIF.ratio()

  make.plot <- function () {
    make.plot1 <- function () {
      pancancer.TCGA.EIF.ratio.long$label <- sub(".*\n", "", 
        pancancer.TCGA.EIF.ratio.long$variable)
      pancancer.TCGA.EIF.ratio.long$label <- factor(
        pancancer.TCGA.EIF.ratio.long$label, 
        levels = c("EIF4G1","EIF4E+EIF4EBP1","EIF4E","EIF4EBP1"))  
      pancancer.TCGA.EIF.ratio.long1 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E",
                                                      "EIF4A1:\nEIF4G1", 
                                                      "EIF4G1:\nEIF4E", 
                                                      "EIF4EBP1:\nEIF4E"), ]
      pancancer.TCGA.EIF.ratio.long1$variable <- factor(
        pancancer.TCGA.EIF.ratio.long1$variable, 
        levels = c("EIF4A1:\nEIF4G1", 
                   "EIF4EBP1:\nEIF4E",
                   "EIF4A1:\nEIF4E",
                   "EIF4G1:\nEIF4E"))
      pancancer.TCGA.EIF.ratio.long1$value <- 2**(pancancer.TCGA.EIF.ratio.long1$value)

      f1 <- factor(pancancer.TCGA.EIF.ratio.long1$primary.disease)
      f.ordered1 <- fct_rev(f1)
      
      p1 <- ggplot(data = pancancer.TCGA.EIF.ratio.long1,
                   aes(x     = f.ordered1,
                       y     = value, 
                       fill  = variable,
                       color = variable)) +
            geom_boxplot(alpha    = .01,
                         outlier.shape = NA,
                         size     = .75,
                         width    = 1,
                         position = position_dodge2(preserve = "single")) + 
            geom_hline(yintercept = 1, linetype = "dashed") +
            scale_y_continuous(
              limits = quantile(pancancer.TCGA.EIF.ratio.long1$value, 
                                c(0, 0.996))) +
            scale_color_manual(values = c("#009E73","#CC79A7","#0072B2","#D55E00")) + #for color-blind palettes
            facet_wrap(~label, scales = "free_x") +
            facet_grid(~label, scales = "free_x", space = "free") +
            #facet_grid_sc(cols = vars(label), shrink = TRUE,
            #              scales = list(y = scales_y), 
            #              space = "free") +
            labs(x = "primary disease",
                 y = "ratio of RNA counts") +
            coord_flip() +
            theme_bw() +
            theme(plot.title           = black_bold_tahoma_12,
                  axis.title.x         = black_bold_tahoma_12,
                  axis.title.y         = element_blank(),
                  axis.text.x          = black_bold_tahoma_12,
                  axis.text.y          = black_bold_tahoma_12,
                  axis.line.x          = element_line(color = "black"),
                  axis.line.y          = element_line(color = "black"),
                  panel.grid           = element_blank(),
                  legend.title         = element_blank(),
                  legend.position      = "top",
                  legend.justification = "left",
                  legend.box           = "horizontal", 
                  legend.text          = black_bold_tahoma_12,
                  strip.background     = element_blank(),
                  strip.text.x         = element_blank())
      print(p1)
      ggsave(path        = "~/Documents/EIF_output/Expression", 
             filename    = "EIFratio.pdf", 
             plot        = p1,
             width       = 9, 
             height      = 9, 
             useDingbats = FALSE)
      }
    make.plot1()
  
    make.plot2 <- function () {
      pancancer.TCGA.EIF.ratio.long$label <- sub(".*\n", "", pancancer.TCGA.EIF.ratio.long$variable)
      pancancer.TCGA.EIF.ratio.long$label <- factor(
        pancancer.TCGA.EIF.ratio.long$label, 
        levels = c("EIF4G1","EIF4E+EIF4EBP1","EIF4E","EIF4EBP1"))
      pancancer.TCGA.EIF.ratio.long2 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      pancancer.TCGA.EIF.ratio.long2$variable <- factor(
        pancancer.TCGA.EIF.ratio.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
                   "EIF4A1:\nEIF4E+EIF4EBP1",
                   "EIF4G1:\nEIF4E",
                   "EIF4A1:\nEIF4E",
                   "EIF4G1:\nEIF4EBP1",
                   "EIF4A1:\nEIF4EBP1"))
      f2 <- factor(pancancer.TCGA.EIF.ratio.long2$primary.disease)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = pancancer.TCGA.EIF.ratio.long2,
        aes(x     = f.ordered2,
            y     = 2**value, 
            fill  = variable,
            color = variable)) + ylim(0, 101)+
    #scale_y_continuous(trans = log2_trans(), labels = label_number_auto()) +
    #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
                     alpha    = .01,
                     #size     = .75,
                     #width    = 1,
                     position = position_dodge(width = 0.9)
                    ) + 
        facet_wrap(~variable) +
        facet_grid(cols = vars(variable)) +
        stat_summary(fun.y = "median",
                     geom = 'line', # alpha = 0.5,
                     aes(group = variable, colour = variable),
                     position = position_dodge(width = 0.9)) +
    #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2",
                                      "#E69F00","#56B4E9","#D55E00")) + #for color-blind palettes
        labs(x = "primary disease",
             y = paste("ratio of RNA counts")) +
        coord_flip() +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
          ) + 
        guides(col = guide_legend(nrow = 1))
      print(p2)
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFsumratio.pdf", 
        plot        = p2,
        width       = 12.5, 
        height      = 9, 
        useDingbats = FALSE)
      }
    make.plot2()
    
    make.plot3 <- function () {
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1")] <- '1'
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4G1:\nEIF4E","EIF4G1:\nEIF4EBP1")] <- '2'
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4EBP1")] <- '3'
      pancancer.TCGA.EIF.ratio.long2 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      pancancer.TCGA.EIF.ratio.long2$variable <- factor(
        pancancer.TCGA.EIF.ratio.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
                   "EIF4A1:\nEIF4E+EIF4EBP1",
                   "EIF4G1:\nEIF4E",
                   "EIF4G1:\nEIF4EBP1",
                   "EIF4A1:\nEIF4E",
                   "EIF4A1:\nEIF4EBP1"))
      f2 <- factor(pancancer.TCGA.EIF.ratio.long2$primary.disease)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = pancancer.TCGA.EIF.ratio.long2,
        aes(
          x     = f.ordered2,
          y     = 2**value, 
          fill  = variable,
          color = variable)) + ylim(0, 101)+
        #scale_y_continuous(trans = log2_trans(), labels = label_number_auto()) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          position = position_dodge(width = 0.9)
        ) + 
        facet_wrap(~label) +
        facet_grid(cols = vars(label)) +
        stat_summary(
          fun.y = "median",
          geom = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.9)) +
        #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2"
          ,"#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip() +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
        ) + 
        guides(col = guide_legend(nrow = 1))
      print(p2)
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFsumratio2.pdf", 
        plot        = p2,
        width       = 12.5, 
        height      = 9, 
        useDingbats = FALSE)
    }
    make.plot3()
  }
  make.plot()
}
plot.boxgraph.EIF.ratio.TCGA(c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"))

plot.boxgraph.EIF.ratio.GTEX <- function (EIF.gene) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE)
    # check the name of the last column
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  
  get.EIFratio.anno.data <- function() {
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.gene, 
      "sample.type",
      "primary.site"),
      drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% "Normal Tissue", ]

    EIF.TCGA.RNAseq.anno.subset$sample.type <- as.factor(as.character(
      EIF.TCGA.RNAseq.anno.subset$sample.type))
    EIF.TCGA.RNAseq.anno.subset$sum <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score <- EIF.TCGA.RNAseq.anno.subset
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4G1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4EBP1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$sum)
    ratio <- c("EIF4A1:\nEIF4E", "EIF4G1:\nEIF4E", "EIF4EBP1:\nEIF4E", 
               "EIF4A1:\nEIF4G1", "EIF4G1:\nEIF4EBP1","EIF4A1:\nEIF4EBP1",
               "EIF4G1:\nEIF4E+EIF4EBP1",
               "EIF4A1:\nEIF4E+EIF4EBP1")
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.score[c(ratio, "sample.type","primary.site")]
    EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score)
    return(EIF.TCGA.GTEX.score.long)
  }
  EIF.TCGA.GTEX.score.long <- get.EIFratio.anno.data()
  
  make.plot <- function () {
    make.plot1 <- function () {
      EIF.TCGA.GTEX.score.long1 <- EIF.TCGA.GTEX.score.long[
      EIF.TCGA.GTEX.score.long$variable %in% c(
        "EIF4A1:\nEIF4E","EIF4G1:\nEIF4E","EIF4EBP1:\nEIF4E","EIF4A1:\nEIF4G1"), ]
    EIF.TCGA.GTEX.score.long1$variable = factor(
      EIF.TCGA.GTEX.score.long1$variable, 
      levels = c("EIF4A1:\nEIF4G1", 
        "EIF4EBP1:\nEIF4E",
        "EIF4A1:\nEIF4E",
        "EIF4G1:\nEIF4E"))
    EIF.TCGA.GTEX.score.long1$label <- sub(".*\n", "", 
      EIF.TCGA.GTEX.score.long1$variable)
    EIF.TCGA.GTEX.score.long1$label <- factor(
      EIF.TCGA.GTEX.score.long1$label, 
      levels = c("EIF4G1","EIF4E"))  
    f1 <- factor(EIF.TCGA.GTEX.score.long1$primary.site)
    f.ordered1 <- fct_rev(f1)
    p1 <- ggplot(data = EIF.TCGA.GTEX.score.long1,
      aes(x     = f.ordered1,  
        #x     = x.ordered, # order primary disease
        y     = 2**value,
        color = variable)) + ylim(0,100)+
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_wrap(~label, scales = "free_x") +
      facet_grid(~label, scales = "free_x", space = "free") +
      stat_n_text(size     = 5, 
        fontface = "bold", 
        hjust    = 0) +
      geom_boxplot(alpha    = .01,
        outlier.shape = NA,
        size     = .75,
        width    = 1,
        position = position_dodge2(preserve = "single")) + 
      scale_color_manual(
        values = c("#009E73","#CC79A7","#0072B2","#D55E00")) + #for color-blind palettes
      labs(x = "primary disease",
        y = paste("normalized RNA counts")) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_12,
        axis.title.x         = black_bold_tahoma_12,
        axis.title.y         = element_blank(),
        axis.text.x          = black_bold_tahoma_12,
        axis.text.y          = black_bold_tahoma_12,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        legend.title         = element_blank(),
        legend.position      = "top",
        legend.justification = "left",
        legend.box           = "horizontal", 
        legend.text          = black_bold_tahoma_12,
        strip.background     = element_blank(),
        strip.text.x         = element_blank())
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/Expression/GTEX", 
      filename    = "EIFratioGTEX.pdf", 
      plot        = p1,
      width       = 12.5, 
      height      = 8, 
      useDingbats = FALSE)
    }
    make.plot1()
    
    make.plot2 <- function () {
      EIF.TCGA.GTEX.score.long2 <- EIF.TCGA.GTEX.score.long[
        EIF.TCGA.GTEX.score.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      EIF.TCGA.GTEX.score.long2$variable <- factor(
        EIF.TCGA.GTEX.score.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4EBP1",
          "EIF4A1:\nEIF4EBP1"))
      EIF.TCGA.GTEX.score.long2$label <- sub(".*\n", "", EIF.TCGA.GTEX.score.long2$variable)
      EIF.TCGA.GTEX.score.long2$label <- factor(
        EIF.TCGA.GTEX.score.long2$label, 
        levels = c("EIF4G1","EIF4E+EIF4EBP1","EIF4E","EIF4EBP1"))
      f2 <- factor(EIF.TCGA.GTEX.score.long2$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = EIF.TCGA.GTEX.score.long2,
        aes(x     = f.ordered2,
          y     = 2**value, 
          fill  = variable,
          color = variable)) + ylim(0, 101)+
        #scale_y_continuous(trans = log2_trans(), labels = label_number_auto()) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          position = position_dodge(width = 0.9)
        ) + 
        facet_wrap(~variable) +
        facet_grid(cols = vars(variable)) +
        stat_summary(fun.y = "median",
          geom = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.9)) +
        #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2",
          "#E69F00","#56B4E9","#D55E00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip() +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
        ) + 
        guides(col = guide_legend(nrow = 1))
      print(p2)
      ggsave(
        path        = "~/Documents/EIF_output/Expression/GTEX", 
        filename    = "EIFsumratio2GTEX.pdf", 
        plot        = p2,
        width       = 12.5, 
        height      = 8, 
        useDingbats = FALSE)
    }
    make.plot2()
    
    make.plot3 <- function () {
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1")] <- '1'
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4G1:\nEIF4E","EIF4G1:\nEIF4EBP1")] <- '2'
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4EBP1")] <- '3'
      EIF.TCGA.GTEX.score.long2 <- EIF.TCGA.GTEX.score.long[
        EIF.TCGA.GTEX.score.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      EIF.TCGA.GTEX.score.long2$variable <- factor(
        EIF.TCGA.GTEX.score.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4G1:\nEIF4E",
          "EIF4G1:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1"))
      f2 <- factor(EIF.TCGA.GTEX.score.long2$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = EIF.TCGA.GTEX.score.long2,
        aes(
          x     = f.ordered2,
          y     = 2**value, 
          fill  = variable,
          color = variable)) + ylim(0, 101)+
        #scale_y_continuous(trans = log2_trans(), labels = label_number_auto()) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          position = position_dodge(width = 0.9)
        ) + 
        facet_wrap(~label) +
        facet_grid(cols = vars(label)) +
        stat_summary(
          fun.y = "median",
          geom = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.9)) +
        #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2"
          ,"#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip() +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
        ) + 
        guides(col = guide_legend(nrow = 1))
      print(p2)
      ggsave(
        path        = "~/Documents/EIF_output/Expression/GTEX", 
        filename    = "EIFsumratioGTEX.pdf", 
        plot        = p2,
        width       = 12.5, 
        height      = 8, 
        useDingbats = FALSE)
    }
    make.plot3()
    }
  make.plot()
  }
plot.boxgraph.EIF.ratio.GTEX (c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"))

plot.boxgraph.EIF.ratio.TCGA.GTEX <- function (EIF.gene) {
  pan.TCGA.gene <- function(){
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.pancancer <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count",
      data.table = FALSE)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
  

    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    subset <- TCGA.sampletype[ ,c("_sample_type", 
                                  "primary disease or tissue", 
                                  "_primary_site", "_study"), 
                                  drop = FALSE]
    row.names(subset) <- row.names(TCGA.sampletype)
    colnames(subset) <- c("sample.type", "primary.disease","primary.site","study")
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    subset,
                                    by    = "row.names",
                                    all.x = TRUE)
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno <- na.omit(TCGA.RNAseq.anno)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    levels(TCGA.RNAseq.anno$primary.disease)
    TCGA.RNAseq.anno$primary.site <- as.factor(TCGA.RNAseq.anno$primary.site)
    levels(TCGA.RNAseq.anno$primary.site)
    row.names(TCGA.RNAseq.anno) <- TCGA.RNAseq.anno$Row.names
    TCGA.RNAseq.anno$Row.names <- NULL
    
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[ ,
      colnames(TCGA.RNAseq.anno) %in% c(EIF.gene, 
                                        "sample.type",
                                        "primary.disease",
                                        "primary.site",
                                        "study")]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    EIF.TCGA.RNAseq.anno.subset$sum <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - 2)
    EIF.TCGA.GTEX.score <- EIF.TCGA.RNAseq.anno.subset
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4G1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4E:\nEIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4E - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    #EIF.TCGA.GTEX.score$`PABPC1:\nEIF4E` <- 
    #  (EIF.TCGA.RNAseq.anno.subset$PABPC1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4A1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4E:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4E - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4EBP1:\nEIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - EIF.TCGA.RNAseq.anno.subset$sum)
    ratio <- c("EIF4A1:\nEIF4E", "EIF4A1:\nEIF4EBP1", "EIF4G1:\nEIF4E", 
      "EIF4G1:\nEIF4EBP1", "EIF4E:\nEIF4EBP1", "EIF4A1:\nEIF4G1",
      "EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1",
      "EIF4E:\nEIF4E+EIF4EBP1", "EIF4EBP1:\nEIF4E+EIF4EBP1")
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.score[ ,c(ratio,
      "sample.type","primary.disease","primary.site","study"), 
      drop = FALSE]
    return(EIF.TCGA.GTEX.score)
  }
  TCGA.RNAseq.anno <- pan.TCGA.gene()
  
  pancancer.TCGA.EIF.ratio <- function () {
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "TCGA", ]
    #TCGA.RNAseq.anno.subset$study <- NULL
    TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno.subset[
      !TCGA.RNAseq.anno.subset$sample.type %in% "Solid Tissue Normal", ]
    EIF.TCGA.GTEX.score.long <- melt(TCGA.RNAseq.anno.subset)
    # reorder bars by explicitly ordering factor levels
    EIF.TCGA.GTEX.score.long$primary.disease <- as.factor(
      EIF.TCGA.GTEX.score.long$primary.disease)
    EIF.TCGA.GTEX.score.long <- droplevels(EIF.TCGA.GTEX.score.long)
    levels(EIF.TCGA.GTEX.score.long$primary.disease)
    return(EIF.TCGA.GTEX.score.long)
  }
  pancancer.TCGA.EIF.ratio.long <- pancancer.TCGA.EIF.ratio()

  tissue.GTEX.EIF.ratio <- function() {
    EIF.TCGA.RNAseq.anno.subset <- TCGA.RNAseq.anno[TCGA.RNAseq.anno$study == "GTEX", ]
    #EIF.TCGA.RNAseq.anno.subset$study <- NULL
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% "Normal Tissue", ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- as.factor(as.character(
      EIF.TCGA.RNAseq.anno.subset$sample.type))
    EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.GTEX.score.long)
  }
  EIF.TCGA.GTEX.score.long <- tissue.GTEX.EIF.ratio()
  
  make.plot <- function () {
    make.plot1 <- function () {
      pancancer.TCGA.EIF.ratio.long$label <- sub(".*\n", "", 
        pancancer.TCGA.EIF.ratio.long$variable)
      pancancer.TCGA.EIF.ratio.long$label <- factor(
        pancancer.TCGA.EIF.ratio.long$label, 
        levels = c("EIF4G1","EIF4EBP1","EIF4E+EIF4EBP1","EIF4E"))  
      pancancer.TCGA.EIF.ratio.long1 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4A1:\nEIF4G1", 
          "EIF4E:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      pancancer.TCGA.EIF.ratio.long1$variable <- factor(
        pancancer.TCGA.EIF.ratio.long1$variable, 
        levels = c(
          "EIF4A1:\nEIF4G1", 
          "EIF4E:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E"))
      pancancer.TCGA.EIF.ratio.long1$hline[pancancer.TCGA.EIF.ratio.long1$label %in% c("EIF4G1","EIF4EBP1")] <- 1
      pancancer.TCGA.EIF.ratio.long1$hline[pancancer.TCGA.EIF.ratio.long1$label %in% c("EIF4E","EIF4E+EIF4EBP1")] <- 4
      #pancancer.TCGA.EIF.ratio.long1$value <- 2**(pancancer.TCGA.EIF.ratio.long1$value)
      pancancer.TCGA.EIF.ratio.long1 <- na.omit(pancancer.TCGA.EIF.ratio.long1)
      
      f1 <- factor(pancancer.TCGA.EIF.ratio.long1$primary.disease)
      f.ordered1 <- fct_rev(f1)
      
      p1 <- ggplot(data = pancancer.TCGA.EIF.ratio.long1,
        aes(x     = f.ordered1,
          y     = 2**value, 
          #fill  = variable,
          color = variable)) +  
        geom_boxplot(alpha    = .1,
          outlier.shape = NA,
          size     = .75,
          width    = 1,
          position = position_dodge2(preserve = "single")) + 
        geom_hline(data = pancancer.TCGA.EIF.ratio.long1, 
                   aes(yintercept = hline), 
                   linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2"
          ,"#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        facet_wrap(~label, scales = "free_x") +
        facet_grid(~label, scales = "free_x", space = "free") +
        #facet_grid_sc(cols = vars(label), shrink = TRUE,
        #              scales = list(y = scales_y), 
        #              space = "free") +
        guides(colour = guide_legend(nrow = 1))+
        labs(x = "primary disease",
             y = "ratio of RNA counts") +
        coord_flip(ylim = c(0, 25)) +
        theme_bw() +
        theme(plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          legend.text          = black_bold_tahoma_12,
          strip.background     = element_blank(),
          strip.text.x         = element_blank())
      print(p1)
      
      EIF.TCGA.GTEX.score.long1 <- EIF.TCGA.GTEX.score.long[
        EIF.TCGA.GTEX.score.long$variable %in% c(
          "EIF4A1:\nEIF4G1", 
          "EIF4E:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      EIF.TCGA.GTEX.score.long1$variable = factor(
        EIF.TCGA.GTEX.score.long1$variable, 
        levels = c(
          "EIF4A1:\nEIF4G1", 
          "EIF4E:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E"))
      EIF.TCGA.GTEX.score.long1$label <- sub(".*\n", "", 
        EIF.TCGA.GTEX.score.long1$variable)
      EIF.TCGA.GTEX.score.long1$label <- factor(
        EIF.TCGA.GTEX.score.long1$label, 
        levels = c("EIF4G1","EIF4EBP1","EIF4E+EIF4EBP1","EIF4E")) 
      EIF.TCGA.GTEX.score.long1$hline[EIF.TCGA.GTEX.score.long1$label %in% c("EIF4G1","EIF4EBP1")] <- 1
      EIF.TCGA.GTEX.score.long1$hline[EIF.TCGA.GTEX.score.long1$label %in% c("EIF4E","EIF4E+EIF4EBP1")] <- 4
      
      
      f2 <- factor(EIF.TCGA.GTEX.score.long1$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = EIF.TCGA.GTEX.score.long1,
        aes(x     = f.ordered2,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) + 
        geom_hline(data = EIF.TCGA.GTEX.score.long1, 
          aes(yintercept = hline), 
          linetype = "dashed") +        facet_wrap(~label, scales = "free_x") +
        facet_grid(~label, scales = "free_x", space = "free") +
        guides(colour = guide_legend(nrow = 1))+
        stat_n_text(size     = 5, 
          fontface = "bold", 
          hjust    = 0) +
        geom_boxplot(alpha    = .01,
          outlier.shape = NA,
          size     = .75,
          width    = 1,
          position = position_dodge2(preserve = "single")) + 
        scale_color_manual(
          values = c("#009E73","#CC79A7","#0072B2",
                     "#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip(ylim = c(0, 25)) +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          legend.text          = black_bold_tahoma_12,
          strip.background     = element_blank(),
          strip.text.x         = element_blank())
      print(p2)
      
  
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      
      ggsave(path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFratiozoom.pdf", 
        plot        = p,
        width       = 15, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot1()
    
    make.plot2 <- function () {
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1")] <- '1'
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4G1:\nEIF4E","EIF4G1:\nEIF4EBP1")] <- '2'
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4EBP1")] <- '3'
      pancancer.TCGA.EIF.ratio.long2 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      pancancer.TCGA.EIF.ratio.long2$variable <- factor(
        pancancer.TCGA.EIF.ratio.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4G1:\nEIF4E",
          "EIF4G1:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1"))
      pancancer.TCGA.EIF.ratio.long2 <- na.omit(pancancer.TCGA.EIF.ratio.long2)
      f1 <- factor(pancancer.TCGA.EIF.ratio.long2$primary.disease)
      f.ordered1 <- fct_rev(f1)
      p1 <- ggplot(data = pancancer.TCGA.EIF.ratio.long2,
        aes(
          x     = f.ordered1,
          y     = 2**value, 
          fill  = variable,
          color = variable)) + #ylim(0, 100)+
        #scale_y_continuous(trans = log2_trans(), labels = label_number_auto()) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.9)
        ) + 
        facet_wrap(~label) +
        facet_grid(cols = vars(label)) +
        stat_summary(
          fun.y = "median",
          geom = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.9)) +
        #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2"
          ,"#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip(ylim = c(0, 25)) +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
        ) + 
        guides(col = guide_legend(nrow = 1))
      print(p1)
      
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1")] <- '1'
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4G1:\nEIF4E","EIF4G1:\nEIF4EBP1")] <- '2'
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4EBP1")] <- '3'
      EIF.TCGA.GTEX.score.long2 <- EIF.TCGA.GTEX.score.long[
        EIF.TCGA.GTEX.score.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      EIF.TCGA.GTEX.score.long2$variable <- factor(
        EIF.TCGA.GTEX.score.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4G1:\nEIF4E",
          "EIF4G1:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1"))
      EIF.TCGA.GTEX.score.long2 <- na.omit(EIF.TCGA.GTEX.score.long2)
      f2 <- factor(EIF.TCGA.GTEX.score.long2$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = EIF.TCGA.GTEX.score.long2,
        aes(
          x     = f.ordered2,
          y     = 2**value, 
          fill  = variable,
          color = variable)) + #ylim(0, 100)+
        #scale_y_continuous(trans = log2_trans(), labels = label_number_auto()) +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.9)
        ) + 
        facet_wrap(~label) +
        facet_grid(cols = vars(label)) +
        stat_summary(
          fun.y = "median",
          geom = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.9)) +
        #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2"
          ,"#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip(ylim = c(0, 25)) +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
        ) + 
        guides(col = guide_legend(nrow = 1))
      print(p2)
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFsumratio2.pdf", 
        plot        = p,
        width       = 12.5, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot2()
    
    make.plot3 <- function () {
      pancancer.TCGA.EIF.ratio.long1 <- pancancer.TCGA.EIF.ratio.long
      pancancer.TCGA.EIF.ratio.long1$label <- sub(".*\n", "", 
        pancancer.TCGA.EIF.ratio.long1$variable)
      pancancer.TCGA.EIF.ratio.long1$label <- factor(
        pancancer.TCGA.EIF.ratio.long1$label, 
        levels = c("EIF4G1","EIF4EBP1","EIF4E"))  
      pancancer.TCGA.EIF.ratio.long1 <- pancancer.TCGA.EIF.ratio.long1[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E",
          "EIF4A1:\nEIF4G1", 
          "EIF4G1:\nEIF4E", 
          "EIF4E:\nEIF4EBP1"), ]
      pancancer.TCGA.EIF.ratio.long1$variable <- factor(
        pancancer.TCGA.EIF.ratio.long1$variable, 
        levels = c("EIF4A1:\nEIF4G1", 
          "EIF4E:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E"))
      #pancancer.TCGA.EIF.ratio.long1$value <- 2**(pancancer.TCGA.EIF.ratio.long1$value)
      pancancer.TCGA.EIF.ratio.long1 <- na.omit(pancancer.TCGA.EIF.ratio.long1)
      
      f1 <- factor(pancancer.TCGA.EIF.ratio.long1$primary.disease)
      f.ordered1 <- fct_rev(f1)
      
      p1 <- ggplot(data = pancancer.TCGA.EIF.ratio.long1,
        aes(x     = f.ordered1,
          y     = 2**value, 
          fill  = variable,
          color = variable)) +  
        geom_boxplot(alpha    = .1,
          outlier.shape = NA,
          size     = .75,
          width    = 1,
          position = position_dodge2(preserve = "single")) + 
        geom_hline(yintercept = 1, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2","#D55E00")) + #for color-blind palettes
        facet_wrap(~label, scales = "free_x") +
        facet_grid(~label, scales = "free_x", space = "free") +
        #facet_grid_sc(cols = vars(label), shrink = TRUE,
        #              scales = list(y = scales_y), 
        #              space = "free") +
        labs(x = "primary disease",
          y = "ratio of RNA counts") +
        coord_flip(ylim = c(0, 100)) +
        theme_bw() +
        theme(plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          legend.text          = black_bold_tahoma_12,
          strip.background     = element_blank(),
          strip.text.x         = element_blank())
      print(p1)
      
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1")] <- '1'
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4G1:\nEIF4E","EIF4G1:\nEIF4EBP1")] <- '2'
      pancancer.TCGA.EIF.ratio.long$label[
        pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4EBP1")] <- '3'
      pancancer.TCGA.EIF.ratio.long2 <- pancancer.TCGA.EIF.ratio.long[
        pancancer.TCGA.EIF.ratio.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      pancancer.TCGA.EIF.ratio.long2$variable <- factor(
        pancancer.TCGA.EIF.ratio.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4G1:\nEIF4E",
          "EIF4G1:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1"))
      pancancer.TCGA.EIF.ratio.long2 <- na.omit(pancancer.TCGA.EIF.ratio.long2)
      f2 <- factor(pancancer.TCGA.EIF.ratio.long2$primary.disease)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = pancancer.TCGA.EIF.ratio.long2,
        aes(
          x     = f.ordered2,
          y     = 2**value, 
          fill  = variable,
          color = variable)) + #ylim(0, 100)+
        geom_hline(yintercept = 1, linetype = "dashed") +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.9)
        ) + 
        facet_wrap(~label) +
        facet_grid(cols = vars(label)) +
        stat_summary(
          fun.y = "median",
          geom = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.9)) +
        #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2"
          ,"#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip(ylim = c(0, 100)) +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
        ) + 
        guides(col = guide_legend(nrow = 1))
      print(p2)
      
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      
      ggsave(path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFratiosumTCGA.pdf", 
        plot        = p,
        width       = 11, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot3()
    
    make.plot4 <- function () {
      EIF.TCGA.GTEX.score.long1 <- EIF.TCGA.GTEX.score.long[
        EIF.TCGA.GTEX.score.long$variable %in% c(
          "EIF4A1:\nEIF4E","EIF4G1:\nEIF4E","EIF4E:\nEIF4EBP1","EIF4A1:\nEIF4G1"), ]
      EIF.TCGA.GTEX.score.long1$variable = factor(
        EIF.TCGA.GTEX.score.long1$variable, 
        levels = c("EIF4A1:\nEIF4G1", 
          "EIF4E:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E"))
      EIF.TCGA.GTEX.score.long1$label <- sub(".*\n", "", 
        EIF.TCGA.GTEX.score.long1$variable)
      EIF.TCGA.GTEX.score.long1$label <- factor(
        EIF.TCGA.GTEX.score.long1$label, 
        levels = c("EIF4G1","EIF4EBP1","EIF4E"))  
      f1 <- factor(EIF.TCGA.GTEX.score.long1$primary.site)
      f.ordered1 <- fct_rev(f1)
      p1 <- ggplot(data = EIF.TCGA.GTEX.score.long1,
        aes(x     = f.ordered1,  
          #x     = x.ordered, # order primary disease
          y     = 2**value,
          color = variable)) + 
        geom_hline(yintercept = 1, linetype = "dashed") +
        facet_wrap(~label, scales = "free_x") +
        facet_grid(~label, scales = "free_x", space = "free") +
        geom_boxplot(alpha    = .01,
          outlier.shape = NA,
          size     = .75,
          width    = 1,
          position = position_dodge2(preserve = "single")) + 
        scale_color_manual(
          values = c("#009E73","#CC79A7","#0072B2","#D55E00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip(ylim = c(0, 100)) +
        theme_bw() +
        theme(
          plot.title           = black_bold_tahoma_12,
          axis.title.x         = black_bold_tahoma_12,
          axis.title.y         = element_blank(),
          axis.text.x          = black_bold_tahoma_12,
          axis.text.y          = black_bold_tahoma_12,
          axis.line.x          = element_line(color = "black"),
          axis.line.y          = element_line(color = "black"),
          panel.grid           = element_blank(),
          legend.title         = element_blank(),
          legend.position      = "top",
          legend.justification = "left",
          legend.box           = "horizontal", 
          legend.text          = black_bold_tahoma_12,
          strip.background     = element_blank(),
          strip.text.x         = element_blank())
      print(p1)
      
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1")] <- '1'
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4G1:\nEIF4E","EIF4G1:\nEIF4EBP1")] <- '2'
      EIF.TCGA.GTEX.score.long$label[
        EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4EBP1")] <- '3'
      EIF.TCGA.GTEX.score.long2 <- EIF.TCGA.GTEX.score.long[
        EIF.TCGA.GTEX.score.long$variable %in% c(
          "EIF4A1:\nEIF4E",
          "EIF4G1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1",
          "EIF4G1:\nEIF4EBP1",
          "EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1"), ]
      EIF.TCGA.GTEX.score.long2$variable <- factor(
        EIF.TCGA.GTEX.score.long2$variable, 
        levels = c("EIF4G1:\nEIF4E+EIF4EBP1",
          "EIF4A1:\nEIF4E+EIF4EBP1",
          "EIF4G1:\nEIF4E",
          "EIF4G1:\nEIF4EBP1",
          "EIF4A1:\nEIF4E",
          "EIF4A1:\nEIF4EBP1"))
      EIF.TCGA.GTEX.score.long2 <- na.omit(EIF.TCGA.GTEX.score.long2)
      f2 <- factor(EIF.TCGA.GTEX.score.long2$primary.site)
      f.ordered2 <- fct_rev(f2)
      p2 <- ggplot(data = EIF.TCGA.GTEX.score.long2,
        aes(
          x     = f.ordered2,
          y     = 2**value, 
          fill  = variable,
          color = variable)) + #ylim(0, 100)+
        geom_hline(yintercept = 1, linetype = "dashed") +
        #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
        geom_boxplot(
          alpha    = .01,
          #size     = .75,
          #width    = 1,
          outlier.shape = NA,
          position = position_dodge(width = 0.9)
        ) + 
        facet_wrap(~label) +
        facet_grid(cols = vars(label)) +
        stat_summary(
          fun.y = "median",
          geom = 'line', # alpha = 0.5,
          aes(group = variable, colour = variable),
          position = position_dodge(width = 0.9)) +
        #geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#009E73","#CC79A7","#0072B2"
          ,"#D55E00","#56B4E9","#E69F00")) + #for color-blind palettes
        labs(x = "primary disease",
          y = paste("ratio of RNA counts")) +
        coord_flip(ylim = c(0, 100)) +
        theme_bw() +   
        theme(
          plot.title       = black_bold_tahoma_12,
          axis.title.x     = black_bold_tahoma_12,
          axis.title.y     = element_blank(),
          axis.text.x      = black_bold_tahoma_12,
          axis.text.y      = black_bold_tahoma_12,
          axis.line.x      = element_line(color = "black"),
          axis.line.y      = element_line(color = "black"),
          panel.grid       = element_blank(),
          legend.title     = element_blank(),
          legend.text      = black_bold_tahoma_12,
          legend.position  = "top",
          strip.background = element_blank(),
          strip.text.x     = element_blank()
          #strip.text           = black_bold_tahoma_12
        ) + 
        guides(col = guide_legend(nrow = 1))
      print(p2)
      g1grob <- ggplotGrob(p1)
      g2grob <- ggplotGrob(p2)
      g2grob$widths <- g1grob$widths
      grid.arrange(g1grob, g2grob)
      p <-  arrangeGrob(g1grob, g2grob) #generates g
      print(p)
      
      ggsave(
        path        = "~/Documents/EIF_output/Expression", 
        filename    = "EIFratiosumGTEX.pdf", 
        plot        = p,
        width       = 10, 
        height      = 16, 
        useDingbats = FALSE)
    }
    make.plot4()
    
  }
  make.plot()
}
plot.boxgraph.EIF.ratio.TCGA.GTEX(c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"))

#################################################################
## violin plot for EIF expression in tumors vs adjacent normal ##
#################################################################
plot.violingraph.EIF.RNAseq.TCGA <- function (EIF.gene) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), 
                               drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.t,
                                       subset,
                                       by    = "row.names",
                                       all.x = TRUE)
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
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno[ ,c(EIF.gene, 
                                                            "sample.type",
                                                            "primary.site"),
                                                        drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Metastatic",
                                                     "Primary Tumor",
                                                     #"Normal Tissue",
                                                     "Solid Tissue Normal"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- as.factor(as.character(
      EIF.TCGA.RNAseq.anno.subset$sample.type))
    EIF.TCGA.RNAseq.anno.subset.long <- melt(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset.long$primary.site <- as.factor(
      EIF.TCGA.RNAseq.anno.subset.long$primary.site)
    return(EIF.TCGA.RNAseq.anno.subset.long)}
  EIF.TCGA.RNAseq.anno.subset.long <- get.subset.data()

  make.plot <- function (){
    EIF.TCGA.RNAseq.anno.subset.long1 <- EIF.TCGA.RNAseq.anno.subset.long[
      EIF.TCGA.RNAseq.anno.subset.long$variable %in% c("EIF4E",
                                                       "EIF4G1",
                                                       "EIF4A1",
                                                       "EIF4EBP1",
                                                       "EIF4E+EIF4EBP1"), ]
    EIF.TCGA.RNAseq.anno.subset.long1$variable <- factor(
      EIF.TCGA.RNAseq.anno.subset.long1$variable, 
      levels = c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1","EIF4E+EIF4EBP1"))
    p1 <- ggplot(data = EIF.TCGA.RNAseq.anno.subset.long1,
                 aes(x     = sample.type,
                     y     = 2**value, 
                     color = sample.type,
                     fill  = sample.type
                    )) +
          stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
          facet_grid(. ~ variable,
                     scales = "free",
                     space  = "free") +
          # facet_wrap(~ variable, ncol = 6) +
          geom_violin(trim = TRUE) + 
          #scale_color_manual(values = c("#D55E00","#0072B2","#CC79A7","#009E73")) + #for color-blind palettes
          geom_boxplot(
            alpha    = .01,
            width    = .25,
            color    = "black",
            position = position_dodge(width = .9)
            ) +
          labs(x = "sample type",
               y = "normalized RNA counts") +
          scale_x_discrete(labels = c("Metastatic Tumor",
                                      "Primary Tumor",
                                      #"Normal Tissue",
                                      "Adjacent Normal")
            ) + 
          scale_y_continuous(trans    = log2_trans(), 
                             labels   = label_comma(),
            breaks = c(1,128,2048,32768),
            #labels = c("1","8","64","512"),
                             position = 'right') +
          scale_color_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          scale_fill_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          theme_bw() +
          theme(
            plot.title         = black_bold_tahoma_16,
            axis.title.x       = element_blank(),
            axis.title.y.right = black_bold_tahoma_16_right,
            axis.text.x        = black_bold_tahoma_16_90,
            axis.text.y        = black_bold_tahoma_16_90,
            axis.line.x        = element_line(color = "black"),
            axis.line.y        = element_line(color = "black"),
            panel.grid         = element_blank(),
            legend.position    = "none",
            strip.text         = black_bold_tahoma_16
          ) +
          stat_compare_means(
            comparisons = list(c("Metastatic",    "Solid Tissue Normal"),
                               c("Primary Tumor", "Solid Tissue Normal"),
                               c("Metastatic",    "Primary Tumor")),
            method      = "t.test", 
            label       = "p.signif", 
            size        = 6)
        print(p1)
        ggsave(path        = "~/Documents/EIF_output/Expression", 
               filename    = "EIFexpressionviolin.pdf", 
               plot        = p1,
               width       = 7.5, 
               height      = 9, 
               useDingbats = FALSE)
    
    EIF.TCGA.RNAseq.anno.subset.long2 <- EIF.TCGA.RNAseq.anno.subset.long[
      EIF.TCGA.RNAseq.anno.subset.long$variable %in% c("EIF4E+EIF4EBP1"), ]
    #EIF.TCGA.RNAseq.anno.subset.long2$variable <- factor(EIF.TCGA.RNAseq.anno.subset.long2$variable, levels = c("EIF4G1","EIF4E+EIF4EBP1"))
    p2 <- ggplot(data = EIF.TCGA.RNAseq.anno.subset.long2,
                 aes(x     = sample.type,
                     y     = 2**value, 
                     color = sample.type,
                     fill  = sample.type
                    )) +
          stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
          geom_violin(trim = TRUE) + 
          #scale_color_manual(values = c("#D55E00","#0072B2","#CC79A7","#009E73")) + #for color-blind palettes
          geom_boxplot(
            alpha    = .01,
            width    = .25,
            color    = "black",
            position = position_dodge(width = .9)
          ) +
          facet_grid(. ~ variable,
            scales = "free",
            space  = "free") +
          labs(x = "sample type",
               y = "normalized RNA counts") +
          scale_x_discrete(labels = c("Metastatic Tumor",
                                      "Primary Tumor",
                                      "Adjacent Normal")
                                    ) + 
          scale_y_continuous(trans    = log2_trans(), 
            labels   = label_comma(),
            breaks = c(1,512,2048,8192,32768),
            #labels = c("1","8","64","512"),
            position = 'right') +
          scale_color_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          scale_fill_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          theme_bw() +
          theme(plot.title         = black_bold_tahoma_16,
                axis.title.x       = element_blank(),
                axis.title.y.right = black_bold_tahoma_16_right,
                axis.text.x        = black_bold_tahoma_16_90,
                axis.text.y        = black_bold_tahoma_16_90,
                axis.line.x        = element_line(color = "black"),
                axis.line.y        = element_line(color = "black"),
                panel.grid         = element_blank(),
                legend.position    = "none",
                strip.text         = black_bold_tahoma_16
            ) +
          stat_compare_means(
            comparisons = list(c("Metastatic",    "Solid Tissue Normal"),
                               c("Primary Tumor", "Solid Tissue Normal"),
                               c("Metastatic",    "Primary Tumor")),
            method = "t.test", label = "p.signif", size = 6)
    print(p2)
    ggsave(path        = "~/Documents/EIF_output/Expression", 
           filename    = "sumviolin.pdf", 
           plot        = p2,
           width       = 1.8, 
           height      = 9, 
           useDingbats = FALSE)
    
    p3 <- ggplot(data = EIF.TCGA.RNAseq.anno.subset.long,
                 aes(x     = variable,
                     y     = value,
                     fill  = variable)) +
          stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
          facet_grid(~ sample.type,
            scales = "free",
            space  = "free") +
          facet_wrap(~ sample.type, ncol = 6) +
          geom_violin(trim = FALSE) +
          geom_boxplot(
            alpha    = .01,
            width    = .25,
            color    = "black",
            position = position_dodge(width = .9)
          ) +
      labs(x = "EIF4F complex",
           y = paste("log2(normalized RNA counts)")) +
      scale_y_continuous(position = 'right') +
      #scale_fill_manual(values = c("#E69F00","#56B4E9","#009E73","#CC79A7")) +
      theme_bw() +
      theme(
        plot.title         = black_bold_tahoma_16,
        axis.title         = black_bold_tahoma_16,
        axis.title.y.right = black_bold_tahoma_16_right,
        axis.text.x        = black_bold_tahoma_16_90,
        axis.text.y        = black_bold_tahoma_16_90,
        axis.line.x        = element_line(color = "black"),
        axis.line.y        = element_line(color = "black"),
        panel.grid         = element_blank(),
        legend.position    = "none",
        strip.text         = black_bold_tahoma_16
      ) +
      stat_compare_means(comparisons = list(c("EIF4A1", "EIF4E"),
                                            c("EIF4G1", "EIF4E"),
                                            c("EIF4EBP1", "EIF4E")),
                         method      = "t.test", 
                         label       = "p.signif", 
                         size        = 6, 
                         hjust       = 0)
    print(p3)
    }
  make.plot()
}
plot.violingraph.EIF.RNAseq.TCGA (c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"))

############################################################
## violin plot for EIF ratio in tumors vs adjacent normal ##
############################################################
plot.violingraph.EIF.ratio.TCGA <- function (EIF.gene) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE)
    # check the name of the last column
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  
  get.EIFratio.anno.data <- function() {
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.gene, 
                                                            "sample.type",
                                                            "primary.site"),
                                                        drop = FALSE]
  EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[!EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
  EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Metastatic",
                                                   "Primary Tumor",
                                                   "Solid Tissue Normal"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- as.factor(as.character(
      EIF.TCGA.RNAseq.anno.subset$sample.type))
    EIF.TCGA.RNAseq.anno.subset$sum <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score <- EIF.TCGA.RNAseq.anno.subset
    EIF.TCGA.GTEX.score$`EIF4A1:EIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4A1:EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4A1:EIF4G1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4G1:EIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4E:EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4E - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4G1:EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4G1:EIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4A1:EIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$sum)
    ratio <- c("EIF4A1:EIF4E", "EIF4G1:EIF4E", "EIF4E:EIF4EBP1", 
               "EIF4A1:EIF4G1", "EIF4G1:EIF4EBP1","EIF4A1:EIF4EBP1",
               "EIF4G1:EIF4E+EIF4EBP1",
               "EIF4A1:EIF4E+EIF4EBP1")
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.score[c(ratio, "sample.type","primary.site")]
    EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score)
    return(EIF.TCGA.GTEX.score.long)
    }
  EIF.TCGA.GTEX.score.long <- get.EIFratio.anno.data()
  
  make.plot <- function () {
    new.label <- c(`EIF4A1:EIF4E`    = "EIF4A1:\nEIF4E",
                   `EIF4G1:EIF4E`    = "EIF4G1:\nEIF4E",
                   `EIF4E:EIF4EBP1`  = "EIF4E:\nEIF4EBP1",
                   `EIF4A1:EIF4G1`   = "EIF4A1:\nEIF4G1",
                   `EIF4G1:EIF4EBP1` = "EIF4G1:\nEIF4EBP1",
                   `EIF4A1:EIF4EBP1` = "EIF4A1:\nEIF4EBP1",
                   `EIF4G1:EIF4E+EIF4EBP1` = "EIF4G1:\nEIF4E+EIF4EBP1",
                   `EIF4A1:EIF4E+EIF4EBP1` = "EIF4A1:\nEIF4E+EIF4EBP1")
    EIF.TCGA.GTEX.score.long1 <- EIF.TCGA.GTEX.score.long[
      EIF.TCGA.GTEX.score.long$variable %in% c(
        "EIF4A1:EIF4E","EIF4G1:EIF4E","EIF4E:EIF4EBP1","EIF4A1:EIF4G1",
        "EIF4A1:EIF4E+EIF4EBP1","EIF4G1:EIF4E+EIF4EBP1"), ]
    EIF.TCGA.GTEX.score.long1$variable = factor(
      EIF.TCGA.GTEX.score.long1$variable, 
      levels = c("EIF4A1:EIF4G1","EIF4E:EIF4EBP1",
                 "EIF4A1:EIF4E","EIF4G1:EIF4E",
                 "EIF4A1:EIF4E+EIF4EBP1","EIF4G1:EIF4E+EIF4EBP1"
                 ))
    p1 <- ggplot(data = EIF.TCGA.GTEX.score.long1,
          aes(x     = sample.type,
              y     = 2**value, 
              color = sample.type,
              fill  = sample.type)) +
          geom_violin(trim = FALSE) +
          geom_boxplot(
            alpha    = 0.01,
            width    = 0.25,
            color    = "black",
            outlier.colour=NA
          ) + 
          stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
          facet_grid(~ variable, 
                     scales = "free",
                     space  = "free") +
          facet_wrap(~ variable, 
                     labeller = labeller(variable = as_labeller(new.label)),
                     ncol = 6) +
          scale_color_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          scale_fill_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          labs(x = "sample type",
               y = "ratio of RNA counts") +
          scale_x_discrete(labels = c("Metastatic Tumor",
                                      "Primary Tumor",
                                      "Adjacent Normal")
                           ) +
          scale_y_continuous(trans = log2_trans(), 
                            #labels = label_comma(),
                             breaks = c(0.125,1,8,64,512),
                             labels = c("0.125","1","8","64","512"),
                             position = 'right') + 
          geom_hline(yintercept = 1, linetype = "dashed") +
          theme_bw() +
          theme(plot.title         = black_bold_tahoma_16,
                axis.title.x       = element_blank(),
                axis.title.y.right = black_bold_tahoma_16_right,
                axis.text.x        = black_bold_tahoma_16_90,
                axis.text.y        = black_bold_tahoma_16_90,
                axis.line.x        = element_line(color = "black"),
                axis.line.y        = element_line(color = "black"),
                panel.grid         = element_blank(),
                legend.position    = "none",
                strip.text         = black_bold_tahoma_16
              ) +
          stat_compare_means(comparisons = list(
            c("Metastatic", "Solid Tissue Normal"),
            c("Primary Tumor", "Solid Tissue Normal"),
            c("Metastatic", "Primary Tumor")),
                             method = "t.test", 
                             label  = "p.signif", 
                             size   = 6, 
                             hjust  = 0)
          print(p1)
          ggsave(path        = "~/Documents/EIF_output/Expression", 
                 filename    = "EIFratioviolin.pdf", 
                 plot        = p1,
                 width       = 9, 
                 height      = 9, 
                 useDingbats = FALSE)
          
          
          EIF.TCGA.GTEX.score.long2 <- EIF.TCGA.GTEX.score.long[
            EIF.TCGA.GTEX.score.long$variable %in% c(
              "EIF4G1:EIF4E",
              "EIF4G1:EIF4EBP1",
              "EIF4G1:EIF4E+EIF4EBP1"), ]
          EIF.TCGA.GTEX.score.long2$variable <- factor(
            EIF.TCGA.GTEX.score.long2$variable, 
            levels = c("EIF4G1:EIF4E",
                       "EIF4G1:EIF4EBP1",
                       "EIF4G1:EIF4E+EIF4EBP1"))
          p2 <- ggplot(data = EIF.TCGA.GTEX.score.long2,
                       aes(x     = sample.type,
                           y     = 2**value, 
                           color = sample.type,
                           fill  = sample.type)) +
                stat_n_text(size     = 6, 
                  fontface = "bold", 
                  angle    = 90, 
                  hjust    = 0) +
                facet_grid(~ variable, 
                  scales = "free",
                  space  = "free") +
                facet_wrap(~ variable, 
                  labeller = labeller(variable = as_labeller(new.label)),
                  ncol = 6) +
                geom_violin(trim = FALSE) +
                geom_boxplot(alpha    = 0.01,
                             width    = 0.25,
                             color    = "black",
                             outlier.colour = NA) + 
                scale_color_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
                scale_fill_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
                labs(x = "sample type",
                     y = "log2(ratio)") +
                scale_x_discrete(labels = c("Metastatic Tumor",
                                            "Primary Tumor",
                                            "Adjacent Normal")) +
                scale_y_continuous(trans = log2_trans(), 
                  #labels = label_comma(),
                  breaks = c(0.125,1,8,64,512),
                  labels = c("0.125","1","8","64","512"),
                  position = 'right') + 
                geom_hline(yintercept = 1, linetype = "dashed") +
                theme_bw() +
                theme(plot.title         = black_bold_tahoma_16,
                  axis.title.x       = element_blank(),
                  axis.title.y.right = black_bold_tahoma_16_right,
                  axis.text.x        = black_bold_tahoma_16_90,
                  axis.text.y        = black_bold_tahoma_16_90,
                  axis.line.x        = element_line(color = "black"),
                  axis.line.y        = element_line(color = "black"),
                  panel.grid         = element_blank(),
                  legend.position    = "none",
                  strip.text         = black_bold_tahoma_16) +
                stat_compare_means(comparisons = list(
                  c("Metastatic", "Solid Tissue Normal"),
                  c("Primary Tumor", "Solid Tissue Normal"),
                  c("Metastatic", "Primary Tumor")),
                  method      = "t.test", 
                  label       = "p.signif", 
                  size        = 6, 
                  hjust       = 0)
          print(p2)
          ggsave(
            path        = "~/Documents/EIF_output/Expression", 
            filename    = "EIFsumratioviolin.pdf", 
            plot        = p2,
            width       = 4.8, 
            height      = 9, 
            useDingbats = FALSE) 
          
    EIF.TCGA.GTEX.score.long2 <- EIF.TCGA.GTEX.score.long[
      EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:EIF4E",
                                               "EIF4A1:EIF4EBP1",
                                               "EIF4A1:EIF4E+EIF4EBP1"
        ), ]
    EIF.TCGA.GTEX.score.long2$variable <- factor(
      EIF.TCGA.GTEX.score.long2$variable, 
      levels = c("EIF4A1:EIF4E",
                 "EIF4A1:EIF4EBP1",
                 "EIF4A1:EIF4E+EIF4EBP1"
        ))
    p3 <- ggplot(data = EIF.TCGA.GTEX.score.long2,
                 aes(x     = sample.type,
                     y     = 2**value, 
                     color = sample.type,
                     fill  = sample.type)) +
          stat_n_text(size     = 6, 
                      fontface = "bold", 
                      angle    = 90, 
                      hjust    = 0) +
          facet_grid(~ variable, 
                     scales = "free",
                     space  = "free") +
          facet_wrap(~ variable, 
                     labeller = labeller(variable = as_labeller(new.label)),
                     ncol = 6) +
          geom_violin(trim = FALSE) +
          geom_boxplot(alpha    = 0.01,
                       width    = 0.25,
                       color    = "black",
                       outlier.colour = NA) + 
          scale_color_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          scale_fill_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
          labs(x = "sample type",
               y = "ratio of RNA counts") +
          scale_x_discrete(labels = c("Metastatic Tumor",
                                      "Primary Tumor",
                                      "Adjacent Normal")
                                    ) +
          scale_y_continuous(trans = log2_trans(), 
            #labels = label_comma(),
            breaks = c(0.125,1,8,64,512),
            labels = c("0.125","1","8","64","512"),
            position = 'right') + 
          geom_hline(yintercept = 1, linetype = "dashed") +
          theme_bw() +
          theme(plot.title         = black_bold_tahoma_16,
                axis.title.x       = element_blank(),
                axis.title.y.right = black_bold_tahoma_16_right,
                axis.text.x        = black_bold_tahoma_16_90,
                axis.text.y        = black_bold_tahoma_16_90,
                axis.line.x        = element_line(color = "black"),
                axis.line.y        = element_line(color = "black"),
                panel.grid         = element_blank(),
                legend.position    = "none",
                strip.text         = black_bold_tahoma_16) +
          stat_compare_means(comparisons = list(
            c("Metastatic", "Solid Tissue Normal"),
            c("Primary Tumor", "Solid Tissue Normal"),
            c("Metastatic", "Primary Tumor")),
                             method      = "t.test", 
                             label       = "p.signif", 
                             size        = 6, 
                             hjust       = 0)
          print(p3)
          ggsave(
            path        = "~/Documents/EIF_output/Expression", 
            filename    = "EIFsumratioviolin2.pdf", 
            plot        = p3,
            width       = 4.8, 
            height      = 9, 
            useDingbats = FALSE) 
  
  }
  make.plot()
}
plot.violingraph.EIF.ratio.TCGA (c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"))


################################################################
## stacked bar plots for eIF4F CNV status across tumor groups ## 
################################################################
plot.bargraph.EIF.CNV.TCGA <- function (EIF) {
  pan.TCGA.CNV <- function(){
    # download https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz
    TCGA.pancancer <- fread(
      "~/Downloads/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", 
      #"~/Downloads/GDC-PANCAN.gistic.tsv",
      data.table = FALSE)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$Sample),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$Sample
    TCGA.pancancer1$Sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    TCGA.sampletype,
                                    by    = "row.names",
                                    all.x = TRUE)
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    sample.type.list <- levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    cancer.type.list <- levels(TCGA.RNAseq.anno$primary.disease)
    return(TCGA.RNAseq.sampletype)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV()
  
  pancancer.TCGA.EIF <- function(){
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno.subset <- TCGA.CNV.anno#[
    #  TCGA.CNV.anno$sample.type %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Recurrent Tumor"), ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[ ,
      colnames(TCGA.CNV.anno.subset) %in% c(EIF, 
                                            "sample.type",
                                            "primary.disease")]
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.CNV.anno.subset.long$primary.disease)
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c("sample.type",
                                                 "primary.disease",
                                                 "variable",
                                                 "CNV")
    
    CNV.sum <- table(EIF.TCGA.CNV.anno.subset.long[,c("CNV","primary.disease")])
    CNV.sum <- as.data.frame(CNV.sum)
   # CNV.sum$TCGAstudy <- str_remove(CNV.sum$TCGAstudy, regex('_.*\n*.*'))
    CNV.sum$primary.disease <- ordered(CNV.sum$primary.disease, levels = rev(levels(factor(CNV.sum$primary.disease))))
    CNV.sum$CNV <- factor(CNV.sum$CNV, levels = c("-2", "-1", "0", "1", "2"))
    return(CNV.sum)
    }
  CNV.sum <- pancancer.TCGA.EIF()
  
  levels(CNV.sum$CNV)
  # reorder bars by explicitly ordering factor levels
  make.plot <- function (EIF) {
    p1 <- ggplot(CNV.sum,
      aes(fill = CNV,  order = as.numeric(CNV),
          y    = Freq, 
          x    = primary.disease)) + 
      geom_bar(stat = "identity", position = "fill") +
      labs(x = "Tumor types (TCGA pan cancer atlas 2018)",
           y = paste0("Percentage with ", EIF, " CNV")) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_12,
        axis.title.x         = black_bold_tahoma_12,
        axis.title.y         = element_blank(),
        axis.text.x          = black_bold_tahoma_12,
        axis.text.y          = black_bold_tahoma_12,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        legend.title         = element_blank(),
        legend.text          = black_bold_tahoma_12,
        legend.position      = "top",
        legend.justification = "left",
        legend.box           = "horizontal", 
        strip.text           = black_bold_tahoma_12) +
      scale_y_continuous(labels = scales::percent_format())+
      guides(fill = guide_legend(reverse = TRUE))+ # Flip ordering of legend without altering ordering in plot
      scale_fill_manual(name   = "Copy number variation",
        breaks = c("-2", "-1", "0", "1", "2"),
        labels = c("Homdel\n 0","Hetlos\n 1","Diploid\n 2","Gain\n 3","Amp\n 3+"),
        values = c('darkblue','blue',
                   'lightgreen','red',
                   'darkred')) 
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/CNV", 
      filename    = paste0(EIF, "pancancerCNV.pdf"), 
      plot        = p1,
      width       = 9, 
      height      = 9, 
      useDingbats = FALSE)}
  make.plot(EIF)
}
lapply(c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","MYC","PTEN"), 
       plot.bargraph.EIF.CNV.TCGA)

plot.bargraph.EIF.CNV.sum <- function (EIF) {
  pan.TCGA.CNV <- function(EIF){
    # download https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
    TCGA.CNV <- fread(
      "~/Downloads/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", 
      data.table = FALSE)
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.CNV1 <- TCGA.CNV[!duplicated(TCGA.CNV$Sample),
      !duplicated(colnames(TCGA.CNV))]
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
                                 all.x = TRUE)
    TCGA.CNV.anno <- as.data.frame(TCGA.CNV.sampletype)
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    sample.type.list <- levels(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno$primary.disease <- as.factor(TCGA.CNV.anno$primary.disease)
    cancer.type.list <- levels(TCGA.CNV.anno$primary.disease)
    return(TCGA.CNV.anno)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV(EIF)
  
  TCGA.CNV.anno.EIF <- function(EIF){
    TCGA.CNV.anno.subset <- TCGA.CNV.anno[
      !TCGA.CNV.anno$sample.type %in% "Solid Tissue Normal", ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[ ,
      colnames(TCGA.CNV.anno.subset) %in% c(EIF, 
        "sample.type",
        "primary.disease")]
    return(EIF.TCGA.CNV.anno.subset)}
  EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.EIF(EIF)
  
  make.CNV.sum.plot <- function (EIF) {
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.CNV.anno.subset.long$primary.disease)
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c("sample.type","primary.disease",
      "variable","CNV")
    CNV.sum <- table(EIF.TCGA.CNV.anno.subset.long[,c("CNV","variable")])
    CNV.sum <- as.data.frame(CNV.sum)
    # CNV.sum$TCGAstudy <- str_remove(CNV.sum$TCGAstudy, regex('_.*\n*.*'))
    CNV.sum$CNV <- factor(CNV.sum$CNV, levels = c("-2", "-1", "0", "1", "2"))
    CNV.sum$variable <- factor(CNV.sum$variable, 
      levels = c("PTEN", "EIF4E", "EIF4A1", "MYC", "EIF4EBP1", "EIF4G1"))
  # reorder bars by explicitly ordering factor levels
    p1 <- ggplot(CNV.sum, aes(fill = CNV, 
                              y    = Freq, 
                              x    = variable)) + 
      geom_bar(stat = "identity", position = "fill") + geom_col() +
      geom_text(aes(label = paste0(Freq/100,"%")), 
        position = position_stack(vjust = 0.5), size = 4) + 
      #scale_y_continuous(labels = scales::percent_format())+
      labs(x = "Tumor types (TCGA pan cancer atlas 2018)",
           y = "All tumors from TCGA)") +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_16,
        axis.title.x         = black_bold_tahoma_16,
        axis.title.y         = element_blank(),
        axis.text.x          = black_bold_tahoma_16,
        axis.text.y          = black_bold_tahoma_16,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        legend.title         = element_blank(),
        legend.text          = black_bold_tahoma_16,
        legend.position      = "top",
        legend.justification = "left",
        legend.box           = "horizontal", 
        strip.text           = black_bold_tahoma_16) +
      guides(fill = guide_legend(reverse = TRUE))+ # Flip ordering of legend without altering ordering in plot
      scale_fill_manual(name   = "Copy number variation",
        breaks = c("-2", "-1", "0", "1", "2"),
        labels = c("Homdel\n 0","Hetlos\n 1","Diploid\n 2","Gain\n 3","Amp\n 3+"),
        values = c('darkblue','blue','lightgreen','red','darkred')) 
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/CNV", 
      filename    = "EIFCNVsum.pdf", 
      plot        = p1,
      width       = 7, 
      height      = 7, 
      useDingbats = FALSE)
    }
  make.CNV.sum.plot(EIF)
  
  }
plot.bargraph.EIF.CNV.sum(c("PTEN", "EIF4A1", "EIF4E", "MYC", "EIF4EBP1", "EIF4G1"))

plot.bargraph.EIF.CNV.corr <- function (EIF) {
  pan.TCGA.CNV <- function(EIF){
    # https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz
    TCGA.CNV <- fread(
      "~/Downloads/Gistic2_CopyNumber_Gistic2_all_data_by_genes", 
      data.table = FALSE)
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.CNV1 <- TCGA.CNV[!duplicated(TCGA.CNV$Sample),
      !duplicated(colnames(TCGA.CNV))]
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
      all.x = TRUE)
    TCGA.CNV.anno <- as.data.frame(TCGA.CNV.sampletype)
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    sample.type.list <- levels(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno$primary.disease <- as.factor(TCGA.CNV.anno$primary.disease)
    cancer.type.list <- levels(TCGA.CNV.anno$primary.disease)
    return(TCGA.CNV.anno)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV(EIF)
  
  TCGA.CNV.anno.EIF <- function(EIF){
    TCGA.CNV.anno.subset <- TCGA.CNV.anno[
      !TCGA.CNV.anno$sample.type %in% "Solid Tissue Normal", ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[ ,
      colnames(TCGA.CNV.anno.subset) %in% c(EIF, 
        "sample.type",
        "primary.disease")]
    return(EIF.TCGA.CNV.anno.subset)}
  EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.EIF(EIF)
  
  plot.EIF.CNV.cor <- function(){
    df1 <- EIF.TCGA.CNV.anno.subset[1:(length(EIF.TCGA.CNV.anno.subset)-2)]
    # correlation plot
    res <- cor(df1,  method = "pearson")
    cor_5 <- rcorr(as.matrix(df1), type = "spearman")
    M <- cor_5$r
    p_mat <- cor_5$P
    pdf(file.path(
      path        = "~/Documents/EIF_output/CNV", 
      filename    = "EIFCNVcormatrix.pdf"), 
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    corrplot(
      res, 
      method      = "color", 
      tl.cex      = 1, 
      number.cex  = 1, 
      addgrid.col = "gray",
      addCoef.col = "black", 
      tl.col      = "black",
      #type        = "upper", 
      order       = "FPC", tl.srt = 45, 
      p.mat       = p_mat, 
      sig.level   = 0.05, #insig = "blank" 
    )
    dev.off()
  }
  plot.EIF.CNV.cor()
}
plot.bargraph.EIF.CNV.corr(c("PTEN", "EIF4A1", "EIF4E", "MYC", "EIF4EBP1", "EIF4G1"))

plot.boxgraph.EIF.CNV.RNAseq <- function (EIF) {
  pan.TCGA.gene <- function(EIF){
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    TCGA.CNV <- fread(
      "~/Downloads/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", 
      data.table = FALSE)
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.RNAseq1 <- TCGA.RNAseq[!duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- "RNAseq"
    
    TCGA.CNV1 <- TCGA.CNV[!duplicated(TCGA.CNV$Sample),
      !duplicated(colnames(TCGA.CNV))]
    row.names(TCGA.CNV1) <- TCGA.CNV1$Sample
    TCGA.CNV1$Sample <- NULL
    TCGA.CNV1 <- TCGA.CNV1[EIF, ]
    TCGA.CNV_transpose <- data.table::transpose(TCGA.CNV1)
    rownames(TCGA.CNV_transpose) <- colnames(TCGA.CNV1)
    colnames(TCGA.CNV_transpose) <- rownames(TCGA.CNV1)
    colnames(TCGA.CNV_transpose) <- "CNV"
    
    TCGA.RNAseq.CNV <- merge(TCGA.RNAseq_transpose,
                             TCGA.CNV_transpose,
                             by    = "row.names",
                             all.x = TRUE)
    
    TCGA.RNAseq.CNV <- as.data.frame(TCGA.RNAseq.CNV)
    TCGA.RNAseq.CNV$CNV <- as.factor(TCGA.RNAseq.CNV$CNV)
    TCGA.RNAseq.CNV$CNV <- factor(TCGA.RNAseq.CNV$CNV, 
                                   levels = c("2", "1", "0", "-1", "-2"))
    TCGA.RNAseq.CNV <- na.omit(TCGA.RNAseq.CNV)
    return(TCGA.RNAseq.CNV)
  }
  TCGA.RNAseq.CNV <- pan.TCGA.gene(EIF)
  make.plot <- function (EIF) {
  p1 <- ggplot(data = TCGA.RNAseq.CNV,
               aes(x     = CNV,
                   y     = 2**RNAseq-1, 
                   color = CNV,
                   fill  = CNV)) +    
               scale_y_continuous(trans = log2_trans(), 
                                  labels = label_comma()) +
        geom_violin(trim = FALSE) +
        geom_boxplot(
          alpha      = .01,
          width      = 0.25,
          color      = "black") +
    #scale_fill_manual(values = c("#0072B2", "#56B4E9", "#009E73", 
    #                             "#CC79A7", "#D55E00")) + #for color-blind palettes
    scale_color_manual(values = c("dark red", "red", "light green", 
                                  "blue", "dark blue"))+
    stat_n_text(size = 6, fontface = "bold") + 
    scale_fill_manual(name   = "Copy number variation",
                      breaks = c("2", "1", "0", "-1", "-2"),
                      # labels = c("Amp\n 3+","Gain\n 3","Diploid\n 2","Hetlos\n 1","Homdel\n 0"),
                      values = c('darkred','red','lightgreen','blue','darkblue')) +
    scale_x_discrete(breaks = c("2", "1", "0", "-1", "-2"),
                     labels = c("Amp\n3+",
                                "Gain\n3",
                                "Diploid\n2",
                                "Hetlos\n1",
                                "Homdel\n0")) +
    labs(x = paste(EIF, "copy number variation"),
         y = paste0(EIF, " RNA counts")) +
    theme_bw() +
    theme(
      axis.title      = black_bold_tahoma_16,
      axis.text       = element_text(size   = 16,
                                     hjust  = 0.5,
                                     face   = "bold",
                                     color  = "black"),
      axis.line.x     = element_line(color  = "black"),
      axis.line.y     = element_line(color  = "black"),
      panel.grid      = element_blank(),
      strip.text      = black_bold_tahoma_16,
      legend.position = "none") +
    stat_compare_means(comparisons = list(c("1","0"),
                                          c("2","0"),
                                          c("-1","0"),
                                          c("-2","0")),
                                          method = "t.test", 
                                          label  = "p.signif",
                                          size   = 6)

    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/CNV", 
      filename    = paste0(EIF,"CNV&RNAseq.pdf"), 
      plot        = p1,
      width       = 6.5, 
      height      = 7.5, 
      useDingbats = FALSE)}
  make.plot(EIF)
}
lapply(c("PTEN", "EIF4A1", "EIF4E", "MYC", "EIF4EBP1", "EIF4G1"), 
       plot.boxgraph.EIF.CNV.RNAseq)
plot.boxgraph.EIF.CNV.RNAseq("EIF4A1")

plot.bargraph.EIF.CNVratio.TCGA <- function (EIF) {
  pan.TCGA.CNV <- function(){
    #download https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz
    TCGA.pancancer <- fread(
      #"~/Downloads/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", 
      "~/Downloads/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena",
      data.table = FALSE)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$sample),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
      TCGA.sampletype,
      by    = "row.names",
      all.x = TRUE)
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    sample.type.list <- levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    cancer.type.list <- levels(TCGA.RNAseq.anno$primary.disease)
    return(TCGA.RNAseq.sampletype)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV()
  
  pancancer.TCGA.EIF <- function(){
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno.subset <- TCGA.CNV.anno#[
      #  TCGA.CNV.anno$sample.type %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Recurrent Tumor"), ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[ ,
      colnames(TCGA.CNV.anno.subset) %in% c(EIF, 
                                            "sample.type",
                                            "primary.disease")]
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
        EIF.TCGA.CNV.anno.subset.long$primary.disease)
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c("sample.type",
                                                 "primary.disease",
                                                 "variable",
                                                 "CNV")
    return(EIF.TCGA.CNV.anno.subset.long)
  }
  EIF.TCGA.CNV.anno.subset.long <- pancancer.TCGA.EIF()
  class(EIF.TCGA.CNV.anno.subset.long$CNV)
  # reorder bars by explicitly ordering factor levels
  make.plot <- function (EIF) {
    sts <- boxplot.stats(EIF.TCGA.CNV.anno.subset.long$CNV)$stats
    
    f1 <- factor(EIF.TCGA.CNV.anno.subset.long$primary.disease)
    f.ordered1 <- fct_rev(f1)
    p1 <- ggplot(data = EIF.TCGA.CNV.anno.subset.long,
                 aes(y     = 2**CNV,
                     x     = f.ordered1, 
                     color = primary.disease)) +    
          ylim(0,3)+
          geom_hline(yintercept = 1, linetype = "dashed") +
          stat_n_text(size     = 5, 
                      fontface = "bold", 
                      hjust    = 0) +
          geom_boxplot(
            alpha    = .01, outlier.colour = NA,
            #size     = .75,
            #width    = 1,
            position = position_dodge(width = .9)
          ) +
          labs(x = "primary disease",
               y = paste(EIF, "CNV ratio", "(tumor/normal)")) +
          coord_cartesian(ylim = c(sts[2]/2, max(sts)*1.05)) +
          coord_flip() +
          theme_bw() +
          theme(
            plot.title           = black_bold_tahoma_12,
            axis.title.x         = black_bold_tahoma_12,
            axis.title.y         = element_blank(),
            axis.text.x          = black_bold_tahoma_12,
            axis.text.y          = black_bold_tahoma_12,
            axis.line.x          = element_line(color = "black"),
            axis.line.y          = element_line(color = "black"),
            panel.grid           = element_blank(),
            legend.title         = element_blank(),
            legend.text          = black_bold_tahoma_12,
            legend.position      = "none",
            legend.justification = "left",
            legend.box           = "horizontal", 
            strip.text           = black_bold_tahoma_12
          )
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/CNV", 
      filename    = paste0(EIF, "pancancerCNVratio.pdf"), 
      plot        = p1,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
  }
    make.plot(EIF)
  }
lapply(c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","MYC","PTEN"), 
       plot.bargraph.EIF.CNVratio.TCGA)

plot.bargraph.EIF.CNV.GDC <- function (EIF) {
  pan.TCGA.CNV <- function(){
    # download https://gdc.xenahubs.net/download/GDC-PANCAN.gistic.tsv.gz
    TCGA.pancancer <- fread(
      "~/Downloads/GDC-PANCAN.gistic.tsv",
      data.table = FALSE)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$V1),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$V1
    TCGA.pancancer1$V1 <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    
    # download https://gdc.xenahubs.net/download/GDC-PANCAN.basic_phenotype.tsv.gz
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/GDC-PANCAN.basic_phenotype.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.sampletype <- TCGA.sampletype[ ,c("sample","sample_type", "project_id")]
    
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
      TCGA.sampletype,
      by    = "row.names",
      all.x = TRUE)
    TCGA.RNAseq.anno <- as.data.frame(TCGA.RNAseq.sampletype)
    TCGA.RNAseq.anno$sample.type <- as.factor(TCGA.RNAseq.anno$sample.type)
    sample.type.list <- levels(TCGA.RNAseq.anno$sample.type)
    TCGA.RNAseq.anno$primary.disease <- as.factor(TCGA.RNAseq.anno$primary.disease)
    cancer.type.list <- levels(TCGA.RNAseq.anno$primary.disease)
    return(TCGA.RNAseq.sampletype)
  }
  TCGA.CNV.anno <- pan.TCGA.CNV()
  
  pancancer.TCGA.EIF <- function(){
    TCGA.CNV.anno$sample.type <- as.factor(TCGA.CNV.anno$sample.type)
    TCGA.CNV.anno.subset <- TCGA.CNV.anno[grep("TCGA", TCGA.CNV.anno$primary.disease), ]
    TCGA.CNV.anno.subset <- TCGA.CNV.anno#[
    #  TCGA.CNV.anno$sample.type %in% c("Primary Blood Derived Cancer - Peripheral Blood", "Recurrent Tumor"), ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[ ,
      colnames(TCGA.CNV.anno.subset) %in% c("ENSG00000161960", 
        "sample.type",
        "primary.disease")]
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.CNV.anno.subset.long$primary.disease)
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c("sample.type",
      "primary.disease",
      "variable",
      "CNV")
    
    CNV.sum <- table(EIF.TCGA.CNV.anno.subset.long[,c("CNV","primary.disease")])
    CNV.sum <- as.data.frame(CNV.sum)
    # CNV.sum$TCGAstudy <- str_remove(CNV.sum$TCGAstudy, regex('_.*\n*.*'))
    CNV.sum$primary.disease <- ordered(CNV.sum$primary.disease, levels = rev(levels(factor(CNV.sum$primary.disease))))
    CNV.sum$CNV <- factor(CNV.sum$CNV, levels = c("-2", "-1", "0", "1", "2"))
    return(CNV.sum)
  }
  CNV.sum <- pancancer.TCGA.EIF()
  
  levels(CNV.sum$CNV)
  # reorder bars by explicitly ordering factor levels
  make.plot <- function (EIF) {
    p1 <- ggplot(CNV.sum,
      aes(fill = CNV,  order = as.numeric(CNV),
        y    = Freq, 
        x    = primary.disease)) + 
      geom_bar(stat = "identity", position = "fill") +
      labs(x = "Tumor types (TCGA pan cancer atlas 2018)",
        y = paste0("Percentage with ", EIF, " CNV")) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_12,
        axis.title.x         = black_bold_tahoma_12,
        axis.title.y         = element_blank(),
        axis.text.x          = black_bold_tahoma_12,
        axis.text.y          = black_bold_tahoma_12,
        axis.line.x          = element_line(color = "black"),
        axis.line.y          = element_line(color = "black"),
        panel.grid           = element_blank(),
        legend.title         = element_blank(),
        legend.text          = black_bold_tahoma_12,
        legend.position      = "top",
        legend.justification = "left",
        legend.box           = "horizontal", 
        strip.text           = black_bold_tahoma_12) +
      scale_y_continuous(labels = scales::percent_format())+
      guides(fill = guide_legend(reverse = TRUE))+ # Flip ordering of legend without altering ordering in plot
      scale_fill_manual(name   = "Copy number variation",
        breaks = c("-2", "-1", "0", "1", "2"),
        labels = c("Homdel\n 0","Hetlos\n 1","Diploid\n 2","Gain\n 3","Amp\n 3+"),
        values = c('darkblue','blue',
          'lightgreen','red',
          'darkred')) 
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/CNV", 
      filename    = paste0(EIF, "pancancerCNV.pdf"), 
      plot        = p1,
      width       = 9, 
      height      = 9, 
      useDingbats = FALSE)}
  make.plot(EIF)
}
lapply(c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","MYC","PTEN"), 
  plot.bargraph.EIF.CNV.GDC)

#################################################################
##  PCA plots on EIF4F RNA-seq data from TCGA and GTEx groups  ##
#################################################################
plot.EIF.TCGA.GTEX.PCA.all <- function (EIF.list) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
                           !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE)
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
  EIF.list <- c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1", 
                "PABPC1","MKNK1","MKNK2")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.list, 
                                                            "sample.type",
                                                            "primary.site"),
                                                            drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site == "Brain", ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Metastatic",
                                                     "Primary Tumor",
                                                     "Normal Tissue"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c("Normal Tissue", 
                 #"Solid Tissue Normal", 
                 "Primary Tumor", 
                 "Metastatic"),
      labels = c("Healthy Tissue (GTEx)", 
                 #"Adjacent Normal Tissue (TCGA)",
                 "Primary Tumor (TCGA)", 
                 "Metastatic Tumor (TCGA)"))
    #EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
    }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)
  ## remove the last two columns 
  df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-1)]
  rownames(df1) <- NULL
  
  plot.PCA.prcomp <- function(){
    # the variables should be scaled to have unit variance 
    PCA <- prcomp(df1[1:(length(df1)-1)], center = TRUE, scale = TRUE)
   # Extract PC axes for plotting
    PCAvalues <- data.frame(Sample.type = EIF.TCGA.RNAseq.anno.subset$sample.type, 
      PCA$x)
    # Extract loadings of the variables
    PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
    # Plot
    plot.PCA.3D <- function(){
      get_colors <- function(groups, group.col = palette()){
      groups <- as.factor(groups)
      ngrps <- length(levels(groups))
      if(ngrps > length(group.col)) 
        group.col <- rep(group.col, ngrps)
      color <- group.col[as.numeric(groups)]
      names(color) <- as.vector(groups)
      return(color)
      }
      cols <- get_colors(PCAvalues$Sample.type, brewer.pal(n=4, name="Dark2"))
      plot3d(PCAvalues[,2:4], col = cols)
      legend3d("topright", legend =  levels(PCAvalues$Sample.type), 
               pch = 16, col = brewer.pal(n = 4, name="Dark2"), 
               cex = 1)}
    plot.PCA.3D()
    p <- ggplot(PCAvalues,
      aes(x      = PC1,
          y      = PC2,
          colour = Sample.type)) +
      #xlim(-5, 8) +
      geom_point(size = 0.2) +
      geom_segment(data = PCAloadings,
        aes(
          x     = 0,
          y     = 0,
          xend  = (PC1*5),
          yend  = (PC2*5)),
          arrow = arrow(length = unit(1/3, "picas")),
          color = "black") +
      annotate("text",
        size     = 6,
        fontface = "bold",
        x        = (PCAloadings$PC1*5),
        y        = (PCAloadings$PC2*5),
        label    = PCAloadings$Variables) +
      # stat_n_text(geom = "label") +
      ggtitle("Principal Component Analysis") +
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(fill   = 'transparent',
                                         color  = 'black',
                                         size   = 1),
        axis.title        = black_bold_tahoma_16,
        axis.text         = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0.3, 0.93),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16,
        legend.key        = element_blank())
    print(p)
  }
  plot.PCA.prcomp()
  
  plot.pca.factomineR <- function(){
    res.pca <- PCA(df1[1:(length(df1)-1)], 
      scale.unit = TRUE, 
      ncp        = 10, 
      graph      = FALSE)
    biplot <- fviz_pca_biplot(res.pca, 
      axes       = c(1, 2),
      labelsize  = 5,
      col.ind    = df1$sample.type, 
      palette    = c("#D55E00","#009E73","#CC79A7","#0072B2"), 
      #palette    = c("#CC79A7","#0072B2"), 
      pointshape = 20,
      pointsize  = 0.75,
      #addEllipses = TRUE,
      title      = "PCA - Biplot (All)",
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      xlim(-7, 8) + ylim (-6, 7.5)+ # for EIF 8
      #xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
      theme_classic() + 
      #scale_alpha_manual(values=c(0.1, 0.1, 0.1, 0.1),guide=F)+
      #scale_x_continuous(breaks = seq(-6, 8, 2), limits=c(-5, 8)) +
      #scale_y_continuous(breaks = seq(-4, 6, 2), limits=c(-4, 7)) +
      theme(
        plot.background  = element_blank(),
        plot.title       = black_bold_tahoma_16,
        panel.background = element_rect(fill   = 'transparent',
          color  = 'black',
          size   = 1),
        axis.title.x     = black_bold_tahoma_16,
        axis.title.y     = black_bold_tahoma_16,
        axis.text.x      = black_bold_tahoma_16,
        axis.text.y      = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0, 0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(biplot)
    ggsave(
      path        = "~/Documents/EIF_output/PCA/TCGA", 
      filename    = "EIFPCAall.pdf", 
      plot        = biplot,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    

    plot.selected.PCA <- function (sample, color){
      test <- df1[df1$sample.type == sample, ]
      sample.type <- levels(df1$sample.type)
      biplot <- fviz_pca_biplot(res.pca, 
        axes       = c(1, 2),
        labelsize  = 5,
        col.ind    = df1$sample.type, 
        #palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"), 
        palette    = color, 
        select.ind = list(name = row.names(test)),
        pointshape = 20,
        pointsize  = 0.75,
        #addEllipses = TRUE,
        title      = "PCA - Biplot (All)",
        label      = "var",
        col.var    = "black", 
        repel      = TRUE) +
        xlim(-7, 8) + ylim (-6, 7.5)+ # for EIF 8
        #xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
        theme_classic() + 
        #scale_alpha_manual(values=c(0.1, 0.1, 0.1, 0.1),guide=F)+
        #scale_x_continuous(breaks = seq(-6, 8, 2), limits=c(-5, 8)) +
        #scale_y_continuous(breaks = seq(-4, 6, 2), limits=c(-4, 7)) +
        theme(
          plot.background  = element_blank(),
          plot.title       = black_bold_tahoma_16,
          panel.background = element_rect(fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x     = black_bold_tahoma_16,
          axis.title.y     = black_bold_tahoma_16,
          axis.text.x      = black_bold_tahoma_16,
          axis.text.y      = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = c(0, 0),
          legend.justification = c(0,0),
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(biplot)
      ggsave(
        path        = "~/Documents/EIF_output/PCA/TCGA", 
        filename    = paste0("EIFPCAall",sample,".pdf"), 
        plot        = biplot,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.selected.PCA ("Healthy Tissue (GTEx)", "#D55E00")
    plot.selected.PCA ("Adjacent Normal Tissue (TCGA)", "#CC79A7")
    plot.selected.PCA ("Primary Tumor (TCGA)", "#009E73")
    plot.selected.PCA ("Metastatic Tumor (TCGA)", "#CC79A7")
    
    indplot <- fviz_pca_ind(res.pca,
      labelsize   = 5,
      col.ind     = EIF.TCGA.RNAseq.anno.subset$sample.type, 
      palette     = "dark1", 
      #pointshape  = 20,
      pointsize   = 0.5,
      addEllipses = TRUE, 
      label       = "var",
      col.var     = "black", 
      repel       = TRUE) +
      theme_classic() + 
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(fill   = 'transparent',
          color  = 'black',
          size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0, 0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(indplot)
    varplot <- fviz_pca_var(res.pca,
      labelsize  = 5,
      col.ind    = EIF.TCGA.RNAseq.anno.subset$sample.type, 
      palette    = "dark1", 
      pointshape = 20,
      #addEllipses = TRUE, 
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      theme_classic() + 
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(fill   = 'transparent',
          color  = 'black',
          size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0, 0),
        legend.justification = c(0,0),
        #legend.position   = c(0.75, 0.93),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(varplot)
    
    eig <- fviz_eig(res.pca, 
                    labelsize = 6,
                    geom      = "bar", 
                    width     = 0.7, 
                    addlabels = TRUE) + 
                    # geom_text(aes(label = res.pca$eig, size = 18)) +
                    theme_classic() +
                    theme(
                      plot.background   = element_blank(),
                      plot.title        = black_bold_tahoma_16,
                      panel.background  = element_rect(
                        fill   = 'transparent',
                        color  = 'black',
                        size   = 1),
                      axis.title.x      = black_bold_tahoma_16,
                      axis.title.y      = black_bold_tahoma_16,
                      axis.text.x       = black_bold_tahoma_16,
                      axis.text.y       = black_bold_tahoma_16)
    print(eig)
    ggsave(
      path        = "~/Documents/EIF_output/PCA/TCGA", 
      filename    = "EIFPCAeig.pdf", 
      plot        = eig,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    
    
    fviz_pca_var(res.pca, col.var="contrib")
    fviz_contrib(res.pca, choice="var", axes = 2, top = 10 )

    
    var <- get_pca_var(res.pca)

    pdf(file.path(path        = "~/Documents/EIF_output/PCA/TCGA", 
                  filename    = "EIFPCAcor.pdf"), 
                  width       = 9, 
                  height      = 9, 
                  useDingbats = FALSE)
    corrplot(var$cos2, #cos2 is better than contribute
      is.corr     = FALSE, 
      tl.cex      = 1.5, 
      number.cex  = 1.5, 
      method      = "color", 
      addgrid.col = "gray",
      addCoef.col = "black", 
      tl.col      = "black")
    dev.off()
    
    #corrplot(var$contrib, is.corr=FALSE)    
    contribplot <- function(x){
      fviz_contrib(res.pca,
        choice = "var",
        axes   = x,
        top    = 10,
        fill   = "lightblue",
        color  = "black") +
        theme_minimal() +
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(fill   = 'transparent',
                                           color  = 'black',
                                           size   = 1),
          axis.title.x      = element_blank(),
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16_45,
          axis.text.y       = black_bold_tahoma_16)
      }
    lapply(c(1,2), contribplot)
  }
  plot.pca.factomineR()
  
  get.EIFsum.TCGA.GTEX <- function(EIF.list) {
    #EIF.list <- c("EIF4E", "EIF4G1", "EIF4G2", "EIF4A1","EIF4EBP1", "PABPC1",
    #              "MKNK1","MKNK2", "MYC","JUN","YY1")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.list,
      "MYC","JUN","YY1",
      "sample.type",
      "primary.site"),
      drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -2 + 1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[c(EIF.list,"MYC","JUN","YY1","sample.type","primary.site")]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4G1 == 0, ]
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site == "Brain", ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Metastatic",
                                                     "Primary Tumor",
                                                     "Normal Tissue"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c("Normal Tissue", 
                 #"Solid Tissue Normal", 
                 "Primary Tumor", 
                 "Metastatic"),
      labels = c("Healthy Tissue (GTEx)", 
                 #"Adjacent Normal Tissue (TCGA)",
                 "Primary Tumor (TCGA)", 
                 "Metastatic Tumor (TCGA)"))
    #EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIFsum.TCGA.RNAseq.anno.subset <- get.EIFsum.TCGA.GTEX(EIF.list)
  ## remove the last two columns 
  df2 <- EIFsum.TCGA.RNAseq.anno.subset[1:(length(EIFsum.TCGA.RNAseq.anno.subset)-2)]
  rownames(df2) <- NULL
  
  #my_data <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.gene, "sum")]
  res <- cor(df2, method = "pearson")
  cor_5 <- rcorr(as.matrix(df2))
  M <- cor_5$r
  p_mat <- cor_5$P
  
  corrplot(
    res, 
    method      = "color", 
    tl.cex      = 1, 
    number.cex  = 1, 
    addgrid.col = "gray",
    addCoef.col = "black", 
    tl.col      = "black",
    #type        = "upper", 
    order       = "FPC", 
    p.mat       = p_mat, 
    sig.level   = 0.0, #insig = "blank" 
  )
  
  plot.pca.factomineR.sum <- function(){
    res.pca <- PCA(df2, 
      scale.unit = TRUE, 
      ncp        = 10, 
      graph      = FALSE)
    
    biplot <- fviz_pca_biplot(res.pca, 
      axes       = c(1, 2),
      labelsize  = 5,
      col.ind    = EIFsum.TCGA.RNAseq.anno.subset$sample.type, 
      palette    = c("#D55E00","#009E73","#CC79A7","#0072B2"), 
      pointshape = 20,
      pointsize  = 0.75,
      #addEllipses = TRUE, 
      title      = "PCA - Biplot (All)",
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      theme_classic() + 
      xlim(-7, 8) + ylim (-6, 7.5)+ # for EIF 8
      #scale_y_continuous(breaks = seq(-4, 6, 2), limits=c(-4, 7)) +
      theme(
        plot.background  = element_blank(),
        plot.title       = black_bold_tahoma_16,
        panel.background = element_rect(fill   = 'transparent',
                                        color  = 'black',
                                        size   = 1),
        axis.title.x     = black_bold_tahoma_16,
        axis.title.y     = black_bold_tahoma_16,
        axis.text.x      = black_bold_tahoma_16,
        axis.text.y      = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0, 0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(biplot)
      ggsave(
      path        = "~/Documents/EIF_output/PCA/TCGA", 
      filename    = "EIFsumPCAall.pdf", 
      plot        = biplot,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
      
      plot.selected.PCA <- function (sample, color){
        test <- df2[df1$sample.type == sample, ]
        sample.type <- levels(df2$sample.type)
        biplot <- fviz_pca_biplot(res.pca, 
          axes       = c(1, 2),
          labelsize  = 5,
          col.ind    = df1$sample.type, 
          #palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"), 
          palette    = color, 
          select.ind = list(name = row.names(test)),
          pointshape = 20,
          pointsize  = 0.75,
          #addEllipses = TRUE,
          title      = "PCA - Biplot (All)",
          label      = "var",
          col.var    = "black", 
          repel      = TRUE) +
          xlim(-7, 8) + ylim (-6, 7.5)+ # for EIF 8
          #xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
          theme_classic() + 
          #scale_alpha_manual(values=c(0.1, 0.1, 0.1, 0.1),guide=F)+
          #scale_x_continuous(breaks = seq(-6, 8, 2), limits=c(-5, 8)) +
          #scale_y_continuous(breaks = seq(-4, 6, 2), limits=c(-4, 7)) +
          theme(
            plot.background  = element_blank(),
            plot.title       = black_bold_tahoma_16,
            panel.background = element_rect(fill   = 'transparent',
              color  = 'black',
              size   = 1),
            axis.title.x     = black_bold_tahoma_16,
            axis.title.y     = black_bold_tahoma_16,
            axis.text.x      = black_bold_tahoma_16,
            axis.text.y      = black_bold_tahoma_16,
            legend.title      = element_blank(),
            legend.position   = c(0, 0),
            legend.justification = c(0,0),
            legend.background = element_blank(),
            legend.text       = black_bold_tahoma_16)
        print(biplot)
        ggsave(
          path        = "~/Documents/EIF_output/PCA/TCGA", 
          filename    = paste0("EIFsumPCAall",sample,".pdf"), 
          plot        = biplot,
          width       = 8, 
          height      = 8, 
          useDingbats = FALSE)
      }
      plot.selected.PCA ("Healthy Tissue (GTEx)", "#D55E00")
      plot.selected.PCA ("Adjacent Normal Tissue (TCGA)", "#CC79A7")
      plot.selected.PCA ("Primary Tumor (TCGA)", "#009E73")
      plot.selected.PCA ("Metastatic Tumor (TCGA)", "#CC79A7")  
      
    indplot <- fviz_pca_ind(res.pca,
      labelsize   = 5,
      col.ind     = EIFsum.TCGA.RNAseq.anno.subset$sample.type, 
      palette     = "dark1", 
      #pointshape  = 20,
      pointsize   = 0.5,
      addEllipses = TRUE, 
      label       = "var",
      col.var     = "black", 
      repel       = TRUE) +
      theme_classic() + 
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(fill   = 'transparent',
                                         color  = 'black',
                                         size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0, 0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(indplot)
    corrplot(var$contrib, is.corr=FALSE)    
    varplot <- fviz_pca_var(res.pca,
      labelsize  = 5,
      col.ind    = EIFsum.TCGA.RNAseq.anno.subset$sample.type, 
      palette    = "dark1", 
      pointshape = 20,
      #addEllipses = TRUE, 
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      theme_classic() + 
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(fill   = 'transparent',
                                         color  = 'black',
                                         size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0, 0),
        legend.justification = c(0,0),
        #legend.position   = c(0.75, 0.93),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(varplot)
    eig <- fviz_eig(res.pca, 
      labelsize = 6,
      geom      = "bar", 
      width     = 0.7, 
      addlabels = TRUE) + 
      # geom_text(aes(label = res.pca$eig, size = 18)) +
      theme_classic() +
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(fill   = 'transparent',
                                         color  = 'black',
                                         size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16)
    print(eig)
    ggsave(
      path        = "~/Documents/EIF_output/PCA/TCGA", 
      filename    = "EIFsumPCAeig.pdf", 
      plot        = eig,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    var <- get_pca_var(res.pca)
    #fviz_pca_var(res.pca, col.var="contrib")
    
    pdf(file.path(
      path        = "~/Documents/EIF_output/PCA/TCGA", 
      filename    = "EIFsumPCAcor.pdf"), 
      width       = 9, 
      height      = 9, 
      useDingbats = FALSE)
    corrplot(var$cos2, #cos2 is better than contribute
             is.corr     = FALSE, 
             tl.cex      = 1.5, 
             number.cex  = 1.5, 
             method      = "color", 
             addgrid.col = "gray",
             addCoef.col = "black", 
             tl.col      = "black")
    dev.off()
    
    contribplot <- function(x){
      fviz_contrib(res.pca,
        choice = "var",
        axes   = x,
        top    = 10,
        fill   = "lightblue",
        color  = "black") +
        theme_minimal() +
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(fill   = 'transparent',
                                           color  = 'black',
                                           size   = 1),
          axis.title.x      = element_blank(),
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16_45,
          axis.text.y       = black_bold_tahoma_16)
    }
    lapply(c(1,2), contribplot)
  }
  plot.pca.factomineR.sum()
}
plot.EIF.TCGA.GTEX.PCA.all(c("EIF4G1", "EIF4A1","EIF4E", "EIF4EBP1", 
                             "PABPC1", "MKNK1","MKNK2"))

plot.EIF.TCGA.GTEX.PCA.each <- function (EIF.list, tissue) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    #TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site", "primary disease or tissue"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site", "primary.disease")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
                           !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE)
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
    #EIF.list <- c("EIF4E", "EIF4G1", "EIF4G2", "EIF4A1","EIF4EBP1", "PABPC1",
    #              "MKNK1","MKNK2", "MTOR", "RPTOR", "RPS6KB1","MYC")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.list, 
      "sample.type",
      "primary.site","primary.disease"),
      drop = FALSE]

    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[c(EIF.list,"sample.type","primary.site","primary.disease")]
    
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Metastatic",
                                                     "Primary Tumor",
                                                     "Normal Tissue",
                                                     "Solid Tissue Normal"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c("Normal Tissue", 
                "Solid Tissue Normal", 
                "Primary Tumor", 
                "Metastatic"),
      labels = c("Healthy Tissue (GTEx)", 
                "Adjacent Normal Tissue (TCGA)",
                "Primary Tumor (TCGA)", 
                "Metastatic Tumor (TCGA)"))
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset$primary.disease <- as.factor(EIF.TCGA.RNAseq.anno.subset$primary.disease)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno <- get.EIF.TCGA.GTEX(EIF.list)

  EIF.PCA.tissue <- function(x){
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno[
      EIF.TCGA.RNAseq.anno$primary.site == x, ]
    print(summary(EIF.TCGA.RNAseq.anno.subset))
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-3)]
    rownames(df1) <- NULL
    plot.pca.factomineR <- function(x){
      res.pca <- PCA(df1, 
        scale.unit = TRUE, 
        ncp        = 10, 
        graph      = FALSE)
      
      biplot <- fviz_pca_biplot(res.pca, 
        axes       = c(1, 2),
        labelsize  = 5,
        col.ind    = EIF.TCGA.RNAseq.anno.subset$sample.type, 
        palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"), 
        pointshape = 20,
        pointsize  = 0.75,
        #addEllipses = TRUE, ellipse.level = 0.9,
        label      = "var",
        col.var    = "black", 
        repel      = TRUE,
        title      = paste0("PCA - Biplot (", x,")") ) +
        #scale_x_continuous(breaks = seq(-6, 12, 2), limits=c(-5, 12)) +
        #scale_y_continuous(breaks = seq(-4, 10, 2), limits=c(-4, 10)) +
        theme_classic() + 
        theme(
          plot.background  = element_blank(),
          plot.title       = black_bold_tahoma_16,
          panel.background = element_rect(fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x     = black_bold_tahoma_16,
          axis.title.y     = black_bold_tahoma_16,
          axis.text.x      = black_bold_tahoma_16,
          axis.text.y      = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = c(0, 0),
          legend.justification = c(0,0),
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(biplot)
      ggsave(
        path        = "~/Documents/EIF_output/PCA/TCGA", 
        filename    = paste0(x,"EIFPCA.pdf"), 
        plot        = biplot,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      eig <- fviz_eig(res.pca, 
        labelsize = 6,
        geom      = "bar", 
        width     = 0.7, 
        addlabels = TRUE) + 
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background  = element_blank(),
          plot.title       = black_bold_tahoma_16,
          panel.background = element_rect(
                                          fill   = 'transparent',
                                          color  = 'black',
                                          size   = 1),
          axis.title.x    = black_bold_tahoma_16,
          axis.title.y    = black_bold_tahoma_16,
          axis.text.x     = black_bold_tahoma_16,
          axis.text.y     = black_bold_tahoma_16)
      print(eig)
      ggsave(
        path        = "~/Documents/EIF_output/PCA/TCGA", 
        filename    = paste0(x,"EIFeig.pdf"), 
        plot        = eig,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      contribplot <- function(x){
        p <- fviz_contrib(res.pca,
          choice = "var",
          axes   = x,
          top    = 10,
          fill   = "lightblue",
          color  = "black") +
          theme_classic() +
          theme(
            plot.background   = element_blank(),
            plot.title        = black_bold_tahoma_16,
            panel.background  = element_rect(fill   = 'transparent',
              color  = 'black',
              size   = 1),
            axis.title.x      = element_blank(),
            axis.title.y      = black_bold_tahoma_16,
            axis.text.x       = black_bold_tahoma_16_45,
            axis.text.y       = black_bold_tahoma_16)
        print(p)
        ggsave(
          path        = "~/Documents/EIF_output/PCA/TCGA", 
          filename    = paste0("EIFcontri",x,".pdf"), 
          plot        = p,
          width       = 8, 
          height      = 4, 
          useDingbats = FALSE)
      }
      lapply(c(1,2), contribplot)

      var <- get_pca_var(res.pca)
      #corrplot(var$contrib, is.corr=FALSE)    
      #fviz_pca_var(res.pca, col.var="contrib")
      pdf(file.path(
        path        = "~/Documents/EIF_output/PCA/TCGA", 
        filename    = paste0(x,"EIFPCAcor.pdf")), 
        width       = 9, 
        height      = 9, 
        useDingbats = FALSE)
      corrplot(var$cos2, #cos2 is better than contribute
        is.corr     = FALSE, 
        tl.cex      = 1.5, 
        number.cex  = 1.5, 
        method      = "color", 
        addgrid.col = "gray",
        addCoef.col = "black", 
        tl.col      = "black")
      dev.off()
    }
    plot.pca.factomineR(x)}
  #lapply(disease.list, EIF.PCA.tissue)
  EIF.PCA.tissue(tissue)
  
  get.EIFsum.TCGA.GTEX <- function(EIF.list) {
    #EIF.list <- c("EIF4E", "EIF4G1", "EIF4G2", "EIF4A1","EIF4EBP1", "PABPC1",
    #              "MKNK1","MKNK2", "MTOR", "RPTOR", "RPS6KB1","MYC")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.list, 
      "sample.type",
      "primary.site"),
      drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -2 +1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[c("EIF4G1", "EIF4A1","EIF4E+EIF4EBP1","PABPC1","MKNK1","MKNK2","MYC","sample.type","primary.site")]
    
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4G1 == 0, ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Metastatic",
        "Primary Tumor",
        "Normal Tissue",
        "Solid Tissue Normal"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c("Normal Tissue", 
                 "Solid Tissue Normal", 
                 "Primary Tumor", 
                 "Metastatic"),
      labels = c("Healthy Tissue (GTEx)", 
                 "Adjacent Normal Tissue (TCGA)",
                 "Primary Tumor (TCGA)", 
                 "Metastatic Tumor (TCGA)"))
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIFsum.TCGA.RNAseq.anno <- get.EIFsum.TCGA.GTEX(EIF.list)
  EIFsum.PCA.tissue <- function(x){
    EIFsum.TCGA.RNAseq.anno <- EIFsum.TCGA.RNAseq.anno[
      EIF.TCGA.RNAseq.anno$primary.site == x, ]
    df1 <- EIFsum.TCGA.RNAseq.anno[1:(length(EIFsum.TCGA.RNAseq.anno)-2)]
    rownames(df1) <- NULL
    plot.pca.factomineR <- function(x){
      res.pca <- PCA(df1, 
        scale.unit = TRUE, 
        ncp        = 10, 
        graph      = FALSE)
      
      biplot <- fviz_pca_biplot(res.pca, 
        axes       = c(1, 2),
        labelsize  = 5,
        col.ind    = EIFsum.TCGA.RNAseq.anno$sample.type, 
        palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"), 
        pointshape = 20,
        pointsize  = 0.75,
        #addEllipses = TRUE, ellipse.level = 0.95,
        label      = "var",
        col.var    = "black", 
        repel      = TRUE,
        title      = paste0("PCA - Biplot (", x,")") ) +
        #scale_x_continuous(breaks = seq(-6, 12, 2), limits=c(-5, 12)) +
        #scale_y_continuous(breaks = seq(-4, 10, 2), limits=c(-4, 10)) +
        theme_classic() + 
        theme(
          plot.background  = element_blank(),
          plot.title       = black_bold_tahoma_16,
          panel.background = element_rect(fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x     = black_bold_tahoma_16,
          axis.title.y     = black_bold_tahoma_16,
          axis.text.x      = black_bold_tahoma_16,
          axis.text.y      = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = c(0, 0),
          legend.justification = c(0,0),
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(biplot)
      ggsave(
        path        = "~/Documents/EIF_output/PCA/TCGA", 
        filename    = paste0(x,"EIFsumPCA.pdf"), 
        plot        = biplot,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      eig <- fviz_eig(res.pca, 
        labelsize = 6,
        geom      = "bar", 
        width     = 0.7, 
        addlabels = TRUE) + 
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background  = element_blank(),
          plot.title       = black_bold_tahoma_16,
          panel.background = element_rect(
            fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x    = black_bold_tahoma_16,
          axis.title.y    = black_bold_tahoma_16,
          axis.text.x     = black_bold_tahoma_16,
          axis.text.y      = black_bold_tahoma_16)
      print(eig)
      ggsave(
        path        = "~/Documents/EIF_output/PCA/TCGA", 
        filename    = paste0(x,"EIFsumeig.pdf"), 
        plot        = eig,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      contribplot <- function(x){
        p <- fviz_contrib(res.pca,
          choice = "var",
          axes   = x,
          top    = 10,
          fill   = "lightblue",
          color  = "black") +
          theme_classic() +
          theme(
            plot.background   = element_blank(),
            plot.title        = black_bold_tahoma_16,
            panel.background  = element_rect(fill   = 'transparent',
              color  = 'black',
              size   = 1),
            axis.title.x      = element_blank(),
            axis.title.y      = black_bold_tahoma_16,
            axis.text.x       = black_bold_tahoma_16_45,
            axis.text.y       = black_bold_tahoma_16)
        print(p)
        ggsave(
          path        = "~/Documents/EIF_output/PCA/TCGA", 
          filename    = paste0("EIFsumcontri",x,".pdf"), 
          plot        = p,
          width       = 8, 
          height      = 4, 
          useDingbats = FALSE)
      }
      lapply(c(1,2), contribplot)

      var <- get_pca_var(res.pca)
      #corrplot(var$contrib, is.corr=FALSE)    
      #fviz_pca_var(res.pca, col.var="contrib")
      pdf(file.path(
        path        = "~/Documents/EIF_output/PCA/TCGA", 
        filename    = paste0(x,"EIFsumPCAcor.pdf")), 
        width       = 9, 
        height      = 9, 
        useDingbats = FALSE)
      corrplot(var$cos2, #cos2 is better than contribute
        is.corr     = FALSE, 
        tl.cex      = 1.5, 
        number.cex  = 1.5, 
        method      = "color", 
        addgrid.col = "gray",
        addCoef.col = "black", 
        tl.col      = "black")
      dev.off()
    }
    plot.pca.factomineR(x)}
  #lapply(disease.list, EIF.PCA.tissue)
  EIFsum.PCA.tissue(tissue)
  }
plot.EIF.TCGA.GTEX.PCA.each(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1", 
                              "PABPC1","MKNK1","MKNK2","MYC"), "Pancreas")

plot.EIF.TCGA.PCA.all <- function (EIF.list) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "primary disease or tissue"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE)
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)] 
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  
  get.EIF.TCGA.GTEX <- function(EIF.list) {
   # EIF.list <- c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","PABPC1","MKNK1","MKNK2", "MYC")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.list, 
      "sample.type",
      "primary.site"),
      drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site %in% c("Brain","Blood","Pancreas"), ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c(
        "Metastatic",
        "Primary Tumor",
        "Normal Tissue",
        "Solid Tissue Normal"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c("Normal Tissue", 
                 "Solid Tissue Normal", 
                 "Primary Tumor", 
                 "Metastatic"),
      labels = c("Healthy Tissue (GTEx)", 
                 "Adjacent Normal Tissue (TCGA)",
                 "Primary Tumor (TCGA)", 
                 "Metastatic Tumor (TCGA)"))
    EIF.TCGA.RNAseq.anno.subset$primary.site <- as.factor(EIF.TCGA.RNAseq.anno.subset$primary.site)
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)
  
  color <- function(){
    n <- 32
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))}
  col_vector <- color()
  
  plot.PCA <- function() {
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type  %in% c("Primary Tumor (TCGA)",
        "Metastatic Tumor (TCGA)"), ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[ ,c(EIF.list,"sample.type","primary.site")]
    EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
  ## remove the last two columns 
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-2)]
    rownames(df1) <- NULL
    plot.pca.factomineR <- function(){
      res.pca <- PCA(df1, 
        scale.unit = TRUE, 
        ncp        = 10, 
        graph      = FALSE)
      
      biplot <- fviz_pca_biplot(res.pca, 
        axes       = c(1, 2),
        labelsize  = 5,
        col.ind    = EIF.TCGA.RNAseq.anno.subset$primary.site, 
        palette    = col_vector, 
        pointshape = 20,
        pointsize  = 0.75,
        #addEllipses = TRUE, 
        label      = "var",
        col.var    = "black", 
        repel      = TRUE) +
        theme_classic() + 
        theme(
          plot.background  = element_blank(),
          plot.title       = black_bold_tahoma_16,
          panel.background = element_rect(fill  = 'transparent',
                                          color = 'black',
                                          size  = 1),
          axis.title.x     = black_bold_tahoma_16,
          axis.title.y     = black_bold_tahoma_16,
          axis.text.x      = black_bold_tahoma_16,
          axis.text.y      = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = "none",
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(biplot)
      ggsave(
        path        = "~/Documents/EIF_output/PCA", 
        filename    = "EIFPCATCGA.pdf", 
        plot        = biplot,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      varplot <- fviz_pca_var(res.pca,
        labelsize  = 5,
        col.ind    = EIF.TCGA.RNAseq.anno.subset$sample.type, 
        palette    = "dark1", 
        pointshape = 20,
        #addEllipses = TRUE, 
        label      = "var",
        col.var    = "black", 
        repel      = TRUE) +
        theme_classic() + 
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(
            fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x      = black_bold_tahoma_16,
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16,
          axis.text.y       = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = c(0.75, 0.93),
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(varplot)
      eig <- fviz_eig(res.pca, 
        labelsize = 6,
        geom      = "bar", 
        width     = 0.7, 
        addlabels = TRUE) + 
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(
            fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x      = black_bold_tahoma_16,
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16,
          axis.text.y       = black_bold_tahoma_16)
      print(eig)
      ggsave(
        path        = "~/Documents/EIF_output/PCA", 
        filename    = "EIFeigTCGA.pdf", 
        plot        = eig,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      var <- get_pca_var(res.pca)
      #fviz_pca_var(res.pca, col.var="contrib")
      
      pdf(file.path(path = "~/Documents/EIF_output/PCA", 
        filename = "EIFTCGAcor.pdf"), 
        width       = 9, 
        height      = 9, 
        useDingbats = FALSE)
      corrplot(var$cos2, #cos2 is better than contribute
        is.corr     = FALSE, 
        tl.cex      = 1.5, 
        number.cex  = 1.5, 
        method      = "color", 
        addgrid.col = "gray",
        addCoef.col = "black", 
        tl.col      = "black")
      dev.off()
      
      contribplot <- function(x){
        fviz_contrib(res.pca,
          choice = "var",
          axes   = x,
          top    = 10,
          fill   = "lightblue",
          color  = "black") +
          theme_minimal() +
          theme(
            plot.background   = element_blank(),
            plot.title        = black_bold_tahoma_16,
            panel.background  = element_rect(fill   = 'transparent',
              color  = 'black',
              size   = 1),
            axis.title.x      = element_blank(),
            axis.title.y      = black_bold_tahoma_16,
            axis.text.x       = black_bold_tahoma_16_45,
            axis.text.y       = black_bold_tahoma_16)
      }
      lapply(c(1,2), contribplot)
    }
    plot.pca.factomineR()}
  plot.PCA()
  
  plot.sum.PCA <- function() {
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
        EIF.TCGA.RNAseq.anno.subset$sample.type  %in% c("Primary Tumor (TCGA)",
          "Metastatic Tumor (TCGA)"), ]
    EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) + 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) -1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[ ,c(EIF.list,"sample.type","primary.site")]  
      ## remove the last two columns 
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-2)]
    rownames(df1) <- NULL
    plot.pca.factomineR <- function(){
      res.pca <- PCA(df1, 
        scale.unit = TRUE, 
        ncp        = 10, 
        graph      = FALSE)
        
      biplot <- fviz_pca_biplot(res.pca, 
          axes       = c(1, 2),
          labelsize  = 5,
          col.ind    = EIF.TCGA.RNAseq.anno.subset$primary.site, 
          palette    = col_vector, 
          pointshape = 20,
          pointsize  = 0.75,
          #addEllipses = TRUE, 
          label      = "var",
          col.var    = "black", 
          repel      = TRUE) +
          theme_classic() + 
          theme(
            plot.background  = element_blank(),
            plot.title       = black_bold_tahoma_16,
            panel.background = element_rect(
              fill  = 'transparent',
              color = 'black',
              size  = 1),
            axis.title.x     = black_bold_tahoma_16,
            axis.title.y     = black_bold_tahoma_16,
            axis.text.x      = black_bold_tahoma_16,
            axis.text.y      = black_bold_tahoma_16,
            legend.title      = element_blank(),
            legend.position   = "none",
            legend.background = element_blank(),
            legend.text       = black_bold_tahoma_16)
      print(biplot)
        ggsave(
          path        = "~/Documents/EIF_output/PCA", 
          filename    = "EIFsumPCATCGA.pdf", 
          plot        = biplot,
          width       = 8, 
          height      = 8, 
          useDingbats = FALSE)
        
      varplot <- fviz_pca_var(res.pca,
          labelsize  = 5,
          col.ind    = EIF.TCGA.RNAseq.anno.subset$sample.type, 
          palette    = "dark1", 
          pointshape = 20,
          #addEllipses = TRUE, 
          label      = "var",
          col.var    = "black", 
          repel      = TRUE) +
          theme_classic() + 
          theme(
            plot.background   = element_blank(),
            plot.title        = black_bold_tahoma_16,
            panel.background  = element_rect(
              fill   = 'transparent',
              color  = 'black',
              size   = 1),
            axis.title.x      = black_bold_tahoma_16,
            axis.title.y      = black_bold_tahoma_16,
            axis.text.x       = black_bold_tahoma_16,
            axis.text.y       = black_bold_tahoma_16,
            legend.title      = element_blank(),
            legend.position   = c(0.75, 0.93),
            legend.background = element_blank(),
            legend.text       = black_bold_tahoma_16)
      print(varplot)
      eig <- fviz_eig(res.pca, 
          labelsize = 6,
          geom      = "bar", 
          width     = 0.7, 
          addlabels = TRUE) + 
          # geom_text(aes(label = res.pca$eig, size = 18)) +
          theme_classic() +
          theme(
            plot.background   = element_blank(),
            plot.title        = black_bold_tahoma_16,
            panel.background  = element_rect(
              fill   = 'transparent',
              color  = 'black',
              size   = 1),
            axis.title.x      = black_bold_tahoma_16,
            axis.title.y      = black_bold_tahoma_16,
            axis.text.x       = black_bold_tahoma_16,
            axis.text.y       = black_bold_tahoma_16)
      print(eig)
      ggsave(
          path        = "~/Documents/EIF_output/PCA", 
          filename    = "EIFsumeigTCGA.pdf", 
          plot        = eig,
          width       = 8, 
          height      = 8, 
          useDingbats = FALSE)
      var <- get_pca_var(res.pca)
        #fviz_pca_var(res.pca, col.var="contrib")
        
      pdf(file.path(path = "~/Documents/EIF_output/PCA", 
          filename = "EIFsumTCGAcor.pdf"), 
          width       = 9, 
          height      = 9, 
          useDingbats = FALSE)
      corrplot(var$cos2, #cos2 is better than contribute
          is.corr     = FALSE, 
          tl.cex      = 1.5, 
          number.cex  = 1.5, 
          method      = "color", 
          addgrid.col = "gray",
          addCoef.col = "black", 
          tl.col      = "black")
      dev.off()
        
      contribplot <- function(x){
          fviz_contrib(res.pca,
            choice = "var",
            axes   = x,
            top    = 10,
            fill   = "lightblue",
            color  = "black") +
            theme_minimal() +
            theme(
              plot.background   = element_blank(),
              plot.title        = black_bold_tahoma_16,
              panel.background  = element_rect(fill   = 'transparent',
                color  = 'black',
                size   = 1),
              axis.title.x      = element_blank(),
              axis.title.y      = black_bold_tahoma_16,
              axis.text.x       = black_bold_tahoma_16_45,
              axis.text.y       = black_bold_tahoma_16)
        }
        lapply(c(1,2), contribplot)
      }
    plot.pca.factomineR()
    }
  plot.sum.PCA()
}
plot.EIF.TCGA.PCA.all(c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
                        "PABPC1", "MKNK1", "MKNK2"))

plot.EIF.GTEX.PCA.all <- function (EIF.list) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE)
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)] 
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
    row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$Row.names
    TCGA.GTEX.sampletype$Row.names <- NULL
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  get.EIF.TCGA.GTEX <- function(EIF.list) {
    #EIF.list <- c("EIF4E", "EIF4G1", "EIF4G2", "EIF4A1","EIF4EBP1", "PABPC1",
    #              "MKNK1","MKNK2", "MTOR", "RPTOR", "RPS6KB1","MYC")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.list, 
                                                        "sample.type",
                                                        "primary.site"),
                                                        drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      !EIF.TCGA.RNAseq.anno.subset$EIF4E == 0, ]
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  !EIF.TCGA.RNAseq.anno.subset$primary.site %in% c("Brain","Blood","Pancreas"), ]
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Metastatic",
        "Primary Tumor",
        "Normal Tissue",
        "Solid Tissue Normal"), ]
    EIF.TCGA.RNAseq.anno.subset$sample.type <- factor(
      EIF.TCGA.RNAseq.anno.subset$sample.type,
      levels = c("Normal Tissue", 
                 "Solid Tissue Normal", 
                 "Primary Tumor", 
                 "Metastatic"),
      labels = c("Healthy Tissue (GTEx)", 
                 "Adjacent Normal Tissue (TCGA)",
                 "Primary Tumor (TCGA)", 
                 "Metastatic Tumor (TCGA)"))
    EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)
  
  color <- function(){
    n <- 32
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))}
  col_vector <- color()

  
  plot.PCA <- function(){
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type == "Healthy Tissue (GTEx)", ]
    EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
    ## remove the last two columns 
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-2)]
    rownames(df1) <- NULL
    
  plot.PCA.prcomp <- function(){
    # the variables should be scaled to have unit variance 
    PCA <- prcomp(df1, scale = TRUE)
    # Extract PC axes for plotting
    PCAvalues <- data.frame(Sample.type = EIF.TCGA.RNAseq.anno.subset$primary.site, 
      PCA$x)
    # Extract loadings of the variables
    PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
    # Plot
    plot.PCA.3D <- function(){
      get_colors <- function(groups, group.col = palette()){
        groups <- as.factor(groups)
        ngrps <- length(levels(groups))
        if(ngrps > length(group.col)) 
          group.col <- rep(group.col, ngrps)
        color <- group.col[as.numeric(groups)]
        names(color) <- as.vector(groups)
        return(color)
      }
      cols <- get_colors(PCAvalues$Sample.type, brewer.pal(n=4, name="Dark2"))
      plot3d(PCAvalues[,2:4], col = cols)
      legend3d("topright", legend =  levels(PCAvalues$Sample.type), 
        pch = 16, col = brewer.pal(n = 4, name="Dark2"), 
        cex = 1)}
    plot.PCA.3D()
    p <- ggplot(PCAvalues,
      aes(x      = PC1,
        y      = PC2,
        colour = Sample.type)) +
      geom_point(size = 0.2) +
      geom_segment(data = PCAloadings,
        aes(
          x     = 0,
          y     = 0,
          xend  = (PC1*5),
          yend  = (PC2*5)),
        arrow = arrow(length = unit(1/3, "picas")),
        color = "black") +
      annotate("text",
        size     = 6,
        fontface = "bold",
        x        = (PCAloadings$PC1*5),
        y        = (PCAloadings$PC2*5),
        label    = PCAloadings$Variables) +
      # stat_n_text(geom = "label") +
      ggtitle("Principal Component Analysis") +
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(
          fill  = 'transparent',
          color = 'black',
          size  = 1),
        axis.title        = black_bold_tahoma_16,
        axis.text         = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0.3, 0.93),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16,
        legend.key        = element_blank())
    print(p)
  }
  plot.PCA.prcomp()
  
  plot.pca.factomineR <- function(){
    res.pca <- PCA(df1, 
      scale.unit = TRUE, 
      ncp        = 10, 
      graph      = FALSE)
    
    biplot <- fviz_pca_biplot(res.pca, 
      axes       = c(1, 2),
      labelsize  = 5,
      col.ind    = EIF.TCGA.RNAseq.anno.subset$primary.site, 
      title      = "PCA - Biplot (all healthy tissues)",
      palette    = col_vector, , 
      pointshape = 20,
      pointsize  = 0.75,
      #addEllipses = TRUE, 
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      theme_classic() + 
      theme(
        plot.background  = element_blank(),
        plot.title       = black_bold_tahoma_16,
        panel.background = element_rect(
          fill  = 'transparent',
          color = 'black',
          size  = 1),
        axis.title.x     = black_bold_tahoma_16,
        axis.title.y     = black_bold_tahoma_16,
        axis.text.x      = black_bold_tahoma_16,
        axis.text.y      = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = "none",
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(biplot)
    ggsave(
      path        = "~/Documents/EIF_output/PCA/GTEX", 
      filename    = "EIFPCAGTEX.pdf", 
      plot        = biplot,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    
    indplot <- fviz_pca_ind(res.pca,
      labelsize   = 5,
      col.ind     = EIF.TCGA.RNAseq.anno.subset$primary.site, 
      palette     = "dark1", 
      #pointshape  = 20,
      pointsize   = 0.95,
      addEllipses = TRUE, 
      # habillage = EIF.TCGA.RNAseq.anno.subset$primary.site,
      label       = "var",
      col.var     = "black", 
      repel       = TRUE) +
      theme_classic() + 
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(fill   = 'transparent',
                                         color  = 'black',
                                         size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0.75, 0.93),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(indplot)
    
    varplot <- fviz_pca_var(res.pca,
      labelsize  = 5,
      col.ind    = EIF.TCGA.RNAseq.anno.subset$sample.type, 
      palette    = "dark1", 
      pointshape = 20,
      #addEllipses = TRUE, 
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      theme_classic() + 
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(
          fill   = 'transparent',
          color  = 'black',
          size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0.75, 0.93),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(varplot)
    eig <- fviz_eig(res.pca, 
      labelsize = 6,
      geom      = "bar", 
      width     = 0.7, 
      addlabels = TRUE) + 
      # geom_text(aes(label = res.pca$eig, size = 18)) +
      theme_classic() +
      theme(
        plot.background   = element_blank(),
        plot.title        = black_bold_tahoma_16,
        panel.background  = element_rect(
          fill   = 'transparent',
          color  = 'black',
          size   = 1),
        axis.title.x      = black_bold_tahoma_16,
        axis.title.y      = black_bold_tahoma_16,
        axis.text.x       = black_bold_tahoma_16,
        axis.text.y       = black_bold_tahoma_16)
    print(eig)
    ggsave(
      path        = "~/Documents/EIF_output/PCA/GTEX", 
      filename    = "EIFPCAeigGTEX.pdf", 
      plot        = eig,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    var <- get_pca_var(res.pca)
    #fviz_pca_var(res.pca, col.var="contrib")
    
    pdf(file.path(path = "~/Documents/EIF_output/PCA/GTEX", 
      filename = "EIFPCAcorGTEX.pdf"), 
      width       = 9, 
      height      = 9, 
      useDingbats = FALSE)
    corrplot(var$cos2, #cos2 is better than contribute
      is.corr     = FALSE, 
      tl.cex      = 1.5, 
      number.cex  = 1.5, 
      method      = "color", 
      addgrid.col = "gray",
      addCoef.col = "black", 
      tl.col      = "black")
    dev.off()
    
    contribplot <- function(x){
      fviz_contrib(res.pca,
        choice = "var",
        axes   = x,
        top    = 10,
        fill   = "lightblue",
        color  = "black") +
        theme_minimal() +
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(
            fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x      = element_blank(),
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16_45,
          axis.text.y       = black_bold_tahoma_16)
    }
    lapply(c(1,2), contribplot)
    }
  plot.pca.factomineR()
  }
  plot.PCA()
  
  plot.sum.PCA <- function(){
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
      EIF.TCGA.RNAseq.anno.subset$sample.type == "Healthy Tissue (GTEx)", ]
    EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**(EIF.TCGA.RNAseq.anno.subset$EIF4E) + 2**(EIF.TCGA.RNAseq.anno.subset$EIF4EBP1) -1)
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[ ,c(EIF.list,"EIF4E+EIF4EBP1","sample.type","primary.site")] 
    ## remove the last two columns 
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-2)]
    rownames(df1) <- NULL
    
    plot.PCA.prcomp <- function(){
      # the variables should be scaled to have unit variance 
      PCA <- prcomp(df1, scale = TRUE)
      # Extract PC axes for plotting
      PCAvalues <- data.frame(Sample.type = EIF.TCGA.RNAseq.anno.subset$primary.site, 
        PCA$x)
      # Extract loadings of the variables
      PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
      # Plot
      plot.PCA.3D <- function(){
        get_colors <- function(groups, group.col = palette()){
          groups <- as.factor(groups)
          ngrps <- length(levels(groups))
          if(ngrps > length(group.col)) 
            group.col <- rep(group.col, ngrps)
          color <- group.col[as.numeric(groups)]
          names(color) <- as.vector(groups)
          return(color)
        }
        cols <- get_colors(PCAvalues$Sample.type, brewer.pal(n=4, name="Dark2"))
        plot3d(PCAvalues[,2:4], col = cols)
        legend3d("topright", legend =  levels(PCAvalues$Sample.type), 
          pch = 16, col = brewer.pal(n = 4, name="Dark2"), 
          cex = 1)}
      plot.PCA.3D()
      p <- ggplot(PCAvalues,
        aes(x      = PC1,
          y      = PC2,
          colour = Sample.type)) +
        geom_point(size = 0.2) +
        geom_segment(data = PCAloadings,
          aes(
            x     = 0,
            y     = 0,
            xend  = (PC1*5),
            yend  = (PC2*5)),
          arrow = arrow(length = unit(1/3, "picas")),
          color = "black") +
        annotate("text",
          size     = 6,
          fontface = "bold",
          x        = (PCAloadings$PC1*5),
          y        = (PCAloadings$PC2*5),
          label    = PCAloadings$Variables) +
        # stat_n_text(geom = "label") +
        ggtitle("Principal Component Analysis") +
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(
            fill  = 'transparent',
            color = 'black',
            size  = 1),
          axis.title        = black_bold_tahoma_16,
          axis.text         = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = c(0.3, 0.93),
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16,
          legend.key        = element_blank())
      print(p)
    }
    plot.PCA.prcomp()
    
    plot.pca.factomineR <- function(){
      res.pca <- PCA(df1, 
        scale.unit = TRUE, 
        ncp        = 10, 
        graph      = FALSE)
      
      biplot <- fviz_pca_biplot(res.pca, 
        axes       = c(1, 2),
        labelsize  = 5,
        col.ind    = EIF.TCGA.RNAseq.anno.subset$primary.site, 
        title      = "PCA - Biplot (all healthy tissues)",
        palette    = col_vector, 
        pointshape = 20,
        pointsize  = 0.75,
        #addEllipses = TRUE, 
        label      = "var",
        col.var    = "black", 
        repel      = TRUE) +
        theme_classic() + 
        theme(
          plot.background  = element_blank(),
          plot.title       = black_bold_tahoma_16,
          panel.background = element_rect(
            fill  = 'transparent',
            color = 'black',
            size  = 1),
          axis.title.x     = black_bold_tahoma_16,
          axis.title.y     = black_bold_tahoma_16,
          axis.text.x      = black_bold_tahoma_16,
          axis.text.y      = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = "none",
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(biplot)
      ggsave(
        path        = "~/Documents/EIF_output/PCA/GTEX", 
        filename    = "EIFsumPCAGTEX.pdf", 
        plot        = biplot,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      indplot <- fviz_pca_ind(res.pca,
        labelsize   = 5,
        col.ind     = EIF.TCGA.RNAseq.anno.subset$primary.site, 
        palette     = "dark1", 
        #pointshape  = 20,
        pointsize   = 0.95,
        addEllipses = TRUE, 
        # habillage = EIF.TCGA.RNAseq.anno.subset$primary.site,
        label       = "var",
        col.var     = "black", 
        repel       = TRUE) +
        theme_classic() + 
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x      = black_bold_tahoma_16,
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16,
          axis.text.y       = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = c(0.75, 0.93),
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(indplot)
      
      varplot <- fviz_pca_var(res.pca,
        labelsize  = 5,
        col.ind    = EIF.TCGA.RNAseq.anno.subset$sample.type, 
        palette    = "dark1", 
        pointshape = 20,
        #addEllipses = TRUE, 
        label      = "var",
        col.var    = "black", 
        repel      = TRUE) +
        theme_classic() + 
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(
            fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x      = black_bold_tahoma_16,
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16,
          axis.text.y       = black_bold_tahoma_16,
          legend.title      = element_blank(),
          legend.position   = c(0.75, 0.93),
          legend.background = element_blank(),
          legend.text       = black_bold_tahoma_16)
      print(varplot)
      eig <- fviz_eig(res.pca, 
        labelsize = 6,
        geom      = "bar", 
        width     = 0.7, 
        addlabels = TRUE) + 
        # geom_text(aes(label = res.pca$eig, size = 18)) +
        theme_classic() +
        theme(
          plot.background   = element_blank(),
          plot.title        = black_bold_tahoma_16,
          panel.background  = element_rect(
            fill   = 'transparent',
            color  = 'black',
            size   = 1),
          axis.title.x      = black_bold_tahoma_16,
          axis.title.y      = black_bold_tahoma_16,
          axis.text.x       = black_bold_tahoma_16,
          axis.text.y       = black_bold_tahoma_16)
      print(eig)
      ggsave(
        path        = "~/Documents/EIF_output/PCA/GTEX", 
        filename    = "EIFsumPCAeigGTEX.pdf", 
        plot        = eig,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      var <- get_pca_var(res.pca)
      #fviz_pca_var(res.pca, col.var="contrib")
      
      pdf(file.path(path = "~/Documents/EIF_output/PCA/GTEX", 
        filename = "EIFsumPCAcorGTEX.pdf"), 
        width       = 9, 
        height      = 9, 
        useDingbats = FALSE)
      corrplot(var$cos2, #cos2 is better than contribute
        is.corr     = FALSE, 
        tl.cex      = 1.5, 
        number.cex  = 1.5, 
        method      = "color", 
        addgrid.col = "gray",
        addCoef.col = "black", 
        tl.col      = "black")
      dev.off()
      
      contribplot <- function(x){
        fviz_contrib(res.pca,
          choice = "var",
          axes   = x,
          top    = 10,
          fill   = "lightblue",
          color  = "black") +
          theme_minimal() +
          theme(
            plot.background   = element_blank(),
            plot.title        = black_bold_tahoma_16,
            panel.background  = element_rect(fill   = 'transparent',
              color  = 'black',
              size   = 1),
            axis.title.x      = element_blank(),
            axis.title.y      = black_bold_tahoma_16,
            axis.text.x       = black_bold_tahoma_16_45,
            axis.text.y       = black_bold_tahoma_16)
      }
      lapply(c(1,2), contribplot)
    }
    plot.pca.factomineR()
  }
  plot.sum.PCA()
}
plot.EIF.GTEX.PCA.all(c("EIF4E", "EIF4G1", "EIF4A1","EIF4EBP1",
                        "PABPC1", "MKNK1","MKNK2"))

plot.EIF.CPTAC.PCA.LUAD <- function(){
  CPTAC.LUAD.Sample <- read_excel(
    "~/Downloads/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx")
  CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample[ ,c("Aliquot (Specimen Label)", "Type")]
  CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample.ID[
    !duplicated(CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`), ]
  row.names(CPTAC.LUAD.Sample.ID) <- CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`
  CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)` <- NULL
  
  CPTAC.LUAD.Proteomics <- fread(
  "~/Downloads/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv",data.table = FALSE)
  EIF.CPTAC.LUAD.Proteomics <- CPTAC.LUAD.Proteomics[CPTAC.LUAD.Proteomics$Gene %in% c("EIF4E", "EIF4G1", "EIF4A1","PABPC1", "MKNK1","MKNK2", "MYC","EIF4EBP1"), ]
  row.names(EIF.CPTAC.LUAD.Proteomics) <- EIF.CPTAC.LUAD.Proteomics$Gene
  EIF.CPTAC.LUAD.Proteomics <- select(EIF.CPTAC.LUAD.Proteomics, -contains("Unshared"))
  
  EIF.CPTAC.LUAD.Proteomics$Gene <- NULL
  EIF.CPTAC.LUAD.Proteomics <- EIF.CPTAC.LUAD.Proteomics[1:(length(EIF.CPTAC.LUAD.Proteomics)-6)]
  EIF.CPTAC.LUAD.Proteomics.t <- data.table::transpose(EIF.CPTAC.LUAD.Proteomics)
  rownames(EIF.CPTAC.LUAD.Proteomics.t) <- colnames(EIF.CPTAC.LUAD.Proteomics)
  colnames(EIF.CPTAC.LUAD.Proteomics.t) <- rownames(EIF.CPTAC.LUAD.Proteomics)
  rownames(EIF.CPTAC.LUAD.Proteomics.t) <- sub(" Log Ratio","",rownames(EIF.CPTAC.LUAD.Proteomics.t)) 
  EIF.CPTAC.LUAD.Proteomics.Sampletype <- merge(EIF.CPTAC.LUAD.Proteomics.t,
                                                CPTAC.LUAD.Sample.ID,
                                                by    = "row.names",
                                                all.x = TRUE)
  rownames(EIF.CPTAC.LUAD.Proteomics.Sampletype) <- EIF.CPTAC.LUAD.Proteomics.Sampletype$Row.names
  EIF.CPTAC.LUAD.Proteomics.Sampletype$Row.names <- NULL
  EIF.CPTAC.LUAD.Proteomics.Sampletype$Type <- factor(
    EIF.CPTAC.LUAD.Proteomics.Sampletype$Type,
    levels = c("Normal", "Tumor"),
    labels = c("Adjacent Normal Tissue (CPTAC)", "Primary Tumor (CPTAC)"))
  EIF.CPTAC.LUAD.Proteomics.Sampletype <- EIF.CPTAC.LUAD.Proteomics.Sampletype[!is.na(EIF.CPTAC.LUAD.Proteomics.Sampletype$Type), ]
  EIF.CPTAC.LUAD.Proteomics.Sampletype <- EIF.CPTAC.LUAD.Proteomics.Sampletype[ , c("EIF4G1", "EIF4A1","EIF4E", "EIF4EBP1","PABPC1", "MKNK1","MKNK2", "MYC","Type")]

 plot.cor.LUAD.tumor <- function(){
  EIF.CPTAC.LUAD.Proteomics.Sampletype.subset <- EIF.CPTAC.LUAD.Proteomics.Sampletype[
    EIF.CPTAC.LUAD.Proteomics.Sampletype$Type == "Primary Tumor (CPTAC)", ]
  EIF.CPTAC.LUAD.Proteomics.Sampletype.subset <- EIF.CPTAC.LUAD.Proteomics.Sampletype.subset[ , c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1","Type")]
  df1 <- EIF.CPTAC.LUAD.Proteomics.Sampletype.subset[1:(length(EIF.CPTAC.LUAD.Proteomics.Sampletype.subset)-1)]
  rownames(df1) <- NULL
  
  cor_5 <- rcorr(as.matrix(df1))
  M <- cor_5$r
  p_mat <- cor_5$P
  
  pdf(file.path(
    path        = "~/Documents/EIF_output/PCA/CPTAC", 
    filename    = "EIFsumLUADcor.pdf"), 
    width       = 6, 
    height      = 6, 
    useDingbats = FALSE)
  
  corrplot(
    M, 
    method      = "color", 
    tl.cex      = 1, tl.srt = 45,
    number.cex  = 1, 
    addgrid.col = "gray",
    addCoef.col = "black", 
    tl.col      = "black",
    #type        = "upper", 
    order       = "FPC", 
    p.mat       = p_mat, 
    sig.level   = 0.05, #insig = "blank" 
  )
  dev.off()
 }
 plot.cor.LUAD.tumor()
  
  
  df1 <- EIF.CPTAC.LUAD.Proteomics.Sampletype[1:(length(EIF.CPTAC.LUAD.Proteomics.Sampletype)-1)]
  rownames(df1) <- NULL
  nb <- missMDA::estim_ncpPCA(df1)
  res.comp <- missMDA::imputePCA(df1,ncp = nb$ncp, nboot = 1000)
  res.pca <- PCA(res.comp$completeObs, 
    scale.unit = TRUE, 
    ncp        = 10, 
    graph = FALSE) 
  
  
  biplot <- fviz_pca_biplot(res.pca, 
    axes       = c(1, 2),
    labelsize  = 5,
    col.ind    = EIF.CPTAC.LUAD.Proteomics.Sampletype$Type, 
    palette    = c("#D55E00","#009E73"), 
    pointshape = 20,
    pointsize  = 0.75,
    title      = "PCA - Biplot (LUAD)",
    label      = "var",
    col.var    = "black", 
    repel      = TRUE) +
    theme_classic() + 
    theme(
      plot.background  = element_blank(),
      plot.title       = black_bold_tahoma_16,
      panel.background = element_rect(
        fill   = 'transparent',
        color  = 'black',
        size   = 1),
      axis.title.x     = black_bold_tahoma_16,
      axis.title.y     = black_bold_tahoma_16,
      axis.text.x      = black_bold_tahoma_16,
      axis.text.y      = black_bold_tahoma_16,
      legend.title      = element_blank(),
      legend.position   = c(0, 0),
      legend.justification = c(0,0),
      legend.background = element_blank(),
      legend.text       = black_bold_tahoma_16)
  print(biplot)
  ggsave(
    path        = "~/Documents/EIF_output/PCA/CPTAC", 
    filename    = "EIFLUADPCA.pdf", 
    plot        = biplot,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  eig <- fviz_eig(res.pca, 
    labelsize = 6,
    geom      = "bar", 
    width     = 0.7, 
    addlabels = TRUE) + 
    # geom_text(aes(label = res.pca$eig, size = 18)) +
    theme_classic() +
    theme(
      plot.background  = element_blank(),
      plot.title       = black_bold_tahoma_16,
      panel.background = element_rect(
        fill   = 'transparent',
        color  = 'black',
        size   = 1),
      axis.title.x    = black_bold_tahoma_16,
      axis.title.y    = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16,
      axis.text.y      = black_bold_tahoma_16)
  print(eig)
  ggsave(
    path        = "~/Documents/EIF_output/PCA/CPTAC", 
    filename    = "EIFLUADEig.pdf", 
    plot        = eig,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  var <- get_pca_var(res.pca)
  #fviz_pca_var(res.pca, col.var="contrib")
  pdf(file.path(
    path        = "~/Documents/EIF_output/PCA/CPTAC", 
    filename    = "EIFLUADcor.pdf"), 
    width       = 9, 
    height      = 9, 
    useDingbats = FALSE)
  corrplot(var$cos2, #cos2 is better than contribute
    is.corr     = FALSE, 
    tl.cex      = 1.5, 
    number.cex  = 1.5, 
    method      = "color", 
    addgrid.col = "gray",
    addCoef.col = "black", 
    tl.col      = "black")
  dev.off()
  }
plot.EIF.CPTAC.PCA.LUAD()

plot.EIF.CPTAC.PCA.BRCA <- function(){
  CPTAC.BRCA.Sample <- read_excel(
    "~/Downloads/S039_Breast_Cancer_Prospective_Collection_Specimens_r1.xlsx")
  CPTAC.BRCA.Sample.ID <- CPTAC.BRCA.Sample[ ,c("Sample Type", "Specimen Label")]
  CPTAC.BRCA.Sample.ID <- CPTAC.BRCA.Sample.ID[
    !duplicated(CPTAC.BRCA.Sample.ID$`Specimen Label`), ]
  CPTAC.BRCA.Sample.ID <- na.omit(CPTAC.BRCA.Sample.ID)
  row.names(CPTAC.BRCA.Sample.ID) <- CPTAC.BRCA.Sample.ID$`Specimen Label`
  CPTAC.BRCA.Sample.ID$`Specimen Label` <- NULL
  
  CPTAC.BRCA.Proteomics <- fread(
    "~/Downloads/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv",
    data.table = FALSE)
  EIF.CPTAC.BRCA.Proteomics <- CPTAC.BRCA.Proteomics[CPTAC.BRCA.Proteomics$Gene %in% c("EIF4E", "EIF4G1", "EIF4A1","EIF4EBP1","PABPC1","MKNK1","MKNK2", "MYC"), ,drop = FALSE]
  EIF.CPTAC.BRCA.Proteomics <- select(EIF.CPTAC.BRCA.Proteomics, -contains("Unshared"))
  row.names(EIF.CPTAC.BRCA.Proteomics) <- EIF.CPTAC.BRCA.Proteomics$Gene
  EIF.CPTAC.BRCA.Proteomics$Gene <- NULL
  EIF.CPTAC.BRCA.Proteomics <- EIF.CPTAC.BRCA.Proteomics[1:(length(EIF.CPTAC.BRCA.Proteomics)-6)]
  EIF.CPTAC.BRCA.Proteomics.t <- data.table::transpose(EIF.CPTAC.BRCA.Proteomics)
  rownames(EIF.CPTAC.BRCA.Proteomics.t) <- colnames(EIF.CPTAC.BRCA.Proteomics)
  colnames(EIF.CPTAC.BRCA.Proteomics.t) <- rownames(EIF.CPTAC.BRCA.Proteomics)
  rownames(EIF.CPTAC.BRCA.Proteomics.t) <- sub(" Log Ratio","",rownames(EIF.CPTAC.BRCA.Proteomics.t)) 
  
  EIF.CPTAC.BRCA.Proteomics.Sampletype <- merge(EIF.CPTAC.BRCA.Proteomics.t,
                                                CPTAC.BRCA.Sample.ID,
                                                by    = "row.names",
                                                all.x = TRUE)
  rownames(EIF.CPTAC.BRCA.Proteomics.Sampletype) <- EIF.CPTAC.BRCA.Proteomics.Sampletype$Row.names
  EIF.CPTAC.BRCA.Proteomics.Sampletype$Row.names <- NULL
  EIF.CPTAC.BRCA.Proteomics.Sampletype$`Sample Type` <- factor(
    EIF.CPTAC.BRCA.Proteomics.Sampletype$`Sample Type`,
    levels = c("Adjacent_Normal", "Tumor"),
    labels = c("Adjacent Normal Tissue (CPTAC)", "Primary Tumor (CPTAC)"))
  EIF.CPTAC.BRCA.Proteomics.Sampletype <- EIF.CPTAC.BRCA.Proteomics.Sampletype[!is.na(EIF.CPTAC.BRCA.Proteomics.Sampletype$`Sample Type`), ]
  EIF.CPTAC.BRCA.Proteomics.Sampletype <- EIF.CPTAC.BRCA.Proteomics.Sampletype[ , c("EIF4G1", "EIF4A1","EIF4E","EIF4EBP1","PABPC1", "MKNK1","MKNK2", "MYC","Sample Type")]
  df1 <- EIF.CPTAC.BRCA.Proteomics.Sampletype[1:(length(EIF.CPTAC.BRCA.Proteomics.Sampletype)-1)]
  rownames(df1) <- NULL
  
  nb <- missMDA::estim_ncpPCA(df1)
  res.comp <- missMDA::imputePCA(df1,ncp = nb$ncp, nboot = 1000)
  res.pca <- PCA(res.comp$completeObs,     
    scale.unit = TRUE, 
    ncp        = 10, 
    graph = FALSE) 
  
  biplot <- fviz_pca_biplot(res.pca, 
    axes       = c(1, 2),
    labelsize  = 5,
    col.ind    = EIF.CPTAC.BRCA.Proteomics.Sampletype$`Sample Type`, 
    palette    = c("#D55E00","#009E73"), 
    pointshape = 20,
    pointsize  = 0.75,
    title      = "PCA - Biplot (BRCA)",
    label      = "var",
    col.var    = "black", 
    repel      = TRUE) +
    theme_classic() + 
    theme(
      plot.background      = element_blank(),
      plot.title           = black_bold_tahoma_16,
      panel.background     = element_rect(
                                          fill   = 'transparent',
                                          color  = 'black',
                                          size   = 1),
      axis.title.x         = black_bold_tahoma_16,
      axis.title.y         = black_bold_tahoma_16,
      axis.text.x          = black_bold_tahoma_16,
      axis.text.y          = black_bold_tahoma_16,
      legend.title         = element_blank(),
      legend.position      = c(0, 0),
      legend.justification = c(0,0),
      legend.background    = element_blank(),
      legend.text          = black_bold_tahoma_16)
  print(biplot)
  ggsave(
    path        = "~/Documents/EIF_output/PCA/CPTAC", 
    filename    = "EIFBRCAPCA.pdf", 
    plot        = biplot,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  eig <- fviz_eig(res.pca, 
    labelsize = 6,
    geom      = "bar", 
    width     = 0.7, 
    addlabels = TRUE) + 
    # geom_text(aes(label = res.pca$eig, size = 18)) +
    theme_classic() +
    theme(
      plot.background  = element_blank(),
      plot.title       = black_bold_tahoma_16,
      panel.background = element_rect(
        fill   = 'transparent',
        color  = 'black',
        size   = 1),
      axis.title.x    = black_bold_tahoma_16,
      axis.title.y    = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16,
      axis.text.y      = black_bold_tahoma_16)
  print(eig)
  ggsave(
    path        = "~/Documents/EIF_output/PCA/CPTAC", 
    filename    = "EIFBRCAEig.pdf", 
    plot        = eig,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  var <- get_pca_var(res.pca)
  #fviz_pca_var(res.pca, col.var="contrib")
  pdf(file.path(
    path        = "~/Documents/EIF_output/PCA/CPTAC", 
    filename    = "EIFBRCAcor.pdf"), 
    width       = 9, 
    height      = 9, 
    useDingbats = FALSE)
  corrplot(var$cos2, #cos2 is better than contribute
    is.corr     = FALSE, 
    tl.cex      = 1.5, 
    number.cex  = 1.5, 
    method      = "color", 
    addgrid.col = "gray",
    addCoef.col = "black", 
    tl.col      = "black")
  dev.off()
  
}
plot.EIF.CPTAC.PCA.BRCA()

################################################################
##  Kaplan-Meier curve with survial and RNASeq data from TCGA ##
################################################################
plot.km.EIF.all.tumors <- function(EIF) {
  pan.TCGA.gene <- function(EIF){
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.RNAseq1 <- TCGA.RNAseq[!duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
   
    ## get OS data ##
    TCGA.OS <- fread(
      "~/Downloads/Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
                        !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE)
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
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
                                       all.x = TRUE)
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease)
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(EIF)
  
  plot.KM <- function(EIF){
    #df <- subset(df, OS.time <= 4000) 
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

    
    KM <- ggplot2::autoplot(
      km, censor = FALSE,
      xlab = "Days",
      ylab = "Survival Probability",
      main = paste0("All TCGA cancer studies (", number, " cases)"),
      #xlim = c(0, 4100),
      color = strata) +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_16,
        axis.title           = black_bold_tahoma_16,
        axis.text            = black_bold_tahoma_16,
        axis.line.x          = element_line(color  = "black"),
        axis.line.y          = element_line(color  = "black"),
        panel.grid           = element_blank(),
        strip.text           = black_bold_tahoma_16,
        legend.text          = black_bold_tahoma_16 ,
        legend.title         = black_bold_tahoma_16 ,
        legend.position      = c(0.98, 0.98),
        legend.justification = c(1, 1)) +
      guides(fill = FALSE) +
      scale_color_manual(
        values = c("red", "blue"),
        name   = paste(EIF, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(paste("Bottom 20%, n =", sub),
                   paste("Top 20%, n =", sub))
      ) +
      #scale_x_continuous(expand = c(0, 0), limits = c(0, 4100)) + 
      #scale_y_continuous(expand = c(0, 0), 
      #                   limits = c(0, 1.05), 
      #                   labels = scales::percent) +
      #geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 10000,
        y        = 0.8,
        label    = paste("log-rank test \n p.val = ", p.val),
        size     = 6.5,
        hjust    = 1,
        fontface = "bold")
    
    print(KM)
    ggsave(
      path        = "~/Documents/EIF_output/KM", 
      filename    = paste(EIF," all tumors KM.pdf"), 
      plot        = KM,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    }
  plot.KM(EIF)
  }
plot.km.EIF.all.tumors("EIF4G1")
lapply(c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1",
         "EIF4EBP1", "PABPC1", "MKNK1","MKNK2", 
         "MTOR", "RPTOR","RPS6KB1", "MYC"), plot.km.EIF.all.tumors)

##
plot.km.EIF.each.tumor <- function(EIF, tumor) {
  pan.TCGA.gene <- function(EIF, tumor){
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.RNAseq1 <- TCGA.RNAseq[!duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    TCGA.OS <- fread(
      "~/Downloads/Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE)
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
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
                                       all.x = TRUE)
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease)
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[TCGA.RNAseq.OS.sampletype$primary.disease %in% tumor, ]
    return(TCGA.RNAseq.OS.sampletype)
  }
  
  plot.KM <- function(EIF, tumor){
    #df <- subset(df, OS.time <= 2000) 
    df <- pan.TCGA.gene(EIF, tumor)
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
  
    KM <- ggplot2::autoplot(
      km, censor = FALSE,
      xlab = "Days",
      ylab = "Survival Probability",
      #xlim = c(0, 2100),
      main = paste0(tumor, " (", number, " cases)")) +
      theme_bw() +
      theme(
        plot.title           = black_bold_tahoma_16,
        axis.title           = black_bold_tahoma_16,
        axis.text            = black_bold_tahoma_16,
        axis.line.x          = element_line(color  = "black"),
        axis.line.y          = element_line(color  = "black"),
        panel.grid           = element_blank(),
        strip.text           = black_bold_tahoma_16,
        legend.text          = black_bold_tahoma_16 ,
        legend.title         = black_bold_tahoma_16 ,
        legend.position      = c(0.98, 0.98),
        legend.justification = c(1, 1)) +
      guides(fill = FALSE) +
      scale_color_manual(
        values = c("red", "blue"),
        name   = paste(EIF, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(paste("Bottom 20%, n =", sub),
                   paste("Top 20%, n =", sub))
      ) +
      #scale_x_continuous(expand = c(0, 0), limits = c(0, 2100)) + 
      #scale_y_continuous(expand = c(0, 0), 
      #                   limits = c(0, 1.05), 
      #                   labels = scales::percent) +
      #geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 7000,
        y        = 0.8,
        label    = paste("log-rank test \n p.val = ", p.val),
        size     = 6.5,
        hjust    = 1,
        fontface = "bold")
    print(KM)
    ggsave(
      path        = "~/Documents/EIF_output/KM", 
      filename    = paste(EIF, tumor,"KM.pdf"), 
      plot        = KM,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
  }
  plot.KM(EIF, tumor)
}
plot.km.EIF.each.tumor(EIF = "EIF4G1", 
                       tumor = c("lung adenocarcinoma"))

lapply(c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1",
         "EIF4EBP1", "PABPC1", "MKNK1","MKNK2", 
         "MTOR", "RPTOR","RPS6KB1", "MYC"), 
         plot.km.EIF.each.tumor, 
         tumor = "lung adenocarcinoma")

## Cox regression model
plot.coxph.EIF.all.tumors <- function(){
  pan.TCGA.gene <- function(EIF){
  ## get TCGA pancancer RNAseq data ##
  # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
  TCGA.RNAseq <- fread(
    "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
    data.table = FALSE)
  # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
  TCGA.RNAseq1 <- TCGA.RNAseq[!duplicated(TCGA.RNAseq$sample),
    !duplicated(colnames(TCGA.RNAseq))]
  row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
  TCGA.RNAseq1$sample <- NULL
  TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
  TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
  rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
  colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
  colnames(TCGA.RNAseq_transpose) <- EIF
  
  ## get OS data ##
  TCGA.OS <- fread(
    "~/Downloads/Survival_SupplementalTable_S1_20171025_xena_sp", 
    data.table = FALSE)
  TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
    !duplicated(colnames(TCGA.OS))]
  row.names(TCGA.OS1) <- TCGA.OS1$sample
  TCGA.OS1$sample <- NULL
  TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
  
  ## get sample type data ##
  TCGA.sampletype <- readr::read_tsv(
    "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
  row.names(TCGA.sampletype) <- TCGA.sampletype$sample
  TCGA.sampletype$sample <- NULL
  TCGA.sampletype$sample_type_id <- NULL
  colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
  
  ## combine OS and sample type data ##
  TCGA.OS.sampletype <- merge(TCGA.OS1,
                              TCGA.sampletype,
                              by    = "row.names",
                              all.x = TRUE)
  TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
  row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
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
                                     all.x = TRUE)
  TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
  row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
  TCGA.RNAseq.OS.sampletype$Row.names <- NULL
  ## remove all rows with NA in the primary disease section
  TCGA.RNAseq.OS.sampletype <-
    TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
  TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
    TCGA.RNAseq.OS.sampletype$primary.disease)
  return(TCGA.RNAseq.OS.sampletype)
}
  df <- pan.TCGA.gene(c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1",
                        "EIF4EBP1","PABPC1", "MKNK1","MKNK2", 
                        "MTOR", "RPTOR","RPS6KB1", "MYC"))
  #df$`EIF4E+EIF4EBP1` <- log2(2**df$EIF4E + 2**df$EIF4EBP1 -2 + 1)
  #df <- df[c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1","EIF4E+EIF4EBP1",
  #           "EIF4EBP1","PABPC1", "MKNK1","MKNK2", 
  #           "MTOR", "RPTOR","RPS6KB1", "MYC","OS","OS.time")]


  # lapply(univ_models, forest_model)
  # Use survivalAnalysis package to draw forest plot of multiple univariate #
  covariate_names <- c(EIF4E = "EIF4E", EIF4G1 = "EIF4G1", 
                      #`EIF4E+EIF4EBP1` = "EIF4E+EIF4EBP1",
                       EIF4G2 = "EIF4G2", EIF4A1 = "EIF4A1",
                       EIF4EBP1 = "EIF4EBP1", PABPC1 = "PABPC1", 
                       MKNK1 = "MKNK1", MKNK2 = "MKNK2", 
                       MTOR = "MTOR", RPTOR = "RPTOR",
                       RPS6KB1 = "RPS6KB1", MYC = "MYC")
  
  plot.univariate <- function(){
    result <- map(vars(EIF4E, EIF4G1, EIF4G2, EIF4A1,#EIF4E+EIF4EBP1,
      EIF4EBP1, PABPC1, MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC), 
      function(by){analyse_multivariate(df,
        vars(OS.time, OS),
        covariates = list(by), # covariates expects a list
        covariate_name_dict = covariate_names)
    })

  HR.table <- function (x) {
    list <- as.data.frame(result[[x]]["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
    return(list)
  }
  a <- lapply(c(1:12), HR.table)
  b <- do.call(rbind.data.frame, a)
  
  # Testing proportional Hazards assumption on univariate analysis
  univ_formulas <- sapply(covariate_names,
    function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
  univ_models <- lapply(univ_formulas, 
                        function(x){cox.zph(coxph(x, data = df))})
  coxassump <- function(x){
    c <- print(univ_models[[x]])
    return(c)}
  univ_results <- lapply(covariate_names,coxassump)
  d <- do.call(rbind.data.frame, univ_results)
  e <- d[-grep("GLOBAL", rownames(d)), ]
  rownames(e) <- gsub("\\..*","",rownames(e))
  f <- e[ ,3, drop=FALSE]
  colnames(f) <- "pinteraction"

  data <- merge(b, f, by.x= "factor.id", by.y="row.names")
  data[ ,4:11] <- round(data[ ,4:11], digits = 3)
  data[ ,4:6] <- round(data[ ,4:6], digits = 2)
  data$np <- 10235
  data <- as.data.frame(data)
  data$HRCI <- paste0(data$HR," (", data$Lower_CI, "-", data$Upper_CI, ")")
  data$p[data$p < 0.001]<-"<0.001"
  data$pinteraction[data$pinteraction < 0.001]<-"<0.001"
  data <-data[order(data$HR),]
  tabletext1 <- cbind(
    c("Gene", data$factor.id), 
    c("No. of\nPatients", data$np), 
    c("Hazard Ratio\n(95% CI)", data$HRCI),
    c("P Value", data$p),
    c("P Value for\nInteraction", data$pinteraction))
  
  
  pdf(file   = '~/Documents/EIF_output/Cox/EIFUniCox.pdf',
      width  = 10, 
      height = 8,onefile=F) 
  p <- forestplot(
    labeltext  =tabletext1, 
    graph.pos  = 3, graphwidth = unit(6, "cm"),
    hrzl_lines =list(
      "1" = gpar(lwd=1, col="black"),
      "2" = gpar(lwd=1, col="black")),
    #  "3.75" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
    #  "3" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
    #  "7" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
    #  "9" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922")),
    mean  = c(NA,data$HR), 
    lower = c(NA,data$Lower_CI), 
    upper = c(NA,data$Upper_CI),
    title = "Univariate Cox proportional-hazards regression analysis (all tumor types)",
    xlab  = "     <---Good prognosis---    ---Poor prognosis--->",
    txt_gp = fpTxtGp(label =gpar(cex=1.2),
                     ticks = gpar(cex = 1.2),
                     xlab  = gpar(cex = 1.2),
                     title = gpar(cex = 1.2)),
    col    = fpColors(box = "black", lines = "black"),
    xticks = c(0.6,0.8,1,1.2,1.4,1.6, 1.8),
    clip  = c(0.6, 1.8),# range of x axis
    zero  = 1, 
    cex   = 1.2, 
    lineheight = "auto", #height of the graph
    boxsize    = 0.2, 
    colgap     = unit(3,"mm"), #the gap between column
    lwd.ci = 2, 
    ci.vertices = FALSE,
    ci.vertices.height = 0.02,
  new_page = getOption("forestplot_new_page", FALSE))
  dev.off()
  print(p)
  
  
  p2 <- map(vars(EIF4E, EIF4G1, EIF4G2, EIF4A1,EIF4E+EIF4EBP1,
                EIF4EBP1, PABPC1, MKNK1, MKNK2, 
                MTOR, RPTOR, RPS6KB1, MYC), function(by)
  {
    analyse_multivariate(df,
      vars(OS.time, OS),
      covariates = list(by), # covariates expects a list
      covariate_name_dict = covariate_names)
  }) %>%
    forest_plot(
      factor_labeller = covariate_names,
      endpoint_labeller = c(OS.time = "OS"),
      orderer           = ~order(HR),
      labels_displayed  = c("factor"),
      values_displayed  = c("HR", "CI", "p"),
      value_headers     = c(HR = "HR", CI = "95%CI", p = "p", n = "N"),
      relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
      label_headers     = c(factor = "Gene"),
      title             = "Univariate Cox proportional-hazards regression analysis\nAll Tumors (N=10295)",
      HR_x_breaks       = seq(0.7, 1.8, 0.2), 
      HR_x_limits       = c(0.7, 1.8),
      ggtheme           = ggplot2::theme_bw(base_size = 7))
  print(p2)
  }
  plot.univariate()

  plot.multivariate <- function(){
    df %>% 
      analyse_multivariate(vars(OS.time, OS), covariates = vars(EIF4E, EIF4A1, EIF4G1, EIF4G2, PABPC1, EIF4EBP1, MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC), covariate_name_dict = covariate_names) -> result1
    data <- as.data.frame(result1["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
  
    # Testing proportional Hazards assumption on univariate analysis
    mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4E + EIF4A1 + EIF4G1 + EIF4EBP1 +       EIF4G2 + PABPC1 + MKNK1 + MKNK2 + 
        MTOR + RPTOR + RPS6KB1 + MYC, data = df)
    test.ph <- cox.zph(mv_fit)
    test <- print(test.ph) 
    f <- test[ ,3, drop=FALSE]
    colnames(f) <- "pinteraction"
    
    data <- merge(data, f, by.x= "factor.id", by.y="row.names")
    data[ ,4:11] <- round(data[ ,4:11], digits = 3)
    data[ ,4:6] <- round(data[ ,4:6], digits = 2)
    data$np <- 10235
    data <- as.data.frame(data)
    data$HRCI <- paste0(data$HR," (", data$Lower_CI, "-", data$Upper_CI, ")")
    data$p[data$p < 0.001]<-"<0.001"
    data$pinteraction[data$pinteraction < 0.001]<-"<0.001"
    data <-data[order(data$HR),]
    tabletext1 <- cbind(
      c("Gene", data$factor.id), 
      c("No. of\nPatients", data$np), 
      c("Hazard Ratio\n(95% CI)", data$HRCI),
      c("P Value", data$p),
      c("P Value for\nInteraction", data$pinteraction))
    
  
    pdf(file    = '~/Documents/EIF_output/Cox/EIFmultiCox.pdf',
        width   = 10, 
        height  = 8,
        onefile =F) 
    p <- forestplot(
      labeltext  =tabletext1, 
      graph.pos  = 3, graphwidth = unit(6, "cm"),
      hrzl_lines =list(
        "1" = gpar(lwd=1, col="black"),
        "2" = gpar(lwd=1, col="black")),
      #  "3.75" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
      #  "3" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
      #  "7" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
      #  "9" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922")),
      mean  = c(NA,data$HR), 
      lower = c(NA,data$Lower_CI), 
      upper = c(NA,data$Upper_CI),
      title = "Multivariate Cox proportional-hazards regression analysis (all tumor types)",
      xlab  = "     <---Good prognosis---    ---Poor prognosis--->",
      txt_gp = fpTxtGp(label=gpar(cex=1.2),
        ticks  = gpar(cex = 1.2),
        xlab   = gpar(cex = 1.2),
        title  = gpar(cex = 1.2)),
      col    = fpColors(box = "black", lines = "black"),
      xticks = c(0.6,0.8,1,1.2,1.4,1.6, 1.8),
      clip  = c(0.6, 1.8),# range of x axis
      zero  = 1, 
      cex   = 1.2, 
      lineheight = "auto", #height of the graph
      boxsize    = 0.2, 
      colgap     = unit(3,"mm"), #the gap between column
      lwd.ci = 2, 
      ci.vertices = FALSE,
      ci.vertices.height = 0.02,
      new_page = getOption("forestplot_new_page", FALSE))
    dev.off()
    print(p)

  
    p2 <- forest_plot(result1,
      factor_labeller   = covariate_names, #label the subgroups.
      endpoint_labeller = c(OS.time="OS"), #label the endpoint.
      orderer           = ~order(HR), #order by hazard ratio
      labels_displayed  = c("factor"),
      values_displayed  = c("HR", "CI", "p"),
      value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
      ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
      relative_widths   = c(0.4, 1.8, 1.2), #more space for the plot, less space for the tables
      label_headers     = c(factor = "Gene", n = "Sample size"),
      title             = "Multivariate Cox proportional-hazards regression analysis\nAll Tumors (N=10295)",
      HR_x_breaks       = seq(0.6, 1.7, 0.2), 
      HR_x_limits       = c(0.6, 1.7)) #more breaks on the X axis
    print(p2)

  df %>% 
    analyse_multivariate(vars(OS.time, OS),
      covariates = vars(EIF4G1, EIF4A1, EIF4E, EIF4EBP1,PABPC1, 
                        MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC),
      covariate_name_dict = covariate_names) -> result2
  p3 <- forest_plot(result2,
    factor_labeller   = covariate_names, #label the subgroups.
    endpoint_labeller = c(OS.time="OS"), #label the endpoint.
    orderer           = ~order(HR), #order by hazard ratio
    labels_displayed  = c("factor"),
    values_displayed  = c("HR", "CI", "p"),
    value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
    ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
    relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
    label_headers     = c(factor = "Gene", n = "Sample size"),
    title             = "Multivariate Cox proportional-hazards regression analysis\nAll Tumors (N=10295)",
    HR_x_breaks       = seq(0.6, 1.7, 0.2), 
    HR_x_limits       = c(0.6, 1.7)) #more breaks on the X axis
  print(p3)
  }
  plot.multivariate()

  # Aalen’s additive regression model #
  # provide detailed information about temporal influence of each of the covariates not available in Cox’s model
  aa_fit <-aareg(Surv(OS.time, OS) ~ EIF4E + EIF4G1 + EIF4G2 + EIF4A1 + EIF4EBP1 + PABPC1 + MKNK1 + MKNK2 + MTOR + RPTOR + RPS6KB1 + MYC, data = df)
  summary(aa_fit)
  p3 <- ggplot2::autoplot(aa_fit)
  print(p3)
  # Lasso with an elastic-net penalty #
  df1 <- na.omit(df)
  df1 <- df1[df1$OS.time != 0, ] ## must remove 0 in OS time
  x <- as.matrix(df1[,1:12])
  y <- as.double(df1$OS.time)
  STATUS <- as.double(df1$OS)
  surv <- Surv(y,STATUS)

  fit.lasso <- glmnet(x, surv, family = "cox", alpha = 1)
  fit.ridge <- glmnet(x, surv, family = "cox", alpha = 0)
  fit.elnet <- glmnet(x, surv, family = "cox", alpha = .5)
  # 10-fold Cross validation for each alpha = 0, 0.1, ... , 0.9, 1.0
  fit.lasso.cv <- cv.glmnet(x, surv, family = "cox", alpha = 1)
  fit.ridge.cv <- cv.glmnet(x, surv, family = "cox", alpha = 0)
  fit.elnet.cv <- cv.glmnet(x, surv, family = "cox", alpha = .5)
  for (i in 0:10) {
    assign(paste("fit", i, sep = ""), 
      cv.glmnet(x, surv, family = "cox", alpha = i/10))
  }
  # par(mfrow=c(3,2))
  # For plotting options, type '?plot.glmnet' in R console
  pdf(file.path(
    path        = "~/Documents/EIF_output/Cox", 
    filename    = "EIFlasso.pdf"), 
    width       = 6, 
    height      = 6, 
    useDingbats = FALSE)
  
  plot_glmnet(fit.lasso, main="LASSO (Alpha = 1)")
  dev.off()
  plot(fit10, main="Lasso penalty\n\n")
  
  plot_glmnet(fit.ridge, main="Ridge")
  plot(fit0, main="Ridge penalty\n\n")
  
  plot_glmnet(fit.elnet, main="Elastic Net")
  plot(fit5, main="Elastic Net")

  # correlation plot
  my_data <- df[, c(1:13)]
  cor_5 <- rcorr(as.matrix(my_data))
  M <- cor_5$r
  p_mat <- cor_5$P
  pdf(file.path(
    path        = "~/Documents/EIF_output/Cox", 
    filename    = "EIFcormatrix.pdf"), 
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  corrplot(
    M, 
    method      = "color", 
    tl.cex      = 1, tl.srt = 45,
    number.cex  = 1, 
    addgrid.col = "gray",
    addCoef.col = "black", 
    tl.col      = "black",
    #type        = "upper", 
    order       = "FPC", 
    p.mat       = p_mat, 
    sig.level   = 0.05, #insig = "blank" 
    )
  dev.off()
  }
plot.coxph.EIF.all.tumors()

plot.coxph.EIF.lung.tumor <- function(tumor){
  pan.TCGA.gene <- function(EIF, tumor){
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.RNAseq1 <- TCGA.RNAseq[!duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    TCGA.OS <- fread(
      "~/Downloads/Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE)
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
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
                                       all.x = TRUE)
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease)
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[TCGA.RNAseq.OS.sampletype$primary.disease %in% tumor, ]
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(c("EIF4E", "EIF4G1", "EIF4G2", "EIF4A1", 
                        "EIF4EBP1", "PABPC1", "MKNK1", "MKNK2", 
                        "MTOR", "RPTOR","RPS6KB1", "MYC"), tumor)
  #df$`EIF4E+EIF4EBP1` <- log2(2**df$EIF4E + 2**df$EIF4EBP1 -2 + 1)
  df <- df[c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1",
             "EIF4EBP1","PABPC1", "MKNK1","MKNK2", 
             "MTOR", "RPTOR","RPS6KB1", "MYC", #"EIF4E+EIF4EBP1",
             "OS","OS.time")]
# Use survivalAnalysis package to draw forest plot of multiple univariate #
  covariate_names <- c(EIF4E = "EIF4E", EIF4G1 = "EIF4G1",
                       EIF4G2 = "EIF4G2", EIF4A1 = "EIF4A1",
                       EIF4EBP1 = "EIF4EBP1", PABPC1 = "PABPC1", 
                       MKNK1 = "MKNK1", MKNK2 = "MKNK2", 
                       #`EIF4E+EIF4EBP1` = "EIF4E+EIF4EBP1",
                       MTOR = "MTOR", RPTOR = "RPTOR",
                       RPS6KB1 = "RPS6KB1", MYC = "MYC")
  
  plot.univariate <- function(){
    result <- map(vars(EIF4E, EIF4G1, EIF4G2, EIF4A1,#EIF4E+EIF4EBP1,
      EIF4EBP1, PABPC1, MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC), 
      function(by){analyse_multivariate(df,
        vars(OS.time, OS),
        covariates = list(by), # covariates expects a list
        covariate_name_dict = covariate_names)
      })
    
    HR.table <- function (x) {
      list <- as.data.frame(result[[x]]["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
      return(list)
    }
    a <- lapply(c(1:12), HR.table)
    b <- do.call(rbind.data.frame, a)
    
    # Testing proportional Hazards assumption on univariate analysis
    univ_formulas <- sapply(covariate_names,
      function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
    univ_models <- lapply(univ_formulas, 
      function(x){cox.zph(coxph(x, data = df))})
    coxassump <- function(x){
      c <- print(univ_models[[x]])
      return(c)}
    univ_results <- lapply(covariate_names,coxassump)
    d <- do.call(rbind.data.frame, univ_results)
    e <- d[-grep("GLOBAL", rownames(d)), ]
    rownames(e) <- gsub("\\..*","",rownames(e))
    f <- e[ ,3, drop=FALSE]
    colnames(f) <- "pinteraction"
    
    data <- merge(b, f, by.x= "factor.id", by.y="row.names")
    data[ ,4:11] <- round(data[ ,4:11], digits = 3)
    data[ ,4:6] <- round(data[ ,4:6], digits = 2)
    data$np <- 517
    data <- as.data.frame(data)
    data$HRCI <- paste0(data$HR," (", data$Lower_CI, "-", data$Upper_CI, ")")
    data$p[data$p < 0.001]<-"<0.001"
    data$pinteraction[data$pinteraction < 0.001]<-"<0.001"
    data <-data[order(data$HR),]
    tabletext1 <- cbind(
      c("Gene", data$factor.id), 
      c("No. of\nPatients", data$np), 
      c("Hazard Ratio\n(95% CI)", data$HRCI),
      c("P Value", data$p),
      c("P Value for\nInteraction", data$pinteraction))
    
    pdf(file    = paste0("~/Documents/EIF_output/Cox/", 
                         tumor, 
                         "EIFuniCox.pdf"),
        width   = 10, 
        height  = 8,
        onefile = F) 
    p <- forestplot(
      labeltext  = tabletext1, 
      graph.pos  = 3, graphwidth = unit(6, "cm"),
      hrzl_lines =list(
        "1" = gpar(lwd=1, col="black"),
        "2" = gpar(lwd=1, col="black")),
      #  "3.75" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
      #  "3" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
      #  "7" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
      #  "9" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922")),
      mean  = c(NA,data$HR), 
      lower = c(NA,data$Lower_CI), 
      upper = c(NA,data$Upper_CI),
      title = "Univariate Cox proportional-hazards regression analysis (lung adenocarcinoma)",
      xlab  = "<---Good prognosis---    ---Poor prognosis--->",
      txt_gp = fpTxtGp(label = gpar(cex = 1.2),
                       ticks = gpar(cex = 1.2),
                       xlab  = gpar(cex = 1.2),
                       title = gpar(cex = 1.2)),
      col    = fpColors(box = "black", lines = "black"),
      xticks = c(0.5,0.9,1.3,1.7,2.1), 
      #xlog = 0,
      clip   = c(0.5, 2.1),# range of x axis
      zero   = 1, 
      cex    = 1.2, 
      lineheight = "auto", #height of the graph
      boxsize    = 0.2, 
      colgap     = unit(3,"mm"), #the gap between column
      lwd.ci = 2, 
      ci.vertices = FALSE,
      ci.vertices.height = 0.02,
      new_page = getOption("forestplot_new_page", FALSE))
    dev.off()
    print(p)
  }
  plot.univariate()
  
  plot.multivariate <- function(){
    df %>% 
      analyse_multivariate(vars(OS.time, OS), covariates = vars(EIF4E, EIF4A1, EIF4G1, EIF4G2, PABPC1, EIF4EBP1, MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC), covariate_name_dict = covariate_names) -> result1
    data <- as.data.frame(result1["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
    
    # Testing proportional Hazards assumption on univariate analysis
    mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4E + EIF4A1 + EIF4G1 + EIF4EBP1 +       EIF4G2 + PABPC1 + MKNK1 + MKNK2 + 
        MTOR + RPTOR + RPS6KB1 + MYC, data = df)
    test.ph <- cox.zph(mv_fit)
    test <- print(test.ph) 
    f <- test[ ,3, drop=FALSE]
    colnames(f) <- "pinteraction"
    
    data <- merge(data, f, by.x= "factor.id", by.y="row.names")
    data[ ,4:11] <- round(data[ ,4:11], digits = 3)
    data[ ,4:6] <- round(data[ ,4:6], digits = 2)
    data$np <- 517
    data <- as.data.frame(data)
    data$HRCI <- paste0(data$HR," (", data$Lower_CI, "-", data$Upper_CI, ")")
    data$p[data$p < 0.001]<-"<0.001"
    data$pinteraction[data$pinteraction < 0.001]<-"<0.001"
    data <-data[order(data$HR),]
    tabletext1 <- cbind(
      c("Gene", data$factor.id), 
      c("No. of\nPatients", data$np), 
      c("Hazard Ratio\n(95% CI)", data$HRCI),
      c("P Value", data$p),
      c("P Value for\nInteraction", data$pinteraction))
    
    
    pdf(
      file   = paste0("~/Documents/EIF_output/Cox/", 
               tumor, 
              "EIFmultiCox.pdf"),
      width  = 10, 
      height = 8,onefile=F) 
    p <- forestplot(
      labeltext  =tabletext1, 
      graph.pos  = 3, graphwidth = unit(6, "cm"),
      hrzl_lines =list(
        "1" = gpar(lwd=1, col="black"),
        "2" = gpar(lwd=1, col="black")),
      #  "3.75" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
      #  "3" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
      #  "7" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
      #  "9" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922")),
      mean  = c(NA,data$HR), 
      lower = c(NA,data$Lower_CI), 
      upper = c(NA,data$Upper_CI),
      title = "Multivariate Cox proportional-hazards regression analysis (lung adenocarcinoma)",
      xlab  = "<---Good prognosis---    ---Poor prognosis--->",
      txt_gp = fpTxtGp(label=gpar(cex=1.2),
      ticks  = gpar(cex = 1.2),
      xlab   = gpar(cex = 1.2),
      title  = gpar(cex = 1.2)),
      col    = fpColors(box = "black", lines = "black"),
      xticks = c(0.5,0.8,1.1,1.4,1.7,2,2.3,2.6),
      clip   = c(0.5, 2.6),# range of x axis
      zero   = 1, 
      cex    = 1.2, 
      lineheight = "auto", #height of the graph
      boxsize    = 0.2, 
      colgap     = unit(3,"mm"), #the gap between column
      lwd.ci = 2, 
      ci.vertices = FALSE,
      ci.vertices.height = 0.02,
      new_page = getOption("forestplot_new_page", FALSE))
    dev.off()
    print(p)
    
    
    p2 <- forest_plot(result1,
      factor_labeller   = covariate_names, #label the subgroups.
      endpoint_labeller = c(OS.time="OS"), #label the endpoint.
      orderer           = ~order(HR), #order by hazard ratio
      labels_displayed  = c("factor"),
      values_displayed  = c("HR", "CI", "p"),
      value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
      ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
      relative_widths   = c(0.4, 1.8, 1.2), #more space for the plot, less space for the tables
      label_headers     = c(factor = "Gene", n = "Sample size"),
      title             = "Multivariate Cox proportional-hazards regression analysis (lung adenocarcinoma)",
      HR_x_breaks       = seq(0.6, 1.7, 0.2), 
      HR_x_limits       = c(0.6, 1.7)) #more breaks on the X axis
    print(p2)
    
    df %>% 
      analyse_multivariate(vars(OS.time, OS),
        covariates = vars(EIF4G1, EIF4A1, EIF4E, EIF4EBP1,PABPC1, 
          MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC),
        covariate_name_dict = covariate_names) -> result2
    p3 <- forest_plot(result2,
      factor_labeller   = covariate_names, #label the subgroups.
      endpoint_labeller = c(OS.time="OS"), #label the endpoint.
      orderer           = ~order(HR), #order by hazard ratio
      labels_displayed  = c("factor"),
      values_displayed  = c("HR", "CI", "p"),
      value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
      ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
      relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
      label_headers     = c(factor = "Gene", n = "Sample size"),
      title             = "Multivariate Cox proportional-hazards regression analysis (lung adenocarcinoma)",
      HR_x_breaks       = seq(0.6, 1.7, 0.2), 
      HR_x_limits       = c(0.6, 1.7)) #more breaks on the X axis
    print(p3)
  }
  plot.multivariate()
  
  
  p <- map(vars(EIF4E, EIF4G1, EIF4G2, EIF4A1,EIF4E, EIF4EBP1,
  EIF4EBP1, PABPC1, MKNK1, MKNK2, 
  MTOR, RPTOR, RPS6KB1, MYC), function(by)
  {
    analyse_multivariate(df,
      vars(OS.time, OS),
      covariates = list(by), # covariates expects a list
      covariate_name_dict = covariate_names)
  }) %>%
  forest_plot(
    factor_labeller   = covariate_names,
    endpoint_labeller = c(OS.time = "OS"),
    orderer           = ~order(HR),
    labels_displayed  = c("factor"),
    values_displayed  = c("HR", "CI", "p"),
    value_headers     = c(HR = "HR", CI = "95%CI", p = "p", n = "N"),
    relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
    label_headers     = c(factor = "Gene"),
    title             = "Univariate Cox proportional-hazards regression analysis\nLung Adenocarcinoma (N=517)",
    HR_x_breaks       = seq(0.5, 2.1, 0.2), 
    HR_x_limits       = c(0.5, 2.1),
    ggtheme           = ggplot2::theme_bw(base_size = 7))
  print(p)

  df %>% 
    analyse_multivariate(vars(OS.time, OS),
      covariates = vars(EIF4G1, EIF4G2, EIF4A1, EIF4E, PABPC1, 
                        MKNK1, MTOR, RPTOR, RPS6KB1),
      covariate_name_dict = covariate_names) -> result1
  p2 <- forest_plot(result1,
    factor_labeller   = covariate_names, #label the subgroups.
    endpoint_labeller = c(OS.time="OS"), #label the endpoint.
    orderer           = ~order(HR), #order by hazard ratio
    labels_displayed  = c("factor"),
    values_displayed  = c("HR", "CI", "p"),
    value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
    ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
    relative_widths   = c(0.4, 1.8, 1.2), #more space for the plot, less space for the tables
    label_headers     = c(factor = "Gene", n = "Sample size"),
    title             = "Multivariate Cox proportional-hazards regression analysis\nLung Adenocarcinoma (N=517)",
    HR_x_breaks       = seq(0.5, 2.7, 0.3), 
    HR_x_limits       = c(0.5, 2.7)) #more breaks on the X axis
  print(p2)


  df %>% 
    analyse_multivariate(vars(OS.time, OS),
      covariates = vars(EIF4G1, EIF4G2, EIF4A1,EIF4E, EIF4EBP1, PABPC1, 
                        MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC),
      covariate_name_dict = covariate_names) -> result2
  p3 <- forest_plot(result2,
    factor_labeller   = covariate_names, #label the subgroups.
    endpoint_labeller = c(OS.time="OS"), #label the endpoint.
    orderer           = ~order(HR), #order by hazard ratio
    labels_displayed  = c("factor"),
    values_displayed  = c("HR", "CI", "p"),
    value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
    ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
    relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
    label_headers     = c(factor = "Gene", n = "Sample size"),
    title             = "Multivariate Cox proportional-hazards regression analysis\nLung Adenocarcinoma (N=517)",
    HR_x_breaks       = seq(0.5, 2.7, 0.3), 
    HR_x_limits       = c(0.5, 2.7)) #more breaks on the X axis
  print(p3)


  # correlation plot
  my_data <- df[, c(1:13)]
  cor_5 <- rcorr(as.matrix(my_data))
  M <- cor_5$r
  p_mat <- cor_5$P
  pdf(file.path(
    path        = "~/Documents/EIF_output/Cox", 
    filename    = paste0(tumor, "EIFcormatrix.pdf")), 
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  corrplot(
    M, 
    method      = "color", 
    tl.cex      = 1, tl.srt = 45,
    number.cex  = 1, 
    addgrid.col = "gray",
    addCoef.col = "black", 
    tl.col      = "black",
    #type        = "upper", 
    order       = "FPC", 
    p.mat       = p_mat, 
    sig.level   = 0.05, #insig = "blank" 
    )
  dev.off()
  }
plot.coxph.EIF.lung.tumor(c("lung adenocarcinoma"))

plot.coxph.EIF.breast.tumor <- function(tumor){
  pan.TCGA.gene <- function(EIF, tumor){
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.RNAseq1 <- TCGA.RNAseq[!duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    TCGA.OS <- fread(
      "~/Downloads/Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
      TCGA.sampletype,
      by    = "row.names",
      all.x = TRUE)
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
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
      all.x = TRUE)
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease)
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[TCGA.RNAseq.OS.sampletype$primary.disease %in% tumor, ]
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(c("EIF4E", "EIF4G1", "EIF4G2", "EIF4A1", 
    "EIF4EBP1", "PABPC1", "MKNK1", "MKNK2", 
    "MTOR", "RPTOR","RPS6KB1", "MYC"), tumor)
  df$`EIF4E+EIF4EBP1` <- log2(2**df$EIF4E + 2**df$EIF4EBP1 -2 + 1)
  df <- df[c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1",
    "EIF4EBP1","PABPC1", "MKNK1","MKNK2", 
    "MTOR", "RPTOR","RPS6KB1", "MYC","EIF4E+EIF4EBP1",
    "OS","OS.time")]
  # Use survivalAnalysis package to draw forest plot of multiple univariate #
  covariate_names <- c(EIF4E = "EIF4E", EIF4G1 = "EIF4G1",
    EIF4G2 = "EIF4G2", EIF4A1 = "EIF4A1",
    EIF4EBP1 = "EIF4EBP1", PABPC1 = "PABPC1", 
    MKNK1 = "MKNK1", MKNK2 = "MKNK2", 
    `EIF4E+EIF4EBP1` = "EIF4E+EIF4EBP1",
    MTOR = "MTOR", RPTOR = "RPTOR",
    RPS6KB1 = "RPS6KB1", MYC = "MYC")
  p <- map(vars(EIF4E, EIF4G1, EIF4G2, EIF4A1,`EIF4E+EIF4EBP1`,
    EIF4EBP1, PABPC1, MKNK1, MKNK2, 
    MTOR, RPTOR, RPS6KB1, MYC), function(by)
    {
      analyse_multivariate(df,
        vars(OS.time, OS),
        covariates = list(by), # covariates expects a list
        covariate_name_dict = covariate_names)
    }) %>%
    forest_plot(
      factor_labeller   = covariate_names,
      endpoint_labeller = c(OS.time = "OS"),
      orderer           = ~order(HR),
      labels_displayed  = c("factor"),
      values_displayed  = c("HR", "CI", "p", "n"),
      value_headers     = c(HR = "HR", CI = "95%CI", p = "p", n = "N"),
      relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
      label_headers     = c(factor = "Gene"),
      title             = "Univariate Cox proportional-hazards regression analysis\nBreast Invasive Carcinoma (N=1101)",
      HR_x_breaks       = seq(0.5, 2.3, 0.2), 
      HR_x_limits       = c(0.5, 2.3),
      ggtheme           = ggplot2::theme_bw(base_size = 7))
  print(p)
  ggsave(
    path        = "~/Documents/EIF_output/Cox", 
    filename    = paste0(tumor, "EIFuniCox.pdf"),
    plot        = p,
    width       = 5, 
    height      = 4, 
    useDingbats = FALSE)
  
  df %>% 
    analyse_multivariate(vars(OS.time, OS),
      covariates = vars(EIF4G1, EIF4G2, EIF4A1, EIF4E, EIF4EBP1,PABPC1, 
        MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC),
      covariate_name_dict = covariate_names) -> result1
  p2 <- forest_plot(result1,
    factor_labeller   = covariate_names, #label the subgroups.
    endpoint_labeller = c(OS.time="OS"), #label the endpoint.
    orderer           = ~order(HR), #order by hazard ratio
    labels_displayed  = c("factor"),
    values_displayed  = c("HR", "CI", "p"),
    value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
    ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
    relative_widths   = c(0.4, 1.8, 1.2), #more space for the plot, less space for the tables
    label_headers     = c(factor = "Gene", n = "Sample size"),
    title             = "Multivariate Cox proportional-hazards regression analysis\nBreast Invasive Carcinoma (N=1101)",
    HR_x_breaks       = seq(0.4, 2.9, 0.3), 
    HR_x_limits       = c(0.4, 2.9)) #more breaks on the X axis
  print(p2)
  ggsave(
    path        = "~/Documents/EIF_output/Cox", 
    filename    = paste0(tumor, "EIFmultiCox.pdf"), 
    plot        = p2,
    width       = 5, 
    height      = 4, 
    useDingbats = FALSE)
  
  df %>% 
    analyse_multivariate(vars(OS.time, OS),
      covariates = vars(EIF4G1, EIF4G2, EIF4A1,`EIF4E+EIF4EBP1`, PABPC1, 
        MKNK1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC),
      covariate_name_dict = covariate_names) -> result2
  p3 <- forest_plot(result2,
    factor_labeller   = covariate_names, #label the subgroups.
    endpoint_labeller = c(OS.time="OS"), #label the endpoint.
    orderer           = ~order(HR), #order by hazard ratio
    labels_displayed  = c("factor"),
    values_displayed  = c("HR", "CI", "p"),
    value_headers     = c(HR = "HR", CI = "95%CI", p = "p"),
    ggtheme           = ggplot2::theme_bw(base_size = 7), #Adjust font size
    relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
    label_headers     = c(factor = "Gene", n = "Sample size"),
    title             = "Multivariate Cox proportional-hazards regression analysis\nBreast Invasive Carcinoma (N=1101)",
    HR_x_breaks       = seq(0.5, 2.7, 0.3), 
    HR_x_limits       = c(0.5, 2.7)) #more breaks on the X axis
  print(p3)
  ggsave(
    path        = "~/Documents/EIF_output/Cox", 
    filename    = paste0(tumor, "EIFsummultiCox.pdf"), 
    plot        = p3,
    width       = 5, 
    height      = 4, 
    useDingbats = FALSE)
  
  # correlation plot
  my_data <- df[, c(1:13)]
  cor_5 <- rcorr(as.matrix(my_data))
  M <- cor_5$r
  p_mat <- cor_5$P
  pdf(file.path(
    path        = "~/Documents/EIF_output/Cox", 
    filename    = paste0(tumor, "EIFcormatrix.pdf")), 
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  corrplot(
    M, 
    method      = "color", 
    tl.cex      = 1, tl.srt = 45,
    number.cex  = 1, 
    addgrid.col = "gray",
    addCoef.col = "black", 
    tl.col      = "black",
    #type        = "upper", 
    order       = "FPC", 
    p.mat       = p_mat, 
    sig.level   = 0.05, #insig = "blank" 
  )
  dev.off()
}
plot.coxph.EIF.breast.tumor(c("breast invasive carcinoma"))

#####################################
## Heatmap of correlation analysis ##
#####################################

### use TCGA-TARGET-GTEX dataset ###
plot.heatmap.total <- function() {
  Lung <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
  Lung <- Lung[!duplicated(Lung$sample), ]
  Lung <- na.omit(Lung)
  row.names(Lung) <- Lung$sample
  Lung$sample <- NULL
  Sample.ID <- row.names(Lung)
  subset <- as.data.frame(Lung$`_sample_type`)
  row.names(subset) <- row.names(Lung)
  colnames(subset) <- "sample.type"
  tissue.GTEX.TCGA.gene <- function(){
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
                           !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.t,
                                       subset,
                                       by    = "row.names",
                                       all.x = TRUE)
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)] 
    TCGA.GTEX.Lung.sampletype <- na.omit(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  gene.name <- names(TCGA.GTEX.sampletype)
  gene.name <- gene.name [! gene.name %in% c("Row.names", "sample.type")]
  
  ### TO BE CONTINUED
  EIF.correlation <- function(y,z){
    TCGA.GTEX.tumor.lung <- TCGA.GTEX.sampletype[
      TCGA.GTEX.sampletype$sample.type %in% y, ]
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.tumor.lung[[x]],
                         TCGA.GTEX.tumor.lung[[y]],
                         method = "pearson")
      res <- data.frame(x,
                        y,
                        result[c("estimate",
                                 "p.value",
                                 "statistic",
                                 "method")],
                        stringsAsFactors=FALSE)
      }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
                          lapply(gene.name,
                                 correlation.coefficient,
                                 y = x))
      rownames(cor.data) <- cor.data[,1]
      #cor.data1 <- cor.data[cor.data[, "p.value"] <= 0.05,]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    plot.pos.Venn <- function(){
      c3 <- cbind(EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
                  EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
                  EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
        B       = a[3, "Counts"],
        C       = a[2, "Counts"],
        "A&B"   = a[7, "Counts"],
        "A&C"   = a[6, "Counts"],
        "B&C"   = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
        #key = TRUE,
        Title = "Tumors",
        lwd = 0,
        fill = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E posCOR",
                                 "EIF4G1 posCOR",
                                 "EIF4A1 posCOR"),
                                 cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all",z,"posVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,,
                  EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,,
                  EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,,
                  EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05,)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts")
      vennDiagram(b)
      pos.Venn2 <- euler(c(
        A       = b[9, "Counts"], #EIF4E
        B       = b[5, "Counts"], #EIF4G1
        C       = b[3, "Counts"], #EIF4A1
        D       = b[2, "Counts"], #EIF4EBP1
        "A&B"   = b[13, "Counts"],
        "A&C"   = b[11, "Counts"],
        "A&D"   = b[10, "Counts"],
        "B&C"   = b[7, "Counts"],
        "B&D"   = b[6, "Counts"],
        "A&B&C" = b[15, "Counts"],
        "A&B&D" = b[14, "Counts"],
        "A&C&D" = b[12, "Counts"],
        "B&C&D" = b[8, "Counts"]))
      p2 <- plot(pos.Venn2,
        #key = TRUE,
        Title = "Tumors",
        lwd = 0,
        fill = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E posCOR",
          "EIF4G1 posCOR",
          "EIF4A1 posCOR",
          "EIF4EBP1 posCOR"),
          cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all", z, "pos4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.Venn()
    plot.neg.Venn <- function(){
      c3 <- cbind(
        EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
        B       = a[3, "Counts"],
        C       = a[2, "Counts"],
        "A&B"   = a[7, "Counts"],
        "A&C"   = a[6, "Counts"],
        "B&C"   = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
        #key = TRUE,
        Title = "Tumors",
        lwd = 0,
        fill = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E negCOR",
          "EIF4G1 negCOR",
          "EIF4A1 negCOR"),
          cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all",z,"negVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(
        EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,,
        EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,,
        EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,,
        EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05,)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts")
      vennDiagram(b)
      pos.Venn2 <- euler(c(
        A       = b[9, "Counts"], #EIF4E
        B       = b[5, "Counts"], #EIF4G1
        C       = b[3, "Counts"], #EIF4A1
        D       = b[2, "Counts"], #EIF4EBP1
        "A&B"   = b[13, "Counts"],
        "A&C"   = b[11, "Counts"],
        "A&D"   = b[10, "Counts"],
        "B&C"   = b[7, "Counts"],
        "B&D"   = b[6, "Counts"],
        "A&B&C" = b[15, "Counts"],
        "A&B&D" = b[14, "Counts"],
        "A&C&D" = b[12, "Counts"],
        "B&C&D" = b[8, "Counts"]))
      p2 <- plot(pos.Venn2,
        #key = TRUE,
        Title = "Tumors",
        lwd = 0,
        fill = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E negCOR",
                                 "EIF4G1 negCOR",
                                 "EIF4A1 negCOR",
                                 "EIF4EBP1 negCOR"),
          cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all", z, "neg4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')),
                      setNames(data.frame(EIF4EBP1.cor[3]), c('EIF4EBP1')))
    return(cor.data)
  }
  all.sample.type <- levels(subset$sample.type)
  all.tumor.type <- all.sample.type [! all.sample.type %in% c("Cell Line", 
                                                              "Normal Tissue", 
                                                              "Solid Tissue Normal")]
  EIF.cor.tumor <- EIF.correlation(y = all.tumor.type, z = "tumor")
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"), z = "normal")
  cor.data <- cbind(setNames(data.frame(EIF.cor.tumor[1:4]),
                    c('EIF4E.tumor',
                      'EIF4G1.tumor',
                      'EIF4A1.tumor', 
                      'EIF4EBP1.tumor')),
                    setNames(data.frame(EIF.cor.normal[1:4]),
                    c('EIF4E.normal',
                      'EIF4G1.normal',
                      'EIF4A1.normal',
                      'EIF4EBP1.normal')))
  DF <- as.matrix(na.omit(cor.data[cor.data$EIF4E.tumor > 0.3 |
                                   cor.data$EIF4E.tumor < -0.3 |
                                   cor.data$EIF4G1.tumor > 0.3 |
                                   cor.data$EIF4G1.tumor < -0.3 |
                                   cor.data$EIF4A1.tumor > 0.3 |
                                   cor.data$EIF4EBP1.tumor < -0.3 |
                                   cor.data$EIF4E.normal > 0.3 |
                                   cor.data$EIF4E.normal < -0.3 |
                                   cor.data$EIF4G1.normal > 0.3 |
                                   cor.data$EIF4G1.normal < -0.3 |
                                   cor.data$EIF4A1.normal > 0.3 |
                                   cor.data$EIF4EBP1.normal < -0.3
    , ]))
  pheatmap::pheatmap(DF,
    # main = "Correlation Coefficient Heatmap",
    # annotation_row = my_sample_col[,"_sample_type",drop=FALSE],
    angle_col     = c("0"),
    fontsize      = 12,
    fontface      = "bold",
    color         = colorRampPalette(rev(brewer.pal(n    = 7, 
                                                    name = "RdYlBu")))(100),
    show_rownames = FALSE,
    show_colnames = FALSE)
  # DF_scaled = t(scale(t(DF)))
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options
  ht1 = Heatmap(DF,
    name     = "Correlation Coefficient Heatmap (all)",
    heatmap_legend_param = list(labels_gp = gpar(font = 15), 
                                legend_width = unit(6, "cm"), 
                                direction = "horizontal"),
    # clustering_distance_rows = function(x, y) 1 - cor(x, y),
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    bottom_annotation    = HeatmapAnnotation(
      annotation_legend_param = list(direction = "horizontal"),
      type      = c("tumor", "tumor","tumor","tumor",
                    "normal","normal","normal","normal"),
      col       = list(type = c("tumor" = "pink", 
                                "normal"  = "royalblue")),
      cn        = anno_text(gsub("\\..*","",colnames(DF)),
                           location = 0,
                           rot      = 0,
                           just     = "center",
                           gp       = gpar(fontsize = 15,
                                           fontface = "bold"))),
     #cluster_rows  = as.dendrogram(hclust(dist(DF))),
     #row_split     = 3,
     row_km     = 3,
     #row_km_repeats = 100,
     row_title    = "cluster_%s", 
     row_title_gp = gpar(fontsize = 15,
                         fontface = "bold"),
     border       = TRUE,
     col          = circlize::colorRamp2(c(-1, 0, 1),
                                         c("blue", "#EEEEEE", "red")))
  ht = draw(ht1,
            merge_legends          = TRUE,
            heatmap_legend_side    = "top", 
            annotation_legend_side = "top")

  pdf(file.path(
    path        = "~/Documents/EIF_output/Heatmap", 
    filename    = "All tumors heatmap.pdf"), 
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  ht = draw(ht1,
            merge_legends       = TRUE,
            heatmap_legend_side = "top")
  dev.off()
  
  ## try to extract clusters from heatmap
  # Saving row names of cluster one  
  plot.cluster.pathway <- function() {
    cluster.geneID.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]],])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez = mapIds(org.Hs.eg.db,
                         keys      = c1$V1,
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")
      # c1 <- c1[!is.na(c1)]
      c1 <- na.omit(c1)
      return(c1$entrez)
    }
    cluster.num <- as.character(c(1:3))
    names(cluster.num) <- paste("cluster", 1:3)
    cluster.data <- lapply(cluster.num, cluster.geneID.list)
    ck.GO <- compareCluster(geneCluster = cluster.data,
                            fun         = "enrichGO",
                            OrgDb       = 'org.Hs.eg.db')
    ck.KEGG <- compareCluster(geneCluster = cluster.data,
                              fun         = "enrichKEGG")
    ck.REACTOME <- compareCluster(geneCluster = cluster.data,
                                  fun         = "enrichPathway")
    p1 <- dotplot(ck.GO,
                  title        = "The Most Enriched GO Pathways",
                  showCategory = 8,
                  font.size    = 18,
                  includeAll   = FALSE) +
          theme_bw() +
          theme(plot.title   = black_bold_tahoma_16,
                axis.title   = black_bold_tahoma_16,
                axis.text.x  = black_bold_tahoma_16,
                axis.text.y  = black_bold_tahoma_16,
                axis.line.x  = element_line(color = "black"),
                axis.line.y  = element_line(color = "black"),
                panel.grid   = element_blank(),
                legend.title = black_bold_tahoma_16,
                legend.text  = black_bold_tahoma_16,
                strip.text   = black_bold_tahoma_16)
    print(p1)    
    ggsave(
          path        = "~/Documents/EIF_output/Heatmap", 
          filename    = paste("all tumors GO.pdf"), 
          plot        = p1,
          width       = 10, 
          height      = 8, 
          useDingbats = FALSE)
    
    p2 <- dotplot(ck.KEGG,
      title        = "The Most Enriched KEGG Pathways",
      showCategory = 8,
      font.size    = 18,
      includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title   = black_bold_tahoma_16,
              axis.title   = black_bold_tahoma_16,
              axis.text.x  = black_bold_tahoma_16,
              axis.text.y  = black_bold_tahoma_16,
              axis.line.x  = element_line(color = "black"),
              axis.line.y  = element_line(color = "black"),
              panel.grid   = element_blank(),
              legend.title = black_bold_tahoma_16,
              legend.text  = black_bold_tahoma_16,
              strip.text   = black_bold_tahoma_16)
    print(p2)    
    ggsave(
      path        = "~/Documents/EIF_output/Heatmap", 
      filename    = paste("all tumors KEGG.pdf"), 
      plot        = p2,
      width       = 12, 
      height      = 8, 
      useDingbats = FALSE)
    
    p3 <- dotplot(ck.REACTOME,
      title        = "The Most Enriched REACTOME Pathways",
      showCategory = 8,
      font.size    = 16,
      includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title   = black_bold_tahoma_16,
              axis.title   = black_bold_tahoma_16,
              axis.text.x  = black_bold_tahoma_16,
              axis.text.y  = black_bold_tahoma_16,
              axis.line.x  = element_line(color = "black"),
              axis.line.y  = element_line(color = "black"),
              panel.grid   = element_blank(),
              legend.title = black_bold_tahoma_16,
              legend.text  = black_bold_tahoma_16,
              strip.text   = black_bold_tahoma_16)
    print(p3)
    ggsave(
      path        = "~/Documents/EIF_output/Heatmap", 
      filename    = paste("all tumors REACTOME.pdf"), 
      plot        = p3,
      width       = 13.5, 
      height      = 8, 
      useDingbats = FALSE)
  }
  plot.cluster.pathway()
}
plot.heatmap.total()

  ### use TCGA-TARGET-GTEX dataset
plot.heatmap.GTEX <- function() {
  get.GTEX.annotation <- function(){
    Sample.type.annotation <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
    Sample.type.annotation <- Sample.type.annotation[
      !duplicated(Sample.type.annotation$sample), ]
    Sample.type.annotation <- as.data.frame(Sample.type.annotation)
    Sample.type.annotation$`_study` <- as.factor(Sample.type.annotation$`_study`)
  # levels(Sample.type.annotation$`_study`), "GTEX"   "TARGET" "TCGA" 
    GTEX.annotation <- Sample.type.annotation[Sample.type.annotation$`_study` == "GTEX",]
    GTEX.annotation$`_study` <- droplevels(GTEX.annotation$`_study`)
    GTEX.annotation <- na.omit(GTEX.annotation)
    row.names(GTEX.annotation) <- GTEX.annotation$sample
    GTEX.annotation$sample <- NULL
    subset <- as.data.frame(GTEX.annotation$`_sample_type`)
    row.names(subset) <- row.names(GTEX.annotation)
    colnames(subset) <- "sample.type"
    return(subset)
    }
  subset <- get.GTEX.annotation()
  Sample.ID <- row.names(subset)
  tissue.GTEX.TCGA.gene <- function(){
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
                           !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE)
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)] 
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  gene.name <- names(TCGA.GTEX.sampletype)
  gene.name <- gene.name [! gene.name %in% c("Row.names", "sample.type")]
  
  ### TO BE CONTINUED
  EIF.correlation <- function(y){
    GTEX.tissue <- TCGA.GTEX.sampletype[
      TCGA.GTEX.sampletype$sample.type %in% y, ]
    correlation.coefficient <- function(x, y) {
      result <- cor.test(GTEX.tissue[[x]],
                         GTEX.tissue[[y]],
                         method = "pearson")
      res <- data.frame(x,
                        y,
                        result[c("estimate",
                                 "p.value",
                                 "statistic",
                                 "method")],
                        stringsAsFactors=FALSE)
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
                          lapply(gene.name,
                                 correlation.coefficient,
                                 y = x))
      rownames(cor.data) <- cor.data[,1]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    plot.pos.Venn <- function(){
      c3 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
                          B       = a[3, "Counts"],
                          C       = a[2, "Counts"],
                          "A&B"   = a[7, "Counts"],
                          "A&C"   = a[6, "Counts"],
                          "B&C"   = a[4, "Counts"],
                          "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
                 #key = TRUE,
                 Title = "Normal Tissues",
                 lwd = 0,
                 fill = c("#999999", "#E69F00", "#56B4E9"),
                 quantities = list(cex = 1.25),
                 labels = list(labels = c("EIF4E posCOR",
                                          "EIF4G1 posCOR",
                                          "EIF4A1 posCOR"),
                               cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all healthy Venn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      c4 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3,
                  EIF4EBP1.cor$estimate > 0.3)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "EIF4EBP1",
                       "Counts")
      vennDiagram(b)
      }
    plot.pos.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')),
                      setNames(data.frame(EIF4EBP1.cor[3]), c('EIF4EBP1')))
    return(cor.data)
  }
  # all.sample.type <- levels(subset$sample.type)
  EIF.cor.normal <- EIF.correlation(y = "Normal Tissue")
  DF <- as.matrix(na.omit(EIF.cor.normal[EIF.cor.normal$EIF4E > 0.3 |
                                         EIF.cor.normal$EIF4E < -0.3 |
                                         EIF.cor.normal$EIF4G1 > 0.3 |
                                         EIF.cor.normal$EIF4G1 < -0.3 |
                                         EIF.cor.normal$EIF4EBP1 > 0.3 |
                                         EIF.cor.normal$EIF4EBP1 < -0.3 |
                                         EIF.cor.normal$EIF4A1 > 0.3 |
                                         EIF.cor.normal$EIF4A1 < -0.3, ]))
  pheatmap::pheatmap(DF,
    main = "Correlation Coefficient Heatmap GTEx",
    angle_col     = c("0"),
    fontsize      = 12,
    fontface      = "bold",
    color         = colorRampPalette(rev(brewer.pal(n    = 7, 
                                                    name = "RdYlBu")))(100),
    show_rownames = FALSE,
    show_colnames = TRUE)
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options
  ht1 = Heatmap(DF,
    name                 = "Correlation Coefficient Heatmap in GTEx",
    heatmap_legend_param = list(direction = "horizontal"),
    # clustering_distance_rows = function(x, y) 1 - cor(x, y),
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    bottom_annotation    = HeatmapAnnotation(
      annotation_legend_param = list(direction = "horizontal"),
      cn       = anno_text(gsub("\\..*","",colnames(DF)),
      location = 0,
      rot      = 0,
      just     = "right",
      gp       = gpar(fontsize = 14,
                      fontface = "bold"))),
    row_km         = 4,
    row_km_repeats = 100,
    row_title      = "cluster_%s",
    row_title_gp   = gpar(fontsize = 14,
                          fontface = "bold"),
    border         = TRUE,
    col            = circlize::colorRamp2(c(-1, 0, 1),
                                          c("blue", "#EEEEEE", "red")))
  ht = draw(ht1,
            merge_legends       = TRUE,
            heatmap_legend_side = "top")
  ## try to extract clusters from heatmap
  # Saving row names of cluster one  
  plot.cluster.pathway <- function() {
    cluster.geneID.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]],])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez = mapIds(org.Hs.eg.db,
        keys      = c1$V1,
        column    = "ENTREZID",
        keytype   = "SYMBOL",
        multiVals = "first")
      # c1 <- c1[!is.na(c1)]
      c1 <- na.omit(c1)
      return(c1$entrez)
    }
    cluster.num <- as.character(c(1:4))
    names(cluster.num) <- paste("cluster", 1:4)
    cluster.data <- lapply(cluster.num, cluster.geneID.list)
    ck.GO <- compareCluster(geneCluster = cluster.data,
      fun         = "enrichGO",
      OrgDb       = 'org.Hs.eg.db')
    ck.KEGG <- compareCluster(geneCluster = cluster.data,
      fun         = "enrichKEGG")
    ck.REACTOME <- compareCluster(geneCluster = cluster.data,
      fun         = "enrichPathway")
    print(dotplot(ck.GO,
      title        = "The Most Enriched GO Pathways",
      showCategory = 8,
      font.size    = 18,
      includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title   = black_bold_tahoma_16,
          axis.title   = black_bold_tahoma_16,
          axis.text.x  = black_bold_tahoma_16,
          axis.text.y  = black_bold_tahoma_16,
          axis.line.x  = element_line(color = "black"),
          axis.line.y  = element_line(color = "black"),
          panel.grid   = element_blank(),
          legend.title = black_bold_tahoma_16,
          legend.text  = black_bold_tahoma_16,
          strip.text   = black_bold_tahoma_16))
    print(dotplot(ck.KEGG,
      title        = "The Most Enriched KEGG Pathways",
      showCategory = 8,
      font.size    = 18,
      includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title   = black_bold_tahoma_16,
          axis.title   = black_bold_tahoma_16,
          axis.text.x  = black_bold_tahoma_16,
          axis.text.y  = black_bold_tahoma_16,
          axis.line.x  = element_line(color = "black"),
          axis.line.y  = element_line(color = "black"),
          panel.grid   = element_blank(),
          legend.title = black_bold_tahoma_16,
          legend.text  = black_bold_tahoma_16,
          strip.text   = black_bold_tahoma_16))
    print(dotplot(ck.REACTOME,
      title        = "The Most Enriched REACTOME Pathways",
      showCategory = 8,
      font.size    = 16,
      includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title   = black_bold_tahoma_16,
          axis.title   = black_bold_tahoma_16,
          axis.text.x  = black_bold_tahoma_16,
          axis.text.y  = black_bold_tahoma_16,
          axis.line.x  = element_line(color = "black"),
          axis.line.y  = element_line(color = "black"),
          panel.grid   = element_blank(),
          legend.title = black_bold_tahoma_16,
          legend.text  = black_bold_tahoma_16,
          strip.text   = black_bold_tahoma_16))
  }
  plot.cluster.pathway()
}
plot.heatmap.GTEX()

### use TCGA-TARGET-GTEX dataset
plot.heatmap.TCGA <- function() {
  get.TCGA.annotation <- function(){
    Sample.type.annotation <- read_tsv(
      "~/Downloads/TcgaTargetGTEX_phenotype.txt")
    Sample.type.annotation <- Sample.type.annotation[
      !duplicated(Sample.type.annotation$sample), ]
    Sample.type.annotation <- as.data.frame(Sample.type.annotation)
    Sample.type.annotation$`_study` <- as.factor(Sample.type.annotation$`_study`)
    # levels(Sample.type.annotation$`_study`), "GTEX"   "TARGET" "TCGA" 
    GTEX.annotation <- Sample.type.annotation[
      Sample.type.annotation$`_study` == "TCGA",]
    GTEX.annotation$`_study` <- droplevels(GTEX.annotation$`_study`)
    GTEX.annotation <- na.omit(GTEX.annotation)
    row.names(GTEX.annotation) <- GTEX.annotation$sample
    GTEX.annotation$sample <- NULL
    subset <- as.data.frame(GTEX.annotation$`_sample_type`)
    row.names(subset) <- row.names(GTEX.annotation)
    colnames(subset) <- "sample.type"
    return(subset)
  }
  subset <- get.TCGA.annotation()
  Sample.ID <- row.names(subset)
  tissue.GTEX.TCGA.gene <- function(){
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
                           !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.sampletype <- merge(TCGA.GTEX.t,
                                  subset,
                                  by    = "row.names",
                                  all.x = TRUE)
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)] 
    TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
    return(TCGA.GTEX.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  gene.name <- names(TCGA.GTEX.sampletype)
  gene.name <- gene.name [! gene.name %in% c("Row.names", "sample.type")]
  
  ### TO BE CONTINUED
  EIF.correlation <- function(y){
    GTEX.tissue <- TCGA.GTEX.sampletype[
      !TCGA.GTEX.sampletype$sample.type %in% y, ]
    correlation.coefficient <- function(x, y) {
      result <- cor.test(GTEX.tissue[[x]],
                         GTEX.tissue[[y]],
                         method = "pearson")
      res <- data.frame(x,
                        y,
                        result[c("estimate",
                                 "p.value",
                                 "statistic",
                                 "method")],
                        stringsAsFactors=FALSE)
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
                          lapply(gene.name,
                                 correlation.coefficient,
                                 y = x))
      rownames(cor.data) <- cor.data[,1]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    
    plot.pos.Venn <- function(){
      c3 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
                          B       = a[3, "Counts"],
                          C       = a[2, "Counts"],
                          "A&B"   = a[7, "Counts"],
                          "A&C"   = a[6, "Counts"],
                          "B&C"   = a[4, "Counts"],
                          "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
                 #key = TRUE,
                 Title = "Tumors",
                 lwd = 0,
                 fill = c("#999999", "#E69F00", "#56B4E9"),
                 quantities = list(cex = 1.25),
                 labels = list(labels = c("EIF4E posCOR",
                                          "EIF4G1 posCOR",
                                          "EIF4A1 posCOR"),
                               cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all tumor Venn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3,
                  EIF4EBP1.cor$estimate > 0.3)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "EIF4EBP1",
                       "Counts")
      vennDiagram(b)
      pos.Venn2 <- euler(c(A       = b[9, "Counts"], #EIF4E
                           B       = b[5, "Counts"], #EIF4G1
                           C       = b[3, "Counts"], #EIF4A1
                           D       = b[2, "Counts"], #EIF4EBP1
                           "A&B"   = b[13, "Counts"],
                           "A&C"   = b[11, "Counts"],
                           "A&D"   = b[10, "Counts"],
                           "B&C"   = b[7, "Counts"],
                           "B&D"   = b[6, "Counts"],
                           "A&B&C" = b[15, "Counts"],
                           "A&B&D" = b[14, "Counts"],
                           "A&C&D" = b[12, "Counts"],
                           "B&C&D" = b[8, "Counts"]))
      p2 <- plot(pos.Venn2,
        #key = TRUE,
        Title = "Tumors",
        lwd = 0,
        fill = c("#999999", "#E69F00", "#56B4E9","#009E73"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E posCOR",
          "EIF4G1 posCOR",
          "EIF4A1 posCOR",
          "EIF4EBP1 posCOR"),
          cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all tumor 4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      }
    plot.pos.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[, c(3,4)]), c('EIF4E','EIF4E.pvalue')),
                      setNames(data.frame(EIF4G1.cor[, c(3,4)]), c('EIF4G1','EIF4G1.pvalue')),
                      setNames(data.frame(EIF4A1.cor[, c(3,4)]), c('EIF4A1','EIF4A1.pvalue')),
                      setNames(data.frame(EIF4EBP1.cor[, c(3,4)]), c('EIF4EBP1','EIF4EBP1.pvalue')))
    return(cor.data)
  }
  # all.sample.type <- levels(subset$sample.type)
  EIF.cor.tumor <- EIF.correlation(y = "Solid Tissue Normal")
  DF <- as.matrix(na.omit(EIF.cor.tumor[
    EIF.cor.tumor$EIF4E > 0.3 & EIF.cor.tumor$EIF4E.pvalue <= 0.05 |
    EIF.cor.tumor$EIF4E < -0.3 & EIF.cor.tumor$EIF4E.pvalue <= 0.05 |
    EIF.cor.tumor$EIF4G1 > 0.3 & EIF.cor.tumor$EIF4G1.pvalue <= 0.05 |
    EIF.cor.tumor$EIF4G1 < -0.3 & EIF.cor.tumor$EIF4G1.pvalue <= 0.05 |
    EIF.cor.tumor$EIF4A1 > 0.3 & EIF.cor.tumor$EIF4A1.pvalue <= 0.05 |
    EIF.cor.tumor$EIF4A1 < -0.3 & EIF.cor.tumor$EIF4A1.pvalue <= 0.05 |
    EIF.cor.tumor$EIF4EBP1 > 0.3 & EIF.cor.tumor$EIF4EBP1.pvalue <= 0.05 |
    EIF.cor.tumor$EIF4EBP1 < -0.3 & EIF.cor.tumor$EIF4EBP1.pvalue <= 0.05, ]))
  pheatmap::pheatmap(DF[, c(1,3,5,7)],
                     main = "Correlation Coefficient Heatmap TCGA",
                     angle_col     = c("0"),
                     fontsize      = 12,
                     fontface      = "bold",
                     color         = colorRampPalette(rev(brewer.pal(n    = 7, 
                                                                     name = "RdYlBu")))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE)
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options
  ht1 = Heatmap(DF[, c(1,3,5,7)],
                name                 = "Correlation Coefficient Heatmap in TCGA",
                heatmap_legend_param = list(direction = "horizontal"),
                # clustering_distance_rows = function(x, y) 1 - cor(x, y),
                show_row_names       = FALSE,
                show_column_names    = FALSE,
                bottom_annotation    = HeatmapAnnotation(
                  annotation_legend_param = list(direction = "horizontal"),
                  cn       = anno_text(gsub("\\..*","",colnames(DF[, c(1,3,5,7)])),
                                       location = 0,
                                       rot      = 0,
                                       just     = "center",
                                       gp       = gpar(fontsize = 14,
                                                       fontface = "bold"))),
                row_km         = 4,
                row_km_repeats = 100,
                row_title      = "cluster_%s",
                row_title_gp   = gpar(fontsize = 14,
                                      fontface = "bold"),
                border         = TRUE,
                col            = circlize::colorRamp2(c(-1, 0, 1),
                                                      c("blue", "#EEEEEE", "red")))
  ht = draw(ht1,
            merge_legends       = TRUE,
            heatmap_legend_side = "top")
  ## try to extract clusters from heatmap
  # Saving row names of cluster one  
  plot.cluster.pathway <- function() {
    cluster.geneID.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]],])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez = mapIds(org.Hs.eg.db,
                         keys      = c1$V1,
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")
      # c1 <- c1[!is.na(c1)]
      c1 <- na.omit(c1)
      return(c1$entrez)
    }
    cluster.num <- as.character(c(1:4))
    names(cluster.num) <- paste("cluster", 1:4)
    cluster.data <- lapply(cluster.num, cluster.geneID.list)
    ck.GO <- compareCluster(geneCluster = cluster.data,
                            fun         = "enrichGO",
                            OrgDb       = 'org.Hs.eg.db')
    ck.KEGG <- compareCluster(geneCluster = cluster.data,
                              fun         = "enrichKEGG")
    ck.REACTOME <- compareCluster(geneCluster = cluster.data,
                                  fun         = "enrichPathway")
    print(dotplot(ck.GO,
                  title        = "The Most Enriched GO Pathways",
                  showCategory = 8,
                  font.size    = 18,
                  includeAll   = FALSE) +
            theme_bw() +
            theme(plot.title   = black_bold_tahoma_16,
                  axis.title   = black_bold_tahoma_16,
                  axis.text.x  = black_bold_tahoma_16,
                  axis.text.y  = black_bold_tahoma_16,
                  axis.line.x  = element_line(color = "black"),
                  axis.line.y  = element_line(color = "black"),
                  panel.grid   = element_blank(),
                  legend.title = black_bold_tahoma_16,
                  legend.text  = black_bold_tahoma_16,
                  strip.text   = black_bold_tahoma_16))
    print(dotplot(ck.KEGG,
                  title        = "The Most Enriched KEGG Pathways",
                  showCategory = 8,
                  font.size    = 18,
                  includeAll   = FALSE) +
            theme_bw() +
            theme(plot.title   = black_bold_tahoma_16,
                  axis.title   = black_bold_tahoma_16,
                  axis.text.x  = black_bold_tahoma_16,
                  axis.text.y  = black_bold_tahoma_16,
                  axis.line.x  = element_line(color = "black"),
                  axis.line.y  = element_line(color = "black"),
                  panel.grid   = element_blank(),
                  legend.title = black_bold_tahoma_16,
                  legend.text  = black_bold_tahoma_16,
                  strip.text   = black_bold_tahoma_16))
    print(dotplot(ck.REACTOME,
                  title        = "The Most Enriched REACTOME Pathways",
                  showCategory = 8,
                  font.size    = 16,
                  includeAll   = FALSE) +
            theme_bw() +
            theme(plot.title   = black_bold_tahoma_16,
                  axis.title   = black_bold_tahoma_16,
                  axis.text.x  = black_bold_tahoma_16,
                  axis.text.y  = black_bold_tahoma_16,
                  axis.line.x  = element_line(color = "black"),
                  axis.line.y  = element_line(color = "black"),
                  panel.grid   = element_blank(),
                  legend.title = black_bold_tahoma_16,
                  legend.text  = black_bold_tahoma_16,
                  strip.text   = black_bold_tahoma_16))
  }
  plot.cluster.pathway()
}
plot.heatmap.TCGA()

    ### use TCGA pan-cancer dataset
plot.heatmap.all.TCGA <- function () {
  pan.TCGA.gene <- function(){
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.pancancer <- fread(
      "~/Downloads/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena", 
      data.table = FALSE)
    # download https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz
    TCGA.sampletype <- read_tsv(
      "~/Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$sample),
                                      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    TCGA.sampletype,
                                    by    = "row.names",
                                    all.x = TRUE)
    TCGA.RNAseq.sampletype <- as.data.frame(TCGA.RNAseq.sampletype)
    return(TCGA.RNAseq.sampletype)
  }
  TCGA.sampletype.all <- pan.TCGA.gene()
  sample.type.list <- levels(as.factor(TCGA.sampletype.all$sample_type))
  EIF.cor.tumor <- function (){
    # Tumors <- sample.type.list[! sample.type.list %in% "Solid Tissue Normal"]
    TCGA.RNAseq.sampletype <- TCGA.sampletype.all[
      !TCGA.sampletype.all$sample_type %in% "Solid Tissue Normal", ]
    EIF.correlation <- function (x, y) {
      result <- cor.test(TCGA.RNAseq.sampletype[[x]],
                         TCGA.RNAseq.sampletype[[y]],
                         method = "pearson")
      res <- data.frame(x,
                        y,
                        result[c("estimate",
                                 "p.value",
                                 "statistic",
                                 "method")],
                        stringsAsFactors=FALSE)
      }
  # find all genes positively correlate with EIF4F expression
  # lapply function gives a large list, need to convert it to a dataframe
    gene.name <- names(TCGA.RNAseq.sampletype)
    gene.name <- gene.name [! gene.name %in% c("Row.names",
                                               "sample_type",
                                               "_primary_disease")]
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
                          lapply(gene.name,
                                 EIF.correlation,
                                 y = x))
      rownames(cor.data) <- cor.data[,1]
        return(cor.data)
      }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    plot.pos.Venn <- function(){
      c3 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
                          B       = a[3, "Counts"],
                          C       = a[2, "Counts"],
                          "A&B"   = a[7, "Counts"],
                          "A&C"   = a[6, "Counts"],
                          "B&C"   = a[4, "Counts"],
                          "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
                 #key = TRUE,
                 lwd        = 0,
                 fill       = c("#999999", "#E69F00", "#56B4E9"),
                 quantities = list(cex = 1.25),
                 main       = "Tumors",
                 labels     = list(labels = c("EIF4E posCOR",
                                              "EIF4G1 posCOR",
                                              "EIF4A1 posCOR"),
                                   cex    = 1.25))
      print(p1)
      c4 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3,
                  EIF4EBP1.cor$estimate > 0.3)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "EIF4EBP1",
                       "Counts")
      vennDiagram(b)
      }
    plot.pos.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')),
                      setNames(data.frame(EIF4EBP1.cor[3]), c('EIF4EBP1')))
    return(cor.data)
    }
  EIF.cor.tumor <- EIF.cor.tumor()
  DF  <- as.matrix(na.omit(EIF.cor.tumor[EIF.cor.tumor$EIF4E  > 0.3 |
                                         EIF.cor.tumor$EIF4G1 > 0.3 |
                                         EIF.cor.tumor$EIF4A1 > 0.3 |
                                         EIF.cor.tumor$EIF4EBP1 > 0.3, ]))
  pheatmap::pheatmap(DF,
                     # main = "Correlation Coefficient Heatmap",
                     angle_col     = c("0"),
                     fontsize      = 12,
                     fontface      = "bold",
                     color         = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                     show_rownames = FALSE,
                     show_colnames = TRUE)
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options)
  ht1 = Heatmap(DF,
               name                 = "Correlation Coefficient Heatmap in TCGA",
               heatmap_legend_param = list(direction = "horizontal"),
               show_row_names       = FALSE,
               show_column_names    = FALSE,
               bottom_annotation    = HeatmapAnnotation(
                cn   = anno_text(colnames(DF),
                                 location = 0,
                                 rot      = 0,
                                 just     = "center",
                                 gp       = gpar(fontsize = 14,
                                                 fontface = "bold"))),
               row_km           = 3,
               row_title        = "cluster_%s",
               row_title_gp     = gpar(fontsize = 14,
                                       fontface = "bold"),
               border           = TRUE,
               col              = colorRamp2(c(-1, 0, 1),
                                             c("blue", "white", "red")))
  ht = draw(ht1, heatmap_legend_side = "top")
  ## try to extract clusters from heatmap
  # Saving row names of cluster one
  plot.cluster.pathway <- function() {
    cluster.gene.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]],])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez = mapIds(org.Hs.eg.db,
                         keys      = c1$V1,
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")
      # c1 <- c1[!is.na(c1)]
      c1 <- na.omit(c1)
      return(c1$entrez)
      }
    cluster.num <- as.character(c(1:3))
    names(cluster.num) <- paste("cluster", 1:3)
    cluster.data <- lapply(cluster.num, cluster.gene.list)
    ck.GO <- compareCluster(geneCluster       = cluster.data,
                            fun               = "enrichGO",
                            OrgDb             = 'org.Hs.eg.db')
    ck.KEGG <- compareCluster(geneCluster     = cluster.data,
                              fun             = "enrichKEGG")
    ck.REACTOME <- compareCluster(geneCluster = cluster.data,
                                  fun         = "enrichPathway")
    print(dotplot(ck.GO,
                  title        = "The Most Enriched GO Pathways",
                  showCategory = 8,
                  font.size    = 18,
                  includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title       = black_bold_tahoma_16,
              axis.title       = black_bold_tahoma_16,
              axis.text.x      = black_bold_tahoma_16,
              axis.text.y      = black_bold_tahoma_16,
              axis.line.x      = element_line(color = "black"),
              axis.line.y      = element_line(color = "black"),
              panel.grid       = element_blank(),
              legend.title     = black_bold_tahoma_16,
              legend.text      = black_bold_tahoma_16,
              strip.text       = black_bold_tahoma_16))
    print(dotplot(ck.KEGG,
                  title        = "The Most Enriched KEGG Pathways",
                  showCategory = 8,
                  font.size    = 18,
                  includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title       = black_bold_tahoma_16,
              axis.title       = black_bold_tahoma_16,
              axis.text.x      = black_bold_tahoma_16,
              axis.text.y      = black_bold_tahoma_16,
              axis.line.x      = element_line(color = "black"),
              axis.line.y      = element_line(color = "black"),
              panel.grid       = element_blank(),
              legend.title     = black_bold_tahoma_16,
              legend.text      = black_bold_tahoma_16,
              strip.text       = black_bold_tahoma_16))
    print(dotplot(ck.REACTOME,
                  title        = "The Most Enriched REACTOME Pathways",
                  showCategory = 8,
                  font.size    = 16,
                  includeAll   = FALSE) +
        theme_bw() +
        theme(plot.title       = black_bold_tahoma_16,
              axis.title       = black_bold_tahoma_16,
              axis.text.x      = black_bold_tahoma_16,
              axis.text.y      = black_bold_tahoma_16,
              axis.line.x      = element_line(color = "black"),
              axis.line.y      = element_line(color = "black"),
              panel.grid       = element_blank(),
              legend.title     = black_bold_tahoma_16,
              legend.text      = black_bold_tahoma_16,
              strip.text       = black_bold_tahoma_16))
  }
  plot.cluster.pathway()
    }
plot.heatmap.all.TCGA()

### use GTEx dataset
plot.heatmap.all.GTEx <- function () {
    # download https://toil.xenahubs.net/download/gtex_RSEM_Hugo_norm_count.gz
    TCGA.pancancer <- fread("~/Downloads/gtex_RSEM_Hugo_norm_count", data.table=FALSE)
    # download https://toil.xenahubs.net/download/GTEX_phenotype.gz
    TCGA.sampletype <- read_tsv("~/Downloads/GTEX_phenotype")
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$sample),
                                      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$sample
    TCGA.pancancer1$sample <- NULL
    # TCGA.pancancer1 <- TCGA.pancancer1[,colnames(TCGA.pancancer1) %in% Sample.ID]
    TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer1)
    rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer1)
    colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer1)
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    TCGA.RNAseq.sampletype <- merge(TCGA.pancancer_transpose,
                                    TCGA.sampletype,
                                    by    = "row.names",
                                    all.x = TRUE)
    TCGA.RNAseq.sampletype <- as.data.frame(TCGA.RNAseq.sampletype)

  EIF.cor.all <- function (){
    EIF.correlation <- function (x, y) {
      result <- cor.test(TCGA.RNAseq.sampletype[[x]],
                         TCGA.RNAseq.sampletype[[y]],
                         method = "pearson")
      res <- data.frame(x,
                        y,
                        result[c("estimate",
                                 "p.value",
                                 "statistic",
                                 "method")],
                        stringsAsFactors=FALSE)
      }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    gene.name <- names(TCGA.RNAseq.sampletype)
    gene.name <- gene.name [! gene.name %in% c("Row.names",
                                               names(TCGA.sampletype))]
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
                          lapply(gene.name,
                                 EIF.correlation,
                                 y = x))
      rownames(cor.data) <- cor.data[,1]
      return(cor.data)
      }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    plot.pos.Venn <- function(){
      c3 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
                          B       = a[3, "Counts"],
                          C       = a[2, "Counts"],
                          "A&B"   = a[7, "Counts"],
                          "A&C"   = a[6, "Counts"],
                          "B&C"   = a[4, "Counts"],
                          "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
        #key = TRUE,
                 lwd = 0,
                 fill = c("#999999", "#E69F00", "#56B4E9"),
                 quantities = list(cex = 1.25),
                 main = "Normal Tissues",
                 labels = list(labels = c("EIF4E posCOR",
                                          "EIF4G1 posCOR",
                                          "EIF4A1 posCOR"),
                               cex    = 1.25))
      print(p1)}
    plot.pos.Venn()
    plot.neg.Venn <- function(){
      c3 <- cbind(EIF4E.cor$estimate < -0.3,
                  EIF4G1.cor$estimate < -0.3,
                  EIF4A1.cor$estimate < -0.3)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
                          B       = a[3, "Counts"],
                          C       = a[2, "Counts"],
                          "A&B"   = a[7, "Counts"],
                          "A&C"   = a[6, "Counts"],
                          "B&C"   = a[4, "Counts"],
                          "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
                 #key = TRUE,
                 lwd = 0,
                 fill = c("#999999", "#E69F00", "#56B4E9"),
                 quantities = list(cex = 1.25),
                 main = "Normal Tissues",
                 labels = list(labels = c("EIF4E negCOR",
                                          "EIF4G1 negCOR",
                                          "EIF4A1 negCOR"),
                               cex    = 1.25))
      print(p1)}
    plot.neg.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')))
    return(cor.data)
  }
  EIF.cor.all <- EIF.cor.all()
  DF  <- as.matrix(na.omit(EIF.cor.all[EIF.cor.all$EIF4E  > 0.3 |
                                       EIF.cor.all$EIF4G1 > 0.3 |
                                       EIF.cor.all$EIF4A1 > 0.3 , ]))
  pheatmap::pheatmap(DF,
    # main = "Correlation Coefficient Heatmap",
    angle_col     = c("0"),
    fontsize      = 12,
    fontface      = "bold",
    color         = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    show_rownames = FALSE,
    show_colnames = TRUE)
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options)
  ht1 = Heatmap(DF,
    name                 = "Correlation Coefficient Heatmap in GTEx",
    heatmap_legend_param = list(direction = "horizontal"),
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    bottom_annotation    = HeatmapAnnotation(
      cn   = anno_text(colnames(DF),
                       location = 0,
                       rot  = 0,
                       just = "center",
                       gp   = gpar(fontsize = 14,
                                   fontface = "bold"))),
    row_km           = 3,
    row_title        = "cluster_%s",
    border           = TRUE,
    col              = colorRamp2(c(-1, 0, 1),
                                  c("blue", "white", "red")))
  ht = draw(ht1, heatmap_legend_side = "top")
  ## try to extract clusters from heatmap
  # Saving row names of cluster one
  plot.cluster.pathway <- function() {
    cluster.gene.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]],])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez = mapIds(org.Hs.eg.db,
                         keys      = c1$V1,
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")
      # c1 <- c1[!is.na(c1)]
      c1 <- na.omit(c1)
      return(c1$entrez)
    }
    cluster.num <- as.character(c(1:3))
    names(cluster.num) <- paste("cluster", 1:3)
    cluster.data <- lapply(cluster.num, cluster.gene.list)
    ck.GO <- compareCluster(geneCluster = cluster.data,
                            fun         = "enrichGO",
                            OrgDb       = 'org.Hs.eg.db')
    ck.KEGG <- compareCluster(geneCluster = cluster.data,
                              fun         = "enrichKEGG")
    ck.REACTOME <- compareCluster(geneCluster = cluster.data,
                                  fun         = "enrichPathway")
    print(dotplot(ck.GO,
                  title         = "The Most Enriched GO Pathways",
                  showCategory  = 8,
                  font.size     = 18,
                  includeAll    = FALSE) +
          theme_bw() +
          theme(plot.title      = black_bold_tahoma_16,
                axis.title      = black_bold_tahoma_16,
                axis.text.x     = black_bold_tahoma_16,
                axis.text.y     = black_bold_tahoma_16,
                axis.line.x     = element_line(color = "black"),
                axis.line.y     = element_line(color = "black"),
                panel.grid      = element_blank(),
                legend.title    = black_bold_tahoma_16,
                legend.text     = black_bold_tahoma_16,
                strip.text      = black_bold_tahoma_16))
    print(dotplot(ck.KEGG,
                  title         = "The Most Enriched KEGG Pathways",
                  showCategory  = 8,
                  font.size     = 18,
                  includeAll    = FALSE) +
          theme_bw() +
          theme(plot.title      = black_bold_tahoma_16,
                axis.title      = black_bold_tahoma_16,
                axis.text.x     = black_bold_tahoma_16,
                axis.text.y     = black_bold_tahoma_16,
                axis.line.x     = element_line(color = "black"),
                axis.line.y     = element_line(color = "black"),
                panel.grid      = element_blank(),
                legend.title    = black_bold_tahoma_16,
                legend.text     = black_bold_tahoma_16,
                strip.text      = black_bold_tahoma_16))
    print(dotplot(ck.REACTOME,
                  title         = "The Most Enriched REACTOME Pathways",
                  showCategory  = 8,
                  font.size     = 16,
                  includeAll    = FALSE) +
          theme_bw() +
          theme(plot.title      = black_bold_tahoma_16,
                axis.title      = black_bold_tahoma_16,
                axis.text.x     = black_bold_tahoma_16,
                axis.text.y     = black_bold_tahoma_16,
                axis.line.x     = element_line(color = "black"),
                axis.line.y     = element_line(color = "black"),
                panel.grid      = element_blank(),
                legend.title    = black_bold_tahoma_16,
                legend.text     = black_bold_tahoma_16,
                strip.text      = black_bold_tahoma_16))
  }
  plot.cluster.pathway()
}
plot.heatmap.all.GTEx()

### use TCGA-TARGET-GTEX dataset
### plot heatmap, Venn diagram and pathway analysis on clusters ###
plot.heatmap.lung <- function(x) {
  # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
  Sampletype <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
  tissue.GTEX.TCGA.gene <- function(x){
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE)
    Lung <- Sampletype[Sampletype$`_primary_site` == x,]
    Lung.ID <- as.vector(Lung$sample)
    Lung.ID <- na.omit(Lung.ID) # NA in the vector
    TCGA.GTEX.Lung <- TCGA.GTEX %>% select("sample", Lung.ID)
    #   TCGA.GTEX.Lung1 <- as.data.frame(TCGA.GTEX.Lung)
    TCGA.GTEX.Lung <- TCGA.GTEX.Lung[!duplicated(TCGA.GTEX.Lung$sample), 
                                     !duplicated(colnames(TCGA.GTEX.Lung))]
    row.names(TCGA.GTEX.Lung) <- TCGA.GTEX.Lung$sample
    TCGA.GTEX.Lung$sample <- NULL
    TCGA.GTEX.Lung.t <- data.table::transpose(TCGA.GTEX.Lung)
    rownames(TCGA.GTEX.Lung.t) <- colnames(TCGA.GTEX.Lung)
    colnames(TCGA.GTEX.Lung.t) <- rownames(TCGA.GTEX.Lung)
    
    Lung <- Lung[!duplicated(Lung$sample), ]
    Lung <- na.omit(Lung)
    row.names(Lung) <- Lung$sample
    Lung$sample <- NULL
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.Lung.t,
                                       Lung,
                                       by    = "row.names",
                                       all.x = TRUE)
    TCGA.GTEX.Lung.sampletype <- as.data.frame(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype.lung <- tissue.GTEX.TCGA.gene(x)
  row.names(TCGA.GTEX.sampletype.lung) <- TCGA.GTEX.sampletype.lung$Row.names
  TCGA.GTEX.sampletype.lung$Row.names <- NULL
  TCGA.GTEX.sampletype.lung$`_sample_type` <- as.factor(
    TCGA.GTEX.sampletype.lung$`_sample_type`)
  # remove Solid Tissue Normal data
  #TCGA.GTEX.sampletype.lung <- TCGA.GTEX.sampletype.lung[
  #  !(TCGA.GTEX.sampletype.lung$`_sample_type` %in% "Solid Tissue Normal"), ]
  
  geneID <- colnames(Sampletype[Sampletype$`_primary_site` == x,])
  TCGA.GTEX.lung <- TCGA.GTEX.sampletype.lung[ ,
    !names(TCGA.GTEX.sampletype.lung) %in% geneID]
  gene.name <- names(TCGA.GTEX.lung)
  EIF.correlation <- function(y){
    TCGA.GTEX.subset.lung <- TCGA.GTEX.sampletype.lung[
      TCGA.GTEX.sampletype.lung$`_sample_type` %in% y, ]
    TCGA.GTEX.subset.lung$`_sample_type` <- droplevels(
      TCGA.GTEX.subset.lung$`_sample_type`)
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.subset.lung[[x]],
                         TCGA.GTEX.subset.lung[[y]],
                         method = "pearson")
      res <- data.frame(x,
                        y,
                        result[c("estimate",
                                 "p.value",
                                 "statistic",
                                 "method")],
                        stringsAsFactors=FALSE)
      }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
                          lapply(gene.name,
                                 correlation.coefficient,
                                 y = x))
      rownames(cor.data) <- cor.data[,1]
      return(cor.data)
      }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    
    plot.pos.Venn <- function(y){
      c3 <- cbind(EIF4E.cor$estimate > 0.3,
                  EIF4G1.cor$estimate > 0.3,
                  EIF4A1.cor$estimate > 0.3)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(A       = a[5, "Counts"],
                          B       = a[3, "Counts"],
                          C       = a[2, "Counts"],
                          "A&B"   = a[7, "Counts"],
                          "A&C"   = a[6, "Counts"],
                          "B&C"   = a[4, "Counts"],
                          "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
                #key = TRUE,
                 lwd        = 0,
                 border     = "black",
                 fill       = c("#999999", "#E69F00", "#56B4E9"),
                 quantities = list(cex = 1.5),
                 # main       = y,
                 labels     = list(labels = c("EIF4E posCOR",
                                              "EIF4G1 posCOR",
                                              "EIF4A1 posCOR"),
                 cex        = 1.5))
      print(p1)
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("Lung Venn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.Venn(y)
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')),
                      setNames(data.frame(EIF4EBP1.cor[3]), c('EIF4EBP1')))
    return(cor.data)
    }
  EIF.cor.tumor <- EIF.correlation(y = c("Primary Tumor", 
                                         "Metastatic", 
                                         "Recurrent Tumor"))
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"))
  cor.data <- cbind(setNames(data.frame(EIF.cor.tumor[1:4]),
                             c('EIF4E.tumor',
                               'EIF4G1.tumor',
                               'EIF4A1.tumor',
                               'EIF4EBP1.tumor')),
                    setNames(data.frame(EIF.cor.normal[1:4]),
                             c('EIF4E.normal',
                               'EIF4G1.normal',
                               'EIF4A1.normal',
                               'EIF4EBP1.normal')))
  DF <- as.matrix(na.omit(cor.data[cor.data$EIF4E.tumor > 0.3 |
                                   cor.data$EIF4E.tumor < -0.3 |
                                   cor.data$EIF4G1.tumor > 0.3 |
                                   cor.data$EIF4G1.tumor < -0.3 |
                                   cor.data$EIF4A1.tumor > 0.3 |
                                   cor.data$EIF4EBP1.tumor < -0.3 |
                                   cor.data$EIF4E.normal > 0.3 |
                                   cor.data$EIF4E.normal < -0.3 |
                                   cor.data$EIF4G1.normal > 0.3 |
                                   cor.data$EIF4G1.normal < -0.3 |
                                   cor.data$EIF4A1.normal > 0.3 |
                                   cor.data$EIF4EBP1.normal < -0.3, ]))
  # DF <- as.matrix(na.omit(cor.data))
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options)
  pheatmap::pheatmap(DF,
    # main = "Correlation Coefficient Heatmap",
    angle_col     = c("0"),
    fontsize      = 12,
    fontface      = "bold",
    color         = colorRampPalette(rev(brewer.pal(n    = 7, 
                                                    name = "RdYlBu")))(100),
    show_rownames = FALSE,
    show_colnames = TRUE)
  ht1 = Heatmap(DF,
    name                 = "Correlation Coefficient Heatmap(Lung)",
    heatmap_legend_param = list(direction = "horizontal",
                                legend_width = unit(6, "cm")),
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    bottom_annotation    = HeatmapAnnotation(
      annotation_legend_param = list(direction = "horizontal"),
      type     = c("tumor", "tumor","tumor","tumor",
                   "normal","normal","normal","normal"),
      col      = list(type = c("tumor"  = "pink", 
                               "normal" = "royalblue")),
      cn       = anno_text(gsub("\\..*","",colnames(DF)),
      location = 0,
      rot      = 0,
      just     = "center",
      gp       = gpar(fontsize = 15,
                      fontface = "bold"))),
    row_km         = 3,
    #row_km_repeats = 100,
    row_title      = "cluster_%s",
    row_title_gp   = gpar(fontsize = 15,
                          fontface = "bold"),
    border         = TRUE,
    col            = circlize::colorRamp2(c(-1, 0, 1),
                                          c("blue", "#EEEEEE", "red")))
  ht = draw(ht1,
            merge_legends       = TRUE,
            heatmap_legend_side = "top")
  pdf(file.path(
    path        = "~/Documents/EIF_output/Heatmap", 
    filename    = "Lung tumors heatmap.pdf"), 
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  ht = draw(ht1,
    merge_legends       = TRUE,
    heatmap_legend_side = "top")
  dev.off()
  ## try to extract clusters from heatmap
  # Saving row names of cluster one  
  plot.cluster.pathway <- function() {
    cluster.geneID.list <- function(x) {
      c1 <- t(t(row.names(DF[row_order(ht1)[[x]],])))
      c1 <- as.data.frame(c1)
      c1$V1 <- as.character(c1$V1)
      c1$entrez = mapIds(org.Hs.eg.db,
                         keys      = c1$V1,
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")
      # c1 <- c1[!is.na(c1)]
      c1 <- na.omit(c1)
      return(c1$entrez)
      }
      cluster.num <- as.character(c(1:3))
      names(cluster.num) <- paste("cluster", 1:3)
    cluster.data <- lapply(cluster.num, cluster.geneID.list)
    ck.GO <- compareCluster(geneCluster = cluster.data,
                            fun         = "enrichGO",
                            OrgDb       = 'org.Hs.eg.db')
    ck.KEGG <- compareCluster(geneCluster = cluster.data,
                              fun         = "enrichKEGG")
    ck.REACTOME <- compareCluster(geneCluster = cluster.data,
                                 fun         = "enrichPathway")
    p1 <- dotplot(ck.GO,
                  title        = "The Most Enriched GO Pathways",
                  showCategory = 8,
                  font.size    = 18,
                  includeAll   = FALSE) +
          theme_bw() +
          theme(plot.title   = black_bold_tahoma_16,
                axis.title   = black_bold_tahoma_16,
                axis.text.x  = black_bold_tahoma_16,
                axis.text.y  = black_bold_tahoma_16,
                axis.line.x  = element_line(color = "black"),
                axis.line.y  = element_line(color = "black"),
                panel.grid   = element_blank(),
                legend.title = black_bold_tahoma_16,
                legend.text  = black_bold_tahoma_16,
                strip.text   = black_bold_tahoma_16)
    print(p1)
    ggsave(
      path        = "~/Documents/EIF_output/Heatmap", 
      filename    = paste("Lung Go.pdf"), 
      plot        = p1,
      width       = 10, 
      height      = 8, 
      useDingbats = FALSE)
    p2 <- dotplot(ck.KEGG,
          title        = "The Most Enriched KEGG Pathways",
          showCategory = 8,
          font.size    = 18,
          includeAll   = FALSE) +
          theme_bw() +
          theme(plot.title   = black_bold_tahoma_16,
                axis.title   = black_bold_tahoma_16,
                axis.text.x  = black_bold_tahoma_16,
                axis.text.y  = black_bold_tahoma_16,
                axis.line.x  = element_line(color = "black"),
                axis.line.y  = element_line(color = "black"),
                panel.grid   = element_blank(),
                legend.title = black_bold_tahoma_16,
                legend.text  = black_bold_tahoma_16,
                strip.text   = black_bold_tahoma_16)
      print(p2)
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("Lung KEGG.pdf"), 
        plot        = p2,
        width       = 12, 
        height      = 8, 
        useDingbats = FALSE)
    p3 <- dotplot(ck.REACTOME,
                  title        = "The Most Enriched REACTOME Pathways",
                  showCategory = 8,
                  font.size    = 16,
                  includeAll   = FALSE) +
          theme_bw() +
          theme(plot.title   = black_bold_tahoma_16,
                axis.title   = black_bold_tahoma_16,
                axis.text.x  = black_bold_tahoma_16,
                axis.text.y  = black_bold_tahoma_16,
                axis.line.x  = element_line(color = "black"),
                axis.line.y  = element_line(color = "black"),
                panel.grid   = element_blank(),
                legend.title = black_bold_tahoma_16,
                legend.text  = black_bold_tahoma_16,
                strip.text   = black_bold_tahoma_16)
    print(p3)
    ggsave(
      path        = "~/Documents/EIF_output/Heatmap", 
      filename    = paste("Lung REACTOME.pdf"), 
      plot        = p3,
      width       = 14, 
      height      = 8, 
      useDingbats = FALSE)
    }
  plot.cluster.pathway()
  
  ### select posCOR for gene expression data ###
  DF.tumor <- as.matrix(na.omit(cor.data[cor.data$EIF4E.tumor > 0.3 |
                                         cor.data$EIF4G1.tumor > 0.3 |
                                         cor.data$EIF4A1.tumor > 0.3, ]))
  gene.list <- row.names(DF.tumor)
  TCGA.GTEX.lung.genelist <- TCGA.GTEX.lung[ , colnames(TCGA.GTEX.lung) %in% gene.list]
  DF2 <- as.matrix(na.omit(TCGA.GTEX.lung.genelist))
  DF2.t <- t(DF2)

  sample.ID <- row.names(DF2)
  sample.ID.type <- Sampletype[Sampletype$sample %in% sample.ID, ]
  
  my_sample_col <- sample.ID.type[ ,
    colnames(sample.ID.type) %in% c("sample", "_sample_type")]
  row.names(my_sample_col) <- my_sample_col$sample
  my_sample_col$sample <- NULL
  my_sample_col$`_sample_type` <- as.factor(my_sample_col$`_sample_type`)
  my_sample_col <- as.data.frame(my_sample_col)
  breaksList = seq(-3, 3, by = 1)
  p4 <- pheatmap::pheatmap(
    DF2.t,
    main           = "Gene expression Heatmap",
    annotation_col =  my_sample_col[,"_sample_type",drop=FALSE],
    scale          = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean", 
    clustering_method        = "complete", # has to be complete
    treeheight_row = 0, treeheight_col = 0,
    cutree_cols = 2,
    angle_col      = c("0"),
    fontsize       = 12,
    fontface       = "bold",
    color          = colorRampPalette(rev(brewer.pal(n    = 7, 
                                                     name = "RdYlBu")))(length(breaksList)),
    breaks = breaksList,
    show_rownames  = FALSE,
    show_colnames  = FALSE)
  print(p4)
  ggsave(
    path        = "~/Documents/EIF_output/Heatmap", 
    filename    = paste("Lung cluster3 heatmap.pdf"), 
    plot        = p4,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  ##############################################
  
}
plot.heatmap.lung(x = "Lung")

### find pathways in the overlapping CORs from all cancer cases
plot.Venn.all <- function() {
  Lung <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
  Lung <- Lung[!duplicated(Lung$sample), ]
  Lung <- na.omit(Lung)
  row.names(Lung) <- Lung$sample
  Lung$sample <- NULL
  Sample.ID <- row.names(Lung)
  subset <- as.data.frame(Lung$`_sample_type`)
  row.names(subset) <- row.names(Lung)
  colnames(subset) <- "sample.type"
  tissue.GTEX.TCGA.gene <- function(){
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE) # data.table = FALSE gives data.frame
    # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
    TCGA.GTEX <- TCGA.GTEX[!duplicated(TCGA.GTEX$sample),
      !duplicated(colnames(TCGA.GTEX))]
    row.names(TCGA.GTEX) <- TCGA.GTEX$sample
    TCGA.GTEX$sample <- NULL
    TCGA.GTEX <- TCGA.GTEX[,colnames(TCGA.GTEX) %in% Sample.ID]
    TCGA.GTEX.t <- data.table::transpose(TCGA.GTEX)
    rownames(TCGA.GTEX.t) <- colnames(TCGA.GTEX)
    colnames(TCGA.GTEX.t) <- rownames(TCGA.GTEX)
    # NA in the vector
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.t,
      subset,
      by    = "row.names",
      all.x = TRUE)
    # check the name of the last column
    # colnames(TCGA.GTEX.Lung.sampletype)[ncol(TCGA.GTEX.Lung.sampletype)] 
    TCGA.GTEX.Lung.sampletype <- na.omit(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype <- tissue.GTEX.TCGA.gene()
  gene.name <- names(TCGA.GTEX.sampletype)
  gene.name <- gene.name [! gene.name %in% c("Row.names", "sample.type")]
  
  ### TO BE CONTINUED
  EIF.correlation <- function(y,z){
    TCGA.GTEX.tumor.lung <- TCGA.GTEX.sampletype[
      TCGA.GTEX.sampletype$sample.type %in% y, ]
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.tumor.lung[[x]],
        TCGA.GTEX.tumor.lung[[y]],
        method = "pearson")
      res <- data.frame(x,
        y,
        result[c("estimate",
          "p.value",
          "statistic",
          "method")],
        stringsAsFactors=FALSE)
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
        lapply(gene.name,
          correlation.coefficient,
          y = x))
      rownames(cor.data) <- cor.data[,1]
      #cor.data1 <- cor.data[cor.data[, "p.value"] <= 0.05,]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    plot.pos.Venn <- function(){
      c3 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(
        A       = a[5, "Counts"],
        B       = a[3, "Counts"],
        C       = a[2, "Counts"],
        "A&B"   = a[7, "Counts"],
        "A&C"   = a[6, "Counts"],
        "B&C"   = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
        #key = TRUE,
        main       = z,
        lwd        = 0,
        fill       = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.25),
        labels     = list(labels = c("EIF4E posCOR",
                                     "EIF4G1 posCOR",
                                     "EIF4A1 posCOR"),
                                     cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all",z,"posVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "EIF4EBP1",
                       "Counts")
      vennDiagram(b)
      pos.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], #EIF4E
        EIF4G1      = b[5, "Counts"], #EIF4G1
        EIF4A1      = b[3, "Counts"], #EIF4A1
        EIF4EBP1    = b[2, "Counts"], #EIF4EBP1
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
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]))
      p2 <- plot(pos.Venn2,
        #key = TRUE,
        main       = paste(z,"posCOR"),
        lwd        = 0,
        fill       = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E",
                                 "EIF4G1",
                                 "EIF4A1",
                                 "EIF4EBP1"),
                                 cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all", z, "pos4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.Venn()
    plot.neg.Venn <- function(){
      c3 <- cbind(
        EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      neg.Venn <- euler(c(A       = a[5, "Counts"],
                          B       = a[3, "Counts"],
                          C       = a[2, "Counts"],
                          "A&B"   = a[7, "Counts"],
                          "A&C"   = a[6, "Counts"],
                          "B&C"   = a[4, "Counts"],
                          "A&B&C" = a[8, "Counts"]))
      p1 <- plot(neg.Venn,
        #key = TRUE,
        main       = z,
        lwd        = 0,
        fill       = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.25),
        labels     = list(labels = c("EIF4E negCOR",
                                     "EIF4G1 negCOR",
                                     "EIF4A1 negCOR"),
                                     cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all",z,"negVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(
        EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "EIF4EBP1",
                       "Counts")
      vennDiagram(b)
      neg.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], #EIF4E
        EIF4G1      = b[5, "Counts"], #EIF4G1
        EIF4A1      = b[3, "Counts"], #EIF4A1
        EIF4EBP1    = b[2, "Counts"], #EIF4EBP1
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
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]))
      p2 <- plot(neg.Venn2,
        #key = TRUE,
        main       = paste(z,"negCOR"),
        lwd        = 0,
        fill       = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels     = list(labels = c("EIF4E",
                                     "EIF4G1",
                                     "EIF4A1",
                                     "EIF4EBP1"),
                                     cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all", z, "neg4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.neg.Venn()
    plot.pos.neg.Venn <- function(){
      c3 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.neg.Venn <- euler(c(
        A       = a[5, "Counts"],
        B       = a[3, "Counts"],
        C       = a[2, "Counts"],
        "A&B"   = a[7, "Counts"],
        "A&C"   = a[6, "Counts"],
        "B&C"   = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.neg.Venn,
        #key = TRUE,
        main       = z,
        lwd        = 0,
        fill       = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.25),
        labels     = list(labels = c("EIF4E posCOR",
                                     "EIF4G1 posCOR",
                                     "EIF4A1 posCOR"),
                                     cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all",z,"posnegVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "EIF4EBP1",
                       "Counts")
      vennDiagram(b)
      pos.neg.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], #EIF4E
        EIF4G1      = b[5, "Counts"], #EIF4G1
        EIF4A1      = b[3, "Counts"], #EIF4A1
        EIF4EBP1    = b[2, "Counts"], #EIF4EBP1
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
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]))
      p2 <- plot(pos.neg.Venn2,
        #key = TRUE,
        main       = z,
        lwd        = 0,
        fill       = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E posCOR",
          "EIF4G1 posCOR",
          "EIF4A1 posCOR",
          "EIF4EBP1 negCOR"),
          cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all", z, "posneg4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.neg.Venn()
    
  }
  all.sample.type <- levels(subset$sample.type)
  all.tumor.type <- all.sample.type [! all.sample.type %in% c("Cell Line", 
    "Normal Tissue", 
    "Solid Tissue Normal")]
  EIF.cor.tumor <- EIF.correlation(y = all.tumor.type, z = "tumor")
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"), z = "normal")
  }
plot.Venn.all()

### find pathways in the overlapping CORs from lung cancer cases
plot.Venn.lung <- function(x) {
  # download https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
  Sampletype <- read_tsv("~/Downloads/TcgaTargetGTEX_phenotype.txt")
  tissue.GTEX.TCGA.gene <- function(x){
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(
      "~/Downloads/TcgaTargetGtex_RSEM_Hugo_norm_count", 
      data.table = FALSE)
    Lung <- Sampletype[Sampletype$`_primary_site` == x,]
    Lung.ID <- as.vector(Lung$sample)
    Lung.ID <- na.omit(Lung.ID) # NA in the vector
    TCGA.GTEX.Lung <- TCGA.GTEX %>% select("sample", Lung.ID)
    #   TCGA.GTEX.Lung1 <- as.data.frame(TCGA.GTEX.Lung)
    TCGA.GTEX.Lung <- TCGA.GTEX.Lung[!duplicated(TCGA.GTEX.Lung$sample), 
      !duplicated(colnames(TCGA.GTEX.Lung))]
    row.names(TCGA.GTEX.Lung) <- TCGA.GTEX.Lung$sample
    TCGA.GTEX.Lung$sample <- NULL
    TCGA.GTEX.Lung.t <- data.table::transpose(TCGA.GTEX.Lung)
    rownames(TCGA.GTEX.Lung.t) <- colnames(TCGA.GTEX.Lung)
    colnames(TCGA.GTEX.Lung.t) <- rownames(TCGA.GTEX.Lung)
    
    Lung <- Lung[!duplicated(Lung$sample), ]
    Lung <- na.omit(Lung)
    row.names(Lung) <- Lung$sample
    Lung$sample <- NULL
    TCGA.GTEX.Lung.sampletype <- merge(TCGA.GTEX.Lung.t,
      Lung,
      by    = "row.names",
      all.x = TRUE)
    TCGA.GTEX.Lung.sampletype <- as.data.frame(TCGA.GTEX.Lung.sampletype)
    return(TCGA.GTEX.Lung.sampletype)
  }
  TCGA.GTEX.sampletype.lung <- tissue.GTEX.TCGA.gene(x)
  row.names(TCGA.GTEX.sampletype.lung) <- TCGA.GTEX.sampletype.lung$Row.names
  TCGA.GTEX.sampletype.lung$Row.names <- NULL
  TCGA.GTEX.sampletype.lung$`_sample_type` <- as.factor(
    TCGA.GTEX.sampletype.lung$`_sample_type`)
  # remove Solid Tissue Normal data
  #TCGA.GTEX.sampletype.lung <- TCGA.GTEX.sampletype.lung[
  #  !(TCGA.GTEX.sampletype.lung$`_sample_type` %in% "Solid Tissue Normal"), ]
  
  geneID <- colnames(Sampletype[Sampletype$`_primary_site` == x,])
  TCGA.GTEX.lung <- TCGA.GTEX.sampletype.lung[ ,
    !names(TCGA.GTEX.sampletype.lung) %in% geneID]
  gene.name <- names(TCGA.GTEX.lung)
  EIF.correlation <- function(y,z){
    TCGA.GTEX.subset.lung <- TCGA.GTEX.sampletype.lung[
      TCGA.GTEX.sampletype.lung$`_sample_type` %in% y, ]
    TCGA.GTEX.subset.lung$`_sample_type` <- droplevels(
      TCGA.GTEX.subset.lung$`_sample_type`)
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.subset.lung[[x]],
        TCGA.GTEX.subset.lung[[y]],
        method = "pearson")
      res <- data.frame(x,
        y,
        result[c("estimate",
          "p.value",
          "statistic",
          "method")],
        stringsAsFactors=FALSE)
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
        lapply(gene.name,
          correlation.coefficient,
          y = x))
      rownames(cor.data) <- cor.data[,1]
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")
    
    plot.pos.Venn <- function(y){
      c3 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.Venn <- euler(c(
        A       = a[5, "Counts"],
        B       = a[3, "Counts"],
        C       = a[2, "Counts"],
        "A&B"   = a[7, "Counts"],
        "A&C"   = a[6, "Counts"],
        "B&C"   = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.Venn,
        #key = TRUE,  
        lwd        = 0,
        # border     = "black",
        fill       = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.5),
        main       = z,
        labels     = list(labels = c("EIF4E posCOR",
                                     "EIF4G1 posCOR",
                                     "EIF4A1 posCOR"),
                                     cex        = 1.5))
      print(p1)
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("Lung",z,"posVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts")
      vennDiagram(b)
      pos.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], #EIF4E
        EIF4G1      = b[5, "Counts"], #EIF4G1
        EIF4A1      = b[3, "Counts"], #EIF4A1
        EIF4EBP1    = b[2, "Counts"], #EIF4EBP1
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
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]))
      p2 <- plot(pos.Venn2,
        #key = TRUE,
        lwd        = 0,
        fill       = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        main       = paste(z,"posCOR"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E",
                                 "EIF4G1",
                                 "EIF4A1",
                                 "EIF4EBP1"),
                                  cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("Lung", z, "pos4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.Venn(y)
    
    plot.neg.Venn <- function(y){
      c3 <- cbind(
        EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      neg.Venn <- euler(c(
        A       = a[5, "Counts"],
        B       = a[3, "Counts"],
        C       = a[2, "Counts"],
        "A&B"   = a[7, "Counts"],
        "A&C"   = a[6, "Counts"],
        "B&C"   = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]))
      p1 <- plot(neg.Venn,
        #key = TRUE,
        lwd        = 0,
        # border     = "black",
        fill       = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.5),
        main       = z,
        labels     = list(labels = c("EIF4E negCOR",
                                     "EIF4G1 negCOR",
                                     "EIF4A1 negCOR"),
                                     cex        = 1.5))
      print(p1)
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("Lung",z,"negVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
                  EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
                  EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
                  EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
                       "EIF4G1",
                       "EIF4A1",
                       "EIF4EBP1",
                       "Counts")
      vennDiagram(b)
      neg.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], #EIF4E
        EIF4G1      = b[5, "Counts"], #EIF4G1
        EIF4A1      = b[3, "Counts"], #EIF4A1
        EIF4EBP1    = b[2, "Counts"], #EIF4EBP1
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
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]))
      p2 <- plot(neg.Venn2,
        #key = TRUE,
        lwd        = 0,
        fill       = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        main       = paste(z,"negCOR"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E",
                                 "EIF4G1",
                                 "EIF4A1",
                                 "EIF4EBP1"),
                                 cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("Lung", z, "neg4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.neg.Venn(y)
    
    plot.pos.neg.Venn <- function(){
      c3 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05)
      a <- vennCounts(c3)
      colnames(a) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "Counts")
      vennDiagram(a)
      ## draw Venn diagram for overlapping genes
      pos.neg.Venn <- euler(c(
        A       = a[5, "Counts"],
        B       = a[3, "Counts"],
        C       = a[2, "Counts"],
        "A&B"   = a[7, "Counts"],
        "A&C"   = a[6, "Counts"],
        "B&C"   = a[4, "Counts"],
        "A&B&C" = a[8, "Counts"]))
      p1 <- plot(pos.neg.Venn,
        #key = TRUE,
        main       = z,
        lwd        = 0,
        fill       = c("#999999", "#E69F00", "#56B4E9"),
        quantities = list(cex = 1.25),
        labels     = list(labels = c("EIF4E posCOR",
          "EIF4G1 posCOR",
          "EIF4A1 posCOR"),
          cex    = 1.25))
      print(p1)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("all",z,"posnegVenn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      c4 <- cbind(
        EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
        EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
        EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
        EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05)
      b <- vennCounts(c4)
      colnames(b) <- c("EIF4E",
        "EIF4G1",
        "EIF4A1",
        "EIF4EBP1",
        "Counts")
      vennDiagram(b)
      pos.neg.Venn2 <- euler(c(
        EIF4E       = b[9, "Counts"], #EIF4E
        EIF4G1      = b[5, "Counts"], #EIF4G1
        EIF4A1      = b[3, "Counts"], #EIF4A1
        EIF4EBP1    = b[2, "Counts"], #EIF4EBP1
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
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]))
      p2 <- plot(pos.neg.Venn2,
        #key = TRUE,
        main       = z,
        lwd        = 0,
        fill       = c("#999999", "#009E73","#56B4E9", "#E69F00"),
        quantities = list(cex = 1.25),
        labels = list(labels = c("EIF4E posCOR",
          "EIF4G1 posCOR",
          "EIF4A1 posCOR",
          "EIF4EBP1 negCOR"),
          cex    = 1.25))
      print(p2)    
      ggsave(
        path        = "~/Documents/EIF_output/Heatmap", 
        filename    = paste("Lung", z, "posneg4Venn.pdf"), 
        plot        = p2,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.neg.Venn()
  }
  EIF.cor.tumor <- EIF.correlation(y = c("Primary Tumor", 
                                         "Metastatic", 
                                         "Recurrent Tumor"),
                                   z = "tumor")
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"), 
                                    z = "normal")
  }
plot.Venn.lung(x = "Lung")

### plot the scatter plot
plot.cor.scatter <- function(x){
  tumor.type <- c("Primary Tumor", "Metastatic")
  TCGA.sampletype.all <- TCGA.sampletype.all[
    TCGA.sampletype.all$sample_type %in% tumor.type, ]
  TCGA.sampletype.all$sample_type <- as.factor(
    TCGA.sampletype.all$sample_type)
  TCGA.sampletype.all$sample_type <- factor(TCGA.sampletype.all$sample_type,
    levels = c("Primary Tumor", "Metastatic"))
  TCGA.sampletype.all$`_primary_disease` <- as.factor(
    TCGA.sampletype.all$`_primary_disease`)
  levels(TCGA.sampletype.all$sample_type)
  levels(TCGA.sampletype.all$`_primary_disease`)
  plotdata <- TCGA.sampletype.all %>% select(x,
    c("EIF4E", "EIF4G1", "EIF4A1"), "sample_type")
  black_bold_tahoma_16 <- element_text(
                                       color  = "black",
                                       face   = "bold",
                                       family = "Tahoma",
                                       size   = 16
                                      )
  p1 <- ggscatter(data = plotdata,
                  x          = x,
                  y          = c("EIF4E", "EIF4G1", "EIF4A1"),
                  # ylim      = c(4, 18),
                  combine    = TRUE,
                  ylab       = "log2(RNA expression)",
                  size       = 1,
                  color      = "sample_type",
                  palette    = "aaas",
                  # facet.by = "sample_type", #scales = "free_x",
                  add        = "reg.line", # Add regression line
                  fullrange  = TRUE, # Extending the regression line
                  font.label = "bold",
                  repel      = TRUE, # repel labels to avoid overlapping
                  shape      = "sample_type", # Change point shape by sample_type
                  rug        = TRUE,    # Add marginal rug
                  conf.int   = TRUE,
                  cor.coef   = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(aes(color       = sample_type),
                                            method      = "pearson",
                                            label.y.npc = "bottom",
                                            size        = 6)) +
                  theme_bw() +
                  theme(
                    plot.title      = black_bold_tahoma_16,
                    axis.title      = black_bold_tahoma_16,
                    axis.text.x     = black_bold_tahoma_16,
                    axis.text.y     = black_bold_tahoma_16,
                    axis.line.x     = element_line(color = "black"),
                    axis.line.y     = element_line(color = "black"),
                    panel.grid      = element_blank(),
                    legend.position = c(0.18, 0.9),
                    legend.title    = element_blank(),
                    legend.text     = black_bold_tahoma_16,
                    strip.text      = black_bold_tahoma_16
                  )
  print(p1)
      }
plot.cor.scatter(x = "CCNB1")
lapply(EIF.gene, plot.cor.scatter, y = "CCNB1")






