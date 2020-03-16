source('R/libraries.R')
source('R/globals.R')
  
source('R/EIFExpressionAcrossTumors.R')
source('R/EIFRatiosAcrossTumors.R')
source('R/EIFExpressionTumorVsAdjacentNormal.R')
source('R/EIFRatioTumorVsAdjacentNormal.R')

plot.boxgraph.EIF.RNAseq.TCGA(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1"))

plot.boxgraph.EIF.ratio.TCGA(c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","PABPC1"))

plot.violingraph.EIF.RNAseq.TCGA (c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"))

plot.violingraph.EIF.ratio.TCGA (c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"))

################################################################
## stacked bar plots for eIF4F CNV status across tumor groups ## 
################################################################
plot.bargraph.EIF.CNV.TCGA <- function (EIF) {
  pan.TCGA.CNV <- function(){
    TCGA.pancancer <- Gistic2CopyNumberAllThresholdedByGenes
    TCGA.sampletype <- TcgaPhenotypeDenseData
    # TCGA.pancancer <- as.data.frame(TCGA.pancancer)
    TCGA.pancancer1 <- TCGA.pancancer[!duplicated(TCGA.pancancer$Sample),
      !duplicated(colnames(TCGA.pancancer))]
    row.names(TCGA.pancancer1) <- TCGA.pancancer1$Sample
    TCGA.pancancer1$Sample <- NULL
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
  TCGA.CNV.anno <- pan.TCGA.CNV()
  pancancer.TCGA.EIF <- function(){
    TCGA.CNV.anno.subset <- TCGA.CNV.anno[
      !TCGA.CNV.anno$sample.type %in% "Solid Tissue Normal", ]
    row.names(TCGA.CNV.anno.subset) <- TCGA.CNV.anno.subset$Row.names
    TCGA.CNV.anno.subset$Row.names <- NULL
    EIF.TCGA.CNV.anno.subset <- TCGA.CNV.anno.subset[ ,
      colnames(TCGA.CNV.anno.subset) %in% c(EIF, 
                                            "sample.type",
                                            "primary.disease")]
    EIF.TCGA.CNV.anno.subset.long <- melt(EIF.TCGA.CNV.anno.subset)
    EIF.TCGA.CNV.anno.subset.long$primary.disease <- as.factor(
      EIF.TCGA.CNV.anno.subset.long$primary.disease)
    colnames(EIF.TCGA.CNV.anno.subset.long) <- c("sample.type","primary.disease",
                                                 "variable","CNV")
    
    CNV.sum <- table(EIF.TCGA.CNV.anno.subset.long[,c("CNV","primary.disease")])
    CNV.sum <- as.data.frame(CNV.sum)
   # CNV.sum$TCGAstudy <- str_remove(CNV.sum$TCGAstudy, regex('_.*\n*.*'))
    CNV.sum$primary.disease <- ordered(CNV.sum$primary.disease, levels = rev(levels(factor(CNV.sum$primary.disease))))
    CNV.sum$CNV <- factor(CNV.sum$CNV, levels = c("-2", "-1", "0", "1", "2"))
    return(CNV.sum)}
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
        values = c('dark blue','blue',
                   'light green','red',
                   'dark red')) 
    print(p1)
    ggsave(
      path        = ".",
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
    TCGA.CNV <- Gistic2CopyNumberAllThresholdedByGenes
    TCGA.sampletype <- TcgaPhenotypeDenseData
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
    df1 <- as.data.frame(df1)
    # correlation plot
    res <- cor(df1)
    cor_5 <- rcorr(as.matrix(df1))
    M <- cor_5$r
    p_mat <- cor_5$P
    pdf(file.path(
      path        = ".",
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
      order       = "FPC", tl.srt = 0, 
      p.mat       = p_mat, 
      sig.level   = 0.05, #insig = "blank" 
    )
    dev.off()
  }
  plot.EIF.CNV.cor()
  
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
    CNV.sum$variable <- factor(CNV.sum$variable, levels = c("PTEN", "EIF4E", "EIF4A1", "MYC", "EIF4EBP1", "EIF4G1"))
  # reorder bars by explicitly ordering factor levels
    p1 <- ggplot(CNV.sum, aes(fill = CNV, 
                              y    = Freq, 
                              x    = variable)) + 
      geom_bar(stat = "identity", position = "fill") +
      labs(x = "Tumor types (TCGA pan cancer atlas 2018)",
           y = paste0("Percentage with CNV (all tumors)")) +
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
      scale_y_continuous(labels = scales::percent_format())+
      guides(fill = guide_legend(reverse = TRUE))+ # Flip ordering of legend without altering ordering in plot
      scale_fill_manual(name   = "Copy number variation",
        breaks = c("-2", "-1", "0", "1", "2"),
        labels = c("Homdel\n 0","Hetlos\n 1","Diploid\n 2","Gain\n 3","Amp\n 3+"),
        values = c('dark blue','blue','light green','red','dark red')) 
    print(p1)
    ggsave(
      path        = ".",
      filename    = "EIFCNVsum.pdf", 
      plot        = p1,
      width       = 7, 
      height      = 7, 
      useDingbats = FALSE)
    }
  make.CNV.sum.plot(EIF)
  
  }
plot.bargraph.EIF.CNV.sum(c("PTEN", "EIF4A1", "EIF4E", "MYC", "EIF4EBP1", "EIF4G1"))

plot.boxgraph.EIF.CNV.RNAseq <- function (EIF) {
  pan.TCGA.gene <- function(EIF){
    TCGA.RNAseq <- EbPanCanIlluminaHiseq
    TCGA.CNV <- Gistic2CopyNumberAllThresholdedByGenes
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
  p1 <- ggplot(
    data = TCGA.RNAseq.CNV,
    aes(x     = CNV,
        y     = RNAseq, color = CNV,
        fill  = CNV)) +    
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
      values = c('dark red','red','light green','blue','dark blue')) +
    scale_x_discrete(breaks = c("2", "1", "0", "-1", "-2"),
      labels = c("Amp\n3+","Gain\n3","Diploid\n2","Hetlos\n1","Homdel\n0")) +
    labs(x = paste(EIF, "copy number variation"),
         y = paste0("log2(", EIF, " RNA counts)")) +
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
      path        = ".",
      filename    = paste0(EIF,"CNV&RNAseq.pdf"), 
      plot        = p1,
      width       = 6.5, 
      height      = 7.5, 
      useDingbats = FALSE)}
  make.plot(EIF)
}
lapply(c("PTEN", "EIF4A1", "EIF4E", "MYC", "EIF4EBP1", "EIF4G1"), 
       plot.boxgraph.EIF.CNV.RNAseq)

#################################################################
##  PCA plots on EIF4F RNA-seq data from TCGA and GTEx groups  ##
#################################################################
plot.EIF.TCGA.GTEX.PCA.all <- function (EIF.list) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <-TcgaTargetGtexPhenoType
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    TCGA.GTEX <- TcgaTargetGtexRsemHugoNormCount
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
                                                          "primary.site"),
                                                          drop = FALSE]
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
    #EIF.TCGA.RNAseq.anno.subset <- na.omit(EIF.TCGA.RNAseq.anno.subset)
    return(EIF.TCGA.RNAseq.anno.subset)
    }
  EIF.TCGA.RNAseq.anno.subset <- get.EIF.TCGA.GTEX(EIF.list)
  ## remove the last two columns 
  df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-2)]
  rownames(df1) <- NULL
  
  plot.PCA.prcomp <- function(){
    # the variables should be scaled to have unit variance 
    PCA <- prcomp(df1, scale = TRUE)
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
      #addEllipses = TRUE, 
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      theme_classic() + 
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
      path        = ".",
      filename    = "EIFPCAall.pdf", 
      plot        = biplot,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
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
      path        = ".",
      filename    = "EIFPCAeig.pdf", 
      plot        = eig,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    var <- get_pca_var(res.pca)
    #fviz_pca_var(res.pca, col.var="contrib")
    
    pdf(file.path(path        = ".",
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
    #              "MKNK1","MKNK2", "MTOR", "RPTOR", "RPS6KB1","MYC")
    EIF.TCGA.RNAseq.anno.subset <- TCGA.GTEX.sampletype[ ,c(EIF.list, 
      "sample.type",
      "primary.site"),
      drop = FALSE]
    EIF.TCGA.RNAseq.anno.subset$`EIF4E+EIF4EBP1` <- log2(2**EIF.TCGA.RNAseq.anno.subset$EIF4E + 2**EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 -2 + 1)
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
      palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"), 
      pointshape = 20,
      pointsize  = 0.75,
      #addEllipses = TRUE, 
      label      = "var",
      col.var    = "black", 
      repel      = TRUE) +
      theme_classic() + 
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
      path        = ".",
      filename    = "EIFsumPCAall.pdf", 
      plot        = biplot,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
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
      path        = ".",
      filename    = "EIFsumPCAeig.pdf", 
      plot        = eig,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    var <- get_pca_var(res.pca)
    #fviz_pca_var(res.pca, col.var="contrib")
    
    pdf(file.path(
      path        = ".",
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
                             "PABPC1", "MKNK1","MKNK2", "MYC"))

plot.EIF.TCGA.GTEX.PCA.each <- function (EIF.list) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <-TcgaTargetGtexPhenoType
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    #TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    TCGA.GTEX <- TcgaTargetGtexRsemHugoNormCount
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
      "primary.site"),
      drop = FALSE]

    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[c(EIF.list,"sample.type","primary.site")]
    
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
    return(EIF.TCGA.RNAseq.anno.subset)
  }
  EIF.TCGA.RNAseq.anno <- get.EIF.TCGA.GTEX(EIF.list)

  get.disease.list <- function (x) {
    x$primary.site <- as.factor(x$primary.site)
    disease.list <- levels(x$primary.site)
    names(disease.list) <- disease.list
    return(disease.list)
  }
  disease.list <- get.disease.list(EIF.TCGA.RNAseq.anno)
  ## remove the last two columns 
  EIF.PCA.tissue <- function(x){
    EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno[
      EIF.TCGA.RNAseq.anno$primary.site == x, ]
    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset)-2)]
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
        #addEllipses = TRUE, 
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
        path        = ".",
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
          axis.text.y      = black_bold_tahoma_16)
      print(eig)
      ggsave(
        path        = ".",
        filename    = paste0(x,"EIFeig.pdf"), 
        plot        = eig,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      var <- get_pca_var(res.pca)
      #fviz_pca_var(res.pca, col.var="contrib")
      pdf(file.path(
        path        = ".",
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
  EIF.PCA.tissue("Lung")
  
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
        #addEllipses = TRUE, 
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
        path        = ".",
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
        path        = ".",
        filename    = paste0(x,"EIFsumeig.pdf"), 
        plot        = eig,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      var <- get_pca_var(res.pca)
      #fviz_pca_var(res.pca, col.var="contrib")
      pdf(file.path(
        path        = ".",
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
  EIFsum.PCA.tissue("Skin")
  }
plot.EIF.TCGA.GTEX.PCA.each(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1", 
                              "PABPC1","MKNK1","MKNK2","MYC"))

plot.EIF.GTEX.PCA.all <- function (EIF.list) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <-TcgaTargetGtexPhenoType
    TCGA.GTEX.anno <- TCGA.GTEX.anno[!duplicated(TCGA.GTEX.anno$sample), ]
    TCGA.GTEX.anno <- na.omit(TCGA.GTEX.anno)
    row.names(TCGA.GTEX.anno) <- TCGA.GTEX.anno$sample
    TCGA.GTEX.anno$sample <- NULL
    Sample.ID <- row.names(TCGA.GTEX.anno)
    TCGA.GTEX.anno <- as.data.frame(TCGA.GTEX.anno) # otherwise lose rownames in the next step, use drop = FALSE to keep the row names 
    subset <- TCGA.GTEX.anno[ ,c("_sample_type", "_primary_site"), drop = FALSE]
    row.names(subset) <- row.names(TCGA.GTEX.anno)
    colnames(subset) <- c("sample.type", "primary.site")
    TCGA.GTEX <- TcgaTargetGtexRsemHugoNormCount
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
  EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    EIF.TCGA.RNAseq.anno.subset$sample.type == "Healthy Tissue (GTEx)", ]
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
    res.pca <- PCA(df1, 
      scale.unit = TRUE, 
      ncp        = 10, 
      graph      = FALSE)
    
    biplot <- fviz_pca_biplot(res.pca, 
      axes       = c(1, 2),
      labelsize  = 5,
      col.ind    = EIF.TCGA.RNAseq.anno.subset$primary.site, 
      #palette    = c("#D55E00","#CC79A7","#009E73","#0072B2"), 
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
        panel.background = element_rect(fill   = 'transparent',
          color  = 'black',
          size   = 1),
        axis.title.x     = black_bold_tahoma_16,
        axis.title.y     = black_bold_tahoma_16,
        axis.text.x      = black_bold_tahoma_16,
        axis.text.y      = black_bold_tahoma_16,
        legend.title      = element_blank(),
        legend.position   = c(0.75, 0.93),
        legend.background = element_blank(),
        legend.text       = black_bold_tahoma_16)
    print(biplot)
    ggsave(
      path        = ".",
      filename    = "EIFPCAall.pdf", 
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
      path        = ".",
      filename    = "EIFPCAeig.pdf", 
      plot        = eig,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
    var <- get_pca_var(res.pca)
    #fviz_pca_var(res.pca, col.var="contrib")
    
    pdf(file.path(path = "PCA", 
      filename = "EIFPCAcor.pdf"), 
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
plot.EIF.TCGA.GTEX.PCA.all(c("EIF4E", "EIF4G1", "EIF4A1","EIF4EBP1",
                             "PABPC1", "MKNK1","MKNK2", "MYC"))

plot.EIF.CPTAC.PCA.LUAD <- function(){
  CPTAC.LUAD.Sample <- read_excel(
    "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx")
  CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample[ ,c("Aliquot (Specimen Label)", "Type")]
  CPTAC.LUAD.Sample.ID <- CPTAC.LUAD.Sample.ID[
    !duplicated(CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`), ]
  row.names(CPTAC.LUAD.Sample.ID) <- CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)`
  CPTAC.LUAD.Sample.ID$`Aliquot (Specimen Label)` <- NULL
  
  CPTAC.LUAD.Proteomics <- fread("CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv",data.table = FALSE)
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
    path        = ".",
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
    path        = ".",
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
    path        = ".",
    filename    = "EIFLUADEig.pdf", 
    plot        = eig,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  var <- get_pca_var(res.pca)
  #fviz_pca_var(res.pca, col.var="contrib")
  pdf(file.path(
    path        = ".",
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
    "S039_Breast_Cancer_Prospective_Collection_Specimens_r1.xlsx")
  CPTAC.BRCA.Sample.ID <- CPTAC.BRCA.Sample[ ,c("Sample Type", "Specimen Label")]
  CPTAC.BRCA.Sample.ID <- CPTAC.BRCA.Sample.ID[
    !duplicated(CPTAC.BRCA.Sample.ID$`Specimen Label`), ]
  CPTAC.BRCA.Sample.ID <- na.omit(CPTAC.BRCA.Sample.ID)
  row.names(CPTAC.BRCA.Sample.ID) <- CPTAC.BRCA.Sample.ID$`Specimen Label`
  CPTAC.BRCA.Sample.ID$`Specimen Label` <- NULL
  
  CPTAC.BRCA.Proteomics <- fread("CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv",
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
    path        = ".",
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
    path        = ".",
    filename    = "EIFBRCAEig.pdf", 
    plot        = eig,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  var <- get_pca_var(res.pca)
  #fviz_pca_var(res.pca, col.var="contrib")
  pdf(file.path(
    path        = ".",
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

#########################################################################
##  Kaplan-Meier curve with clinic and EIF RNASeq data all tumor group ##
#########################################################################
plot.km.EIF.all.tumors <- function(EIF) {
  pan.TCGA.gene <- function(EIF){
    ## get TCGA pancancer RNAseq data ##
    TCGA.RNAseq <- EbPanCanIlluminaHiseq
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
    TCGA.OS <- fread("Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
                        !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- TcgaPhenotypeDenseData
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
      path        = ".",
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
    TCGA.RNAseq <- EbPanCanIlluminaHiseq
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
    TCGA.OS <- fread("Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- TcgaPhenotypeDenseData
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
      TCGA.RNAseq.OS.sampletype[TCGA.RNAseq.OS.sampletype$primary.disease == tumor, ]
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
      path        = ".",
      filename    = paste(EIF, tumor,"KM.pdf"), 
      plot        = KM,
      width       = 8, 
      height      = 8, 
      useDingbats = FALSE)
  }
  plot.KM(EIF, tumor)
}
plot.km.EIF.each.tumor(EIF = "EIF4G1", tumor = "lung adenocarcinoma")

lapply(c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1",
         "EIF4EBP1", "PABPC1", "MKNK1","MKNK2", 
         "MTOR", "RPTOR","RPS6KB1", "MYC"), 
         plot.km.EIF.each.tumor, 
         tumor = "lung adenocarcinoma")

## Cox regression model
plot.coxph.EIF.all.tumors <- function(){
  pan.TCGA.gene <- function(EIF){
  ## get TCGA pancancer RNAseq data ##
  TCGA.RNAseq <- EbPanCanIlluminaHiseq
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
  TCGA.OS <- fread("Survival_SupplementalTable_S1_20171025_xena_sp", 
    data.table = FALSE)
  TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
    !duplicated(colnames(TCGA.OS))]
  row.names(TCGA.OS1) <- TCGA.OS1$sample
  TCGA.OS1$sample <- NULL
  TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
  
  ## get sample type data ##
  TCGA.sampletype <- TcgaPhenotypeDenseData
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
  df$`EIF4E+EIF4EBP1` <- log2(2**df$EIF4E + 2**df$EIF4EBP1 -2 + 1)
  df <- df[c("EIF4E", "EIF4G1","EIF4G2", "EIF4A1","EIF4E+EIF4EBP1",
             "EIF4EBP1","PABPC1", "MKNK1","MKNK2", 
             "MTOR", "RPTOR","RPS6KB1", "MYC","OS","OS.time")]


  # lapply(univ_models, forest_model)
  # Use survivalAnalysis package to draw forest plot of multiple univariate #
  covariate_names <- c(EIF4E = "EIF4E", EIF4G1 = "EIF4G1", 
                      `EIF4E+EIF4EBP1` = "EIF4E+EIF4EBP1",
                       EIF4G2 = "EIF4G2", EIF4A1 = "EIF4A1",
                       EIF4EBP1 = "EIF4EBP1", PABPC1 = "PABPC1", 
                       MKNK1 = "MKNK1", MKNK2 = "MKNK2", 
                       MTOR = "MTOR", RPTOR = "RPTOR",
                       RPS6KB1 = "RPS6KB1", MYC = "MYC")
  p <- map(vars(EIF4E, EIF4G1, EIF4G2, EIF4A1,EIF4E+EIF4EBP1,
                EIF4EBP1, PABPC1, MKNK1, MKNK2, 
                MTOR, RPTOR, RPS6KB1, MYC), function(by)
  {
    analyse_multivariate(df,
      vars(OS.time, OS),
      covariates = list(by), # covariates expects a list
      covariate_name_dict = covariate_names)
  }) %>%
    forest_plot(factor_labeller = covariate_names,
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
  print(p)
  ggsave(
    path        = ".",
    filename    = "EIFuniCox.pdf", 
    plot        = p,
    width       = 5, 
    height      = 4, 
    useDingbats = FALSE)

  # Testing proportional Hazards assumption on univariate analysis
  # graphs of the scaled Schoenfeld residuals against the transformed time 
  # 1
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4E, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "EIF4EuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
    plot.title           = black_bold_tahoma_7,
    axis.title           = black_bold_tahoma_7,
    axis.text            = black_bold_tahoma_7,
    axis.line.x          = element_line(color  = "black"),
    axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
    panel.grid           = element_blank(),
    strip.text           = black_bold_tahoma_7))
  dev.off()

  # 2
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4G1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "EIF4G1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 3
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4G2, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "EIF4G2uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 4
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4A1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "EIF4A1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 5
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4EBP1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "EIF4EBP1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 6
  mv_fit <- coxph(Surv(OS.time, OS) ~ PABPC1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "PABPC1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 7
  mv_fit <- coxph(Surv(OS.time, OS) ~ MKNK1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "MKNK1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 8
  mv_fit <- coxph(Surv(OS.time, OS) ~ MKNK2, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "MKNK2uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 9
  mv_fit <- coxph(Surv(OS.time, OS) ~ MTOR, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "MTORuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 10
  mv_fit <- coxph(Surv(OS.time, OS) ~ RPTOR, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "RPTORuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 11
  mv_fit <- coxph(Surv(OS.time, OS) ~ RPS6KB1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "RPS6KB1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 12
  mv_fit <- coxph(Surv(OS.time, OS) ~ MYC, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "MYCuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title       = black_bold_tahoma_7,
      axis.title       = black_bold_tahoma_7,
      axis.text        = black_bold_tahoma_7,
      axis.line.x      = element_line(color  = "black"),
      axis.line.y      = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid       = element_blank(),
      strip.text       = black_bold_tahoma_7))
  dev.off()
  
  
  
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4G1 +  MKNK2 + MTOR + RPTOR + RPS6KB1 + MYC, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  ggcoxzph(test.ph)
  


  df %>% 
  analyse_multivariate(vars(OS.time, OS),
    covariates = vars(EIF4G1, MKNK2, MTOR, RPTOR, RPS6KB1, MYC),
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
    title             = "Multivariate Cox proportional-hazards regression analysis\nAll Tumors (N=10295)",
    HR_x_breaks       = seq(0.6, 1.7, 0.2), 
    HR_x_limits       = c(0.6, 1.7)) #more breaks on the X axis
  print(p2)
  ggsave(
    path        = ".",
    filename    = "EIFmultiCox.pdf", 
    plot        = p2,
    width       = 5, 
    height      = 4, 
    useDingbats = FALSE)
  
  df %>% 
    analyse_multivariate(vars(OS.time, OS),
      covariates = vars(EIF4G1, EIF4G2, EIF4A1, `EIF4E+EIF4EBP1`,PABPC1, 
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
  ggsave(
    path        = ".",
    filename    = "EIFsummultiCox.pdf", 
    plot        = p3,
    width       = 5, 
    height      = 4, 
    useDingbats = FALSE)
  
  ggcoxdiagnostics(mv_fit, 
    type               = "dfbeta",
    linear.predictions = FALSE, 
    ggtheme            = theme_bw())
  cz <- cox.zph(mv_fit)
  ggcoxzph(cz)
  
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
    path        = ".",
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
    path        = ".",
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
    TCGA.RNAseq <- EbPanCanIlluminaHiseq
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
    TCGA.OS <- fread("Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- TcgaPhenotypeDenseData
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
    values_displayed  = c("HR", "CI", "p"),
    value_headers     = c(HR = "HR", CI = "95%CI", p = "p", n = "N"),
    relative_widths   = c(0.6, 1.8, 1.2), #more space for the plot, less space for the tables
    label_headers     = c(factor = "Gene"),
    title             = "Univariate Cox proportional-hazards regression analysis\nLung Adenocarcinoma (N=517)",
    HR_x_breaks       = seq(0.5, 2.1, 0.2), 
    HR_x_limits       = c(0.5, 2.1),
    ggtheme           = ggplot2::theme_bw(base_size = 7))
print(p)
ggsave(
  path        = ".",
  filename    = paste0(tumor, "EIFuniCox.pdf"),
  plot        = p,
  width       = 5, 
  height      = 4, 
  useDingbats = FALSE)

plot.cox.assumption <- function() {
# 1
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4E, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungEIF4EuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()

  # 2
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4G1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungEIF4G1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 3
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4G2, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungEIF4G2uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 4
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4A1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungEIF4A1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 5
  mv_fit <- coxph(Surv(OS.time, OS) ~ EIF4EBP1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungEIF4EBP1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 6
  mv_fit <- coxph(Surv(OS.time, OS) ~ PABPC1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungPABPC1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 7
  mv_fit <- coxph(Surv(OS.time, OS) ~ MKNK1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungMKNK1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 8
  mv_fit <- coxph(Surv(OS.time, OS) ~ MKNK2, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungMKNK2uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 9
  mv_fit <- coxph(Surv(OS.time, OS) ~ MTOR, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungMTORuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 10
  mv_fit <- coxph(Surv(OS.time, OS) ~ RPTOR, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungRPTORuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 11
  mv_fit <- coxph(Surv(OS.time, OS) ~ RPS6KB1, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungRPS6KB1uniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title           = black_bold_tahoma_7,
      axis.title           = black_bold_tahoma_7,
      axis.text            = black_bold_tahoma_7,
      axis.line.x          = element_line(color  = "black"),
      axis.line.y          = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid           = element_blank(),
      strip.text           = black_bold_tahoma_7))
  dev.off()
  
  # 12
  mv_fit <- coxph(Surv(OS.time, OS) ~ MYC, data = df)
  test.ph <- cox.zph(mv_fit)
  test.ph
  pdf(file.path(
    path        = ".",
    filename    = "LungMYCuniCox.pdf"), 
    width       = 3, 
    height      = 3, 
    useDingbats = FALSE)
  ggcoxzph(test.ph, point.size = 0.1,  
    ggtheme              =    theme(
      plot.title       = black_bold_tahoma_7,
      axis.title       = black_bold_tahoma_7,
      axis.text        = black_bold_tahoma_7,
      axis.line.x      = element_line(color  = "black"),
      axis.line.y      = element_line(color  = "black"),
      panel.background = element_blank(),
      panel.grid       = element_blank(),
      strip.text       = black_bold_tahoma_7))
  dev.off()
}
plot.cox.assumption()

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
ggsave(
  path        = ".",
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
  title             = "Multivariate Cox proportional-hazards regression analysis\nLung Adenocarcinoma (N=517)",
  HR_x_breaks       = seq(0.5, 2.7, 0.3), 
  HR_x_limits       = c(0.5, 2.7)) #more breaks on the X axis
print(p3)
ggsave(
  path        = ".",
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
    path        = ".",
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
    TCGA.RNAseq <- EbPanCanIlluminaHiseq
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
    TCGA.OS <- fread("Survival_SupplementalTable_S1_20171025_xena_sp", 
      data.table = FALSE)
    TCGA.OS1 <- TCGA.OS[!duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[ ,c("OS","OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- TcgaPhenotypeDenseData
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
    path        = ".",
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
    path        = ".",
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
    path        = ".",
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
    path        = ".",
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
####################################
### use TCGA-TARGET-GTEX dataset ###
plot.heatmap.total <- function() {
  Lung <-TcgaTargetGtexPhenoType
  Lung <- Lung[!duplicated(Lung$sample), ]
  Lung <- na.omit(Lung)
  row.names(Lung) <- Lung$sample
  Lung$sample <- NULL
  Sample.ID <- row.names(Lung)
  subset <- as.data.frame(Lung$`_sample_type`)
  row.names(subset) <- row.names(Lung)
  colnames(subset) <- "sample.type"
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX <- TcgaTargetGtexRsemHugoNormCount
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
  EIF.correlation <- function(y){
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
      return(cor.data)
    }
    EIF4E.cor <- EIF.cor.list("EIF4E")
    EIF4G1.cor <- EIF.cor.list("EIF4G1")
    EIF4A1.cor <- EIF.cor.list("EIF4A1")
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')))
    return(cor.data)
  }
  all.sample.type <- levels(subset$sample.type)
  all.tumor.type <- all.sample.type [! all.sample.type %in% c("Cell Line", 
                                                              "Normal Tissue", 
                                                              "Solid Tissue Normal")]
  EIF.cor.tumor <- EIF.correlation(y = all.tumor.type)
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"))
  cor.data <- cbind(setNames(data.frame(EIF.cor.tumor[1:3]),
                    c('EIF4E.tumor',
                      'EIF4G1.tumor',
                      'EIF4A1.tumor')),
                    setNames(data.frame(EIF.cor.normal[1:3]),
                    c('EIF4E.normal',
                      'EIF4G1.normal',
                      'EIF4A1.normal')))
  DF <- as.matrix(na.omit(cor.data[cor.data$EIF4E.tumor > 0.3 |
                                   #cor.data$EIF4E.tumor < -0.3 |
                                   cor.data$EIF4G1.tumor > 0.3 |
                                   #cor.data$EIF4G1.tumor < -0.3 |
                                   cor.data$EIF4A1.tumor > 0.3 |
                                   #cor.data$EIF4A1.tumor < -0.3 |
                                   cor.data$EIF4E.normal > 0.3 |
                                   #cor.data$EIF4E.normal < -0.3 |
                                   cor.data$EIF4G1.normal > 0.3 |
                                   #cor.data$EIF4G1.normal < -0.3 |
                                   cor.data$EIF4A1.normal > 0.3 
                                   #cor.data$EIF4A1.normal < -0.3
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
      type      = c("tumor", "tumor","tumor",
                    "normal","normal","normal"),
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
    path        = ".",
    filename    = "All tumors heatmap.pdf"), 
    width       = 6, 
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
          path        = ".",
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
      path        = ".",
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
      path        = ".",
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
    Sample.type.annotation <-TcgaTargetGtexPhenoType
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
    TCGA.GTEX <- TcgaTargetGtexRsemHugoNormCount
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
        path        = ".",
        filename    = paste("all healthy Venn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      }
    plot.pos.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')))
    return(cor.data)
  }
  # all.sample.type <- levels(subset$sample.type)
  EIF.cor.normal <- EIF.correlation(y = "Normal Tissue")
  DF <- as.matrix(na.omit(EIF.cor.normal[EIF.cor.normal$EIF4E > 0.3 |
                                         EIF.cor.normal$EIF4G1 > 0.3 |
                                         EIF.cor.normal$EIF4A1 > 0.3, ]))
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
  get.GTEX.annotation <- function(){
    Sample.type.annotation <-TcgaTargetGtexPhenoType
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
  subset <- get.GTEX.annotation()
  Sample.ID <- row.names(subset)
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX <- TcgaTargetGtexRsemHugoNormCount
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
        path        = ".",
        filename    = paste("all tumor Venn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
      
      }
    plot.pos.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')))
    return(cor.data)
  }
  # all.sample.type <- levels(subset$sample.type)
  EIF.cor.normal <- EIF.correlation(y = "Solid Tissue Normal")
  DF <- as.matrix(na.omit(EIF.cor.normal[EIF.cor.normal$EIF4E > 0.3 |
                                         EIF.cor.normal$EIF4G1 > 0.3 |
                                         EIF.cor.normal$EIF4A1 > 0.3, ]))
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
                name                 = "Correlation Coefficient Heatmap in TCGA",
                heatmap_legend_param = list(direction = "horizontal"),
                # clustering_distance_rows = function(x, y) 1 - cor(x, y),
                show_row_names       = FALSE,
                show_column_names    = FALSE,
                bottom_annotation    = HeatmapAnnotation(
                  annotation_legend_param = list(direction = "horizontal"),
                  cn       = anno_text(gsub("\\..*","",colnames(DF)),
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
    TCGA.pancancer <- EbPanCanIlluminaHiseq
    TCGA.sampletype <- TcgaPhenotypeDenseData
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
      }
    plot.pos.Venn()
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')))
    return(cor.data)
    }
  EIF.cor.tumor <- EIF.cor.tumor()
  DF  <- as.matrix(na.omit(EIF.cor.tumor[EIF.cor.tumor$EIF4E  > 0.3 |
                                         EIF.cor.tumor$EIF4G1 > 0.3 |
                                         EIF.cor.tumor$EIF4A1 > 0.3 , ]))
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
    TCGA.pancancer <- fread("gtex_RSEM_Hugo_norm_count", data.table=FALSE)
    # download https://toil.xenahubs.net/download/GTEX_phenotype.gz
    TCGA.sampletype <- read_tsv("GTEX_phenotype")
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
  Sampletype <-TcgaTargetGtexPhenoType
  tissue.GTEX.TCGA.gene <- function(x){
    TCGA.GTEX <- TcgaTargetGtexRsemHugoNormCount
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
        path        = ".",
        filename    = paste("Lung Venn.pdf"), 
        plot        = p1,
        width       = 8, 
        height      = 8, 
        useDingbats = FALSE)
    }
    plot.pos.Venn(y)
    cor.data <- cbind(setNames(data.frame(EIF4E.cor[3]), c('EIF4E')),
                      setNames(data.frame(EIF4G1.cor[3]), c('EIF4G1')),
                      setNames(data.frame(EIF4A1.cor[3]), c('EIF4A1')))
    return(cor.data)
    }
  EIF.cor.tumor <- EIF.correlation(y = c("Primary Tumor", 
                                         "Metastatic", 
                                         "Recurrent Tumor"))
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"))
  cor.data <- cbind(setNames(data.frame(EIF.cor.tumor[1:3]),
                             c('EIF4E.tumor','EIF4G1.tumor','EIF4A1.tumor')),
                    setNames(data.frame(EIF.cor.normal[1:3]),
                             c('EIF4E.normal','EIF4G1.normal','EIF4A1.normal')))
  DF <- as.matrix(na.omit(cor.data[cor.data$EIF4E.tumor > 0.3 |
                                   cor.data$EIF4G1.tumor > 0.3 |
                                   cor.data$EIF4A1.tumor > 0.3 |
                                   cor.data$EIF4E.normal > 0.3 |
                                   cor.data$EIF4G1.normal > 0.3 |
                                   cor.data$EIF4A1.normal > 0.3 ,  ]))
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
      type     = c("tumor", "tumor","tumor",
                   "normal","normal","normal"),
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
    path        = ".",
    filename    = "Lung tumors heatmap.pdf"), 
    width       = 6, 
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
      path        = ".",
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
        path        = ".",
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
      path        = ".",
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
    path        = ".",
    filename    = paste("Lung cluster3 heatmap.pdf"), 
    plot        = p4,
    width       = 8, 
    height      = 8, 
    useDingbats = FALSE)
  ##############################################
  
}
plot.heatmap.lung(x = "Lung")

### find pathways in the overlapping CORs from all cancer cases
plot.EIF.cor.pathway.all <- function() {
  TCGA.RNAseq.sampletype <- TCGA.sampletype.all[
    TCGA.sampletype.all$sample_type != "Solid Tissue Normal", ]
  EIF.correlation <- function(x, y) {
    result <- cor.test(TCGA.RNAseq.sampletype[[x]],
                       TCGA.RNAseq.sampletype[[y]],
                       method = "pearson")
    res <- data.frame(x,
                      y,
                      result[c("estimate",
                               "p.value",
                               "statistic",
                               "method")],
                      stringsAsFactors = FALSE)
    }
# find all genes positively correlate with EIF4F expression
# lapply function gives a large list, need to convert it to a dataframe
  TCGA.RNAseq.sampletype <- TCGA.RNAseq.sampletype[ ,
    !names(TCGA.RNAseq.sampletype) %in% c("Row.names",
                                          "sample_type",
                                          "_primary_disease")]
  gene.name <- names(TCGA.RNAseq.sampletype)
  EIF.cor.list <- function(x) {
    cor.data <- do.call(rbind.data.frame,
                        lapply(gene.name,
                               EIF.correlation,
                               y = x))}
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
               labels = list(labels = c("EIF4E posCOR",
                                        "EIF4G1 posCOR",
                                        "EIF4A1 posCOR"),
                             cex = 1.25))
    print(p1)}
  plot.pos.Venn()
  plot.neg.Venn <- function(){
    c4 <- cbind(EIF4E.cor$estimate < -0.3,
                EIF4G1.cor$estimate < -0.3,
                EIF4A1.cor$estimate < -0.3)
    b <- vennCounts(c4)
    colnames(b) <- c("EIF4E",
                     "EIF4G1",
                     "EIF4A1",
                     "Counts")
    vennDiagram(b)
    neg.Venn <- euler(c(A       = b[5, "Counts"],
                        B       = b[3, "Counts"],
                        C       = b[2, "Counts"],
                        "A&B"   = b[7, "Counts"],
                        "A&C"   = b[6, "Counts"],
                        "B&C"   = b[4, "Counts"],
                        "A&B&C" = b[8, "Counts"]))
    p2 <- plot(neg.Venn,
      #key = TRUE,
               lwd = 0,
               fill = c("#999999", "#E69F00", "#56B4E9"),
               quantities = list(cex = 1.25),
               labels = list(labels = c("EIF4E negCOR",
                                        "EIF4G1 negCOR",
                                        "EIF4A1 negCOR"),
                             cex    = 1.25))
  print(p2)}
  plot.neg.Venn()
# perform pathway analysis on overlapping genes
  plot.pos.pathway <- function(){
    EIF4E.pos <- na.omit(EIF4E.cor[EIF4E.cor$estimate > 0.3, ])
    EIF4G1.pos <- na.omit(EIF4G1.cor[EIF4G1.cor$estimate > 0.3, ])
    EIF4A1.pos <- na.omit(EIF4A1.cor[EIF4A1.cor$estimate > 0.3, ])
    pos.overlap <- merge(EIF4G1.pos, EIF4A1.pos, by.x = "x", by.y = "x")
  # overlap <- merge(EIF4G1.neg, EIF4A1.neg, by.x = "x", by.y = "x")
    pos.overlap$entrez = mapIds(org.Hs.eg.db,
                                keys      = pos.overlap$x,
                                column    = "ENTREZID",
                                keytype   = "SYMBOL",
                                multiVals = "first")
    reac <- enrichPathway(gene = pos.overlap$entrez, readable = T)
    p3 <- barplot(reac, showCategory = 8)
    p4 <- dotplot(reac,
                  x = "count",
                  font.size = 18,
                  title = "posCOR EIF4A1 & EIF4G1")
    p5 <- emapplot(reac, font.size = 18)
    p6 <- cnetplot(reac)
    print(p3)
    print(p4)
    print(p5)
    print(p6)

    allthree <- Reduce(function(x, y) merge(x, y, by = "x", all=TRUE),
      list(EIF4E.pos, EIF4G1.pos, EIF4A1.pos))
    allthree <- na.omit(allthree)
    allthree$entrez = mapIds(org.Hs.eg.db,
                             keys      = allthree$x,
                             column    = "ENTREZID",
                             keytype   = "SYMBOL",
                             multiVals = "first")
    y <- enrichPathway(gene = allthree$entrez, readable = T)
    p7 <- barplot(y, showCategory = 8)
    p8 <- emapplot(y, font.size = 18)
    p9 <- dotplot(y, x = "count", font.size = 18,
                  title = "posCOR EIF4E & EIF4A1 & EIF4G1")
    p10 <- cnetplot(y, circular = TRUE, colorEdge = TRUE)  +
           theme(legend.position="none")
    print(p7)
    print(p8)
    print(p9)
    print(p10)
    }
  plot.pos.pathway()
  plot.neg.pathway <- function(){
    EIF4E.neg <- na.omit(EIF4E.cor[EIF4E.cor$estimate < -0.3, ])
    EIF4G1.neg <- na.omit(EIF4G1.cor[EIF4G1.cor$estimate < -0.3, ])
    EIF4A1.neg <- na.omit(EIF4A1.cor[EIF4A1.cor$estimate < -0.3, ])
    neg.overlap <- merge(EIF4G1.neg, EIF4A1.neg, by.x = "x", by.y = "x")
    neg.overlap$entrez = mapIds(org.Hs.eg.db,
                                keys      = neg.overlap$x,
                                column    = "ENTREZID",
                                keytype   = "SYMBOL",
                                multiVals = "first")
    x <- enrichPathway(gene = neg.overlap$entrez, readable = T)
    p3 <- barplot(x, showCategory = 8)
    p4 <- dotplot(x, font.size = 18, title = "negCOR with EIF4A1 & EIF4G1")
    # p5 <- emapplot(x, font.size = 18)
    p6 <- cnetplot(x)
    print(p3)
    print(p4)
    print(p5)
    print(p6)

    allthree <- Reduce(function(x, y) merge(x, y, by = "x", all=TRUE),
      list(EIF4E.neg, EIF4G1.neg, EIF4A1.neg))
    allthree <- na.omit(allthree)
    allthree$entrez = mapIds(org.Hs.eg.db,
                             keys      = allthree$x,
                             column    = "ENTREZID",
                             keytype   = "SYMBOL",
                             multiVals = "first")
    y <- enrichPathway(gene = allthree$entrez, readable = T)
    p7 <- barplot(y, showCategory = 8)
    # p8 <- emapplot(y, font.size = 18)
    p9 <- dotplot(y,
                  font.size = 18,
                  title = "negCOR with EIF4E & EIF4A1 & EIF4G1")
    p10 <- cnetplot(y)
    print(p7)
    print(p8)
    print(p9)
    print(p10)
    }
  plot.neg.pathway()
  }
plot.EIF.cor.pathway.all()

### find pathways in the overlapping CORs from lung cancer cases
plot.EIF.cor.pathway.lung <- function() {
  TCGA.RNAseq.sampletype <- TCGA.GTEX.sampletype.lung[
    TCGA.GTEX.sampletype.lung$`_sample_type` == "Normal Tissue", ]
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
  Sampletype <-TcgaTargetGtexPhenoType
  Lung <- Sampletype[Sampletype$`_primary_site` == "Lung",]
  geneID <- colnames(Lung)
  TCGA.RNAseq.sampletype <- TCGA.RNAseq.sampletype[ ,
    !names(TCGA.RNAseq.sampletype) %in% c("Row.names", geneID)]
  gene.name <- names(TCGA.RNAseq.sampletype)
  EIF.cor.list <- function(x) {
    cor.data <- do.call(rbind.data.frame,
                        lapply(gene.name,
                               EIF.correlation,
                               y = x))}
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
               border = "black",
               fill = c("#999999", "#E69F00", "#56B4E9"),
               quantities = list(cex = 1.5),
               labels = list(labels = c("EIF4E posCOR",
                                        "EIF4G1 posCOR",
                                        "EIF4A1 posCOR"),
                             cex    = 1.5))
    print(p1)
    }
  plot.pos.Venn()

  plot.neg.Venn <- function(){
    c4 <- cbind(EIF4E.cor$estimate < -0.3,
                EIF4G1.cor$estimate < -0.3,
                EIF4A1.cor$estimate < -0.3)
    b <- vennCounts(c4)
    colnames(b) <- c("EIF4E",
                     "EIF4G1",
                     "EIF4A1",
                     "Counts")
    vennDiagram(b)
    neg.Venn <- euler(c(A       = b[5, "Counts"],
                        B       = b[3, "Counts"],
                        C       = b[2, "Counts"],
                        "A&B"   = b[7, "Counts"],
                        "A&C"   = b[6, "Counts"],
                        "B&C"   = b[4, "Counts"],
                        "A&B&C" = b[8, "Counts"]))
    p2 <- plot(neg.Venn,
      #key = TRUE,
               lwd = 0,
               fill = c("#999999", "#E69F00", "#56B4E9"),
               quantities = list(cex = 1.25),
               labels = list(labels = c("EIF4E negCOR",
                                        "EIF4G1 negCOR",
                                        "EIF4A1 negCOR"),
                             cex    = 1.25))
    print(p2)}
  plot.neg.Venn()
  # perform pathway analysis on overlapping corr genes
  plot.pos.pathway <- function(){
    EIF4E.pos <- na.omit(EIF4E.cor[EIF4E.cor$estimate > 0.3, ])
    EIF4G1.pos <- na.omit(EIF4G1.cor[EIF4G1.cor$estimate > 0.3, ])
    EIF4A1.pos <- na.omit(EIF4A1.cor[EIF4A1.cor$estimate > 0.3, ])
    pos.overlap <- merge(EIF4G1.pos, EIF4A1.pos, by.x = "x", by.y = "x")
    # overlap <- merge(EIF4G1.neg, EIF4A1.neg, by.x = "x", by.y = "x")
    pos.overlap$entrez = mapIds(org.Hs.eg.db,
                                keys      = pos.overlap$x,
                                column    = "ENTREZID",
                                keytype   = "SYMBOL",
                                multiVals = "first")
    reac <- enrichPathway(gene = pos.overlap$entrez, readable = T)
    p3 <- barplot(reac, showCategory = 8)
    p4 <- dotplot(reac,
                  x         = "count",
                  font.size = 18,
                  title     = "posCOR EIF4A1 & EIF4G1")
    p5 <- emapplot(reac, font.size = 18)
    p6 <- cnetplot(reac)
    print(p3)
    print(p4)
    print(p5)
    #print(p6)

    allthree <- Reduce(function(x, y) merge(x, y, by = "x", all=TRUE),
      list(EIF4E.pos, EIF4G1.pos, EIF4A1.pos))
    allthree <- na.omit(allthree)
    allthree$entrez <- mapIds(org.Hs.eg.db,
                              keys      = allthree$x,
                              column    = "ENTREZID",
                              keytype   = "SYMBOL",
                              multiVals = "first")
    y <- enrichPathway(gene = allthree$entrez, readable = T)
    p7 <- barplot(y, showCategory = 8)
    p8 <- emapplot(y, font.size = 18)
    p9 <- dotplot(y,
                  x         = "count",
                  font.size = 18,
                  title     = "posCOR EIF4E & EIF4A1 & EIF4G1")
    p10 <- cnetplot(y)
    print(p7)
    print(p8)
    print(p9)
    # print(p10)
  }
  plot.pos.pathway()
  plot.neg.pathway <- function(){
    EIF4E.neg <- na.omit(EIF4E.cor[EIF4E.cor$estimate < -0.3, ])
    EIF4G1.neg <- na.omit(EIF4G1.cor[EIF4G1.cor$estimate < -0.3, ])
    EIF4A1.neg <- na.omit(EIF4A1.cor[EIF4A1.cor$estimate < -0.3, ])
    neg.overlap <- merge(EIF4G1.neg, EIF4A1.neg, by.x = "x", by.y = "x")
    neg.overlap$entrez = mapIds(org.Hs.eg.db,
                                keys      = neg.overlap$x,
                                column    = "ENTREZID",
                                keytype   = "SYMBOL",
                                multiVals = "first")
    x <- enrichPathway(gene = neg.overlap$entrez, readable = T)
    p3 <- barplot(x, showCategory = 8)
    p4 <- dotplot(x, font.size = 16, title = "negCOR EIF4A1 & EIF4G1")
    p5 <- emapplot(x, font.size = 18)
    p6 <- cnetplot(x)
    print(p3)
    print(p4)
    print(p5)
    #print(p6)

    allthree <- Reduce(function(x, y) merge(x, y, by = "x", all=TRUE),
      list(EIF4E.neg, EIF4G1.neg, EIF4A1.neg))
    allthree <- na.omit(allthree)
    allthree$entrez = mapIds(org.Hs.eg.db,
                             keys      = allthree$x,
                             column    = "ENTREZID",
                             keytype   = "SYMBOL",
                             multiVals = "first")
    y <- enrichPathway(gene = allthree$entrez, readable = T)
    p7 <- barplot(y, showCategory = 8)
    p8 <- emapplot(y,
      font.size = 18)
    p9 <- dotplot(y,
      font.size = 18,
      title = "negCOR EIF4E & EIF4A1 & EIF4G1")
    p10 <- cnetplot(y)
    print(p7)
    print(p8)
    print(p9)
    #print(p10)
    }
  plot.neg.pathway()
  }
plot.EIF.cor.pathway.lung()

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







###################################################################
###################################################################

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
  EIF.TCGA.GTEX <- read.csv(file.path("project-data","EIFTCGAGTEX.csv"),
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
  EIF.TCGA.GTEX <- read.csv(file.path("project-data","EIFTCGAGTEX.csv"),
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
  EIF.TCGA.GTEX <- read.csv(file.path("project-data","EIFTCGAGTEX.csv"),
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
  EIF.TCGA.GTEX <- read.csv(file.path("project-data","EIFTCGAGTEX.csv"),
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
  EIF.TCGA.GTEX <- read.csv(file.path("project-data","EIFTCGAGTEX.csv"),
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
plot.EIF.correlation <- function(x, y){
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
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA', ]
  cor.test(EIF.TCGA.GTEX$EIF4A1, EIF.TCGA.GTEX$MYC, method = "pearson")
  black_bold_tahoma_16 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 16
  )
  p1 <- ggscatter(EIF.TCGA.GTEX,
            x        = x,
            y        = y,
            size     = 0.3,
            color    = "sample_type",
            # palette  = "jco",
            facet.by = "sample_type", #scales = "free_x",
            add      = "reg.line",
            conf.int = TRUE) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_16,
      axis.title      = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16,
      axis.text.y     = black_bold_tahoma_16,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_16
    ) +
    stat_cor(#aes(color   = "sample_type"),
                 method  = "pearson",
                 label.y = 6)
  p2 <- ggscatter(EIF.TCGA,
            x        = x,
            y        = y,
            size     = 0.3,
            color    = "primary_disease_or_tissue",
            #palette  = "nrc",
            facet.by = "primary_disease_or_tissue", #scales = "free_x",
            add      = "reg.line",
            conf.int = TRUE) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_16,
      axis.title      = black_bold_tahoma_16,
      axis.text.x     = black_bold_tahoma_16,
      axis.text.y     = black_bold_tahoma_16,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_16
    ) +
    # discrete_scale("fill", "manual", palette_Dark2)+
    stat_cor(#aes(color   = "primary_disease_or_tissue"),
                 method  = "pearson",
                 label.y = 6)
  print(p1)
  print(p2)
  }

##
plotEIF.RNAseq.TCGA.GTEX <- function (x) {
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
plotEIF.RNAseq.TCGA.GTEX.tissue <- function (tissue) {
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
plot.box.EIF.RNAseq.TCGA <- function (x) {
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
plot.box.EIF.RNAseq.TCGA (get.EIF.TCGA.RNAseq.long(EIF.gene))

##
plotEIF.RNAseq.TCGA <- function (x) {
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
plot.box.EIF.score.TCGA <- function (x) {
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
plot.box.EIF.score.TCGA (get.EIF.TCGA.score.long(EIF.gene))
  
##
plotEIF.score.TCGA <- function (x) {
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
plotEIF.score.TCGA (get.EIF.TCGA.score.long(EIF.gene))

##
plotEIF.score.TCGA.GTEX.tissue <- function (tissue) {
  get.EIF.TCGA.GTEX.score.tissue <- function (tissue) {
    EIF.TCGA.GTEX <- read.csv(file.path("project-data","EIFTCGAGTEX.csv"),
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
plotEIF.GTEX <- function (x) {
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
  EIF.TCGA.GTEX <- read_csv("project-data/EIFTCGAGTEX.csv",
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
  EIF.TCGA.GTEX <- read_csv("project-data/EIFTCGAGTEX.csv",
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

######################
## plot correlation ##
######################

plot.EIF.correlation.pathway()
#################################################################
##  PCA plots on EIF4F RNA-seq data from TCGA and GTEx groups  ##
#################################################################
## the following script perform PCA on RNA-Seq data of seven EIF genes
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
  EIF.proteomics <- read.csv(file.path("project-data","proteomics.csv"),
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
plot.EIF.correlation(x = "EIF4A1", y = "MYC")

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

