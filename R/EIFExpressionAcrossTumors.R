##############################################
## boxplot for EIF expression across tumors ##
##############################################
plot.boxgraph.EIF.RNAseq.TCGA <- function (EIF.gene) {
  pan.TCGA.gene <- function(){
    TCGA.pancancer <- EbPanCanIlluminaHiseq
    TCGA.sampletype <- TcgaTargetGtexPhenoType
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
      path        = ".",
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
      path        = ".",
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
                 aes(
                   x     = f.ordered1,  
                   #x     = x.ordered, # order primary disease
                   y     = 2**value,
                   color = variable)) +
      scale_y_continuous(trans = log2_trans(), labels = label_comma()) +
      stat_n_text(size = 5, fontface = "bold", hjust = 0) +
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
      path        = ".",
      filename    = "EIFexpression.pdf", 
      plot        = p1,
      width       = 9, 
      height      = 9, 
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
      path        = ".",
      filename    = "EIFsumexpression.pdf", 
      plot        = p2,
      width       = 7, 
      height      = 9, 
      useDingbats = FALSE)
  }
  make.plot()
  
}
