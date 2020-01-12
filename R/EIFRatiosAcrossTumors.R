##########################################
## boxplot for EIF ratios across tumors ##
##########################################
plot.boxgraph.EIF.ratio.TCGA <- function (EIF.gene) {
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
    ratio <- c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4EBP1", "EIF4G1:\nEIF4E", 
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
  
  pancancer.TCGA.EIF.ratio.long$label <- sub(".*\n", "", pancancer.TCGA.EIF.ratio.long$variable)
  pancancer.TCGA.EIF.ratio.long$label <- factor(pancancer.TCGA.EIF.ratio.long$label, levels = c("EIF4E+EIF4EBP1","EIF4E", "EIF4EBP1"))
  make.plot <- function () {
    mean <- within(pancancer.TCGA.EIF.ratio.long[
      pancancer.TCGA.EIF.ratio.long$variable == "EIF4EBP1:\nEIF4E", ], # TCGAstudy is one column in df2
      primary.disease <- reorder(primary.disease, value, median))
    mean$primary.disease <- as.factor(mean$primary.disease)
    neworder <- levels(mean$primary.disease)
    x.ordered <- factor(pancancer.TCGA.EIF.ratio.long$primary.disease, 
                        levels = neworder)
    levels(pancancer.TCGA.EIF.ratio.long$variable)
    
    pancancer.TCGA.EIF.ratio.long1 <- pancancer.TCGA.EIF.ratio.long[
      pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E","EIF4A1:\nEIF4G1", "EIF4G1:\nEIF4E", "EIF4EBP1:\nEIF4E"), ]
    pancancer.TCGA.EIF.ratio.long1$variable <- factor(pancancer.TCGA.EIF.ratio.long1$variable, levels = c("EIF4A1:\nEIF4G1", "EIF4EBP1:\nEIF4E","EIF4A1:\nEIF4E","EIF4G1:\nEIF4E"))
    f1 <- factor(pancancer.TCGA.EIF.ratio.long1$primary.disease)
    f.ordered1 <- fct_rev(f1)
    p1 <- ggplot(data = pancancer.TCGA.EIF.ratio.long1,
                 aes(
                   x     = f.ordered1,
                   y     = 2**value, 
                   fill  = variable,
                   color = variable)) + # ylim(0,100)+
      scale_y_continuous(trans = log2_trans(), 
                         #labels = label_comma(),
                         breaks = c(0.125,1,8,64,512),
                         labels = c("0.125","1","8","64","512")) +
      #stat_n_text(size = 5, fontface = "bold", hjust = 0) +
      geom_boxplot(
        alpha    = .01,
        #size     = .75,
        #width    = 1,
        position = position_dodge(width = .9)
      ) + 
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_color_manual(values = c("#0072B2","#CC79A7","#009E73","#D55E00")) + #for color-blind palettes
      labs(x = "primary disease",
           y = "ratio of RNA counts") +
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
        strip.text           = black_bold_tahoma_12
      )
    print(p1)
    ggsave(
      path        = ".",
      filename    = "EIFratio.pdf", 
      plot        = p1,
      width       = 9, 
      height      = 9, 
      useDingbats = FALSE)
    
    pancancer.TCGA.EIF.ratio.long2 <- pancancer.TCGA.EIF.ratio.long[
      pancancer.TCGA.EIF.ratio.long$variable %in% c("EIF4A1:\nEIF4E","EIF4G1:\nEIF4E","EIF4A1:\nEIF4EBP1","EIF4G1:\nEIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1","EIF4A1:\nEIF4E+EIF4EBP1"), ]
    pancancer.TCGA.EIF.ratio.long2$variable = factor(pancancer.TCGA.EIF.ratio.long2$variable, levels = c("EIF4G1:\nEIF4E+EIF4EBP1","EIF4A1:\nEIF4E+EIF4EBP1","EIF4G1:\nEIF4E","EIF4A1:\nEIF4E","EIF4G1:\nEIF4EBP1","EIF4A1:\nEIF4EBP1"))
    f2 <- factor(pancancer.TCGA.EIF.ratio.long2$primary.disease)
    f.ordered2 <- fct_rev(f2)
    p2 <- ggplot(data = pancancer.TCGA.EIF.ratio.long2,
                 aes(
                   x     = f.ordered2,
                   y     = 2**value, 
                   fill  = variable,
                   color = variable)) + #ylim(0, 101)+
      scale_y_continuous(trans = log2_trans(), labels = label_number_auto()) +
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
      path        = ".",
      filename    = "EIFsumratio.pdf", 
      plot        = p2,
      width       = 12.5, 
      height      = 9, 
      useDingbats = FALSE)
  }
  make.plot()
}
