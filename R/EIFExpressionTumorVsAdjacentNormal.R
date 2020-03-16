#################################################################
## violin plot for EIF expression in tumors vs adjacent normal ##
#################################################################
plot.violingraph.EIF.RNAseq.TCGA <- function (EIF.gene) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- TcgaTargetGtexPhenoType
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
      EIF.TCGA.RNAseq.anno.subset.long$variable %in% c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1"), ]
    EIF.TCGA.RNAseq.anno.subset.long1$variable <- factor(EIF.TCGA.RNAseq.anno.subset.long1$variable, levels = c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1"))
    p1 <- ggplot(data = EIF.TCGA.RNAseq.anno.subset.long1,
                 aes(
                   x     = sample.type,
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
    ggsave(path        = ".",
           filename    = "EIFexpressionviolin.pdf", 
           plot        = p1,
           width       = 6, 
           height      = 9, 
           useDingbats = FALSE)
    
    EIF.TCGA.RNAseq.anno.subset.long2 <- EIF.TCGA.RNAseq.anno.subset.long[
      EIF.TCGA.RNAseq.anno.subset.long$variable %in% c("EIF4G1","EIF4E+EIF4EBP1"), ]
    EIF.TCGA.RNAseq.anno.subset.long2$variable <- factor(EIF.TCGA.RNAseq.anno.subset.long2$variable, levels = c("EIF4G1","EIF4E+EIF4EBP1"))
    p2 <- ggplot(data = EIF.TCGA.RNAseq.anno.subset.long2,
                 aes(
                   x     = sample.type,
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
        method = "t.test", label = "p.signif", size = 6)
    print(p2)
    ggsave(path        = ".",
           filename    = "EIF4G1sumviolin.pdf", 
           plot        = p2,
           width       = 3.2, 
           height      = 9, 
           useDingbats = FALSE)
    
    p3 <- ggplot(data = EIF.TCGA.RNAseq.anno.subset.long,
                 aes(
                   x     = variable,
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
      labs(
        x = "EIF4F complex",
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

