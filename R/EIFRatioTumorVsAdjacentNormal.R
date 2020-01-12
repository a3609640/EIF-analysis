############################################################
## violin plot for EIF ratio in tumors vs adjacent normal ##
############################################################
plot.violingraph.EIF.ratio.TCGA <- function (EIF.gene) {
  tissue.GTEX.TCGA.gene <- function(){
    TCGA.GTEX.anno <- TcgaTargetGtexPhenoType
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
    EIF.TCGA.GTEX.score$`EIF4A1:EIF4G1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$EIF4G1)
    EIF.TCGA.GTEX.score$`EIF4G1:EIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4EBP1:EIF4E` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4EBP1 - EIF.TCGA.RNAseq.anno.subset$EIF4E)
    EIF.TCGA.GTEX.score$`EIF4G1:EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$EIF4EBP1)
    EIF.TCGA.GTEX.score$`EIF4G1:EIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4G1 - EIF.TCGA.RNAseq.anno.subset$sum)
    EIF.TCGA.GTEX.score$`EIF4A1:EIF4E+EIF4EBP1` <- 
      (EIF.TCGA.RNAseq.anno.subset$EIF4A1 - EIF.TCGA.RNAseq.anno.subset$sum)
    ratio <- c("EIF4A1:EIF4E", "EIF4G1:EIF4E", "EIF4EBP1:EIF4E", "EIF4A1:EIF4G1",
               "EIF4G1:EIF4EBP1","EIF4G1:EIF4E+EIF4EBP1","EIF4A1:EIF4E+EIF4EBP1")
    EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX.score[c(ratio, "sample.type","primary.site")]
    EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score)
    return(EIF.TCGA.GTEX.score.long)
  }
  EIF.TCGA.GTEX.score.long <- get.EIFratio.anno.data()
  make.plot <- function () {
    new.label <- c(`EIF4A1:EIF4E`    = "EIF4A1:\nEIF4E",
                   `EIF4G1:EIF4E`    = "EIF4G1:\nEIF4E",
                   `EIF4EBP1:EIF4E`  = "EIF4EBP1:\n EIF4E",
                   `EIF4A1:EIF4G1`   = "EIF4A1:\nEIF4G1",
                   `EIF4G1:EIF4EBP1` = "EIF4G1:\n EIF4EBP1",
                   `EIF4G1:EIF4E+EIF4EBP1` = "EIF4G1:\n EIF4E+EIF4EBP1",
                   `EIF4A1:EIF4E+EIF4EBP1` = "EIF4A1:\n EIF4E+EIF4EBP1")
    EIF.TCGA.GTEX.score.long1 <- EIF.TCGA.GTEX.score.long[
      EIF.TCGA.GTEX.score.long$variable %in% c("EIF4A1:EIF4E","EIF4G1:EIF4E","EIF4EBP1:EIF4E","EIF4A1:EIF4G1"), ]
    EIF.TCGA.GTEX.score.long1$variable = factor(EIF.TCGA.GTEX.score.long1$variable, levels = c("EIF4A1:EIF4E","EIF4G1:EIF4E","EIF4EBP1:EIF4E","EIF4A1:EIF4G1"))
    p1 <- ggplot(data = EIF.TCGA.GTEX.score.long1,
                 aes(
                   x     = sample.type,
                   y     = 2**value, 
                   color = sample.type,
                   fill  = sample.type)) +
      stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
      facet_grid(~ variable, 
                 scales = "free",
                 space  = "free") +
      facet_wrap(~ variable, 
                 labeller = labeller(variable = as_labeller(new.label)),
                 ncol = 6) +
      geom_violin(trim = FALSE) +
      geom_boxplot(
        alpha    = 0.01,
        width    = 0.25,
        color    = "black",
        outlier.colour=NA
      ) + 
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
      stat_compare_means(comparisons = list(
        c("Metastatic", "Solid Tissue Normal"),
        c("Primary Tumor", "Solid Tissue Normal"),
        c("Metastatic", "Primary Tumor")
      ),
      method = "t.test", label = "p.signif", size = 6, hjust = 0)
    print(p1)
    ggsave(
      path        = ".",
      filename    = "EIFratioviolin.pdf", 
      plot        = p1,
      width       = 6, 
      height      = 9, 
      useDingbats = FALSE)
    
    EIF.TCGA.GTEX.score.long2 <- EIF.TCGA.GTEX.score.long[
      EIF.TCGA.GTEX.score.long$variable %in% c("EIF4G1:EIF4E+EIF4EBP1",
                                               "EIF4A1:EIF4E+EIF4EBP1"), ]
    EIF.TCGA.GTEX.score.long2$variable = factor(EIF.TCGA.GTEX.score.long2$variable, levels = c("EIF4G1:EIF4E+EIF4EBP1","EIF4A1:EIF4E+EIF4EBP1"))
    p2 <- ggplot(data = EIF.TCGA.GTEX.score.long2,
                 aes(
                   x     = sample.type,
                   y     = 2**value, 
                   color = sample.type,
                   fill  = sample.type)) +
      stat_n_text(size = 6, fontface = "bold", angle = 90, hjust = 0) +
      facet_grid(~ variable, 
                 scales = "free",
                 space  = "free") +
      facet_wrap(~ variable, 
                 labeller = labeller(variable = as_labeller(new.label)),
                 ncol = 6) +
      geom_violin(trim = FALSE) +
      geom_boxplot(
        alpha    = 0.01,
        width    = 0.25,
        color    = "black",
        outlier.colour=NA
      ) + 
      scale_color_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
      scale_fill_manual(values = c( "#56B4E9", "#009E73", "#D55E00")) + #for color-blind palettes
      labs(x = "sample type",
           y = "log2(ratio)") +
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
      stat_compare_means(comparisons = list(
        c("Metastatic", "Solid Tissue Normal"),
        c("Primary Tumor", "Solid Tissue Normal"),
        c("Metastatic", "Primary Tumor")
      ),
      method = "t.test", label = "p.signif", size = 6, hjust = 0)
    print(p2)
    ggsave(
      path        = ".",
      filename    = "EIFsumratioviolin.pdf", 
      plot        = p2,
      width       = 3.2, 
      height      = 9, 
      useDingbats = FALSE) 
    
    
    
    p3 <- ggplot(data = EIF.TCGA.GTEX.score.long,
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
        width    = 0.25,
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
        strip.text      = black_bold_tahoma_16
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
    print(p3)
  }
  make.plot()
}

