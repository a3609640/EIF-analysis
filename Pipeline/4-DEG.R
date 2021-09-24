# prepare RNA-seq related dataset from TCGA and GTEx----------------------------
get.TCGA.GTEX.RNAseq <- function() {
  TCGA.pancancer <- data.table::fread(
    file.path(data.file.directory, 
              "TcgaTargetGtex_RSEM_Hugo_norm_count"),
    data.table = FALSE
  )  %>% 
    as.data.frame(.) %>% 
    distinct(., sample, .keep_all = TRUE) %>% 
    na.omit(.) %>%
    remove_rownames(.) %>%
    column_to_rownames(var = 'sample')
  
  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and colnames in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return (TCGA.pancancer_transpose)
}
TCGA.GTEX.RNAseq <- get.TCGA.GTEX.RNAseq()   

TCGA.GTEX.sampletype <- readr::read_tsv(
  file.path(data.file.directory, 
            "TcgaTargetGTEX_phenotype.txt")) %>% {
              as.data.frame(.) %>% 
                distinct(., sample, .keep_all = TRUE) %>% 
                remove_rownames() %>%
                column_to_rownames(var = 'sample') %>%
                select("_sample_type",
                       "primary disease or tissue",
                       "_primary_site",
                       "_study") %>%
                rename("sample.type" = "_sample_type", 
                       "primary.disease" = "primary disease or tissue",
                       "primary.site" = "_primary_site",
                       "study" = "_study")}

TCGA.GTEX.RNAseq.sampletype <- merge(TCGA.GTEX.RNAseq,
                                     TCGA.GTEX.sampletype,
                                     by    = "row.names",
                                     all.x = TRUE) %>% {
                                       remove_rownames(.) %>%
                                         column_to_rownames(var = 'Row.names')}


# Differential expression analysis and plotting --------------------------------
RNAseq.all.gene <- function (df){
  df <- df %>% 
    filter(study == "TCGA" & sample.type != "Solid Tissue Normal") 
  #order the EIF genes according to the order of the expression means
  order <- df %>%
    #filter out category equal to 'Lung Adenocarcinoma'
    filter(primary.disease == 'Lung Adenocarcinoma') %>%
    #use the same groups as in the ggplot
    group_by(variable) %>%
    #calculate the means
    summarise(mean_RNAseq = mean(RNAseq)) %>%
    mutate(variable = fct_reorder(variable, mean_RNAseq)) %>%
    .$variable  %>%
    levels()
  
  df <- df %>%
    mutate(variable = factor(variable, levels = rev(order)))
}
RNAseq.grouped.boxplot <- function(df) {
  p1 <- ggplot(data = df,
               aes(
                 x = primary.disease,
                 y = 2**RNAseq)
  ) +
    scale_y_continuous(
      trans = log2_trans(),
      #limits = c(2**4, 2**17),
      labels = label_comma()
    ) +
    stat_n_text(
      size = 5,
      fontface = "bold",
      angle = 90,
      hjust = 0
    ) +
    geom_boxplot(aes(
      colour = factor(variable),
      #fill = factor(variable)
    ),
    outlier.shape = NA,
    position = position_dodge(width = 1)
    ) +
    labs(
      x = "primary disease",
      y = paste("Normalized expression (RNA-Seq counts)")
      ) +
    theme_bw() +
    theme(
      plot.title = black_bold_12(),
      axis.title.x = element_blank(),
      axis.title.y = black_bold_12(),
      axis.text.x = black_bold_12_45(),
      axis.text.y = black_bold_12(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_12(),
      legend.position = c(0, 0.95),
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


RNAseq.ind.gene <- function (df, x) {
  df <- df %>% 
    filter(study == "TCGA") %>% 
    droplevels()  %>% 
    mutate(sample.type = case_when(sample.type != "Solid Tissue Normal" ~ "Tumor", 
                                   sample.type == "Solid Tissue Normal" ~ "Normal")) %>%
    filter(variable == x) %>% 
    mutate(primary.disease = forcats::fct_rev(primary.disease))
  output <- list(df, x) 
  return(output)
}
RNAseq.boxplot <- function(df) {
  p1 <- ggplot(
    data = df[[1]],
    aes(
      x = primary.disease,
      y = 2**RNAseq, 
      color = sample.type)
  ) +
    scale_y_continuous(
      trans = log2_trans(),
      #limits = c(2**11, 2**17),# for 4g
      #limits = c(2**7, 2**14),# for eif4E
      labels = label_comma()
    ) +
    stat_n_text(
      size = 5,
      fontface = "bold",
      hjust = 0
    ) +
    geom_boxplot(
      #alpha = .1,
      #fill = sample.type,
      outlier.shape = NA,
      #size = .75,
      #width = 1,
      position = "dodge"
    ) +
    scale_color_manual(values = c("Tumor" = "#CC79A7", 
                                  "Normal" = "#0072B2"),
                       breaks = c("Tumor", "Normal"),
                       labels = c("Tumor\n", "Normal\n")
    ) +
    labs(
      x = "primary disease",
      y = paste(df[[2]], "expression (RNA-Seq counts)")
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title = black_bold_12(),
      axis.title.x = black_bold_12(),
      axis.title.y = element_blank(),
      axis.text.x = black_bold_12(),
      axis.text.y = black_bold_12(),
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
    filename = paste0(df[[2]], "tumorvsnormal.pdf"),
    plot = p1,
    width = 7.5,
    height = 9,
    useDingbats = FALSE
  )
}


RNAseq.tumortype <- function (df) {
  TCGA.GTEX.RNAseq.sampletype.subset <- df %>%
    filter(study == "TCGA") %>% 
    mutate_if(is.character, as.factor) %>% 
    filter(sample.type %in% c(
      "Metastatic",
      "Primary Tumor",
      "Solid Tissue Normal"))
}
RNAratio.tumortype <- function(df, x){
  RNAratio.data <- df %>% 
    filter(sample.type %in% c("Metastatic",
                              "Primary Tumor",
                              "Solid Tissue Normal")) %>% 
    droplevels() %>% 
    select(all_of(x),
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>% 
    melt(.,id = c("sample.type",
                  "primary.disease",
                  "primary.site",
                  "study"), 
         value.name = "RNAseq") %>%
    mutate_if(is.character, as.factor)  %>%
    mutate(primary.disease = forcats::fct_rev(primary.disease))
}
violinplot <- function(df) {
  p1 <- ggplot(
    data = df,
    #data = EIF.TCGA.RNAseq.anno.subset.long,
    aes(
      x = sample.type,
      y = 2**RNAseq,
      color = sample.type,
      fill = sample.type
    )
  ) +
    stat_n_text(size = 6, 
                fontface = "bold", 
                angle = 90, 
                hjust = 0) +
    ggplot2::facet_grid(. ~ variable,
                        scales = "free",
                        space  = "free"
    ) +
    # facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = TRUE) +
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


RNAratio.EIF.gene <- function(x){
  TCGA.RNAratio.sampletype.subset <- TCGA.GTEX.RNAseq.sampletype %>%
    select(all_of(x),       
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>% 
    as.data.frame(.) %>%      
    filter(EIF4E != 0 & !is.na(primary.site)) %>% 
    mutate(sum = log2(2**EIF4E + 2**EIF4EBP1 - 1),
           #minus = log2(2**EIF4E - 2**EIF4EBP1 + 1),
           `EIF4A1:\nEIF4E` = EIF4A1 - EIF4E,
           `EIF4A1:\nEIF4E2` = EIF4A1 - EIF4E2,
           `EIF4A2:\nEIF4E` = EIF4A2 - EIF4E,
           `EIF4A2:\nEIF4E2` = EIF4A2 - EIF4E2,
           `EIF4G1:\nEIF4E` = EIF4G1 - EIF4E,
           `EIF4G1:\nEIF4E2` = EIF4G1 - EIF4E2,
           `EIF4G3:\nEIF4E` = EIF4G3 - EIF4E,
           `EIF4G3:\nEIF4E2` = EIF4G3 - EIF4E2,
           `EIF4A1:\nEIF4G1` = EIF4A1 - EIF4G1,
           `EIF4A2:\nEIF4G1` = EIF4A2 - EIF4G1,
           `EIF4A1:\nEIF4G2` = EIF4A1 - EIF4G2,
           `EIF4A2:\nEIF4G2` = EIF4A2 - EIF4G2,
           `EIF4E:\nEIF4EBP1` = EIF4E - EIF4EBP1,
           `EIF4E2:\nEIF4E` = EIF4E2 - EIF4E,
           `EIF4G2:\nEIF4G1` = EIF4G2 - EIF4G1,
           `EIF4G1:\nEIF4G3` = EIF4G1 - EIF4G3,
           `EIF4A1:\nEIF4A2` = EIF4A1 - EIF4A2,
           `EIF4G1:\nEIF4E+EIF4EBP1` = EIF4G1 - sum,
           `EIF4A1:\nEIF4E+EIF4EBP1` = EIF4G1 - sum)  %>% 
    select("EIF4A1:\nEIF4E", "EIF4A1:\nEIF4E2", 
           "EIF4A2:\nEIF4E", "EIF4A2:\nEIF4E2",          
           "EIF4G1:\nEIF4E", "EIF4G1:\nEIF4E2", 
           "EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2", 
           "EIF4A1:\nEIF4G1", "EIF4A2:\nEIF4G1", 
           "EIF4A1:\nEIF4G2", "EIF4A2:\nEIF4G2",
           "EIF4E:\nEIF4EBP1", "EIF4E2:\nEIF4E", 
           "EIF4G2:\nEIF4G1", "EIF4G1:\nEIF4G3", 
           "EIF4A1:\nEIF4A2",          
           "EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1",
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>%        
    mutate_if(is.character, as.factor)%>% 
    na.omit(.)
}
RNAratio.selected <- function(df, x){
  RNAratio.data <- df %>% 
    filter(study == "TCGA") %>% 
    droplevels() %>% 
    mutate(sample.type = case_when(sample.type != "Solid Tissue Normal" ~ "Tumor", 
                                   sample.type == "Solid Tissue Normal" ~ "NAT")) %>% 
    select(all_of(x),
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>% 
    melt(.,id = c("sample.type",
                  "primary.disease",
                  "primary.site",
                  "study"), 
         value.name = "RNAseq") %>%
    mutate_if(is.character, as.factor)  %>%
    mutate(primary.disease = forcats::fct_rev(primary.disease))
}
RNAratio.boxplot <- function(df, dashline, ylimit, filename) {        
  p1 <- ggplot(
    data = df,
    aes(
      x = primary.disease,
      #x = f.ordered1,
      y = 2**RNAseq,
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
      #data = df,
      aes(yintercept = dashline),
      linetype = "dashed"
    ) +
    scale_color_manual(values = c("Tumor" = "#CC79A7", "NAT" = "#0072B2"),
                       breaks = c("Tumor", "NAT")) +
    # for color-blind palettes
    facet_wrap(~variable, scales = "free_x") +
    ggplot2::facet_grid(~variable, scales = "free_x", space = "free") +
    guides(colour = guide_legend(nrow = 1)) +
    labs(x = "primary disease",
         y = "Ratio of RNA counts") +
    coord_flip(ylim = ylimit) +
    theme_bw() +
    theme(
      plot.title = black_bold_12(),
      axis.title.x = black_bold_12(),
      axis.title.y = element_blank(),
      axis.text.x = black_bold_12(),
      axis.text.y = black_bold_12(),
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
    filename = filename,
    plot = p1,
    width = 18,
    height = 8,
    useDingbats = FALSE)}


# master functions to call DEG gene analysis and plotting ----------------------
plot.boxgraph.RNAseq.TCGA <- function(EIF) {
  TCGA.GTEX.RNAseq.sampletype.subset <- TCGA.GTEX.RNAseq.sampletype %>%
    select(all_of(EIF),       
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>% 
    as.data.frame(.) %>% 
    melt(.,id = c("sample.type",
                  "primary.disease", 
                  "primary.site", 
                  "study"), 
         value.name = "RNAseq") %>% 
    na.omit(.$primary.site) %>%
    #filter(RNAseq != 0) %>% 
    mutate_if(is.character, as.factor)

  # boxplot to compare relative abundance of genes across tumors
  RNAseq.all.gene(TCGA.GTEX.RNAseq.sampletype.subset) %>% 
    RNAseq.grouped.boxplot ()
  
  # boxplot to compare RNA-seq of one gene in tumor vs adjacent normal
  RNAseq.ind.gene.df <- lapply(EIF, 
                               RNAseq.ind.gene, 
                               df = TCGA.GTEX.RNAseq.sampletype.subset)
  lapply(RNAseq.ind.gene.df, RNAseq.boxplot)
  
  # violin plot to compare  expression in primary, metastatic tumors vs NATs
  RNAseq.tumortype (TCGA.GTEX.RNAseq.sampletype.subset) %>% violinplot()
}

plot.RNAratio.TCGA <- function(EIF) {
  RNAratio.data <- RNAratio.EIF.gene(EIF)
  
  RNAratio.selected(RNAratio.data, c("EIF4G1:\nEIF4E", "EIF4A1:\nEIF4E",
                                     "EIF4A2:\nEIF4E", "EIF4G3:\nEIF4E", 
                                     "EIF4G3:\nEIF4E2", "EIF4G1:\nEIF4G3")) %>%
    RNAratio.boxplot(df = ., 
                     dashline = 1, 
                     ylimit = c(0,25),
                     filename = "RNAratio1.pdf")
  
  RNAratio.selected(RNAratio.data, c("EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", 
                                     "EIF4A1:\nEIF4A2", "EIF4E:\nEIF4EBP1",
                                     "EIF4G1:\nEIF4E+EIF4EBP1", 
                                     "EIF4A1:\nEIF4E+EIF4EBP1")) %>%
    RNAratio.boxplot(df = ., 
                     dashline = 4, 
                     ylimit = c(0,25),
                     filename = "RNAratio2.pdf") 
  
  
  RNAratio.selected(RNAratio.data, c("EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2",
                                     "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", 
                                     "EIF4A1:\nEIF4A2", "EIF4E:\nEIF4EBP1")) %>%
    RNAratio.boxplot(df = ., 
                     dashline = 1, 
                     ylimit = c(0,5),
                     filename = "RNAratio3.pdf") 
  
  RNAratio.tumortype (RNAratio.data, c("EIF4G1:\nEIF4E", "EIF4A1:\nEIF4E", 
                                       "EIF4A2:\nEIF4E", "EIF4G3:\nEIF4E",
                                       "EIF4G3:\nEIF4E2", "EIF4G1:\nEIF4G3",
                                       "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", 
                                       "EIF4A1:\nEIF4A2", "EIF4E:\nEIF4EBP1",
                                       "EIF4G1:\nEIF4E+EIF4EBP1", 
                                       "EIF4A1:\nEIF4E+EIF4EBP1")) %>% 
    violinplot()
}

## this function was not reported in the paper.
plot.cormatrix.RNAseq <- function (EIF) {
  TCGA.GTEX.RNAseq.sampletype.subset <- TCGA.GTEX.RNAseq.sampletype %>%
    select(all_of(EIF),       
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>% 
    as.data.frame(.)     
  
  plot.cor <- function(x){
    EIF.TCGA.RNAseq.subset <- TCGA.GTEX.RNAseq.sampletype.subset %>% 
      filter(sample.type != "Solid Tissue Normal" & study == x) %>%
      select(all_of(EIF))
    
    M = cor(EIF.TCGA.RNAseq.subset)
    testRes = corrplot::cor.mtest(EIF.TCGA.RNAseq.subset, conf.level = 0.95)
    
    #pdf(file.path(output.directory, "CNV", "EIFCNVcormatrix.pdf"),
    #    width = 9,
    #    height = 9,
    #    useDingbats = FALSE
    #)
    corrplot(
      M,
      method      = "color",
      title = paste("correlation (", x, ")"),
      tl.cex      = 1,
      number.cex  = 1,
      addgrid.col = "gray",
      addCoef.col = "black",
      tl.col      = "black",
      #type        = "lower",
      order       = "FPC", 
      #order = 'hclust', 
      tl.srt = 45,
      p.mat       = testRes$p,
      sig.level   = 0.05, # insig = "blank"
    )
  }
  plot.cor("GTEX")
  plot.cor("TCGA")
}



# Run master functions ---------------------------------------------------------
plot.boxgraph.RNAseq.TCGA(c("EIF4G1","EIF4G2","EIF4G3","PABPC1",
                            "EIF4A1","EIF4A2","EIF4B","EIF4H",
                            "EIF4E","EIF4E2","EIF4E3",
                            "EIF4EBP1","EIF3D"))

plot.RNAratio.TCGA(c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
                     "EIF4G1","EIF4G2","EIF4G3","EIF3D",#
                     "EIF4A1","EIF4A2"))

plot.cormatrix.RNAseq(c("EIF4G1", "EIF4G2","EIF4G3",
                        "EIF4A1","EIF4A2", 
                        "EIF4E", "EIF4E2", "EIF4E3", 
                        "EIF4EBP1", "EIF4EBP2","MTOR",
                        "EIF3C","EIF3D","EIF3E","PABPC1", 
                        "MKNK1", "MKNK2",
                        "TP53","MYC"))

