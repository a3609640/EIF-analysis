## prepare RNA-seq related dataset
TCGA.GTEX.RNAseq <- function() {
  TCGA.pancancer <- fread(
    file.path(data.file.directory, 
              "TcgaTargetGtex_RSEM_Hugo_norm_count"),
    data.table = FALSE
  )  %>% 
    as.data.frame(.) %>% 
    distinct(., sample, .keep_all = TRUE) %>% 
    na.omit(.) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'sample')
  
  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and colnames in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return (TCGA.pancancer_transpose)
}
TCGA.GTEX.RNAseq <- TCGA.GTEX.RNAseq()   

TCGA.GTEX.sampletype <- readr::read_tsv(
  file.path(data.file.directory, 
            "TcgaTargetGTEX_phenotype.txt")) %>% 
  as.data.frame(.) %>% 
  distinct(., sample, .keep_all = TRUE) %>% 
  #na.omit(.) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'sample') %>%
  select("_sample_type",
         "primary disease or tissue",
         "_primary_site",
         "_study") %>%
  rename("sample.type" = "_sample_type", 
         "primary.disease" = "primary disease or tissue",
         "primary.site" = "_primary_site",
         "study" = "_study")

TCGA.GTEX.RNAseq.sampletype <- merge(TCGA.GTEX.RNAseq,
                                     TCGA.GTEX.sampletype,
                                     by    = "row.names",
                                     all.x = TRUE) %>% 
  remove_rownames() %>%
  column_to_rownames(var = 'Row.names')

###########################################################
## boxplot for RNA-seq in TCGA tumors vs adjacent normal ##
###########################################################
plot.boxgraph.RNAseq.TCGA <- function(EIF) {
  ## 
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
    filter(RNAseq != 0) %>% 
    mutate_if(is.character, as.factor)

  # boxplot to compare RNA-seq of one gene in tumor vs adjacent normal
  make.box.plot.1 <- function(x) {
    df <- TCGA.GTEX.RNAseq.sampletype.subset %>% 
      filter(study == "TCGA") %>% 
      droplevels()  %>% 
      mutate(sample.type = case_when(sample.type != "Solid Tissue Normal" ~ "Tumor", 
                                     sample.type == "Solid Tissue Normal" ~ "Normal")) %>%
      filter(variable == x) %>% 
      mutate(primary.disease = forcats::fct_rev(primary.disease))
    
    p1 <- ggplot(
      data = df,
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
        y = paste(x, "expression (RNA-Seq counts)")
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
      filename = paste0(x, "tumorvsnormal.pdf"),
      plot = p1,
      width = 7.5,
      height = 9,
      useDingbats = FALSE
    )
  }
  lapply(EIF, make.box.plot.1)
  
  # boxplot to compare relative abundance of genes across tumors
  make.box.plot.2 <- function() {
    df <- TCGA.GTEX.RNAseq.sampletype.subset %>% 
      filter(study == "TCGA" & sample.type != "Solid Tissue Normal") #%>% 

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
      #scale_color_brewer(palette="Dark2") +
      #scale_fill_brewer(palette = "Set2") +
      #scale_color_manual(
      #  values = c(
      #    "EIF4EBP1" = "#B997C7",
      #    "EIF4E3" = "#824D99",
      #    "EIF4E2" = "#4E78C4",
      #    "EIF4E" = "#57A2AC", 
      #    "EIF4G3" = "#7EB875",
      #    "EIF4G2" = "#D0B541",
      #    "EIF4G1" = "#E67F33",
      #    "EIF4A2" = "#CE2220", 
      #    "EIF4A1" = "#521A13"
      #  ),
      #  breaks = rev(c("EIF4EBP1","EIF4E3","EIF4E2","EIF4E",
      #                 "EIF4G3", "EIF4G2","EIF4G1","EIF4A2","EIF4A1"))
      #) +
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
  make.box.plot.2()
}

##################################################################
## violin plot for EIF expression in tumors vs adjacent normals ##
##################################################################
plot.violingraph.RNAseq.TCGA <- function(EIF) {
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
    filter(RNAseq != 0 & study == "TCGA") %>% 
    mutate_if(is.character, as.factor) %>% 
    filter(sample.type %in% c(
      "Metastatic",
      "Primary Tumor",
      "Solid Tissue Normal"
    ))
  
  make.plot <- function() {
    #EIF.TCGA.RNAseq.anno.subset.long1 <- EIF.TCGA.RNAseq.anno.subset.long[
    #  EIF.TCGA.RNAseq.anno.subset.long$variable %in% c(EIF),]
    #EIF.TCGA.RNAseq.anno.subset.long1$variable <- factor(
    #  EIF.TCGA.RNAseq.anno.subset.long1$variable,
    #  levels = c("EIF4G1", "EIF4G2","EIF4G3","EIF4A1","EIF4A2", "PABPC1","PAIP1","PAIP2",
    #             "EIF4E", "EIF4E2", "EIF4E3", "EIF4EBP1" ,"EIF4EBP2", "EIF3D")
    #)
    p1 <- ggplot(
      data = TCGA.GTEX.RNAseq.sampletype.subset,
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
  make.plot()
}


########################################################
## boxplot and violin plots for RNA ratios across tumors ##
########################################################
plot.RNAratio.TCGA <- function(EIF) {
  ## 
  TCGA.RNAratio.sampletype.subset <- TCGA.GTEX.RNAseq.sampletype %>%
    select(all_of(EIF),       
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
           `EIF4A1:\nEIF4E+EIF4EBP1` = EIF4G1 - sum#,
           #`EIF4G1:\nEIF4E-EIF4EBP1` = EIF4G1 - minus,
           #`EIF4A1:\nEIF4E-EIF4EBP1` = EIF4A1 - minus
    )  %>% 
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
           #"EIF4G1:\nEIF4E-EIF4EBP1", "EIF4A1:\nEIF4E-EIF4EBP1",
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>%        
    mutate_if(is.character, as.factor)%>% 
    na.omit(.)
  
  make.boxplot <- function() {
    TCGA.RNAratio.subset <- TCGA.RNAratio.sampletype.subset %>% 
      filter(study == "TCGA") %>% 
      droplevels() %>% 
      mutate(sample.type = case_when(sample.type != "Solid Tissue Normal" ~ "Tumor", 
                                     sample.type == "Solid Tissue Normal" ~ "NAT"))
    
    p1 <- TCGA.RNAratio.subset %>%
      select("EIF4G1:\nEIF4E","EIF4A1:\nEIF4E","EIF4A2:\nEIF4E",
             "EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2" ,"EIF4G1:\nEIF4G3",
             "sample.type",
             "primary.disease",
             "primary.site",
             "study") %>% 
      melt(.,id = c("sample.type",
                    "primary.disease",
                    "primary.site",
                    "study")) %>%
      mutate_if(is.character, as.factor)  %>%
      mutate(primary.disease = forcats::fct_rev(primary.disease)) %>%
      # wrap the plotting code in {...}
      {ggplot(
        data = .,
        aes(
          x = primary.disease,
          #x = f.ordered1,
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
        scale_color_manual(values = c("Tumor" = "#CC79A7", "NAT" = "#0072B2"),
                           breaks = c("Tumor", "NAT")) +
        # for color-blind palettes
        facet_wrap(~variable, scales = "free_x") +
        ggplot2::facet_grid(~variable, scales = "free_x", space = "free") +
        guides(colour = guide_legend(nrow = 1)) +
        labs(x = "primary disease",
             y = "Ratio of RNA counts") +
        coord_flip(ylim = c(0,25)) +
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
          )}
    print(p1)
    ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorratio.pdf",
        plot = p1,
        width = 18,
        height = 8,
        useDingbats = FALSE)
    
    p2 <- TCGA.RNAratio.subset %>%
      select("EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", 
             "EIF4A1:\nEIF4A2", "EIF4E:\nEIF4EBP1",
             "EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1",
             "sample.type",
             "primary.disease",
             "primary.site",
             "study") %>% 
      melt(.,id = c("sample.type",
                    "primary.disease",
                    "primary.site",
                    "study")) %>%
      mutate_if(is.character, as.factor)  %>%
      mutate(primary.disease = forcats::fct_rev(primary.disease)) %>%
      {ggplot(
      data = .,
      aes(
          x = primary.disease,
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
        data = .,
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
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.box = "horizontal",
        legend.text = black_bold_12(),
        strip.text.x = black_bold_12()
        )}
      print(p2)
      
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorratio2.pdf",
        plot = p2,
        width = 18,
        height = 8,
        useDingbats = FALSE)
      
      p3 <- TCGA.RNAratio.subset %>%
        select("EIF4G3:\nEIF4E", "EIF4G3:\nEIF4E2",
               "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", 
               "EIF4A1:\nEIF4A2", "EIF4E:\nEIF4EBP1",
               "sample.type",
               "primary.disease",
               "primary.site",
               "study") %>% 
        melt(.,id = c("sample.type",
                      "primary.disease",
                      "primary.site",
                      "study")) %>%
        mutate_if(is.character, as.factor) %>%
        mutate(primary.disease = forcats::fct_rev(primary.disease)) %>%
        {ggplot(
        data = .,
        aes(
          x = primary.disease,
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
          data = .,
          aes(yintercept = 1),
          linetype = "dashed"
        ) +
        # geom_text(data = pancancer.TCGA.EIF.ratio.long1, aes(f.ordered1, hline, label = hline), vjust = 2, hjust = -0.2) +
        scale_color_manual(values = c("Tumor" = "#CC79A7", 
                                      "NAT" = "#0072B2"),
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
        )}
      print(p3)
      ggplot2::ggsave(
        path = file.path(output.directory, "Expression"),
        filename = "tumorratio3.pdf",
        plot = p3,
        width = 18,
        height = 8,
        useDingbats = FALSE)
    }
  make.boxplot()
  
  make.violinplot <- function() {
    p1 <- TCGA.RNAratio.sampletype.subset %>% 
      filter(sample.type %in% c("Metastatic",
                                "Primary Tumor",
                                "Solid Tissue Normal")) %>% 
      droplevels() %>% 
      select("EIF4G1:\nEIF4E","EIF4A1:\nEIF4E", "EIF4A2:\nEIF4E", 
             "EIF4G3:\nEIF4E","EIF4G3:\nEIF4E2","EIF4G1:\nEIF4G3",
             "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E", "EIF4A1:\nEIF4A2",
             "EIF4E:\nEIF4EBP1","EIF4G1:\nEIF4E+EIF4EBP1", "EIF4A1:\nEIF4E+EIF4EBP1",
             "sample.type",
             "primary.disease",
             "primary.site",
             "study") %>% 
      melt(.,id = c("sample.type",
                    "primary.disease",
                    "primary.site",
                    "study")) %>%
      mutate_if(is.character, as.factor)  %>%
      mutate(primary.disease = forcats::fct_rev(primary.disease)) %>%
      ggplot(
        data = .,
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
  make.violinplot()
}


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
    
    cor_5 <- rcorr(as.matrix(EIF.TCGA.RNAseq.subset), type = "pearson")
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
      p.mat       = p_mat,
      sig.level   = 0.05, # insig = "blank"
    )
  }
  plot.cor("GTEX")
  plot.cor("TCGA")
}



plot.boxgraph.RNAseq.TCGA(c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
                            "EIF4G1","EIF4G2","EIF4G3","PABPC1",
                            #"MKNK1", "MKNK2",
                            "EIF4A1","EIF4A2"#,
                            #"EIF3D","TP53","MYC"
                            ))

plot.violingraph.RNAseq.TCGA(c("EIF4G1", "EIF4G2","EIF4G3",
                               "EIF4A1","EIF4A2", 
                               "EIF4E", "EIF4E2", "EIF4E3", 
                               "EIF4EBP1", "EIF3D",
                               "PABPC1", "MKNK1", "MKNK2"))

plot.RNAratio.TCGA(c("EIF4E","EIF4E2","EIF4E3","EIF4EBP1",
                     "EIF4G1","EIF4G2","EIF4G3","EIF3D",#"PABPC1",
                     "EIF4A1","EIF4A2"))

plot.cormatrix.RNAseq(c("EIF4G1", "EIF4G2","EIF4G3",
                        "EIF4A1","EIF4A2", 
                        "EIF4E", "EIF4E2", "EIF4E3", 
                        "EIF4EBP1", "EIF4EBP2","MTOR",
                        "EIF3C","EIF3D","EIF3E","PABPC1", 
                        "MKNK1", "MKNK2",
                        "TP53","MYC"))

