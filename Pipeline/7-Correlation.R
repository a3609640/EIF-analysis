## Figure 5 ##
######################################
#### RNA protein correlation CCLE ####
######################################
plot.EIF.cor.CCLE <- function(x, y) {
  get_EIF_CCLE_RNA <- function(x, y){
    CCLE_RNA <- fread(
      file.path(data.file.directory, 
                "CCLE_expression_full.csv"),
      data.table = FALSE
      )
    
    # find column names contain EIF gene names
    # colnames(CCLE_RNA[grepl('PABPC1', colnames(CCLE_RNA))])
    
    
    EIF_CCLE_RNA <- CCLE_RNA[c(y,
                               #'EIF4G1 (ENSG00000114867)', 
                               #'EIF4A1 (ENSG00000161960)',
                               #'EIF4E (ENSG00000151247)',
                               #'EIF4EBP1 (ENSG00000187840)', 
                               #'PABPC1 (ENSG00000070756)',
                               "V1")]
    
    CCLE_Anno <- fread(
      file.path(data.file.directory, 
                "sample_info.csv"),
      data.table = FALSE
      )
    CCLE_Anno <- CCLE_Anno[, 1:2]
    CCLE_RNA_Anno <- merge(EIF_CCLE_RNA, CCLE_Anno, 
                           by.x = "V1", by.y = "DepMap_ID", 
                           all = T)
    CCLE_RNA_Anno <- na.omit(CCLE_RNA_Anno)
    CCLE_RNA_Anno$V1 <- NULL
    colnames(CCLE_RNA_Anno) <- c(x, #"EIF4G1", "EIF4A1", "EIF4E", 
                                 #"EIF4EBP1", #"PABPC1",
                                 "Celline")
    return(CCLE_RNA_Anno)
    
  }
  EIF_CCLE_RNA <- get_EIF_CCLE_RNA(x, y)
  
  get_EIF_CCLE_Pro <- function(x){
    CCLE_PRO <- fread(
      file.path(data.file.directory, "protein_quant_current_normalized.csv"),
      data.table = FALSE
    )
    EIF_CCLE_PRO <- CCLE_PRO[CCLE_PRO$Gene_Symbol %in% c(#"EIF4E", 
                                                         #"EIF4G1"#, 
                                                         #"EIF4A1", 
                                                         #"EIF4EBP1"#,
                                                         #"PABPC1"
                                                         x), ]
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
  EIF_CCLE_Pro <- get_EIF_CCLE_Pro(x)
  
  CCLE_RNA_Pro <- merge(EIF_CCLE_RNA, 
                        EIF_CCLE_Pro, 
                        by = "Celline",  
                        suffixes = c("RNA","Pro"),
                        all = T)
  CCLE_RNA_Pro <- na.omit(CCLE_RNA_Pro)
  
  Scatter.plot <- function(x){
    p1 <- ggscatter(CCLE_RNA_Pro, 
                    x = paste0(x, "Pro"), 
                    y = paste0(x, "RNA"), #color = "Type", 
                    add = "reg.line", #conf.int = TRUE, 
                    cor.coef = TRUE, 
                    cor.method = "pearson", 
                    title = x,
                    xlab = "Protein expresion", 
                    ylab = "RNA expression")+
      stat_cor(aes(color = Type),     
               label.x.npc = 0.5, 
               label.y.npc = 0.5, 
               hjust = 0
      )+           # Add correlation coefficient
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12(),
        axis.text.y = black_bold_12(),
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
                    y = paste0(x, "RNA"), 
                    color = "Type", 
                    add = "reg.line", #conf.int = TRUE, 
                    cor.coef = TRUE, 
                    cor.method = "pearson", #title = x,
                    #xlab = "Protein expresion", 
                    #ylab = "RNA expression"
    )+
      stat_cor(aes(color = Type),     
               label.x.npc = 0.5, 
               label.y.npc = 0.5, 
               hjust = 0
      )+           # Add correlation coefficient
      theme_bw() +
      theme(
        plot.title = black_bold_12(),
        axis.title.x = black_bold_12(),
        axis.title.y = black_bold_12(),
        axis.text.x = black_bold_12(),
        axis.text.y = black_bold_12(),
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
  Scatter.plot(x)
}

######################################
#### RNA protein correlation LUAD ####
######################################
plot.EIF.cor.LUAD <- function(EIF) {
  get.proteomics.data <- function(EIF_list) {
    LUAD.Proteomics <- read_excel(
      file.path(data.file.directory, "Protein.xlsx"), 
      col_names = FALSE
    )
    #EIF_list <- c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1","Type","Sample")
    EIF.LUAD.Proteomics <- LUAD.Proteomics[LUAD.Proteomics$...1 %in% c(EIF_list, 
                                                                       "Type","Sample"), ]
    EIF.LUAD.Proteomics.T <- t(EIF.LUAD.Proteomics)
    colnames(EIF.LUAD.Proteomics.T) <- EIF.LUAD.Proteomics.T[1, ]
    EIF.LUAD.Proteomics.T <- EIF.LUAD.Proteomics.T[-1, ] 
    EIF.LUAD.Proteomics.T <- as.data.frame(EIF.LUAD.Proteomics.T, 
                                           stringsAsFactors = FALSE)
    EIF.LUAD.Proteomics.T[EIF_list] <- sapply(EIF.LUAD.Proteomics.T[EIF_list], 
                                              as.numeric)
    return(EIF.LUAD.Proteomics.T)
  }
  EIF.LUAD.Proteomics <- get.proteomics.data(EIF)
  
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
  EIF.LUAD.RNA <- get.rna.data(EIF)
  

    #EIF.LUAD.Proteomics <- get.proteomics.data(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1"))
  EIF.LUAD.Proteomics.Tumor <- EIF.LUAD.Proteomics[EIF.LUAD.Proteomics$Type == "Tumor", ] 
  EIF.LUAD.Proteomics.NAT <- EIF.LUAD.Proteomics[EIF.LUAD.Proteomics$Type == "NAT", ] 
  EIF.LUAD.RNA.Tumor <- EIF.LUAD.RNA[EIF.LUAD.RNA$Type == "Tumor", ] 
  #EIF.LUAD.RNA.NAT <- EIF.LUAD.RNA[EIF.LUAD.RNA$Type == "NAT", ] 
  tumor.total <- merge(EIF.LUAD.Proteomics.Tumor, 
                       EIF.LUAD.RNA.Tumor, 
                       by = c("Sample","Type"), 
                       suffixes = c(".pro",".rna"))
  #NAT.total <- merge(EIF.LUAD.Proteomics.NAT, 
  #                   EIF.LUAD.RNA.NAT, 
  #                   by = c("Sample","Type"), 
  #                   suffixes = c(".pro",".rna"))
    
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
  lapply(EIF, Scatter.plot)
  }

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
  
  p1 <- ggplot(data = EIF.cor,
               aes(x = gene,
                   y = posCORs,
                   fill = label), color = label) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = posCORs),
              position = position_dodge(width = 0.9),
              size = 3.5
    ) +
    scale_fill_manual(values = c("#CC79A7", "#0072B2", "#E69F00", 
                                 "#009E73", "#D55E00")) + # for color-blind palettes
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
    scale_fill_manual(values = c("#CC79A7", "#0072B2", "#E69F00", 
                                 "#009E73", "#D55E00")) + # for color-blind palettes
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


## prepare RNA-seq related dataset
TCGA.GTEX.RNAseq <- function() {
  TCGA.pancancer <- fread(
    file.path(data.file.directory, 
              "TcgaTargetGtex_RSEM_Hugo_norm_count"),
    data.table = FALSE) %>% 
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

### find posCOR and negCOR in the overlapping CORs from all cancer cases
plot.Venn.all <- function() {
  TCGA.GTEX.sampletype <- function () {
    TCGA.GTEX.sampletype <- read_tsv(
      file.path(data.file.directory, 
                "TcgaTargetGTEX_phenotype.txt")
      ) %>% 
      # TCGA.GTEX.sampletype <- as.data.frame(TCGA.GTEX.sampletype)
      as.data.frame(.) %>% 
      # TCGA.GTEX.sampletype <- TCGA.GTEX.sampletype[!duplicated(TCGA.GTEX.sampletype$sample), ]
      distinct(., sample, .keep_all = TRUE) %>% 
      # TCGA.GTEX.sampletype <- na.omit(TCGA.GTEX.sampletype)
      na.omit(.) %>%
      #row.names(TCGA.GTEX.sampletype) <- TCGA.GTEX.sampletype$sample
      #TCGA.GTEX.sampletype$sample <- NULL
      remove_rownames() %>%
      column_to_rownames(var = 'sample') %>%
      #TCGA.GTEX.sampletype <- TCGA.GTEX.sampletype[,"_sample_type", drop = FALSE]
      select("_sample_type", "_primary_site") %>%
      #colnames(TCGA.GTEX.sampletype) <- "sample.type"
      rename("sample.type" = "_sample_type", "primary.site" = "_primary_site") %>%
      #TCGA.GTEX.sampletype$sample.type <- as.factor(TCGA.GTEX.sampletype$sample.type)
      mutate_if(is.character, as.factor)
    return(TCGA.GTEX.sampletype)
    #Sample.ID <- row.names(TCGA.GTEX.sampletype)
    }
  #TCGA.GTEX.sampletype1 <- TCGA.GTEX.sampletype()
    
  TCGA.GTEX.RNAseq <- function() {
    TCGA.GTEX.RNAseq <- fread(file.path(data.file.directory, 
                                        "TcgaTargetGtex_RSEM_Hugo_norm_count"),
                              # data.table = FALSE gives data.frame
                              data.table = FALSE)  %>% 
      #select(c(1:50)) %>% 
      as.data.frame(.) %>% 
      distinct(., sample, .keep_all = TRUE) %>% 
      remove_rownames() %>%
      column_to_rownames(var = 'sample')

    #TCGA.GTEX.RNAseq <- as.data.frame(TCGA.GTEX.RNAseq)
    #TCGA.GTEX.RNAseq <- TCGA.GTEX.RNAseq[!duplicated(TCGA.GTEX.RNAseq$sample),
    #                                     !duplicated(colnames(TCGA.GTEX.RNAseq))]
    #row.names(TCGA.GTEX.RNAseq) <- TCGA.GTEX.RNAseq$sample
    #TCGA.GTEX.RNAseq$sample <- NULL

    TCGA.GTEX.RNAseq.t <- data.table::transpose(TCGA.GTEX.RNAseq)
    rownames(TCGA.GTEX.RNAseq.t) <- colnames(TCGA.GTEX.RNAseq)
    colnames(TCGA.GTEX.RNAseq.t) <- rownames(TCGA.GTEX.RNAseq)
    return (TCGA.GTEX.RNAseq.t)
  }
  #TCGA.GTEX.RNAseq1 <- TCGA.GTEX.RNAseq ()
  
  #TCGA.GTEX.RNAseq.sampletype1 <- merge(TCGA.GTEX.RNAseq1,
  #                                     TCGA.GTEX.sampletype1,
  #                                     by    = "row.names",
  #                                     all.x = TRUE) %>% 
  #  na.omit(.)  
  #all.sample.type <- levels(TCGA.GTEX.sampletype$sample.type)
  #all.tumor.type <- all.sample.type [!all.sample.type %in% c(
  #  "Cell Line","Normal Tissue","Solid Tissue Normal")]
  
  TCGA.GTEX.RNAseq.sampletype <- TCGA.GTEX.RNAseq.sampletype %>% 
    na.omit(.) %>%
    #mutate_if(is.character, as.factor) %>% 
    mutate_at(c("sample.type",
                "primary.disease",
                "primary.site",
                "study"), factor)
  
  all.tumor.type <- TCGA.GTEX.RNAseq.sampletype %>% 
    select(sample.type) %>% 
    mutate_if(is.character, as.factor) %>% 
    {levels(.$sample.type)}%>% 
    .[!.%in% c("Cell Line","Normal Tissue","Solid Tissue Normal")]

  
  #all.tumor.type <- levels(TCGA.GTEX.sampletype1$sample.type1) %>% 
  #  .[!.%in% c("Cell Line","Normal Tissue","Solid Tissue Normal")]
  
  
  EIF.correlation <- function(y, z) {
    TCGA.GTEX.tumor <- TCGA.GTEX.RNAseq.sampletype[
      TCGA.GTEX.RNAseq.sampletype$sample.type %in% y,
    ] %>% na.omit(.)
    Gene.ID <- names(TCGA.GTEX.RNAseq.sampletype) %>% 
      .[!.%in% c("Row.names", 
                 "sample.type",
                 "primary.disease",
                 "primary.site",
                 "study")]
   #Gene.ID <- Gene.ID [!Gene.ID %in% c("Row.names", 
  #                                     "sample.type",
  #                                     "primary.site")]
    correlation.coefficient <- function(x, y) {
      result <- cor.test(TCGA.GTEX.tumor[[x]],
                         TCGA.GTEX.tumor[[y]],
                         method = "pearson"
      )
      res <- data.frame(x,
                        y,
                        result[c("estimate",
                                 "p.value",
                                 "statistic",
                                 "method")],
                        stringsAsFactors = FALSE
      )
    }
    # find all genes positively correlate with EIF4F expression
    # lapply function gives a large list, need to convert it to a dataframe
    EIF.cor.list <- function(x) {
      cor.data <- do.call(rbind.data.frame,
                          lapply(Gene.ID,
                                 correlation.coefficient,
                                 y = x))
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
    
    plot.pos.Venn <- function() {
      c4 <- cbind(EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
                  EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
                  EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
                  EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05)
      summary(c4)
      b <- limma::vennCounts(c4)
      colnames(b) <- c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","Counts")
      vennDiagram(b)
      ## eulerr generates area-proportional Euler diagrams that display set 
      ## relationships (intersections, unions, and disjoints) with circles or ellipses.
      pos.Venn2 <- eulerr::euler(c(
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
        "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]),
        #shape = "ellipse"
        )
      p2 <- plot(pos.Venn2,
                 # key = TRUE,
                 main = paste(z, "posCOR"),
                 lwd = 0,
                 fill = c("#999999", "#009E73", 
                          "#56B4E9", "#E69F00"),
                 quantities = list(cex = 1.25),
                 labels = list(labels = c("EIF4E","EIF4G1",
                                          "EIF4A1","EIF4EBP1"),
                               cex = 1.25)
                 )
      print(p2)
      ggplot2::ggsave(
        path = file.path(output.directory, "Heatmap"),
        filename = paste("all", z, "pos4Venn.pdf"),
        plot = p2,
        width = 8,
        height = 8,
        useDingbats = FALSE)
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
    
    get.CORs <- function (){
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
      return(df2)}
    get.CORs <- get.CORs()
    ## output two variables: cor.data for the heatmap function and get.CORs for 
    ## bargraph function
    output <- list(cor.data, get.CORs) 
    return(output)
    }

  EIF.cor.tumor <- EIF.correlation(y = all.tumor.type, 
                                   z = "tumor")
  EIF.cor.normal <- EIF.correlation(y = c("Normal Tissue"), 
                                    z = "normal")
  
  plot.bargraph.CORs(
    EIF.cor.tumor = EIF.cor.tumor[[2]],
    EIF.cor.normal = EIF.cor.normal[[2]],
    tumor.label = "tumor",
    normal.label = "normal",
    gene = gene,
    posCORs = posCORs,
    negCORs = negCORs,
    output.poscor.filename = paste("all posCORs.pdf"),
    output.negcor.filename = paste("all negCORs.pdf"),
    coord_flip.ylim = 14000)
  
  plot.heatmap.CORs <- function (){
    cor.data <- cbind(
      setNames(
        data.frame(EIF.cor.tumor[[1]][1:8]),
        c(
          "EIF4E.tumor", "EIF4E.p.tumor",
          "EIF4G1.tumor", "EIF4G1.p.tumor",
          "EIF4A1.tumor", "EIF4A1.p.tumor",
          "EIF4EBP1.tumor", "EIF4EBP1.p.tumor"
        )
      ),
      setNames(
        data.frame(EIF.cor.normal[[1]][1:8]),
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
                     type = c("tumor", "tumor", "tumor", "tumor",
                              "normal", "normal", "normal", "normal"),
                     col = list(type = c(
                       "normal" = "royalblue",
                       "tumor" = "pink"
                     )),
                     cn = anno_text(gsub("\\..*", "", colnames(DF)),
                                    location = 0,
                                    rot = 0,
                                    just = "center",
                                    gp = gpar(fontsize = 15,
                                              fontface = "bold")
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
  plot.heatmap.CORs()
}

### find posCOR and negCOR in the overlapping CORs from lung cancer cases
plot.Venn.lung <- function(x) {
  TCGA.GTEX.RNAseq.sampletype.lung <- TCGA.GTEX.RNAseq.sampletype %>% 
    filter(primary.site == x) %>%
    mutate_at(c("sample.type",
                "primary.disease",
                "primary.site",
                "study"), factor)
  
  Gene.ID <- names(TCGA.GTEX.RNAseq.sampletype.lung) %>% 
    .[!.%in% c("sample.type",
               "primary.disease",
               "primary.site",
               "study")]

  all.tumor.type <- TCGA.GTEX.RNAseq.sampletype.lung %>% 
    select(sample.type) %>% 
    mutate_if(is.character, as.factor) %>% 
    {levels(.$sample.type)}%>% 
    .[!.%in% c("Cell Line","Normal Tissue","Solid Tissue Normal")]
  
  EIF.correlation <- function(y, z) {
    #TCGA.GTEX.subset.lung <- TCGA.GTEX.sampletype.lung[
    #  TCGA.GTEX.sampletype.lung$`_sample_type` %in% y,
    #]
    #TCGA.GTEX.subset.lung$`_sample_type` <- droplevels(
    #  TCGA.GTEX.subset.lung$`_sample_type`
    #)
    
    TCGA.GTEX.subset.lung <- TCGA.GTEX.RNAseq.sampletype.lung[
      TCGA.GTEX.RNAseq.sampletype.lung$sample.type %in% y, ] %>% 
      droplevels()
    
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
        lapply(Gene.ID,
               #gene.name,
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
    
    cor.data <- cbind(
      setNames(data.frame(EIF4E.cor[, c(3, 4)]), c("EIF4E", "EIF4E.p")),
      setNames(data.frame(EIF4G1.cor[, c(3, 4)]), c("EIF4G1", "EIF4G1.p")),
      setNames(data.frame(EIF4A1.cor[, c(3, 4)]), c("EIF4A1", "EIF4A1.p")),
      setNames(data.frame(EIF4EBP1.cor[, c(3, 4)]), c("EIF4EBP1", "EIF4EBP1.p"))
    )
    
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
    
    get.CORs <- function (){
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
      return(df2)}
    get.CORs <- get.CORs()
    ## output two variables: cor.data for the heatmap function and get.CORs for 
    ## bargraph function
    ## output two variables: cor.data for the heatmap function and get.CORs for 
    ## bargraph function
    output <- list(cor.data, get.CORs) 
    return(output)
    }
  EIF.cor.tumor <- EIF.correlation(
    y = all.tumor.type,
    #y = c("Primary Tumor",
    #      "Metastatic",
    #      "Recurrent Tumor"),
    z = "tumor"
  )
  EIF.cor.normal <- EIF.correlation(
    y = c("Normal Tissue"),
    z = "normal"
  )
  plot.bargraph.CORs(
    EIF.cor.tumor = EIF.cor.tumor[[2]],
    EIF.cor.normal = EIF.cor.normal[[2]],
    tumor.label = paste(x, "tumor"),
    normal.label = paste("Normal", x),
    gene = gene,
    posCORs = posCORs,
    negCORs = negCORs,
    output.poscor.filename = paste(x, "posCORs barplot.pdf"),
    output.negcor.filename = paste(x, "negCORs barplot.pdf"),
    coord_flip.ylim = 15000)
  
  plot.heatmap.CORs <- function (){
    cor.data <- cbind(
      setNames(
        data.frame(EIF.cor.tumor[[1]][1:8]),
        c(
          "EIF4E.tumor", "EIF4E.p.tumor",
          "EIF4G1.tumor", "EIF4G1.p.tumor",
          "EIF4A1.tumor", "EIF4A1.p.tumor",
          "EIF4EBP1.tumor", "EIF4EBP1.p.tumor"
        )
      ),
      setNames(
        data.frame(EIF.cor.normal[[1]][1:8]),
        c("EIF4E.normal", "EIF4E.p.normal",
          "EIF4G1.normal", "EIF4G1.p.normal",
          "EIF4A1.normal", "EIF4A1.p.normal",
          "EIF4EBP1.normal", "EIF4EBP1.p.normal")
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
                     type = c("tumor", "tumor", "tumor", "tumor",
                              "normal", "normal", "normal", "normal"),
                     col = list(type = c("normal" = "royalblue",
                                         "tumor" = "pink")),
                     cn = anno_text(gsub("\\..*", "", colnames(DF)),
                                    location = 0,
                                    rot = 0,
                                    just = "center",
                                    gp = gpar(fontsize = 15,
                                              fontface = "bold"))
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
  }
  plot.heatmap.CORs()
}

### plot heatmapand pathway analysis on clusters from all cancer cases ###
plot.heatmap.total <- function() {
  Data <- read_tsv(file.path(data.file.directory, 
                             "TcgaTargetGTEX_phenotype.txt"))
  Data <- Data[!duplicated(Data$sample), ]
  Data <- na.omit(Data)
  Data <- as.data.frame(Data)
  row.names(Data) <- Data$sample
  Data$sample <- NULL
  Sample.ID <- row.names(Data)
  subset <- as.data.frame(Data$`_sample_type`)
  row.names(subset) <- row.names(Data)
  colnames(subset) <- "sample.type"
  tissue.GTEX.TCGA.gene <- function() {
    TCGA.GTEX <- fread(file.path(data.file.directory, 
                                 "TcgaTargetGtex_RSEM_Hugo_norm_count"),
                       data.table = FALSE) # data.table = FALSE gives data.frame
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
                   type = c("tumor", "tumor", "tumor", "tumor",
                            "normal", "normal", "normal", "normal"),
                   col = list(type = c(
                     "normal" = "royalblue",
                     "tumor" = "pink"
                   )),
                   cn = anno_text(gsub("\\..*", "", colnames(DF)),
                                  location = 0,
                                  rot = 0,
                                  just = "center",
                                  gp = gpar(fontsize = 15,
                                            fontface = "bold")
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
  Sampletype <- read_tsv(file.path(data.file.directory, 
                                   "TcgaTargetGTEX_phenotype.txt"))
  tissue.GTEX.TCGA.gene <- function(x) {
    # download https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz
    TCGA.GTEX <- fread(file.path(data.file.directory, 
                                 "TcgaTargetGtex_RSEM_Hugo_norm_count"),
                       data.table = FALSE)
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
      c("EIF4E.normal", "EIF4E.p.normal",
        "EIF4G1.normal", "EIF4G1.p.normal",
        "EIF4A1.normal", "EIF4A1.p.normal",
        "EIF4EBP1.normal", "EIF4EBP1.p.normal")
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
                   type = c("tumor", "tumor", "tumor", "tumor",
                            "normal", "normal", "normal", "normal"),
                   col = list(type = c("normal" = "royalblue",
                                       "tumor" = "pink")),
                   cn = anno_text(gsub("\\..*", "", colnames(DF)),
                                  location = 0,
                                  rot = 0,
                                  just = "center",
                                  gp = gpar(fontsize = 15,
                                            fontface = "bold"))
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


# Figure 5
plot.EIF.cor.CCLE("EIF4G1", 'EIF4G1 (ENSG00000114867)')
plot.EIF.cor.CCLE("EIF4A1", 'EIF4A1 (ENSG00000161960)')
plot.EIF.cor.CCLE("EIF4E", 'EIF4E (ENSG00000151247)')
plot.EIF.cor.CCLE("EIF4EBP1", 'EIF4EBP1 (ENSG00000187840)')
plot.EIF.cor.CCLE("PABPC1", 'PABPC1 (ENSG00000070756)')

plot.EIF.cor.LUAD(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1","PABPC1"))


plot.Venn.all()
plot.Venn.lung(x = "Lung")
#plot.heatmap.total()
#plot.heatmap.lung(x = "Lung")
