# prepare CNV related datasets from TCGA ---------------------------------------
get.TCGA.CNV <- function() {
  TCGA.pancancer <- data.table::fread(
    file.path(data.file.directory, 
              "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"),
    data.table = FALSE
  )  %>% 
    as.data.frame(.) %>% 
    distinct(., Sample, .keep_all = TRUE) %>% 
    na.omit(.) %>%
    remove_rownames(.) %>%
    column_to_rownames(var = 'Sample')
  
  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and column names in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return (TCGA.pancancer_transpose)
}
TCGA.CNV <- get.TCGA.CNV()  # CNV data with threshold 

get.TCGA.CNV.value <- function() {
  TCGA.pancancer <- fread(
    file.path(data.file.directory, 
              "Gistic2_CopyNumber_Gistic2_all_data_by_genes"),
    data.table = FALSE
  )  %>% 
    as.data.frame(.) %>% 
    distinct(., Sample, .keep_all = TRUE) %>% 
    na.omit(.) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'Sample')
  
  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and colnames in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return (TCGA.pancancer_transpose)
}
TCGA.CNV.value <- get.TCGA.CNV.value() # CNV data with numeric value

get.TCGA.CNV.ratio <- function() {
  TCGA.pancancer <- fread(
    file.path(data.file.directory, 
              "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena"),
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
TCGA.CNV.ratio <- get.TCGA.CNV.ratio() # CNV ration between tumor and NATs

TCGA.sampletype <- readr::read_tsv(
  file.path(data.file.directory, 
            "TCGA_phenotype_denseDataOnlyDownload.tsv")) %>% {
              as.data.frame(.) %>% 
                distinct(., sample, .keep_all = TRUE) %>% 
                na.omit(.) %>%
                remove_rownames(.) %>%
                column_to_rownames(., var = 'sample') %>%
                select("sample_type", "_primary_disease") %>%
                dplyr::rename("sample.type" = "sample_type", 
                              "primary.disease" = "_primary_disease")}

TCGA.CNV.sampletype <- merge(TCGA.CNV,
                             TCGA.sampletype,
                             by    = "row.names",
                             all.x = TRUE) %>% 
  filter(sample.type != "Solid Tissue Normal") %>% {
    remove_rownames(.) %>%
      column_to_rownames(., var = 'Row.names')}

TCGA.CNVratio.sampletype <- merge(TCGA.CNV.ratio,
                                  TCGA.sampletype,
                                  by    = "row.names",
                                  all.x = TRUE) %>% {
                                    remove_rownames(.) %>%
                                      column_to_rownames(var = 'Row.names')}


# CNV data analysis and plotting -----------------------------------------------
CNV.all.cancer <- function (df) {
  TCGA.CNV.anno.subset.long <- melt(df, 
                                    id = c("sample.type",
                                           "primary.disease"),
                                    value.name = "CNV") %>% 
    mutate_if(is.character, as.factor)
  
  CNV.sum <- table(TCGA.CNV.anno.subset.long[, c("CNV", "variable")]) %>% 
    as.data.frame(.) %>% 
    mutate(CNV = factor(CNV, levels = c("-2", "-1", "0", "1", "2")))
  
  # reorder stack bars by the frequency of duplication.
  Freq.sum <- dcast(CNV.sum, variable~CNV, mean)
  CNV.sum$variable <- factor(CNV.sum$variable, 
                             levels = Freq.sum[order(Freq.sum$`1`),]$variable)
  return (CNV.sum)} 
CNV.sum.barplot <- function(data) {
  p1 <- ggplot(data, 
               aes(fill = CNV,
                   y = Freq,
                   x = variable
               )) +
    geom_bar(stat = "identity", 
             position = "fill") +
    geom_col() +
    geom_text(aes(label = paste0(Freq / 100, "%")),
              position = position_stack(vjust = 0.5), 
              size = 4
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
      labels = c("Deep del\n 0", "Shallow del\n 1", 
                 "Diploid\n 2", "Gain\n 3", "Amp\n 3+"),
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

CNV.ind.cancer <- function (df, x){
  TCGA.CNV.anno.subset.long <- df %>% 
    select(all_of(x), 
           "sample.type", 
           "primary.disease") %>% 
    melt(.,id = c("sample.type",
                  "primary.disease"), 
         value.name = "CNV") %>% 
    mutate_if(is.character, as.factor)
  
  CNV.sum <- table(TCGA.CNV.anno.subset.long[, c("CNV", "primary.disease")]) %>% 
    as.data.frame(.) %>% 
    mutate(CNV = factor(CNV, levels = c("-2", "-1", "0", "1", "2"))) %>% 
    mutate(primary.disease = forcats::fct_rev(primary.disease))
  
  output <- list(CNV.sum, x) 
  return(output)
  }
CNV.barplot <- function(df) {

  p1 <- ggplot(
    df[[1]],
    aes(
      fill = CNV, 
      order = as.numeric(CNV),
      y = Freq,
      x = primary.disease
    )
  ) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      x = "Tumor types (TCGA pan cancer atlas 2018)",
      y = paste0("Percentages of ", df[[2]], " CNVs")
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
    filename = paste0(df[[2]], "pancancerCNV.pdf"),
    plot = p1,
    width = 7.5,
    height = 9,
    useDingbats = FALSE
  )
}

matrix.plot <- function(df) {
  # correlation plot
  M = cor(df)
  testRes = corrplot::cor.mtest(df, conf.level = 0.95)
  
  p1 <-corrplot::corrplot(
    M,
    method      = "color",  
    cl.pos      = "n", # remove color legend
    tl.cex      = 1,
    number.cex  = 1,
    addgrid.col = "gray",
    addCoef.col = "black",
    tl.col      = "black",
    type        = "lower",
    order       = "FPC", 
    tl.srt      = 45,
    p.mat       = testRes$p,
    sig.level   = 0.05, # insig = "blank"
  )
  print(p1) # print correlation matrix on the screen
  # save correlation plot as a pdf file
  pdf(file.path(output.directory, "CNV", "EIFCNVcormatrix.pdf"),
      width = 9,
      height = 9,
      useDingbats = FALSE
  )
  corrplot::corrplot(
    M,
    method      = "color",  
    cl.pos      = "n", # remove color legend
    tl.cex      = 1,
    number.cex  = 1,
    addgrid.col = "gray",
    addCoef.col = "black",
    tl.col      = "black",
    type        = "lower",
    order       = "FPC", tl.srt = 45,
    p.mat       = testRes$p,
    sig.level   = 0.05, # insig = "blank"
  )
  dev.off()
}

CNVratio.tumor <- function(df, x) {
  CNVratio.data <- df %>% 
    select(all_of(x), "sample.type", "primary.disease") %>% 
    melt(.,id = c("sample.type","primary.disease"), value.name = "CNV") %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate(primary.disease = forcats::fct_rev(primary.disease))
  output <- list(CNVratio.data, x) 
  return(output)}
CNVratio.boxplot <- function(df) {
  p1 <- ggplot(
    data = df[[1]],
    aes(
      y = 2**CNV,
      x = primary.disease,
      #x = f.ordered1,
      color = primary.disease
    )
  ) +
    #ylim(0, 3) +
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
      y = paste(df[[2]], "CNV ratio", "(tumor/normal)")
    ) +
    # scale_color_manual(values = col_vector) +
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
      legend.position = "none",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_12()
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output.directory, "CNV"),
    filename = paste0(df[[2]], "pancancerCNVratio.pdf"),
    plot = p1,
    width = 7,
    height = 9,
    useDingbats = FALSE
  )
}


# master function to call CNV data analysis and plotting -----------------------
# stacked bar plots for grouped CNV status in TCGA tumors 
plot.bargraph.CNV.TCGA <- function(EIF) {
  # combine CNV data with annotation data
  TCGA.CNV.sampletype.subset <- TCGA.CNV.sampletype %>%
    select(all_of(EIF), "sample.type", "primary.disease")
      
  #stacked bar plots for CNV status in combined TCGA tumors
  CNV.all.cancer(TCGA.CNV.sampletype.subset) %>% 
    CNV.sum.barplot()
  
  #stacked bar plots for CNV status in each TCGA tumor type
  EIF.CNV.ind.cancer <- lapply(EIF, 
                               CNV.ind.cancer, 
                               df = TCGA.CNV.sampletype.subset)
  
  lapply(EIF.CNV.ind.cancer, CNV.barplot)

}

# correlation matrix for CNV values in combined TCGA tumors
plot.matrix.CNVcorr.TCGA <- function(EIF) {
  TCGA.CNVvalue.subset <- TCGA.CNV.value %>%
    select(all_of(EIF)) %>%
    matrix.plot()
}

# boxplot for EIF ratios in tumors vs adjacent normals 
plot.boxgraph.CNVratio.TCGA <- function(EIF) {
  TCGA.CNVratio.sampletype.subset <- TCGA.CNVratio.sampletype %>%
    select(all_of(EIF), 
           "sample.type", 
           "primary.disease")
  
  EIF.CNVratio.ind.cancer <- lapply(EIF, 
                                    CNVratio.tumor,
                                    df = TCGA.CNVratio.sampletype.subset)
      
  lapply(EIF.CNVratio.ind.cancer, CNVratio.boxplot)
  }



# run master functions  --------------------------------------------------------
plot.bargraph.CNV.TCGA(c("TP53", "EIF4A1","EIF4A2", 
                         "EIF4E", "EIF4E2", "EIF4E3", 
                         "MYC", "EIF3D", "EIF4EBP1", 
                         "EIF4G1","EIF4G2", "EIF4G3",
                         "PABPC1", "MKNK1", "MKNK2"))

plot.matrix.CNVcorr.TCGA(c("TP53", "EIF4A1","EIF4A2", 
                           "EIF4E", "EIF4E2", "EIF4E3", 
                           "MYC", "EIF3D", "EIF4EBP1", 
                           "EIF4G1", "EIF4G2", "EIF4G3",
                           "PABPC1", "MKNK1", "MKNK2"))


plot.boxgraph.CNVratio.TCGA(c("TP53", "EIF4A1","EIF4A2", 
                              "EIF4E", "EIF4E2", "EIF4E3", 
                              "MYC", "EIF3D", "EIF4EBP1", 
                              "EIF4G1","EIF4G2", "EIF4G3",
                              "PABPC1", "MKNK1", "MKNK2"))

