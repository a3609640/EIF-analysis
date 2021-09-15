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
                       "study" = "_study")}

TCGA.GTEX.RNAseq.sampletype <- merge(TCGA.GTEX.RNAseq,
                                     TCGA.GTEX.sampletype,
                                     by    = "row.names",
                                     all.x = TRUE) %>% {
                                       remove_rownames(.) %>%
                                         column_to_rownames(var = 'Row.names')}


### CPTAC data
CPTAC.LUAD.Proteomics <- function() {
  CPTAC.LUAD.Proteomics <- fread(
    file.path(data.file.directory, 
              "CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv"),
    data.table = FALSE) %>%
    #filter(Gene %in% EIF.list) %>% 
    remove_rownames() %>%
    column_to_rownames(var = 'Gene') %>%
    select(-contains("Unshared")) %>%
    select(c(1:(length(.)-6)))
  CPTAC.LUAD.Proteomics.t <- data.table::transpose(CPTAC.LUAD.Proteomics)
  rownames(CPTAC.LUAD.Proteomics.t) <- colnames(CPTAC.LUAD.Proteomics)
  colnames(CPTAC.LUAD.Proteomics.t) <- rownames(CPTAC.LUAD.Proteomics)
  rownames(CPTAC.LUAD.Proteomics.t) <- sub(" Log Ratio", 
                                           "", 
                                           rownames(CPTAC.LUAD.Proteomics.t))
  return (CPTAC.LUAD.Proteomics.t)
}
CPTAC.LUAD.Proteomics <- CPTAC.LUAD.Proteomics()


CPTAC.LUAD.sampletype <- read_excel(
  file.path(data.file.directory, 
            "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx")
) %>%
  as.data.frame(.) %>% 
  select("Aliquot (Specimen Label)", "Type") %>%
  distinct(., `Aliquot (Specimen Label)`, .keep_all = TRUE) %>% 
  remove_rownames() %>%
  column_to_rownames(var = 'Aliquot (Specimen Label)')%>%
  mutate(Type = factor(Type, 
                       levels = c("Normal", 
                                  "Tumor"),
                       labels = c("Adjacent Normal Tissue (CPTAC)", 
                                  "Primary Tumor (CPTAC)")))


CPTAC.LUAD.Proteomics.sampletype <- merge(CPTAC.LUAD.Proteomics,
                                      CPTAC.LUAD.sampletype,
                                      by    = "row.names",
                                      all.x = TRUE) %>% 
  remove_rownames() %>%
  column_to_rownames(var = 'Row.names') 


### functions for analyses and ploting

standardPCA <- function (df) {
  res.pca <- PCA(df %>% select_if(is.numeric), # remove column with characters
                 #df[1:(length(df)-4)],
                 scale.unit = TRUE,
                 ncp = 10,
                 graph = FALSE)
  return(res.pca)
}

imputePCA <- function (df){ 
  #Impute the missing values of a dataset with the Principal Components Analysis model
  nb <- missMDA::estim_ncpPCA(
    df %>% select_if(is.numeric),
    #CPTAC.LUAD.Proteomics.Sample[1:(length(CPTAC.LUAD.Proteomics.Sample) - 1)],
    ncp.max = 5)# estimate the number of dimensions to impute
  res.comp <- missMDA::imputePCA(
    df %>% select_if(is.numeric),
    #CPTAC.LUAD.Proteomics.Sample[1:(length(CPTAC.LUAD.Proteomics.Sample) - 1)], 
    ncp = nb$ncp)
  res.pca <- PCA(res.comp$completeObs,
                 scale.unit = TRUE,
                 ncp = 10,
                 graph = FALSE)
  return(res.pca)
}

biplot <- function(res.pca, df, x, y, color, folder) {
  biplot <- fviz_pca_biplot(res.pca,
                            axes = c(1, 2),
                            labelsize = 5,
                            col.ind = df[,y],
                            #title = paste("PCA - Biplot (", x, ")"),
                            palette = color,
                            pointshape = 20,
                            pointsize = 0.75,
                            label = "var",
                            col.var = "black",
                            repel = TRUE
                            ) +
    {if(x == "All") labs(title = "PCA - Biplot (Healthy Tissues + Tumors)") 
      else labs(title = paste0("PCA - Biplot (", x, ")"))
    } +
    #xlim(-7, 8) + 
    #ylim(-6, 7.5) + # for EIF 8
    theme_classic() +
    theme(
      plot.background = element_blank(),
      plot.title = black_bold_16(),
      panel.background = element_rect(
        fill = "transparent",
        color = "black",
        size = 1
        ),
      axis.title.x = black_bold_16(),
      axis.title.y = black_bold_16(),
      axis.text.x = black_bold_16(),
      axis.text.y = black_bold_16(),
      legend.title = element_blank(),
      legend.position = c(0, 0),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.text = black_bold_16()
      )
  print(biplot)
  ggplot2::ggsave(
    path = file.path(output.directory, "PCA", folder),
    filename = paste("PCA", x, ".pdf"),
      #filename = "EIFPCAprimary.pdf",
    plot = biplot,
    width = 8,
    height = 8,
    useDingbats = FALSE
    )
  
  eig <- fviz_eig(res.pca,
                  labelsize = 6,
                  geom = "bar",
                  width = 0.7,
                  addlabels = TRUE
                  ) +
    {if(x == "All") labs(title = "Scree plot (Healthy Tissues + Tumors)") 
      else labs(title = paste0("Scree plot (", x, ")"))
    } +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      plot.title = black_bold_16(),
      panel.background = element_rect(
        fill = "transparent",
        color = "black",
        size = 1
        ),
      axis.title.x = black_bold_16(),
      axis.title.y = black_bold_16(),
      axis.text.x = black_bold_16(),
      axis.text.y = black_bold_16()
      )
  print(eig)
  
  ggplot2::ggsave(
    path = file.path(output.directory, "PCA", folder),
    filename = paste("eig", x, ".pdf"),
    plot = eig,
    width = 8,
    height = 8,
    useDingbats = FALSE
    )

  corrtitle <-  function (x) {
    if(x == "All") title = "PC (Healthy Tissues + Tumors)" 
    else title = paste0("PC (", x, ")")
    return (title)
  }
  
  title <- corrtitle(x)  
  
  var <- get_pca_var(res.pca)
  pdf(file.path(
    path = file.path(output.directory, "PCA", folder),
    filename = paste("matrix", x, ".pdf")),
    width = 9,
    height = 9,
    useDingbats = FALSE)
  corrplot::corrplot(var$cos2, # cos2 is better than contribute
           #title = paste("PCA (", x, ")"),
           method     = "color",
           title       = title,
           #is.corr     = FALSE,
           tl.cex      = 1.5,
           number.cex  = 1.5,
           addgrid.col = "gray",
           addCoef.col = "black",
           tl.col      = "black"
           )
  dev.off()
  corrplot(var$cos2, # cos2 is better than contribute
           #title = paste("PCA (", x, ")"),
           title = title,
           #is.corr = FALSE,
           tl.cex = 1.5,
           number.cex = 1.5,
           method = "color",
           addgrid.col = "gray",
           addCoef.col = "black",
           tl.col = "black"
    )

  
  contribplot <- function(x) {
      fviz_contrib(res.pca,
                   choice = "var",
                   axes = x,
                   top = 10,
                   fill = "lightblue",
                   color = "black"
      ) +
        theme_minimal() +
        theme(
          plot.background = element_blank(),
          plot.title = black_bold_16(),
          panel.background = element_rect(
            fill = "transparent",
            color = "black",
            size = 1
          ),
          axis.title.x = element_blank(),
          axis.title.y = black_bold_16(),
          axis.text.x = black_bold_16_45(),
          axis.text.y = black_bold_16()
        )
    }
    lapply(c(1, 2), contribplot)
    
}

selected.biplot <- function(res.pca, df, x, y, color) {
  selected.samples <- df %>%  
    filter(sample.type == x) 
  biplot <- fviz_pca_biplot(res.pca,
                            axes = c(1, 2),
                            labelsize = 5,
                            #col.ind = TCGA.GTEX.sampletype.subset$sample.type,
                            col.ind = df[,y],
                            palette = color,
                            select.ind = list(name = row.names(selected.samples)),
                            pointshape = 20,
                            pointsize = 0.75,
                            # addEllipses = TRUE,
                            title = "PCA - Biplot (Healthy Tissues + Tumors)",
                            label = "var",
                            col.var = "black",
                            repel = TRUE
  ) +
    xlim(-7, 8) + ylim(-6, 7.5) + # for EIF 8
    # xlim(-6, 6) + ylim (-7, 7)+ # for EIF 4
    theme_classic() +
    theme(
      plot.background = element_blank(),
      plot.title = black_bold_16(),
      panel.background = element_rect(
        fill = "transparent",
        color = "black",
        size = 1
      ),
      axis.title.x = black_bold_16(),
      axis.title.y = black_bold_16(),
      axis.text.x = black_bold_16(),
      axis.text.y = black_bold_16(),
      legend.title = element_blank(),
      legend.position = c(0, 0),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.text = black_bold_16()
    )
  print(biplot)
  ggplot2::ggsave(
    path = file.path(output.directory, "PCA", "All"),
    filename = paste0("EIFPCAall", x, ".pdf"),
    plot = biplot,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
}

#################################################################
##  PCA plots on EIF4F RNA-seq data from TCGA and GTEx groups  ##
#################################################################
plot.PCA.TCGA.GTEX <- function(EIF.list) {
  TCGA.GTEX.sampletype.subset <- TCGA.GTEX.RNAseq.sampletype %>%
    select(all_of(EIF.list),       
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>% 
    as.data.frame(.) %>%      
    filter(EIF4E != 0 & 
             study %in% c("TCGA", "GTEX") & 
             sample.type %in% c("Metastatic",
                                "Primary Tumor",
                                "Normal Tissue",
                                "Solid Tissue Normal")) %>%  
    mutate_if(is.character, as.factor) %>%
    mutate(sample.type = factor(sample.type, 
                                levels = c("Normal Tissue",
                                           "Primary Tumor",
                                           "Metastatic",
                                           "Solid Tissue Normal"),
                                labels = c("Healthy Tissue (GTEx)",
                                           "Primary Tumor (TCGA)",
                                           "Metastatic Tumor (TCGA)",
                                           "Adjacent Normal Tissue (TCGA)")))
  
  sampletype.subset <- function(x) {
    TCGA.GTEX.sampletype.subset %>%
      filter(if (x != "All") sample.type == x else TRUE) %>%  
      droplevels() 
    } 
  
  df.subset <-sampletype.subset ("Primary Tumor (TCGA)") 
  res.pca <- standardPCA (df.subset)
  biplot(res.pca = res.pca,
         df = df.subset, 
         x = "Primary Tumor (TCGA)", 
         y = "primary.disease", 
         color = col_vector,
         folder = "TCGA")
  
  df.subset <-sampletype.subset ("Metastatic Tumor (TCGA)")
  res.pca <- standardPCA (df.subset)
  biplot(res.pca = res.pca,
         df = df.subset, 
         x = "Metastatic Tumor (TCGA)", 
         y = "primary.disease", 
         color = col_vector,
         folder = "TCGA")

  df.subset <-sampletype.subset ("Healthy Tissue (GTEx)")
  res.pca <- standardPCA (df.subset)
  biplot(res.pca = res.pca,
         df = df.subset, 
         x = "Healthy Tissue (GTEx)", 
         y = "primary.site", 
         color = col_vector,
         folder = "GTEX")
    
  
  df.subset <-sampletype.subset ("All")
  res.pca <- standardPCA (df.subset)
  biplot(res.pca = res.pca,
         df = df.subset, 
         x = "All", 
         y = "sample.type", 
         color = c("#D55E00", "#009E73", "#CC79A7", "#0072B2"),
         folder = "All")
  
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Healthy Tissue (GTEx)", 
                  y = "sample.type", 
                  color = "#D55E00")
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Healthy Tissue (GTEx)", 
                  y = "primary.site", 
                  color = col_vector)
    
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Primary Tumor (TCGA)", 
                  y = "sample.type", 
                  color = "#009E73")
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Primary Tumor (TCGA)", 
                  y = "primary.disease", 
                  color = col_vector)
    
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Metastatic Tumor (TCGA)", 
                  y = "sample.type", 
                  color = "#CC79A7")
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Metastatic Tumor (TCGA)", 
                  y = "primary.disease", 
                  color = col_vector)
  
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Adjacent Normal Tissue (TCGA)", 
                  y = "sample.type", 
                  color = "#0072B2")
  selected.biplot(res.pca = res.pca,
                  df = df.subset,
                  x = "Adjacent Normal Tissue (TCGA)", 
                  y = "primary.disease", 
                  color = col_vector)   
    
}

# PCA for Figure S9
plot.PCA.TCGA.GTEX.tumor <- function(EIF.list, tissue) {
  TCGA.GTEX.sampletype.subset <- TCGA.GTEX.RNAseq.sampletype %>%
    select(all_of(EIF.list),       
           "sample.type",
           "primary.disease",
           "primary.site",
           "study") %>% 
    as.data.frame(.) %>%      
    filter(EIF4E != 0 & 
             study %in% c("TCGA", "GTEX") & 
             sample.type %in% c("Metastatic",
                                "Primary Tumor",
                                "Normal Tissue",
                                "Solid Tissue Normal")) %>%  
    mutate_if(is.character, as.factor) %>%
    mutate(sample.type = factor(sample.type, 
                                levels = c("Normal Tissue",
                                           "Primary Tumor",
                                           "Metastatic",
                                           "Solid Tissue Normal"),
                                labels = c("Healthy Tissue (GTEx)",
                                           "Primary Tumor (TCGA)",
                                           "Metastatic Tumor (TCGA)",
                                           "Adjacent Normal Tissue (TCGA)"))) %>%
    filter(primary.site == tissue)
  res.pca <- standardPCA (TCGA.GTEX.sampletype.subset)
  biplot(res.pca = res.pca,
              df = TCGA.GTEX.sampletype.subset, 
              x = "All", 
              y = "sample.type", 
              color = c("#D55E00", "#009E73", "#CC79A7", "#0072B2"),
              folder = "Lung")  
}

plot.PCA.CPTAC.LUAD <- function(EIF.list) {
  CPTAC.LUAD.Proteomics.Sample.subset <- CPTAC.LUAD.Proteomics.sampletype %>%
    select(all_of(EIF.list),       
           "Type") %>% 
    as.data.frame(.) %>%      
    mutate_if(is.character, as.factor)  %>%
    filter(!is.na(Type))%>%
    remove_rownames()
  
  res.pca <- imputePCA (CPTAC.LUAD.Proteomics.Sample.subset)
  
  biplot(res.pca = res.pca,
              df = CPTAC.LUAD.Proteomics.Sample.subset, 
              x = "LUAD(CPTAC)", 
              y = "Type", 
              color = c("#D55E00", "#009E73"),
              folder = "Lung")
}





plot.PCA.TCGA.GTEX(c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
                     "PABPC1", "MKNK1", "MKNK2"))

plot.PCA.TCGA.GTEX(c(
  "EIF4G1","EIF4A1","EIF4E", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2",
  "EIF4B", "EIF4H", 
  "MYC", "JUN"
))

plot.PCA.TCGA.GTEX(c("EIF4G1", "EIF4G2","EIF4G3",
                     "EIF4A1","EIF4A2", 
                     "EIF4E", "EIF4E2", "EIF4E3", 
                     "EIF4EBP1", "EIF4EBP2","MTOR",
                     "EIF3C","EIF3D","EIF3E","PABPC1", 
                     "MKNK1", "MKNK2",
                     "TP53","MYC"))

plot.PCA.TCGA.GTEX.tumor(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
                           "PABPC1", "MKNK1", "MKNK2"), 
                         "Lung")

plot.PCA.CPTAC.LUAD(c("EIF4E", "EIF4G1", "EIF4A1", "PABPC1", 
                      "MKNK1", "MKNK2", "EIF4EBP1"))
