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
                                           "Solid Tissue Normal",
                                           "Primary Tumor",
                                           "Metastatic"),
                                labels = c("Healthy Tissue (GTEx)",
                                           "Adjacent Normal Tissue (TCGA)",
                                           "Primary Tumor (TCGA)",
                                           "Metastatic Tumor (TCGA)")))
  
  plot.PCA <- function(x, y) {
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[
    #  EIF.TCGA.RNAseq.anno.subset$sample.type %in% c("Primary Tumor (TCGA)"),
    #]
    #EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno.subset[, c(EIF.list, 
    #                                                               "sample.type", 
    #                                                               "primary.site")]
    #EIF.TCGA.RNAseq.anno.subset <- droplevels(EIF.TCGA.RNAseq.anno.subset)
    ## remove the last two columns
    #df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset) - 2)]
    #rownames(df1) <- NULL
    TCGA.sampletype.subset <- TCGA.GTEX.sampletype.subset %>% 
      filter(sample.type == x) %>% 
      droplevels() %>%
      remove_rownames() 
    
    plot.pca.factomineR <- function() {
      res.pca <- PCA(TCGA.sampletype.subset[1:(length(TCGA.sampletype.subset)-4)],
                     scale.unit = TRUE,
                     ncp = 10,
                     graph = FALSE)
      #res.pca <- PCA(df1,
      #               scale.unit = TRUE,
      #               ncp = 10,
      #               graph = FALSE
      #)
      
      biplot <- fviz_pca_biplot(res.pca,
                                axes = c(1, 2),
                                labelsize = 5,
                                col.ind = TCGA.sampletype.subset[,y],
                                #col.ind = TCGA.sampletype.subset$primary.disease,
                                title = paste("PCA - Biplot (", x, ")"),
                                palette = col_vector,
                                pointshape = 20,
                                pointsize = 0.75,
                                # addEllipses = TRUE,
                                label = "var",
                                col.var = "black",
                                repel = TRUE
      ) +
        xlim(-7, 8) + 
        ylim(-6, 7.5) + # for EIF 8
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
          legend.position = c(0, .625),
          legend.background = element_blank(),
          legend.text = black_bold_16()
        )
      print(biplot)
      ggplot2::ggsave(
        path = file.path(output.directory, "PCA", "TCGA"),
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
        # geom_text(aes(label = res.pca$eig, size = 18)) +
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
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = paste("eig", x, ".pdf"),
        plot = eig,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      
      var <- get_pca_var(res.pca)
      pdf(file.path(
        path = file.path(output.directory, "PCA", "TCGA"),
        filename = paste("matrix", x, ".pdf")),
        #filename = "EIFcorprimary.pdf"),
        width = 9,
        height = 9,
        useDingbats = FALSE)
      corrplot(var$cos2, # cos2 is better than contribute
               title = paste("PCA (", x, ")"),
               #title = "PCA (Primary tumors)",
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      dev.off()
      corrplot(var$cos2, # cos2 is better than contribute
               title = paste("PCA (", x, ")"),
               #title = "PCA (Primary tumors)",
               is.corr = FALSE,
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
    plot.pca.factomineR()
  }
  plot.PCA("Primary Tumor (TCGA)", "primary.disease")
  plot.PCA("Metastatic Tumor (TCGA)", "primary.disease")
  plot.PCA("Healthy Tissue (GTEx)", "primary.site")
  
}

plot.PCA.TCGA.GTEX.all <- function(EIF.list) {
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
                                "Solid Tissue Normal",
                                "Normal Tissue")) %>%  
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
  #plot.pca.factomineR <- function() {
    #df1 <- EIF.TCGA.RNAseq.anno.subset
    #res.pca <- PCA(df1[1:(length(df1) - 3)],
    #               scale.unit = TRUE,
    #               ncp = 14,
    #               graph = FALSE
    #)
    res.pca <- PCA(TCGA.GTEX.sampletype.subset[1:(length(TCGA.GTEX.sampletype.subset)-4)],
                   scale.unit = TRUE,
                   ncp = 10,
                   graph = FALSE)
    biplot <- fviz_pca_biplot(res.pca,
                              axes = c(1, 2),
                              labelsize = 5,
                              col.ind = TCGA.GTEX.sampletype.subset$sample.type,
                              palette = c("#D55E00", "#009E73", "#CC79A7", "#0072B2"),
                              # palette    = c("#CC79A7","#0072B2"),
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
      filename = "EIFPCAall.pdf",
      plot = biplot,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )
    
    plot.selected.PCA <- function(x, y, color) {
      selected.samples <- TCGA.GTEX.sampletype.subset %>%  
        filter(sample.type == x) 
      biplot <- fviz_pca_biplot(res.pca,
                                axes = c(1, 2),
                                labelsize = 5,
                                #col.ind = TCGA.GTEX.sampletype.subset$sample.type,
                                col.ind = TCGA.GTEX.sampletype.subset[,y],
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
    plot.selected.PCA("Healthy Tissue (GTEx)", "sample.type", "#D55E00")
    plot.selected.PCA("Healthy Tissue (GTEx)", "primary.site", col_vector)
    
    plot.selected.PCA("Primary Tumor (TCGA)", "sample.type", "#009E73")
    plot.selected.PCA("Primary Tumor (TCGA)", "primary.site", col_vector)
    
    plot.selected.PCA("Metastatic Tumor (TCGA)", "sample.type", "#CC79A7")
    plot.selected.PCA("Metastatic Tumor (TCGA)", "primary.site", col_vector)
    
    plot.selected.PCA("Adjacent Normal Tissue (TCGA)", "sample.type", "#0072B2")
    plot.selected.PCA("Adjacent Normal Tissue (TCGA)", "primary.site", col_vector)    
    
    eig <- fviz_eig(res.pca,
                    labelsize = 6,
                    geom = "bar",
                    width = 0.7,
                    addlabels = TRUE
    ) +
      # geom_text(aes(label = res.pca$eig, size = 18)) +
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
      path = file.path(output.directory, "PCA", "All"),
      filename = "EIFPCAeig.pdf",
      plot = eig,
      width = 8,
      height = 8,
      useDingbats = FALSE
    )
    
    var <- get_pca_var(res.pca)
    corrplot(var$contrib, is.corr = FALSE)
    
    pdf(file.path(
      path = file.path(output.directory, "PCA", "All"),
      filename = "EIFPCAcor.pdf"
    ),
    width = 9,
    height = 9,
    useDingbats = FALSE
    )
    corrplot(var$cos2, # cos2 is better than contribute
             is.corr = FALSE,
             tl.cex = 1.5,
             number.cex = 1.5,
             method = "color",
             addgrid.col = "gray",
             addCoef.col = "black",
             tl.col = "black"
    )
    dev.off()
    corrplot(var$cos2, # cos2 is better than contribute
             is.corr = FALSE,
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
  #}
  #plot.pca.factomineR()
}

# PCA for Figure S9
plot.PCA.TCGA.GTEX.each.tumor <- function(EIF.list, tissue) {
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
                                           "Solid Tissue Normal",
                                           "Primary Tumor",
                                           "Metastatic"),
                                labels = c("Healthy Tissue (GTEx)",
                                           "Adjacent Normal Tissue (TCGA)",
                                           "Primary Tumor (TCGA)",
                                           "Metastatic Tumor (TCGA)"))) %>%
    filter(primary.site == tissue)
  

  #EIF.PCA.tissue <- function(x) {
  #  EIF.TCGA.RNAseq.anno.subset <- EIF.TCGA.RNAseq.anno[
  #    EIF.TCGA.RNAseq.anno$primary.site == x,
  #  ]
  #  print(summary(EIF.TCGA.RNAseq.anno.subset))
#    df1 <- EIF.TCGA.RNAseq.anno.subset[1:(length(EIF.TCGA.RNAseq.anno.subset) - 3)]
  #  rownames(df1) <- NULL
    
    res.pca <- PCA(TCGA.GTEX.sampletype.subset[1:(length(TCGA.GTEX.sampletype.subset)-4)],
                   scale.unit = TRUE,
                   ncp = 10,
                   graph = FALSE)
    
    #plot.pca.factomineR <- function(x) {
    #  res.pca <- PCA(df1,
    #                 scale.unit = TRUE,
    #                 ncp = 10,
    #                 graph = FALSE
    #  )
      
      biplot <- fviz_pca_biplot(res.pca,
                                axes = c(1, 2),
                                labelsize = 5,
                                col.ind = TCGA.GTEX.sampletype.subset$sample.type,
                                #col.ind = EIF.TCGA.RNAseq.anno.subset$sample.type,
                                palette = c("#D55E00", "#CC79A7", "#009E73", "#0072B2"),
                                pointshape = 20,
                                pointsize = 0.75,
                                # addEllipses = TRUE, ellipse.level = 0.9,
                                label = "var",
                                col.var = "black",
                                repel = TRUE,
                                title = paste0("PCA - Biplot (", tissue, ")")
      ) +
        # scale_x_continuous(breaks = seq(-6, 12, 2), limits=c(-5, 12)) +
        # scale_y_continuous(breaks = seq(-4, 10, 2), limits=c(-4, 10)) +
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
        path = file.path(output.directory, "PCA", "Lung"),
        filename = paste0(tissue, "EIFPCA.pdf"),
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
        # geom_text(aes(label = res.pca$eig, size = 18)) +
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
        path = file.path(output.directory, "PCA", "Lung"),
        filename = paste0(tissue, "EIFeig.pdf"),
        plot = eig,
        width = 8,
        height = 8,
        useDingbats = FALSE
      )
      
      contribplot <- function(x) {
        p <- fviz_contrib(res.pca,
                          choice = "var",
                          axes = x,
                          top = 10,
                          fill = "lightblue",
                          color = "black"
        ) +
          theme_classic() +
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
        print(p)
        ggplot2::ggsave(
          path = file.path(output.directory, "PCA", "Lung"),
          filename = paste0("EIFcontri", tissue, ".pdf"),
          plot = p,
          width = 8,
          height = 4,
          useDingbats = FALSE
        )
      }
      lapply(c(1, 2), contribplot)
      
      var <- get_pca_var(res.pca)
      # corrplot(var$contrib, is.corr=FALSE)
      # fviz_pca_var(res.pca, col.var="contrib")
      pdf(file.path(
        path = file.path(output.directory, "PCA", "Lung"),
        filename = paste0(tissue, "EIFPCAcor.pdf")
      ),
      width = 9,
      height = 9,
      useDingbats = FALSE
      )
      corrplot(var$cos2, # cos2 is better than contribute
               is.corr = FALSE,
               tl.cex = 1.5,
               number.cex = 1.5,
               method = "color",
               addgrid.col = "gray",
               addCoef.col = "black",
               tl.col = "black"
      )
      dev.off()
    #}
    #plot.pca.factomineR(x)
  #}
  # lapply(disease.list, EIF.PCA.tissue)
  #EIF.PCA.tissue(tissue)
}

plot.PCA.CPTAC.LUAD <- function(EIF.list) {
  CPTAC.LUAD.Sample <- read_excel(
    file.path(data.file.directory, 
              "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx")
  ) %>%
    as.data.frame(.) %>% 
    select("Aliquot (Specimen Label)", "Type") %>%
    distinct(., `Aliquot (Specimen Label)`, .keep_all = TRUE) %>% 
    remove_rownames() %>%
    column_to_rownames(var = 'Aliquot (Specimen Label)')

  CPTAC.LUAD.Proteomics <- function() {
    CPTAC.LUAD.Proteomics <- fread(
      file.path(data.file.directory, "CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv"),
      data.table = FALSE) %>%
      filter(Gene %in% EIF.list) %>% 
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
  
  CPTAC.LUAD.Proteomics.Sample <- merge(CPTAC.LUAD.Proteomics,
                                        CPTAC.LUAD.Sample,
                                        by    = "row.names",
                                        all.x = TRUE) %>% 
    remove_rownames() %>%
    column_to_rownames(var = 'Row.names') %>%  
    mutate_if(is.character, as.factor) %>%
    mutate(Type = factor(Type, 
                         levels = c("Normal", 
                                    "Tumor"),
                         labels = c("Adjacent Normal Tissue (CPTAC)", 
                                    "Primary Tumor (CPTAC)"))) %>%
    filter(!is.na(Type))%>%
    remove_rownames()

  #Impute the missing values of a dataset with the Principal Components Analysis model
  nb <- missMDA::estim_ncpPCA(
    CPTAC.LUAD.Proteomics.Sample[1:(length(CPTAC.LUAD.Proteomics.Sample) - 1)],
    ncp.max = 5)# estimate the number of dimensions to impute
  res.comp <- missMDA::imputePCA(
    CPTAC.LUAD.Proteomics.Sample[1:(length(CPTAC.LUAD.Proteomics.Sample) - 1)], 
    ncp = nb$ncp)
  
  res.pca <- PCA(res.comp$completeObs,
                 scale.unit = TRUE,
                 ncp = 10,
                 graph = FALSE
  )
  
  biplot <- fviz_pca_biplot(res.pca,
                            axes = c(1, 2),
                            labelsize = 5,
                            col.ind = CPTAC.LUAD.Proteomics.Sample$Type,
                            #col.ind = EIF.CPTAC.LUAD.Proteomics.Sampletype$Type,
                            palette = c("#D55E00", "#009E73"),
                            pointshape = 20,
                            pointsize = 0.75,
                            title = "PCA - Biplot (LUAD)",
                            label = "var",
                            col.var = "black",
                            repel = TRUE
  ) +
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
    path = file.path(output.directory, "PCA", "Lung"),
    filename = "EIFLUADPCA.pdf",
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
    # geom_text(aes(label = res.pca$eig, size = 18)) +
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
    path = file.path(output.directory, "PCA", "Lung"),
    filename = "EIFLUADEig.pdf",
    plot = eig,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
  var <- get_pca_var(res.pca)
  # fviz_pca_var(res.pca, col.var="contrib")
  pdf(file.path(
    path = file.path(output.directory, "PCA", "Lung"),
    filename = "EIFLUADcor.pdf"
  ),
  width = 9,
  height = 9,
  useDingbats = FALSE
  )
  corrplot(var$cos2, # cos2 is better than contribute
           is.corr = FALSE,
           tl.cex = 1.5,
           number.cex = 1.5,
           method = "color",
           addgrid.col = "gray",
           addCoef.col = "black",
           tl.col = "black"
  )
  dev.off()
  corrplot(var$cos2, # cos2 is better than contribute
           is.corr = FALSE,
           tl.cex = 1.5,
           number.cex = 1.5,
           method = "color",
           addgrid.col = "gray",
           addCoef.col = "black",
           tl.col = "black"
  )
}


plot.PCA.TCGA.GTEX(c(
  "EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))

plot.PCA.TCGA.GTEX(c("EIF4G1", "EIF4G2","EIF4G3",
                     "EIF4A1","EIF4A2", 
                     "EIF4E", "EIF4E2", "EIF4E3", 
                     "EIF4EBP1", "EIF4EBP2","MTOR",
                     "EIF3C","EIF3D","EIF3E","PABPC1", 
                     "MKNK1", "MKNK2",
                     "TP53","MYC"))

plot.PCA.TCGA.GTEX.all(c(
  "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))

plot.PCA.TCGA.GTEX.all(c(
  "EIF4G1",#"EIF4G2", "EIF4G3",
  "EIF4A1",#"EIF4A2", 
  "EIF4E", #"EIF4E2", "EIF4E3", 
  #"EIF3D", 
  "EIF4EBP1", "PABPC1", "MKNK1", "MKNK2",
  "EIF4B", "EIF4H", 
  "MYC", "JUN"
))

plot.PCA.TCGA.GTEX.each.tumor(
  c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
    "PABPC1", "MKNK1", "MKNK2"#, "MYC"
  ),
  "Lung"
)

plot.PCA.CPTAC.LUAD(c("EIF4E", "EIF4G1", "EIF4A1", "PABPC1", 
                      "MKNK1", "MKNK2", "EIF4EBP1"))
