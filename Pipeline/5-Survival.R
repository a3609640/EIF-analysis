# prepare RNA-seq related dataset from TCGA and GTEx----------------------------
## prepare TCGA RNA-seq dataset
get.TCGA.RNAseq <- function() {
  TCGA.pancancer <- fread(
    file.path(data.file.directory, 
              "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
    data.table = FALSE) 
  TCGA.RNAseq <- TCGA.pancancer[
    !duplicated(TCGA.pancancer$sample),
    !duplicated(colnames(TCGA.pancancer))]%>% 
    as.data.frame(.) %>% 
    #na.omit(.) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'sample')
  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq)
  # get row and colnames in order
  #rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq)
  colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq)
  #TCGA.RNAseq_transpose$rn <- rownames(TCGA.RNAseq_transpose)
  TCGA.RNAseq_transpose$rn <- colnames(TCGA.RNAseq)
  return (TCGA.RNAseq_transpose)
}
TCGA.RNAseq <- get.TCGA.RNAseq()   

## get OS data ##
TCGA.OS <- data.table::fread(
  file.path(data.file.directory, 
            "Survival_SupplementalTable_S1_20171025_xena_sp"),
  data.table = FALSE) %>% {
  distinct(., sample, .keep_all = TRUE)  %>%
  #remove_rownames() %>%
  #column_to_rownames(var = 'sample') %>%
  select("sample","OS", "OS.time") %>% 
  rename(rn = sample)}

## get sample type data ##
TCGA.sampletype <- readr::read_tsv(
  file.path(data.file.directory, 
            "TCGA_phenotype_denseDataOnlyDownload.tsv")) %>% {
  as.data.frame(.) %>% 
  distinct(., sample, .keep_all = TRUE) %>% 
  select("sample",
         "sample_type",
         "_primary_disease") %>%
  rename(rn = sample,
         sample.type = sample_type, 
         primary.disease = `_primary_disease`)}

## combine OS, sample type and RNAseq data ##    
TCGA.RNAseq.OS.sampletype <- list(TCGA.RNAseq, TCGA.OS, TCGA.sampletype) %>% 
  reduce(full_join, by = "rn") %>%  
  remove_rownames(.) %>%
  column_to_rownames(var = 'rn')


# Survival analysis and plotting -----------------------------------------------
##  KM survival analyses 
KM.curve <- function(gene, data, cutoff, tumor) {
  # df <- subset(df, OS.time <= 4000)
  #number <- nrow(df)
  #sub <- round(number*y, digits = 0)
  # bottom.label <- paste("Bottom 20%, n = ", sub)
  # top.label <- paste("Top 20%, n = ", sub)
  #df$Group[df[[EIF]] < quantile(df[[EIF]], prob = y)] <- "Bottom %"
  #df$Group[df[[EIF]] > quantile(df[[EIF]], prob = (1-y))] <- "Top %"
  # censoring status 0 = censored, 1 = dead
  #df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  #df <- na.omit(df)
  km <- survival::survfit(SurvObj ~ data$Group, data = data, conf.type = "log-log")
  stats <- survival::survdiff(SurvObj ~ data$Group, data = data, rho = 0) # rho = 0 log-rank
  p.val <- (1 - pchisq(stats$chisq, length(stats$n) - 1)) %>%
    signif(., 3)
  
  KM <- ggplot2::autoplot(
    km,
    censor = FALSE,
    xlab = "Days",
    ylab = "Survival Probability",
    #main = paste0("All TCGA cancer studies (", nrow(df), " cases)"),
    # xlim = c(0, 4100),
    color = strata) +
    {if(tumor == "All") labs(title = paste0("All TCGA cancer studies (", nrow(data), " cases)")) 
    else labs(title = paste0(tumor, "studies (", nrow(data), " cases)"))
    }+
    theme_bw() +
    theme(
      plot.title = black_bold_16(),
      axis.title = black_bold_16(),
      axis.text = black_bold_16(),
      #axis.line.x = element_line(color = "black"),
      #axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      strip.text = black_bold_16(),
      legend.text = black_bold_16(),
      legend.title = black_bold_16(),
      legend.position = c(0.9, 0.98),
      legend.justification = c(1, 1)
    ) +
    guides(fill = "none") +
    scale_color_manual(
      values = c("red", "blue"),
      name = paste(gene, "mRNA expression"),
      breaks = c("Bottom %", "Top %"),
      labels = c(
        paste("Bottom ",percent(cutoff),", n =", round((nrow(data))*cutoff, digits = 0)),
        paste("Top ", percent(cutoff),", n =", round((nrow(data))*cutoff, digits = 0))
      )
    ) +
    annotate(
      "text",
      x        = 10000,
      y        = 0.75,
      label    = paste("log-rank test \n p.val = ", p.val),
      size     = 6.5,
      hjust    = 1,
      fontface = "bold"
    )
  
  print(KM)
  ggplot2::ggsave(
    path = file.path(output.directory, "KM"),
    filename = paste(gene, " all tumors KM.pdf"),
    plot = KM,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )
}

## Cox regression model and forest plot
forest.graph <- function (data, output.file, plot.title, x.tics, x.range){
  tabletext1 <- cbind(
    c("Gene", data$factor.id),
    c("No. of\nPatients", data$np),
    c("Hazard Ratio\n(95% CI)", data$HRCI),
    c("P Value", data$p),
    c("P Value for\nInteraction", data$pinteraction)
  )
  
  # print the forest plot on screen
  p <- forestplot::forestplot(
    labeltext = tabletext1,
    graph.pos = 3, graphwidth = unit(12, "cm"),
    hrzl_lines = list(
      "1" = gpar(lwd = 1, col = "black"),
      "2" = gpar(lwd = 1, col = "black")
    ),
    mean = c(NA, data$HR),
    lower = c(NA, data$Lower_CI),
    upper = c(NA, data$Upper_CI),
    title = plot.title,
    xlab = "<---Good prognosis---    ---Poor prognosis--->",
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = fpColors(box = "black", lines = "black"),
    xticks = x.tics,
    # xlog = 0,
    clip = x.range,
    zero = 1,
    cex = 1.2,
    lineheight = "auto", # height of the graph
    boxsize = 0.2,
    colgap = unit(6, "mm"), # the gap between column
    lwd.ci = 2,
    ci.vertices = FALSE,
    ci.vertices.height = 0.02,
    new_page = getOption("forestplot_new_page", TRUE)
  )
  print(p)
  # print the forest plot as a pdf file
  pdf(
    file = output.file,
    width = 14,
    height = 12,
    onefile = F
  )
  p <- forestplot::forestplot(
    labeltext = tabletext1,
    graph.pos = 3, graphwidth = unit(12, "cm"),
    hrzl_lines = list(
      "1" = gpar(lwd = 1, col = "black"),
      "2" = gpar(lwd = 1, col = "black")
    ),
    mean = c(NA, data$HR),
    lower = c(NA, data$Lower_CI),
    upper = c(NA, data$Upper_CI),
    title = plot.title,
    xlab = "<---Good prognosis---    ---Poor prognosis--->",
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = fpColors(box = "black", lines = "black"),
    xticks = x.tics,
    clip = x.range,
    zero = 1,
    cex = 1.2,
    lineheight = "auto", # height of the graph
    boxsize = 0.2,
    colgap = unit(6, "mm"), # the gap between column
    lwd.ci = 2,
    ci.vertices = FALSE,
    ci.vertices.height = 0.02,
    new_page = getOption("forestplot_new_page", TRUE)
  )
  print(p)
  dev.off()
  }

univariable.analysis <- function(df, covariate_names) {
  # Multiple Univariate Analyses
  res.cox <- map(covariate_names, function(gene) {
    survivalAnalysis::analyse_multivariate(df,
                                           vars(OS.time, OS),
                                           covariates = list(gene)) %>% 
      with(summaryAsFrame)
    }) %>% 
    bind_rows()  
    
# To test for the proportional-hazards (PH) assumption
  test.ph <- map(covariate_names, function (x) {
    coxph(as.formula(paste("Surv(OS.time, OS)~", x)), 
          data = df) %>%
      cox.zph() %>%
      print() %>% 
      as.data.frame(.) %>% 
      slice(1)}) %>% 
    bind_rows()  %>% 
    select("p") %>%
    rename("pinteraction" = "p") %>%
    rownames_to_column()
  
  data1 <- full_join(res.cox, test.ph, by = c("factor.id" = "rowname")) %>%
    as.data.frame() %>%
    mutate(across(7:11, round, 3)) %>%
    mutate(across(4:6, round, 2)) %>%
    mutate(np = nrow(df)) %>%
    mutate(HRCI = paste0(HR, " (", Lower_CI, "-", Upper_CI, ")")) %>%
    mutate(p = case_when(p < 0.001 ~ "<0.001",
                         #p > 0.05 ~ paste(p,"*"),
                         TRUE ~ as.character(p))) %>%
    mutate(pinteraction = case_when(pinteraction < 0.001 ~ "<0.001",
                                    pinteraction > 0.05 ~ paste(pinteraction,"*"),
                         TRUE ~ as.character(pinteraction))) %>%
    arrange(desc(HR))
  
  return(data1)
}

multivariable.analysis <- function(df, covariate_names) {
  res.cox <- survivalAnalysis::analyse_multivariate(
    df,
    vars(OS.time, OS),
    covariates = covariate_names) %>% 
    with(summaryAsFrame) #  to extract an element from a list 

  # To test for the proportional-hazards (PH) assumption
  test.ph <- coxph(as.formula(paste("Surv(OS.time, OS)~", 
                                    paste(covariate_names,collapse="+"))), 
                   data = df) %>%
    cox.zph() %>%
    print() %>% 
    as.data.frame(.) %>% 
    select("p") %>%
    rename("pinteraction" = "p") %>%
    rownames_to_column() %>%
    filter(rowname != "GLOBAL") # remove the global test result for graph
  
  data1 <- full_join(res.cox, test.ph, by = c("factor.id" = "rowname")) %>%
    as.data.frame() %>%
    mutate(across(7:11, round, 3)) %>%
    mutate(across(4:6, round, 2)) %>%
    mutate(np = nrow(df)) %>%
    mutate(HRCI = paste0(HR, " (", Lower_CI, "-", Upper_CI, ")")) %>%
    mutate(p = case_when(p < 0.001 ~ "<0.001",
                         #p > 0.05 ~ paste(p,"*"),
                         TRUE ~ as.character(p))) %>%
    mutate(pinteraction = case_when(pinteraction < 0.001 ~ "<0.001",
                                    pinteraction > 0.05 ~ paste(pinteraction,"*"),
                                    TRUE ~ as.character(pinteraction))) %>%
    arrange(desc(HR)) 
  
  return(data1)
}


# master functions to call Survival analysis and plotting ----------------------
plot.km.EIF.tumor <- function(EIF, cutoff, tumor) {
  df <- TCGA.RNAseq.OS.sampletype %>%
    filter(sample.type != "Solid Tissue Normal") %>%
    select(all_of(EIF),
           "OS", 
           "OS.time",       
           "sample.type",
           "primary.disease") %>% 
    #drop_na(EIF) %>% 
    filter(if (tumor != "All") primary.disease == tumor else TRUE) %>%  
    drop_na(.) %>% 
    #na.omit(.) %>%
    as.data.frame(.) %>% 
    rename(RNAseq = EIF) %>%
    mutate(Group = case_when(
      RNAseq < quantile(RNAseq, cutoff) ~ "Bottom %", 
      RNAseq > quantile(RNAseq, (1-cutoff)) ~ "Top %")) %>%
    mutate(SurvObj = Surv(OS.time, OS == 1)) 
  KM.curve(gene = EIF, data = df, cutoff = cutoff, tumor = tumor) 
}

plot.coxph.EIF.tumor <- function(EIF, tumor) {
  df1 <- TCGA.RNAseq.OS.sampletype %>%
    filter(sample.type != "Solid Tissue Normal") %>%
    select(all_of(EIF),
           "OS", 
           "OS.time",       
           "sample.type",
           "primary.disease") %>% 
    #drop_na(EIF) %>% 
    filter(if (tumor != "All") primary.disease == tumor else TRUE) %>%  
    drop_na(.) %>% 
    #na.omit(.) %>%
    as.data.frame(.)  
  
  univariable.result <- univariable.analysis(df = df1,
                                           covariate_names = EIF)
  forest.graph (data = univariable.result, 
                output.file = if(tumor == "All") file.path(output.directory, "Cox", "EIFUniCox.pdf")
                else paste0(output.directory, "/Cox/", tumor, "EIFUniCox.pdf"),
                plot.title = if(tumor == "All") "Univariable Cox proportional-hazards regression analysis (all tumor types)" 
                else paste("Univariable Cox proportional-hazards regression analysis", tumor),
                x.tics = if(tumor == "All") c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8) 
                else c(0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8),
                x.range = if(tumor == "All") c(0.6, 1.8) 
                else c(0.4, 2.8))
  
  
  multivariable.result <- multivariable.analysis(df = df1, 
                                               covariate_names = EIF)
  forest.graph (data = multivariable.result,
                output.file = if(tumor == "All") file.path(output.directory, "Cox", "EIFmultiCox.pdf")
                else paste0(output.directory, "/Cox/", tumor, "EIFmultiCox.pdf"),
                plot.title = if(tumor == "All") "Multivariable Cox proportional-hazards regression analysis (all tumor types)" 
                else paste("Multivariable Cox proportional-hazards regression analysis", tumor),
                x.tics = if(tumor == "All") c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8) 
                else c(0.4, 0.8, 1.2, 1.6, 2.4, 2.8, 3.2),
                x.range = if(tumor == "All") c(0.6, 1.8) 
                else c(0.4, 3.2))
  
  

  
}


# Run master functions ---------------------------------------------------------
lapply(c("EIF4G1","EIF4G2", "EIF4G3",
         "EIF4A1","EIF4A2", 
         "EIF4E", "EIF4E2", "EIF4E3", 
         "EIF3D", "EIF3E",
         "EIF4EBP1", "EIF4EBP2", 
         "EIF4H", "EIF4B", "MYC",
         "PABPC1", "MKNK1", "MKNK2"), 
       plot.km.EIF.tumor, cutoff = 0.2, tumor = "All")

lapply(c("EIF4G1","EIF4G2", "EIF4G3",
         "EIF4A1","EIF4A2", 
         "EIF4E", "EIF4E2", "EIF4E3", 
         "EIF3D", "EIF3E","EIF4EBP1", "EIF4EBP2", 
         "EIF4H", "EIF4B", "MYC",
         "PABPC1", "MKNK1", "MKNK2"),
       plot.km.EIF.tumor, cutoff =  0.2,
       tumor = "lung adenocarcinoma"
)

plot.km.EIF.tumor(EIF = "EIF4E", 
                  cutoff =  0.3, 
                  tumor = "lung adenocarcinoma")

plot.km.EIF.tumor(EIF = "EIF4E", 
                  cutoff =  0.2, 
                  tumor = "skin cutaneous melanoma")

plot.km.EIF.tumor(EIF = "EIF4E", 
                  cutoff =  0.3, 
                  tumor = "All")

plot.coxph.EIF.tumor(c(
  "EIF4E", "EIF4E2", "EIF4E3", 
  "EIF4G1", "EIF4G2", "EIF4G3", 
  "EIF4A1", "EIF4A2", "EIF3D",
  "EIF3E", "EIF4EBP1", "EIF4EBP2", #"PABPC1", 
  "MKNK1", "MKNK2", "EIF4B", "EIF4H",
  "MTOR", #"RPS6KB1", 
  "MYC"
), "All")

plot.coxph.EIF.tumor(c(
  "EIF4E", "EIF4E2", "EIF4E3", 
  "EIF4G1", "EIF4G2", "EIF4G3", 
  "EIF4A1", "EIF4A2", "EIF3D",
  "EIF3E", "EIF4EBP1", "EIF4EBP2", #"PABPC1", 
  "MKNK1", "MKNK2", "EIF4B", "EIF4H",
  "MTOR", #"RPS6KB1", 
  "MYC"
), "lung adenocarcinoma")


