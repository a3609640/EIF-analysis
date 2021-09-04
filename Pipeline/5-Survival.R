## prepare TCGA RNA-seq dataset
TCGA.RNAseq <- function() {
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
TCGA.RNAseq <- TCGA.RNAseq()   

## get OS data ##
TCGA.OS <- fread(
  file.path(data.file.directory, 
            "Survival_SupplementalTable_S1_20171025_xena_sp"),
  data.table = FALSE
  ) %>% 
  distinct(., sample, .keep_all = TRUE)  %>%
  #remove_rownames() %>%
  #column_to_rownames(var = 'sample') %>%
  select("sample","OS", "OS.time") %>% 
  rename(rn = sample)


  ## get sample type data ##
TCGA.sampletype <- readr::read_tsv(
  file.path(data.file.directory, 
            "TCGA_phenotype_denseDataOnlyDownload.tsv")
  ) %>% 
  as.data.frame(.) %>% 
  distinct(., sample, .keep_all = TRUE) %>% 
  #na.omit(.) %>%
  #remove_rownames() %>%
  #column_to_rownames(var = 'sample') %>%
  select("sample",
         "sample_type",
         "_primary_disease") %>%
  rename(rn = sample,
         sample.type = sample_type, 
         primary.disease = `_primary_disease`)

## combine OS, sample type and RNAseq data ##    
TCGA.RNAseq.OS.sampletype <- list(TCGA.RNAseq,TCGA.OS, TCGA.sampletype) %>% 
  reduce(full_join, by = "rn") %>%  
  remove_rownames() %>%
  column_to_rownames(var = 'rn')


plot.KM <- function(gene, data, cutoff, tumor) {
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
  km <- survfit(SurvObj ~ data$Group, data = data, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ data$Group, data = data, rho = 0) # rho = 0 log-rank
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
    guides(fill = FALSE) +
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

#x <- setdiff (df1, df)
###########################
##  KM survival analyses ##
###########################
plot.km.EIF.all.tumors <- function(EIF, cutoff, tumor) {
  df <- TCGA.RNAseq.OS.sampletype %>%
    select(all_of(EIF),
           "OS", 
           "OS.time",       
           "sample.type",
           "primary.disease") %>% 
    #drop_na(EIF) %>% 
    drop_na() %>% 
    #na.omit(.) %>%
    as.data.frame(.) %>%
    rename(RNAseq = EIF) %>%
    mutate(Group = case_when(
      RNAseq < quantile(RNAseq, cutoff) ~ "Bottom %", 
      RNAseq > quantile(RNAseq, (1-cutoff)) ~ "Top %")) %>%
    mutate(SurvObj = Surv(OS.time, OS == 1))
  
  # test for normal distribution by Anderson–Darling test
  plot(density(df$RNAseq), xlab = paste(EIF, "log2(TPM)"), 
       main = "Density plot of RNASeq")
  test <- ad.test(df$RNAseq)
  shapiro.test(df$RNAseq[0:5000])
  print(paste(EIF, test))
  plot.KM(gene = EIF, data = df, cutoff = cutoff, tumor = tumor) 
  }

plot.km.EIF.each.tumor <- function(EIF, cutoff, tumor) {
  df <- TCGA.RNAseq.OS.sampletype %>%
    select(all_of(EIF),
           "OS", 
           "OS.time",       
           "sample.type",
           "primary.disease") %>% 
    #drop_na(EIF) %>% 
    filter(if (tumor != "All") primary.disease == tumor else TRUE) %>%  
    #{if (tumor != "All") filter(primary.disease == tumor) } %>% 
    #filter(primary.disease == tumor) %>% 
    drop_na() %>% 
    #na.omit(.) %>%
    as.data.frame(.) %>%
    rename(RNAseq = EIF) %>%
    mutate(Group = case_when(
      RNAseq < quantile(RNAseq, cutoff) ~ "Bottom %", 
      RNAseq > quantile(RNAseq, (1-cutoff)) ~ "Top %")) %>%
    mutate(SurvObj = Surv(OS.time, OS == 1)) 
  
  plot.KM(gene = EIF, data = df, cutoff = cutoff, tumor = tumor) 
}

##########################
## Cox regression model ##
##########################
plot.univariate <- function(df, covariate_names, data.np, output.file, plot.title, x.tics, x.range) {
  result <- map(
    vars(
      EIF4E, EIF4E2, EIF4E3, 
      EIF4G1, EIF4G2, EIF4G3, 
      EIF4A1, EIF4A2, EIF3D, 
      EIF3E, EIF4EBP1, EIF4EBP2, #PABPC1, 
      MKNK1, MKNK2, EIF4B, EIF4H, 
      MTOR, #RPS6KB1, 
      MYC
    ),
    function(by) {
      analyse_multivariate(df,
                           vars(OS.time, OS),
                           covariates = list(by), # covariates expects a list
                           covariate_name_dict = covariate_names
      )
    }
  )
  
  HR.table <- function(x) {
    list <- as.data.frame(result[[x]]["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
    return(list)
  }
  a <- lapply(c(1:18), HR.table)
  b <- do.call(rbind.data.frame, a)
  
  # Testing proportional Hazards assumption on univariate analysis
  univ_formulas <- sapply(
    covariate_names,
    function(x) as.formula(paste("Surv(OS.time, OS)~", x))
  )
  univ_models <- lapply(
    univ_formulas,
    function(x) {
      cox.zph(coxph(x, data = df))
    }
  )
  coxassump <- function(x) {
    c <- print(univ_models[[x]])
    return(c)
  }
  univ_results <- lapply(covariate_names, coxassump)
  d <- do.call(rbind.data.frame, univ_results)
  e <- d[-grep("GLOBAL", rownames(d)), ]
  rownames(e) <- gsub("\\..*", "", rownames(e))
  f <- e[, 3, drop = FALSE]
  colnames(f) <- "pinteraction"
  
  data <- merge(b, f, by.x = "factor.id", by.y = "row.names")
  data[, 4:11] <- round(data[, 4:11], digits = 3)
  data[, 4:6] <- round(data[, 4:6], digits = 2)
  data$np <- data.np
  data <- as.data.frame(data)
  data$HRCI <- paste0(data$HR, " (", data$Lower_CI, "-", data$Upper_CI, ")")
  data$p[data$p < 0.001] <- "<0.001"
  data$pinteraction[data$pinteraction < 0.001] <- "<0.001"
  data <- data[order(data$HR, decreasing = TRUE), ]
  tabletext1 <- cbind(
    c("Gene", data$factor.id),
    c("No. of\nPatients", data$np),
    c("Hazard Ratio\n(95% CI)", data$HRCI),
    c("P Value", data$p),
    c("P Value for\nInteraction", data$pinteraction)
  )
  
  pdf(
    file = output.file,
    width = 14,
    height = 12,
    onefile = F
  )
  p <- forestplot(
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
    new_page = getOption("forestplot_new_page", FALSE)
  )
  dev.off()
}

plot.multivariate <- function(df, covariate_names, data.np, output.file, plot.title, x.tics, x.range) {
  df %>%
    analyse_multivariate(vars(OS.time, OS), 
                         covariates = vars(EIF4E, EIF4E2, EIF4E3, 
                                           EIF4G1, EIF4G2, EIF4G3, 
                                           EIF4A1, EIF4A2, EIF3D, 
                                           EIF3E, EIF4EBP1, EIF4EBP2, #PABPC1, 
                                           MKNK1, MKNK2, EIF4B, EIF4H,
                                           MTOR, #RPS6KB1, 
                                           MYC), 
                         covariate_name_dict = covariate_names) -> result1
  data <- as.data.frame(result1["summaryAsFrame"], col.names = NULL) # remove summaryASFrame in colnames
  
  # Testing proportional Hazards assumption on univariate analysis
  mv_fit <- coxph(Surv(OS.time, OS) ~  EIF4E + EIF4E2 + EIF4E3 + 
                    EIF4G1 + EIF4G2 + EIF4G3 + 
                    EIF4A1 + EIF4A2 + EIF3D + 
                    EIF3E + EIF4EBP1 + EIF4EBP2 + #PABPC1 + 
                    MKNK1 + MKNK2 + EIF4B + EIF4H +
                    MTOR + #RPS6KB1 + 
                    MYC, data = df)
  test.ph <- cox.zph(mv_fit)
  test <- print(test.ph)
  f <- test[, 3, drop = FALSE]
  colnames(f) <- "pinteraction"
  
  data <- merge(data, f, by.x = "factor.id", by.y = "row.names")
  data[, 4:11] <- round(data[, 4:11], digits = 3)
  data[, 4:6] <- round(data[, 4:6], digits = 2)
  data$np <- data.np
  data <- as.data.frame(data)
  data$HRCI <- paste0(data$HR, " (", data$Lower_CI, "-", data$Upper_CI, ")")
  data$p[data$p < 0.001] <- "<0.001"
  data$pinteraction[data$pinteraction < 0.001] <- "<0.001"
  data <- data[order(data$HR, decreasing = TRUE), ]
  tabletext1 <- cbind(
    c("Gene", data$factor.id),
    c("No. of\nPatients", data$np),
    c("Hazard Ratio\n(95% CI)", data$HRCI),
    c("P Value", data$p),
    c("P Value for\nInteraction", data$pinteraction)
  )
  
  pdf(
    file = output.file,
    width = 14,
    height = 12,
    onefile = F
  )
  p <- forestplot(
    labeltext = tabletext1,
    graph.pos = 3, graphwidth = unit(12, "cm"),
    hrzl_lines = list(
      "1" = gpar(lwd = 1, col = "black"),
      "2" = gpar(lwd = 1, col = "black")
    ),
    #  "3.75" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
    #  "3" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
    #  "7" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922"),
    #  "9" = gpar(lwd=6, lineend="butt", columns=c(2:6), col="#99999922")),
    mean = c(NA, data$HR),
    lower = c(NA, data$Lower_CI),
    upper = c(NA, data$Upper_CI),
    title = plot.title,
    xlab = "     <---Good prognosis---    ---Poor prognosis--->",
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
    new_page = getOption("forestplot_new_page", FALSE)
  )
  dev.off()
  print(p)
}

plot.coxph.EIF.all.tumors <- function() {
  pan.TCGA.gene <- function(EIF) {
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      file.path(data.file.directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
      data.table = FALSE
    )
    TCGA.RNAseq1 <- TCGA.RNAseq[
      !duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))
    ]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1 <- as.data.frame(TCGA.RNAseq1)
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    TCGA.OS <- fread(
      file.path(data.file.directory, "Survival_SupplementalTable_S1_20171025_xena_sp"),
      data.table = FALSE
    )
    TCGA.OS1 <- TCGA.OS[
      !duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))
    ]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1 <- as.data.frame(TCGA.OS1)
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[, c("OS", "OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE
    )
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    TCGA.OS.sampletype$Row.names <- NULL
    TCGA.OS.sampletype$sample.type <- as.factor(TCGA.OS.sampletype$sample.type)
    # remove "solid tissue normal from dataset "
    TCGA.OS.sampletype <-
      TCGA.OS.sampletype[TCGA.OS.sampletype$sample.type != "Solid Tissue Normal", ]
    TCGA.OS.sampletype$sample.type <- droplevels(TCGA.OS.sampletype$sample.type)
    levels(TCGA.OS.sampletype$sample.type)
    
    ## combine OS, sample type and RNAseq data ##
    TCGA.RNAseq.OS.sampletype <- merge(TCGA.RNAseq_transpose,
                                       TCGA.OS.sampletype,
                                       by    = "row.names",
                                       all.x = TRUE
    )
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease
    )
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(c(
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF4G1", "EIF4G2", "EIF4G3", 
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2", #"PABPC1", 
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", #"RPS6KB1", 
    "MYC"
  ))
  
  covariate_names <- c(
    EIF4E = "EIF4E", EIF4E2 = "EIF4E2", EIF4E3 = "EIF4E3", 
    EIF4G1 = "EIF4G1", EIF4G2 = "EIF4G2", EIF4G3 = "EIF4G3",
    EIF4A1 = "EIF4A1", EIF4A2 = "EIF4A2", EIF3D = "EIF3D", 
    EIF3E = "EIF3E", EIF4EBP1 = "EIF4EBP1", EIF4EBP2 = "EIF4EBP2",#PABPC1 = "PABPC1", 
    MKNK1 = "MKNK1", MKNK2 = "MKNK2", EIF4B = "EIF4B", EIF4H = "EIF4H",
    MTOR = "MTOR", #RPS6KB1 = "RPS6KB1", 
    MYC = "MYC"
  )
  
  plot.univariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 10235,
    output.file = file.path(output.directory, "Cox", "EIFUniCox.pdf"),
    plot.title = "Univariate Cox proportional-hazards regression analysis (all tumor types)",
    x.tics = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8),
    x.range = c(0.6, 1.8))
  
  plot.multivariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 10235,
    output.file = file.path(output.directory, "Cox", "EIFmultiCox.pdf"),
    plot.title = "Multivariate Cox proportional-hazards regression analysis (all tumor types)",
    x.tics = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8),
    x.range = c(0.6, 1.8))
}

plot.coxph.EIF.each.tumor <- function(tumor) {
  pan.TCGA.gene <- function(EIF, tumor) {
    ## get TCGA pancancer RNAseq data ##
    # download https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
    TCGA.RNAseq <- fread(
      file.path(data.file.directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"),
      data.table = FALSE
    )
    TCGA.RNAseq1 <- TCGA.RNAseq[
      !duplicated(TCGA.RNAseq$sample),
      !duplicated(colnames(TCGA.RNAseq))
    ]
    row.names(TCGA.RNAseq1) <- TCGA.RNAseq1$sample
    TCGA.RNAseq1 <- as.data.frame(TCGA.RNAseq1)
    TCGA.RNAseq1$sample <- NULL
    TCGA.RNAseq1 <- TCGA.RNAseq1[EIF, ]
    TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq1)
    rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq1)
    colnames(TCGA.RNAseq_transpose) <- EIF
    
    ## get OS data ##
    TCGA.OS <- fread(
      file.path(data.file.directory, "Survival_SupplementalTable_S1_20171025_xena_sp"),
      data.table = FALSE
    )
    TCGA.OS1 <- TCGA.OS[
      !duplicated(TCGA.OS$sample),
      !duplicated(colnames(TCGA.OS))
    ]
    row.names(TCGA.OS1) <- TCGA.OS1$sample
    TCGA.OS1 <- as.data.frame(TCGA.OS1)
    TCGA.OS1$sample <- NULL
    TCGA.OS1 <- TCGA.OS1[, c("OS", "OS.time")]
    
    ## get sample type data ##
    TCGA.sampletype <- readr::read_tsv(
      file.path(data.file.directory, "TCGA_phenotype_denseDataOnlyDownload.tsv")
    )
    row.names(TCGA.sampletype) <- TCGA.sampletype$sample
    TCGA.sampletype <- as.data.frame(TCGA.sampletype)
    TCGA.sampletype$sample <- NULL
    TCGA.sampletype$sample_type_id <- NULL
    colnames(TCGA.sampletype) <- c("sample.type", "primary.disease")
    
    ## combine OS and sample type data ##
    TCGA.OS.sampletype <- merge(TCGA.OS1,
                                TCGA.sampletype,
                                by    = "row.names",
                                all.x = TRUE
    )
    row.names(TCGA.OS.sampletype) <- TCGA.OS.sampletype$Row.names
    TCGA.OS.sampletype <- as.data.frame(TCGA.OS.sampletype)
    TCGA.OS.sampletype$Row.names <- NULL
    TCGA.OS.sampletype$sample.type <- as.factor(TCGA.OS.sampletype$sample.type)
    # remove "solid tissue normal from dataset "
    TCGA.OS.sampletype <-
      TCGA.OS.sampletype[TCGA.OS.sampletype$sample.type != "Solid Tissue Normal", ]
    TCGA.OS.sampletype$sample.type <- droplevels(TCGA.OS.sampletype$sample.type)
    levels(TCGA.OS.sampletype$sample.type)
    
    ## combine OS, sample type and RNAseq data ##
    TCGA.RNAseq.OS.sampletype <- merge(TCGA.RNAseq_transpose,
                                       TCGA.OS.sampletype,
                                       by    = "row.names",
                                       all.x = TRUE
    )
    row.names(TCGA.RNAseq.OS.sampletype) <- TCGA.RNAseq.OS.sampletype$Row.names
    TCGA.RNAseq.OS.sampletype <- as.data.frame(TCGA.RNAseq.OS.sampletype)
    TCGA.RNAseq.OS.sampletype$Row.names <- NULL
    ## remove all rows with NA in the primary disease section
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[!is.na(TCGA.RNAseq.OS.sampletype$primary.disease), ]
    TCGA.RNAseq.OS.sampletype$primary.disease <- as.factor(
      TCGA.RNAseq.OS.sampletype$primary.disease
    )
    TCGA.RNAseq.OS.sampletype <-
      TCGA.RNAseq.OS.sampletype[TCGA.RNAseq.OS.sampletype$primary.disease %in% tumor, ]
    return(TCGA.RNAseq.OS.sampletype)
  }
  df <- pan.TCGA.gene(c(
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF4G1", "EIF4G2", "EIF4G3", 
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2", #"PABPC1", 
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", #"RPS6KB1", 
    "MYC"
  ), tumor)
  # df$`EIF4E+EIF4EBP1` <- log2(2**df$EIF4E + 2**df$EIF4EBP1 -2 + 1)
  df <- df[c(
    "EIF4E", "EIF4E2", "EIF4E3", 
    "EIF4G1", "EIF4G2", "EIF4G3", 
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2",#"PABPC1", 
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", #"RPS6KB1", 
    "MYC", # "EIF4E+EIF4EBP1",
    "OS", "OS.time"
  )]
  # Use survivalAnalysis package to draw forest plot of multiple univariate #
  covariate_names <- c(
    EIF4E = "EIF4E", EIF4E2 = "EIF4E2", EIF4E3 = "EIF4E3", 
    EIF4G1 = "EIF4G1", EIF4G2 = "EIF4G2", EIF4G3 = "EIF4G3",
    EIF4A1 = "EIF4A1", EIF4A2 = "EIF4A2", EIF3D = "EIF3D", 
    EIF3E = "EIF3E", EIF4EBP1 = "EIF4EBP1", EIF4EBP2 = "EIF4EBP2", #PABPC1 = "PABPC1", 
    MKNK1 = "MKNK1", MKNK2 = "MKNK2", EIF4B = "EIF4B", EIF4H = "EIF4H",
    MTOR = "MTOR", #RPS6KB1 = "RPS6KB1", 
    MYC = "MYC"
  )
  
  # TODO: Use file.path below (noting that paste0() has been convenient
  #       for filename construction in the current implementation.)
  plot.univariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 517,
    output.file = paste0(output.directory, "/Cox/", tumor, "EIFUniCox.pdf"),
    plot.title = paste("Univariate Cox proportional-hazards regression analysis", tumor),
    x.tics = c(0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8),
    x.range = c(0.4, 2.8))
  
  plot.multivariate(
    df = df,
    covariate_names = covariate_names,
    data.np = 517,
    output.file = paste0(output.directory, "/Cox/", tumor, "EIFmultiCox.pdf"),
    plot.title = paste("Multivariate Cox proportional-hazards regression analysis", tumor),
    x.tics = c(0.4, 0.8, 1.2, 1.6, 2.4, 2.8, 3.2),
    x.range = c(0.4, 3.2))
}


plot.km.EIF.all.tumors(EIF = "EIF4E", 
                       cutoff =  0.3,
                       tumor = "All")

lapply(c("EIF4G1","EIF4G2", "EIF4G3",
         "EIF4A1","EIF4A2", 
         "EIF4E", "EIF4E2", "EIF4E3", 
         "EIF3D", "EIF3E",
         "EIF4EBP1", "EIF4EBP2", 
         "EIF4H", "EIF4B", "MYC",
         "PABPC1", "MKNK1", "MKNK2"), 
       plot.km.EIF.each.tumor, cutoff = 0.2, tumor = "All")

lapply(c("EIF4G1","EIF4G2", "EIF4G3",
         "EIF4A1","EIF4A2", 
         "EIF4E", "EIF4E2", "EIF4E3", 
         "EIF3D", "EIF3E","EIF4EBP1", "EIF4EBP2", 
         "EIF4H", "EIF4B", "MYC",
         "PABPC1", "MKNK1", "MKNK2"),
       plot.km.EIF.each.tumor, cutoff =  0.2,
       tumor = "lung adenocarcinoma"
)

plot.km.EIF.each.tumor(EIF = "EIF4E", 
                       cutoff =  0.3, 
                       tumor = "lung adenocarcinoma")

plot.km.EIF.each.tumor(EIF = "EIF4E", 
                       cutoff =  0.3, 
                       tumor = "All")

plot.coxph.EIF.all.tumors()

plot.coxph.EIF.each.tumor(c("lung adenocarcinoma"))
