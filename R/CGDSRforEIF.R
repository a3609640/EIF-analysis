# the following script generates the plot and statistics
# for EIF expression and mutations data from TCGA provisional dataset.
# install.packages("cgdsr")

library(car)
library(cgdsr)
library(corrplot)
library(cowplot)
library(data.table)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(httr)
library(plyr)
library(plotly)
library(reshape2)
library(stringr)
library(survival)
library(survMisc)
library(survminer)
library(tidyr)


# Create CGDS object
# mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
mycgds <- CGDS("http://www.cbioportal.org/")
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)
# Get cases from TCGA provisional studies only
EIF.gene <- c("EIF4A1","EIF4E","EIF4G1",
              "EIF4EBP1","EIF4EBP2","EIF4EBP3", 
              "MYC","RPS6KB1","MTOR","RPTOR","PABP")
names(EIF.gene) <- EIF.gene

ONCO <- c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","TP53","KRAS","BRAF")

####################################################################
## plot EIF RNASeq data from TCGA provisional cancer study groups ##
####################################################################
plot.EIF.provisional.tcga <- function(EIF){
  ### Get EIF RNAseq data from all TCGA study groups
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  ### "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  ### use lappy to pull out all the caselists within tcga.study.list
  ### because we named each elements in tcga.study.list,
  ### lappy will return a large list, each element (with a cancer study name)
  ### in that list is a data-table
  tcga.pro.caselist <- lapply(tcga.study.list, caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  ### for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  ### in laml study group.
  ### to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  ### we use the following tcag_provisional_caselist[[1][8,1]
  ### a <- tcga.pro.caselist[[1]][grep("tcga_rna_seq_v2_mrna",
  ### tcga.pro.caselist[[1]]$"tcga_rna_seq_v2_mrna"), 
  ### ][1,1]
  ### b <- tcga.pro.geneticprofile[[1]][
  ### grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  ### tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  ### how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.RNAseq <- function(x) {
    tcga.pro.caselist[[x]][
      grep("tcga_rna_seq_v2_mrna",
           tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pro.geneticprofile[[x]][
  ### double backslash \\ suppress the special meaning of ( )
  ### in regular expression
      grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
           tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  ### test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  ### caselist.RNAseq.1 = caselist.RNAseq ('acc_tcga')
  ### geneticprofile.RNAseq.1 = geneticprofile.RNAseq ('acc_tcga')
  ### Wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  ### within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
  }
  EIF.tcga.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            geneticprofile.RNAseq(y),
                            caselist.RNAseq(y))
  }
  ### EIF.tcga.RNAseq('SCD','acc_tcga')
  EIF.RNAseq.tcga.all <- function(x) {
    test <- lapply(tcga.study.list, 
      function(y) mapply(EIF.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "EIFgene", "TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.RNAseq.tcga.all(EIF)
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  ### plot EIF gene expression across all TCGA groups ##
  m <- paste0(EIF, ".", EIF)
  mean <- within(df2[df2$EIFgene == m,], ### TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  a <- levels(mean$TCGAstudy)
  data.long <- mean
  data.long$names <- rownames(data.long)
  data.wide <- dcast(data.long,  
                     EIFgene + names ~ TCGAstudy, 
                     value.var = "RNAseq")
  ### write.csv(data.wide, file = paste(EIF, ".csv", sep = ""))
  ### with highlight on skin cancer
  colors <- ifelse(a == "meso_tcga", "red", "black")
  print(
    ggplot(mean,
           aes(x        = TCGAstudy,
               y        = log2(RNAseq),
               color    = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      coord_flip() +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(", EIF, " RNA counts)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       hjust  = 1, ### 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, ### 1 means right-justified
                                       face   = "bold",
                                       color  = colors),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
}
plot.EIF.provisional.tcga("EIF4A1")
sapply(EIF.gene, plot.EIF.provisional.tcga)

###########################################################
## plot EIF RNASeq data from TCGA pan cancer study groups ##
############################################################
plot.EIF.pan.tcga <- function(EIF){
  # Get EIF RNAseq data from all TCGA study groups
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list,
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  # for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pro.caselist[[1]][
  # grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[1]]$case_list_id),
  # ][1,1]
  # b <- tcga.pro.geneticprofile[[1]][
  # grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  # tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  # how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.RNAseq <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      # double backslash \\ suppress the special meaning of ( )
      # in regular expression
      grep("mRNA Expression, RSEM", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  # test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  # caselist.RNAseq = caselist.RNAseq ('acc_tcga')
  # geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
  # Wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  # within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.tcga.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
      geneticprofile.RNAseq(y),
      caselist.RNAseq(y))
  }
  EIF.RNAseq.tcga.all <- function(x) {
    test <- lapply(tcga.study.list,
      function(y) mapply(EIF.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "EIFgene", "TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.RNAseq.tcga.all(EIF)
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  ### plot EIF gene expression across all TCGA groups ##
  m <- paste0(EIF, ".", EIF)
  mean <- within(df2[df2$EIFgene == m,], # TCGAstudy is one column in df2
    TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  a <- levels(mean$TCGAstudy)
  ### write.csv(mean, file = paste(EIF, ".mean", sep=""))
  ### with highlight on skin cancer
  colors <- ifelse(a == "hnsc_tcga_pan_can_atlas_2018"
    | a == "cesc_tcga_pan_can_atlas_2018", "red", "black")
  print(
    ggplot(mean,
      aes(x        = TCGAstudy,
          y        = log2(RNAseq),
          color    = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      coord_flip() +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(", EIF, " RNA counts)")) +
      theme(
        axis.title  = element_text(face   = "bold",
                                   size   = 9,
                                   color  = "black"),
        axis.text.x = element_text(size   = 9,
                                   hjust  = 1, # 1 means right-justified
                                   face   = "bold",
                                   color  = "black"),
        axis.text.y = element_text(size   = 9,
                                   angle  = 0,
                                   hjust  = 1, # 1 means right-justified
                                   face   = "bold",
                                   color  = colors),
        axis.line.x = element_line(color  = "black"),
        axis.line.y = element_line(color  = "black"),
        panel.grid  = element_blank(),
        strip.text  = element_text(face   = "bold",
                                   size   = 9,
                                   color  = "black"),
        legend.position = "none"))
}
plot.EIF.pan.tcga("EIF4A1")
sapply(EIF.gene, plot.EIF.pan.tcga)

##############################################################################
## Kaplan-Meier curve with EIF RNAseq data from all TCGA provisional groups ##
##############################################################################
plot.km.all.pro.tcga <- function(EIF) {
  #  mycgds <- CGDS("http://www.cbioportal.org/")
  #  test(mycgds)
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)",
      getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  tcga.pro.caselist <- lapply(tcga.study.list,
    caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list,
    geneticprofile)
  caselist.RNAseq <- function(x) {
    tcga.pro.caselist[[x]][
      grep("tcga_rna_seq_v2_mrna",
        tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pro.geneticprofile[[x]][
      ### double backslash \\ suppress the special meaning of ( )
      ### in regular expression
      grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
        tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  tcga.profiledata.RNAseq <- function(genename,
    geneticprofile,
    caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  tcga.gene.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
      geneticprofile.RNAseq(y),
      caselist.RNAseq(y))
  }
  tcga.EIF.RNAseq <- function(y) {
    tcga.gene.RNAseq(x = EIF, y)
  }
  ## test1 <- SCD.tcga.RNAseq(y = "skcm_tcga")
  ## test ## try to keep patient ID in the rowname
  test2 <- lapply(tcga.study.list, tcga.EIF.RNAseq)
  for(x in 1:32)
  {
    test2[[x]]$case.id <- rownames(test2[[x]])
    message("test2 = ", x)
  }
  df1 <- melt(test2)
  colnames(df1) <- c("case.id",
    "EIFgene",
    "RNAseq",
    "TCGAstudy")
  all.tcga.EIF.RNAseq <- data.frame(df1)
  all.tcga.EIF.RNAseq$TCGAstudy <- as.factor(all.tcga.EIF.RNAseq$TCGAstudy)
  message("RNAseq data retrieved")
  ##### retrieve clinic data from all tcga groups #####
  tcga.clinic.data <- function(x) {
    print(x)
    url <- function(x){
      url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="
      url <- paste0(url, x, "_all")
      return(url)
    }
    # testurl <- url("acc_tcga")
    # tesereq <- GET(url("acc_tcga"))
    req <- function(x) {GET(url(x))}
    # req <- req("acc_tcga")
    clinical_data <- function(x) {content(req(x),
      type      = 'text/tab-separated-values',
      col_names = T,
      col_types = NULL)}
    data <- clinical_data(x)
    data <- data[c("OS_MONTHS",
      "OS_STATUS",
      "CASE_ID")]
  }
  # three datasets donot have OS data and cause bugs remove them
  bug.data.set <- names(tcga.study.list) %in% c("meso_tcga", 
    "pcpg_tcga", 
    "ucs_tcga")
  tcga.study.list <- tcga.study.list[!bug.data.set]
  all.tcga.clinic.data <- lapply(tcga.study.list, tcga.clinic.data)
  all.tcga.clinic.data <- melt(all.tcga.clinic.data)
  all.tcga.clinic.data <- all.tcga.clinic.data[c("OS_STATUS",
    "CASE_ID",
    "value",
    "L1")]
  colnames(all.tcga.clinic.data) <- c("OS_STATUS",
    "case.id",
    "OS_MONTHS",
    "TCGAstudy")
  all.tcga.clinic.data$case.id <- str_replace_all(all.tcga.clinic.data$case.id,
    '-',
    '.')
  message("clinical data retrieved")
  df <- join_all(list(all.tcga.clinic.data[c("OS_MONTHS",
    "OS_STATUS",
    "case.id")],
    all.tcga.EIF.RNAseq[c("case.id",
      "RNAseq",
      "TCGAstudy")]),
    by   = "case.id",
    type = "full")
  df <- na.omit(df)
  message("clinical and RNAseq data combined")
  df$Group[df$RNAseq < quantile(df$RNAseq, prob = 0.2)] = "Bottom 20%"
  df$Group[df$RNAseq > quantile(df$RNAseq, prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  #  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0)
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
    size   = 12,
    colour = "black")
  print(
    ggplot2::autoplot(km,
      xlab = "Months",
      ylab = "Survival Probability",
      main = paste("Kaplan-Meier plot", EIF, 
        "RNA expression in all TCGA provisional groups")) +
      theme(axis.title           = black.bold.12pt,
        axis.text            = black.bold.12pt,
        axis.line.x          = element_line(color  = "black"),
        axis.line.y          = element_line(color  = "black"),
        panel.grid           = element_blank(),
        strip.text           = black.bold.12pt,
        legend.text          = black.bold.12pt ,
        legend.title         = black.bold.12pt ,
        legend.justification = c(1,1)) +
      guides(fill = FALSE) +
      scale_color_manual(values = c("red", "blue"),
        name   = paste(EIF, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c("Bottom 20%, n = 1859",
          "Top 20%, n = 1859")) +
      geom_point(size = 0.25) +
      annotate("text",
        x     = 300,
        y     = 0.85,
        label = paste("log-rank test, p.val = ", p.val),
        size  = 4.5,
        hjust = 1,
        fontface = "bold")
  )
  # rho = 1 the Gehan-Wilcoxon test
  #  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 1)
  print(EIF)
  print(stats)
}
plot.km.all.pro.tcga("PABP")
sapply(EIF.gene, plot.km.all.pro.tcga)

############################################################################
## Kaplan-Meier curve with EIF RNAseq data from all TCGA pancancer groups ##
############################################################################
plot.km.all.pan.tcga <- function(EIF) {
  #  mycgds <- CGDS("http://www.cbioportal.org/")
  #  test(mycgds)
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  pan.tcga.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)",
      getCancerStudies(mycgds)$name), ]
  pan.tcga.study.list <- pan.tcga.studies$cancer_study_id
  names(pan.tcga.study.list) <- pan.tcga.study.list
  pan.tcga.caselist <- lapply(pan.tcga.study.list,
    caselist)
  pan.tcga.geneticprofile <- lapply(pan.tcga.study.list,
    geneticprofile)
  pan.caselist.RNAseq <- function(x) {
    pan.tcga.caselist[[x]][grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna", 
      pan.tcga.caselist[[x]]$case_list_id), ][1, 1]
  } # pancancer group does not contain OS data
  pan.geneticprofile.RNAseq <- function(x) {
    pan.tcga.geneticprofile[[x]][
      grep("mRNA Expression, RSEM",
        pan.tcga.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  tcga.profiledata.RNAseq <- function(genename,
    geneticprofile,
    caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  pan.tcga.gene.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
      pan.geneticprofile.RNAseq(y),
      pan.caselist.RNAseq(y))
  }
  pan.tcga.EIF.RNAseq <- function(y) {
    pan.tcga.gene.RNAseq(x = EIF, y)
  }
  ## test1 <- SCD.tcga.RNAseq(y = "skcm_tcga")
  ## test ## try to keep patient ID in the rowname
  test2 <- lapply(pan.tcga.study.list, pan.tcga.EIF.RNAseq)
  for(x in 1:32)
  {
    test2[[x]]$case.id <- rownames(test2[[x]])
    message("test2 = ", x)
  }
  df1 <- melt(test2)
  colnames(df1) <- c("case.id",
    "EIFgene",
    "RNAseq",
    "TCGAstudy")
  all.tcga.EIF.RNAseq <- data.frame(df1)
  all.tcga.EIF.RNAseq$TCGAstudy <- as.factor(all.tcga.EIF.RNAseq$TCGAstudy)
  message("RNAseq data retrieved")
  ##### retrieve clinic data from all tcga groups #####
  tcga.clinic.data <- function(x) {
    print(x)
    url <- function(x){
      url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="
      url <- paste0(url, x, "_all")
      return(url)
    }
    # testurl <- url("acc_tcga")
    # tesereq <- GET(url("acc_tcga"))
    req <- function(x) {GET(url(x))}
    # req <- req("acc_tcga")
    clinical_data <- function(x) {content(req(x),
      type      = 'text/tab-separated-values',
      col_names = T,
      col_types = NULL)}
    data <- clinical_data(x)
    data <- data[c("OS_MONTHS",
      "OS_STATUS",
      "CASE_ID")]
  }
  pro.tcga.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "pro.tcga.study.list" contains all the tcga provisional cancer studies
  pro.tcga.study.list <- pro.tcga.studies$cancer_study_id
  names(pro.tcga.study.list) <- pro.tcga.study.list
  # three datasets donot have OS data and cause bugs remove them
  bug.data.set <- names(pro.tcga.study.list) %in% c("meso_tcga", 
    "pcpg_tcga", 
    "ucs_tcga")
  pro.tcga.study.list <- pro.tcga.study.list[!bug.data.set]
  all.tcga.clinic.data <- lapply(pro.tcga.study.list, tcga.clinic.data)
  all.tcga.clinic.data <- melt(all.tcga.clinic.data)
  all.tcga.clinic.data <- all.tcga.clinic.data[c("OS_STATUS",
    "CASE_ID",
    "value",
    "L1")]
  colnames(all.tcga.clinic.data) <- c("OS_STATUS",
    "case.id",
    "OS_MONTHS",
    "TCGAstudy")
  all.tcga.clinic.data$case.id <- str_replace_all(all.tcga.clinic.data$case.id,
    '-',
    '.')
  message("clinical data retrieved")
  df <- join_all(list(all.tcga.clinic.data[c("OS_MONTHS",
    "OS_STATUS",
    "case.id")],
    all.tcga.EIF.RNAseq[c("case.id",
      "RNAseq",
      "TCGAstudy")]),
    by   = "case.id",
    type = "full")
  df <- na.omit(df)
  message("clinical and RNAseq data combined")
  number <- round(nrow(df)/5)
  df$Group[df$RNAseq < quantile(df$RNAseq, prob = 0.2)] = "Bottom 20%"
  df$Group[df$RNAseq > quantile(df$RNAseq, prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  #  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0)
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
    size   = 12,
    colour = "black")
  print(
    ggplot2::autoplot(
      km,
      xlab = "Months",
      ylab = "Survival Probability",
      main = paste("Kaplan-Meier plot", EIF, "RNA expression"),
      xlim = c(0, 250)) +
      theme(axis.title           = black.bold.12pt,
        axis.text            = black.bold.12pt,
        axis.line.x          = element_line(color  = "black"),
        axis.line.y          = element_line(color  = "black"),
        panel.grid           = element_blank(),
        strip.text           = black.bold.12pt,
        legend.text          = black.bold.12pt ,
        legend.title         = black.bold.12pt ,
        legend.justification = c(1,1),
        legend.position      = c(1,1))+
      guides(fill = FALSE) +
      scale_color_manual(values = c("red", "blue"),
        name   = paste(EIF, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(paste("Bottom 20%, n =", number),
          paste("Top 20%, n =", number))) +
      geom_point(size = 0.25) +
      annotate("text",
        x     = 250,
        y     = 0.80,
        label = paste("log-rank test, p.val = ", p.val),
        size  = 4.5,
        hjust = 1,
        fontface = "bold"))
  
  # rho = 1 the Gehan-Wilcoxon test
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 1)
  print(EIF)
  print(stats)
}
plot.km.all.pan.tcga("PABP")
sapply(EIF.gene, plot.km.all.pan.tcga)

##########################################################
## plot RNAseq data of EIF complex in TCGA study groups ##
##########################################################
EIF.RNAseq.data <- getProfileData(mycgds,
                                  EIF.gene,
                                  "laml_tcga_rna_seq_v2_mrna",
                                  "laml_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), 
        main="EIF RNAseq in Acute Myeloid Leukemia")


EIF.RNAseq.data <- getProfileData(mycgds,
                                  EIF.gene,
                                  "acc_tcga_rna_seq_v2_mrna",
                                  "acc_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), 
        main="EIF RNAseq data in Head-Neck Squamous Cell Carcinoma")


EIF.RNAseq.data <- getProfileData(mycgds,
                                  EIF.gene,
                                  "skcm_tcga_rna_seq_v2_mrna",
                                  "skcm_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), 
        main="EIF RNAseq data in skin cutaneous melanoma")


EIF.RNAseq.data <- getProfileData(mycgds,
                                  EIF.gene,
                                  "dlbc_tcga_rna_seq_v2_mrna",
                                  "dlbc_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), 
        main="EIF RNAseq data in diffuse large B-cell lymphoma")


EIF.RNAseq.data <- getProfileData(mycgds,
                                  EIF.gene,
                                  "esca_tcga_rna_seq_v2_mrna",
                                  "esca_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), 
        main="EIF RNAseq data in Esophageal cancer")


##################################################################
## density plot of EIF complex RNA-seq in all TCGA study groups ##
##################################################################
EIF.RNAseq.all.tumor <- function(){
  EIF.RNAseq <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.RNAseq <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
      }
    geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("mRNA Expression, RSEM", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
    tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
    EIF.tcga.RNAseq <- function(x, y) {
    EIF.tcga.RNAseq <- tcga.profiledata.RNAseq(x, geneticprofile.RNAseq(y),caselist.RNAseq(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.RNAseq, keep.rownames = TRUE)[]
    return(EIF.tcga.RNAseq)
  }
    EIF.RNAseq.tcga.all <- function(x) {
    ## test[] to keep row.names by lapply function 
    test <- lapply(tcga.study.list, EIF.tcga.RNAseq, x = x)
    df2 <- melt(test)
    colnames(df2) <- c("SampleID", "EIFgene", EIF,"TCGAstudy")
    df2 <- data.frame(df2)
  }
    df2 <- EIF.RNAseq.tcga.all(EIF)
    df2$EIFgene <- NULL
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  # df2 <- na.omit(df2)
    return(df2)
    }
# EIF.RNAseq.data <- EIF.RNAseq(c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1"))
  df3 <- lapply(c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1"), EIF.RNAseq)
## use cbind to convert df3 list into a data frame
  df4 <- do.call(cbind.data.frame, df3)
## remove duplicated columns (TCGAstudy)
  df4 <- df4[, !duplicated(colnames(df4))]
  }

plot.density.RNAseq.all.tumor <- function(x , y ){
  EIF.RNAseq.data <- EIF.RNAseq.all.tumor()
  ggplot(EIF.RNAseq.data, 
    aes(x  = log2(EIF.RNAseq.data[,x]), 
        y  = log2(EIF.RNAseq.data[,y]))) + 
    labs(x = paste0("log2(", x, "RNA-seq counts)"), 
         y = paste0("log2(", y, "RNA-seq counts)"))+
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
}
plot.density.RNAseq.all.tumor("EIF4E","EIF4A1")

## 3D denisty plot by plotly
plot_ly(x = log2(EIF.RNAseq.data$EIF4A1), 
        y = log2(EIF.RNAseq.data$EIF4E), 
        z = log2(EIF.RNAseq.data$EIF4G1),
        type   = "scatter3d", 
        mode   = "markers",
        marker = list(size = 1)) %>%
    layout(scene = list(xaxis = list(title = 'EIF4A1'),
                        yaxis = list(title = 'EIF4E'),
                        zaxis = list(title = 'EIF4G1')))


################################################################
## stacked bar plots for eIF4F CNV status across tumor groups ## 
################################################################
## Get oncogene CNV data from all cancer group ##
plot.CNV.sum <- function(EIF){
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  caselist.CNV <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_cnaseq",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.CNV <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("tcga_pan_can_atlas_2018_gistic", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
  }
  tcga.profiledata.CNV <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.tcga.CNV <- function(x, y) {
    EIF.tcga.CNV <- tcga.profiledata.CNV(x, geneticprofile.CNV(y), caselist.CNV(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.CNV, keep.rownames = TRUE)[]
    return(EIF.tcga.CNV)
  }
  EIF.tcga.CNV.all <- function(x) {
    test <- lapply(tcga.study.list, EIF.tcga.CNV, x = x)
    df2 <- melt(test)
    colnames(df2) <- c("SampleID","Oncogene","CNV","TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.tcga.CNV.all(EIF)
  df2$Oncogene <- as.factor(df2$Oncogene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  CNV.sum <- table(df2[,c("CNV","TCGAstudy")])
  CNV.sum <- as.data.frame(CNV.sum)
  CNV.sum$TCGAstudy <- str_remove(CNV.sum$TCGAstudy, regex('_.*\n*.*'))
  CNV.sum$CNV <- ordered(CNV.sum$CNV, levels = c("2", "1", "0", "-1", "-2"))
  # return(CNV.sum)
  p <- ggplot(CNV.sum, aes(fill = CNV, 
                           y    = Freq, 
                           x    = TCGAstudy)) + 
       geom_bar(stat = "identity", position = "stack") +
       labs(x = "Tumor types (TCGA pan cancer atlas 2018)",
            y = paste0("Tumors with ", EIF, " CNV")) +
       theme(axis.title  = element_text(face   = "bold",
                                        size   = 9,
                                        color  = "black"),
             axis.text.x = element_text(size   = 9,
                                        angle  = 45,
                                        hjust  = 1, # 1 means right-justified
                                        face   = "bold",
                                        color  = "black"),
             axis.text.y = element_text(size   = 9,
                                        angle  = 0,
                                        hjust  = 1, # 1 means right-justified
                                        face   = "bold",
                                        color  = "black"),
             axis.line.x = element_line(color  = "black"),
             axis.line.y = element_line(color  = "black"),
             panel.grid  = element_blank(),
             strip.text  = element_text(face   = "bold",
                                        size   = 9,
                                        color  = "black"),
             legend.title = element_text(colour = "black", 
                                         size   = 9, 
                                         face   = "bold"),
             legend.text = element_text(colour  = "black", 
                                        size    = 9, 
                                        face    = "bold"),
             legend.position = c(0.8, 0.8)) +
        scale_fill_manual(name   = "Copy number variation",
                          breaks = c("2", "1", "0", "-1", "-2"),
                          labels = c("Amp","Gain","Diploid","Hetlos", "Homdel"),
                          values = c('dark red','red',
                                     'light green','blue',
                                     'dark blue')) 
    print(p)
}
plot.CNV.sum("PABP")


######################################################################
## correlation analysis of CNV from all EIF4F subunits cancer types ## 
######################################################################
get.CNV.sum <- function(EIF){
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  caselist.CNV <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_cnaseq",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.CNV <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("tcga_pan_can_atlas_2018_gistic", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
  }
  tcga.profiledata.CNV <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.tcga.CNV <- function(x, y) {
    EIF.tcga.CNV <- tcga.profiledata.CNV(x, geneticprofile.CNV(y), caselist.CNV(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.CNV, keep.rownames = TRUE)[]
    return(EIF.tcga.CNV)
  }
  EIF.tcga.CNV.all <- function(x) {
    test <- lapply(tcga.study.list, EIF.tcga.CNV, x = x)
    df2 <- melt(test)
    colnames(df2) <- c("SampleID","Oncogene","CNV","TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.tcga.CNV.all(EIF)
  df2$Oncogene <- as.factor(df2$Oncogene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  CNV.sum <- table(df2[,c("CNV")])
  CNV.sum <- as.data.frame(CNV.sum)
  colnames(CNV.sum)<- c("CNV", EIF)
  return(CNV.sum)
  }
  
  EIF.CNV.sum <- lapply(c("EIF4A1","EIF4E","EIF4G1","SOX2","EIF4EBP1","MYC","MITF","FASN"), get.CNV.sum)
  EIF.CNV.sum.2 <- do.call(cbind.data.frame, EIF.CNV.sum)
  EIF.CNV.sum.2 <- EIF.CNV.sum.2[, !duplicated(colnames(EIF.CNV.sum.2))]
  EIF.CNV.sum.3 <- EIF.CNV.sum.2[,-1]
  rownames(EIF.CNV.sum.3) <- EIF.CNV.sum.2[,1]
  M <- cor(EIF.CNV.sum.3)
  p.mat <- cor_pmat(EIF.CNV.sum.3)
  corrplot(M, method = "circle")
  corrplot(M, method = "number", type = "upper",  
           order="hclust", tl.col="black", tl.srt = 0)
  
  library(ggcorrplot)

  ggcorrplot(cor(EIF.CNV.sum.3), p.mat = cor_pmat(EIF.CNV.sum.3), hc.order=TRUE, type='lower')
  ggcorrplot(round(cor(EIF.CNV.sum.3), 1), hc.order = TRUE, type = "lower",
    lab = TRUE)
  # Leave blank on no significant coefficient
  ggcorrplot(M, method = "circle", p.mat = p.mat, hc.order = TRUE,
    type = "lower", insig = "blank")




###################################################################
## correlation analysis of RNA-seq and CNV from all cancer types ## 
###################################################################
## Get DNFA gene expression from all cancer groups ##
EIF.RNAseq <- function(EIF){
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  caselist.RNAseq <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
  geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("mRNA Expression, RSEM", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
    }
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
    }
  EIF.tcga.RNAseq <- function(x, y) {
    EIF.tcga.RNAseq <- tcga.profiledata.RNAseq(x, geneticprofile.RNAseq(y),caselist.RNAseq(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.RNAseq, keep.rownames = TRUE)[]
    return(EIF.tcga.RNAseq)
  }
  EIF.RNAseq.tcga.all <- function(x) {
    ## test[] to keep row.names by lapply function 
    test <- lapply(tcga.study.list, EIF.tcga.RNAseq, x = x)
    df2 <- melt(test)
    colnames(df2) <- c("SampleID", "EIFgene", "RNAseq","TCGAstudy")
    df2 <- data.frame(df2)
    }
  df2 <- EIF.RNAseq.tcga.all(EIF)
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  return(df2)
  }
EIF.RNAseq.data <- EIF.RNAseq(EIF.gene)

## Get oncogene CNV data from all cancer groups ##
Onco.CNV <- function(EIF){
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  caselist.CNV <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_cnaseq",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.CNV <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("tcga_pan_can_atlas_2018_gistic", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
  }
  tcga.profiledata.CNV <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.tcga.CNV <- function(x, y) {
    EIF.tcga.CNV <- tcga.profiledata.CNV(x, geneticprofile.CNV(y), caselist.CNV(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.CNV, keep.rownames = TRUE)[]
    return(EIF.tcga.CNV)
    }
  EIF.tcga.CNV.all <- function(x) {
    test <- lapply(tcga.study.list, EIF.tcga.CNV, x = x)
    df2 <- melt(test)
    colnames(df2) <- c("SampleID","Oncogene","CNV","TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.tcga.CNV.all(EIF)
  df2$Oncogene <- as.factor(df2$Oncogene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  return(df2)
  }
CNV.data <- Onco.CNV("EIF4G1")
CNV.data <- lapply(EIF.gene, Onco.CNV)

plot.CNV.RNAseq.all.tumor <- function(ONCO, EIF) {
  EIF.RNAseq <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.RNAseq <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
    geneticprofile.RNAseq <- function(x) {
      tcga.pan.geneticprofile[[x]][
        grep("mRNA Expression, RSEM", 
          tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
    }
    tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
      getProfileData(mycgds,
        genename,
        geneticprofile,
        caselist)
    }
    EIF.tcga.RNAseq <- function(x, y) {
      EIF.tcga.RNAseq <- tcga.profiledata.RNAseq(x, geneticprofile.RNAseq(y),caselist.RNAseq(y))
      ## change the rownames into the first column
      setDT(EIF.tcga.RNAseq, keep.rownames = TRUE)[]
      return(EIF.tcga.RNAseq)
    }
    EIF.RNAseq.tcga.all <- function(x) {
      ## test[] to keep row.names by lapply function 
      test <- lapply(tcga.study.list, EIF.tcga.RNAseq, x = x)
      df2 <- melt(test)
      colnames(df2) <- c("SampleID", "EIFgene", "RNAseq","TCGAstudy")
      df2 <- data.frame(df2)
    }
    df2 <- EIF.RNAseq.tcga.all(EIF)
    df2$EIFgene <- as.factor(df2$EIFgene)
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
    df2 <- na.omit(df2)
    return(df2)
  }
  EIF.RNAseq.data <- EIF.RNAseq(EIF)
  Onco.CNV <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.CNV <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_cnaseq",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
    geneticprofile.CNV <- function(x) {
      tcga.pan.geneticprofile[[x]][
        grep("tcga_pan_can_atlas_2018_gistic", 
          tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
    }
    tcga.profiledata.CNV <- function(genename, geneticprofile, caselist) {
      getProfileData(mycgds,
        genename,
        geneticprofile,
        caselist)
    }
    EIF.tcga.CNV <- function(x, y) {
      EIF.tcga.CNV <- tcga.profiledata.CNV(x, geneticprofile.CNV(y), caselist.CNV(y))
      ## change the rownames into the first column
      setDT(EIF.tcga.CNV, keep.rownames = TRUE)[]
      return(EIF.tcga.CNV)
    }
    EIF.tcga.CNV.all <- function(x) {
      test <- lapply(tcga.study.list, EIF.tcga.CNV, x = x)
      df2 <- melt(test)
      colnames(df2) <- c("SampleID","Oncogene","CNV","TCGAstudy")
      df2 <- data.frame(df2)
    }
    df2 <- EIF.tcga.CNV.all(EIF)
    df2$Oncogene <- as.factor(df2$Oncogene)
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
    df2 <- na.omit(df2)
    return(df2)
  }
  CNV.data <- Onco.CNV(ONCO)
  CNV.DNFA.RNAseq <- merge(CNV.data, EIF.RNAseq.data, by = "SampleID", all = T)
  # na.omit cannot eleminate NaN here!
  CNV.DNFA.RNAseq <- na.omit(CNV.DNFA.RNAseq)
  CNV.DNFA.RNAseq$CNV <- as.factor(CNV.DNFA.RNAseq$CNV)
  Homdel.number <- 
    nrow(CNV.DNFA.RNAseq[CNV.DNFA.RNAseq$CNV == "-2", ])
  Hetlos.number <- 
    nrow(CNV.DNFA.RNAseq[CNV.DNFA.RNAseq$CNV == "-1", ])
  Diploid.number <-
    nrow(CNV.DNFA.RNAseq[CNV.DNFA.RNAseq$CNV == "0", ])
  Gain.number <-
    nrow(CNV.DNFA.RNAseq[CNV.DNFA.RNAseq$CNV == "1", ])
  Amp.number <-
    nrow(CNV.DNFA.RNAseq[CNV.DNFA.RNAseq$CNV == "2", ])
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  print(
    ggplot(data = CNV.DNFA.RNAseq,
      aes(x     = CNV,
          y     = log2(RNAseq),
          color = CNV)) +    geom_violin(trim = FALSE) +
      geom_boxplot(alpha      = .01,
                   width      = .5) +
      
      scale_x_discrete(limits = c("-2", "-1", "0", "1", "2"), # skip NaN data
                       labels = c(
                                  "-2" = paste("Homdel \n n= ", Homdel.number),
                                  "-1" = paste("Hetlos \n n= ", Hetlos.number),
                                  "0"  = paste("Diploid \n n= ", Diploid.number),
                                  "1"  = paste("Gain \n n= ", Gain.number),
                                  "2"  = paste("Amp \n n= ", Amp.number)
                                  )) +
      labs(x = paste(ONCO, "copy number variation"),
           y = paste("log2(", EIF, "RNA counts)")) +
      theme(axis.title      = element_text(face   = "bold",
                                           size   = 9,
                                           color  = "black"),
            axis.text       = element_text(size   = 9,
                                           hjust  = 1,
                                           face   = "bold",
                                           color  = "black"),
            axis.line.x     = element_line(color  = "black"),
            axis.line.y     = element_line(color  = "black"),
            panel.grid      = element_blank(),
            strip.text      = element_text(face   = "bold",
                                           size   = 9,
                                           colour = "black"),
            legend.position = "none") +
    stat_compare_means(comparisons = list(
      c("1","0"),c("2","0"),c("-1","0"),c("-2","0")),
      method = "t.test")
  )
  a <- leveneTest(RNAseq ~ CNV, CNV.DNFA.RNAseq)
  a$data.name <- paste(EIF,'expression by', ONCO, 'status')
  b <- fligner.test(RNAseq ~ CNV, CNV.DNFA.RNAseq)
  b$data.name <- paste(EIF,'expression by', ONCO, 'status')
  print(b)
  print(a)
}
plot.CNV.RNAseq.all.tumor("EIF4A1", "EIF4E")

###############################################################
## Distribution of CNV in EIF subunits from all cancer types ## 
###############################################################
Onco.CNV <- function(EIF){
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  caselist.CNV <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_cnaseq",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.CNV <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("tcga_pan_can_atlas_2018_gistic", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
  }
  tcga.profiledata.CNV <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.tcga.CNV <- function(x, y) {
    EIF.tcga.CNV <- tcga.profiledata.CNV(x, geneticprofile.CNV(y), caselist.CNV(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.CNV, keep.rownames = TRUE)[]
    return(EIF.tcga.CNV)
  }
  EIF.tcga.CNV.all <- function(x) {
    test <- lapply(tcga.study.list, EIF.tcga.CNV, x = x)
    df2 <- melt(test)
    colnames(df2) <- c("SampleID","Oncogene","CNV","TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.tcga.CNV.all(EIF)
#  df2$Oncogene <- as.factor(df2$Oncogene)
#  df2$CNV <- as.factor(df2$CNV)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2$Oncogene <- NULL
  names(df2)[names(df2) == "CNV"] <- EIF
  # df2 <- na.omit(df2)
  return(df2)
}

plot.bubble.CNV.all.tumor <- function(x , y ){
  Onco.CNV <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.CNV <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_cnaseq",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
    geneticprofile.CNV <- function(x) {
      tcga.pan.geneticprofile[[x]][
        grep("tcga_pan_can_atlas_2018_gistic", 
          tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
    }
    tcga.profiledata.CNV <- function(genename, geneticprofile, caselist) {
      getProfileData(mycgds,
        genename,
        geneticprofile,
        caselist)
    }
    EIF.tcga.CNV <- function(x, y) {
      EIF.tcga.CNV <- tcga.profiledata.CNV(x, geneticprofile.CNV(y), caselist.CNV(y))
      ## change the rownames into the first column
      setDT(EIF.tcga.CNV, keep.rownames = TRUE)[]
      return(EIF.tcga.CNV)
    }
    EIF.tcga.CNV.all <- function(x) {
      test <- lapply(tcga.study.list, EIF.tcga.CNV, x = x)
      df2 <- melt(test)
      colnames(df2) <- c("SampleID","Oncogene","CNV","TCGAstudy")
      df2 <- data.frame(df2)
    }
    df2 <- EIF.tcga.CNV.all(EIF)
    df2$Oncogene <- as.factor(df2$Oncogene)
    #  df2$CNV <- as.factor(df2$CNV)
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
    df2$Oncogene <- NULL
    names(df2)[names(df2) == "CNV"] <- EIF
    # df2 <- na.omit(df2)
    return(df2)
  }
  df3 <- lapply(c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1", "MYC"), Onco.CNV)
  ## use cbind to convert df3 list into a data frame
  df4 <- do.call(cbind.data.frame, df3)
  ## remove duplicated columns (TCGAstudy)
  df4 <- df4[, !duplicated(colnames(df4))]
  # cnt <- with(df4, table(EIF4A1, EIF4G1))
  cnt <- with(df4, table(df4[,x], df4[,y]))
  cnt <- as.data.frame(cnt)
  cnt$radius <- sqrt(cnt$Freq / pi)
  p <- ggplot(cnt, aes(x = Var1, y = Var2)) 
  p + geom_point(aes(size   = radius),
                     shape  = 21, 
                     colour = "black", 
                     fill   = "skyblue")+
      geom_text(aes(label = Freq),size = 4)+
      theme(panel.background = element_blank(), 
            panel.border     = element_rect(colour = "blue", 
                                            fill   = NA, 
                                            size   = 1),
            legend.position  = "none")+
      scale_size_area(max_size = 20)+
      scale_x_discrete(limits = c("-2", "-1", "0", "1", "2"), # skip NaN data
                       labels = c(
                                  "-2" = "Homdel",
                                  "-1" = "Hetlos",
                                  "0"  = "Diploid",
                                  "1"  = "Gain",
                                  "2"  = "Amp")) +
      scale_y_discrete(limits = c("-2", "-1", "0", "1", "2"), # skip NaN data
                       labels = c(
                                  "-2" = "Homdel",
                                  "-1" = "Hetlos",
                                  "0"  = "Diploid",
                                  "1"  = "Gain",
                                  "2"  = "Amp")) +
  #Add labels to axes
  labs(x = paste(x, "copy number variation"), 
       y = paste(y, "copy number variation"))
}
plot.bubble.CNV.all.tumor ("EIF4G1", "MYC")

### bugs in function plot.scatter.CNV.all.tumor
library(cowplot)
plot.scatter.CNV.all.tumor <- function(x, y, z){
  Onco.CNV <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.CNV <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_cnaseq",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
    geneticprofile.CNV <- function(x) {
      tcga.pan.geneticprofile[[x]][
        grep("tcga_pan_can_atlas_2018_gistic", 
          tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
    }
    tcga.profiledata.CNV <- function(genename, geneticprofile, caselist) {
      getProfileData(mycgds,
        genename,
        geneticprofile,
        caselist)
    }
    EIF.tcga.CNV <- function(x, y) {
      EIF.tcga.CNV <- tcga.profiledata.CNV(x, geneticprofile.CNV(y), caselist.CNV(y))
      ## change the rownames into the first column
      setDT(EIF.tcga.CNV, keep.rownames = TRUE)[]
      return(EIF.tcga.CNV)
    }
    EIF.tcga.CNV.all <- function(x) {
      test <- lapply(tcga.study.list, EIF.tcga.CNV, x = x)
      df2 <- melt(test)
      colnames(df2) <- c("SampleID","Oncogene","CNV","TCGAstudy")
      df2 <- data.frame(df2)
    }
    df2 <- EIF.tcga.CNV.all(EIF)
    df2$Oncogene <- as.factor(df2$Oncogene)
    # df2$CNV <- as.factor(df2$CNV)
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
    df2$Oncogene <- NULL
    # df2 <- na.omit(df2)
    return(df2)
  }
  CNV.data <- Onco.CNV(z)
  EIF.RNAseq.data <- EIF.RNAseq.all.tumor()
  EIF.CNV.RNAseq <- merge(CNV.data, 
                          EIF.RNAseq.data, 
                          by  = "SampleID", 
                          all = T)
  ggplot(data = EIF.CNV.RNAseq, 
    aes(x    = log2(EIF.CNV.RNAseq[,x]), 
        y    = log2(EIF.CNV.RNAseq[,y]), 
        fill = EIF.CNV.RNAseq[,z])) +
    geom_point(size = 3, 
               alpha = 0.7, 
               shape = 21)+
    labs(x = paste0("log2(", x, "RNAseq)"), 
         y = paste0("log2(", y, "RNAseq)"))
  }
plot.scatter.CNV.all.tumor ("EIF4G1", "EIF4E1","EIF4A1")

########################################################################
## correlation analysis of RNA-seq and mutation from all cancer types ## 
########################################################################
## Get DNFA gene expression from all cancer groups ##
EIF.RNAseq <- function(EIF){
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  caselist.RNAseq <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("mRNA Expression, RSEM", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.tcga.RNAseq <- function(x, y) {
    EIF.tcga.RNAseq <- tcga.profiledata.RNAseq(x, geneticprofile.RNAseq(y),caselist.RNAseq(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.RNAseq, keep.rownames = TRUE)[]
    return(EIF.tcga.RNAseq)
  }
  EIF.RNAseq.tcga.all <- function(x) {
    ## test[] to keep row.names by lapply function 
    test <- lapply(tcga.study.list, EIF.tcga.RNAseq, x = x)
    df2 <- melt(test)
    colnames(df2) <- c("SampleID", "EIFgene", "RNAseq","TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.RNAseq.tcga.all(EIF)
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  return(df2)
}
EIF.RNAseq.data <- EIF.RNAseq("EIF4G1")

## Get oncogene CNV data from SKCM group ##
Onco.Mut <- function(EIF){
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  caselist.Mut <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_sequenced",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.Mut <- function(x) {
    tcga.pan.geneticprofile[[x]][
      grep("tcga_pan_can_atlas_2018_mutations", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
  }
  tcga.profiledata.Mut <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.tcga.Mut <- function(x, y) {
    EIF.tcga.Mut <- tcga.profiledata.Mut(x, geneticprofile.Mut(y), caselist.Mut(y))
    ## change the rownames into the first column
    setDT(EIF.tcga.Mut, keep.rownames = TRUE)[]
    return(EIF.tcga.Mut)
  }
  EIF.tcga.Mut.all <- function(x) {
    test <- lapply(tcga.study.list, EIF.tcga.Mut, x = x)
    df2 <- melt(test)
    drops <- c("variable","value")
    df2 <- df2[ , !(names(df2) %in% drops)]
    colnames(df2) <- c("SampleID",x,"TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.tcga.Mut.all(EIF)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  return (df2)
}
Mut.data <- Onco.Mut("EIF4G1")
Mut.data <- na.omit(Mut.data)
write.csv(Mut.data, file = "EIF4G1Mut.csv")

plot.Mut.RNAseq.all.tumor <- function(ONCO, EIF) {
  EIF.RNAseq <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.RNAseq <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
    geneticprofile.RNAseq <- function(x) {
      tcga.pan.geneticprofile[[x]][
        grep("mRNA Expression, RSEM", 
          tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
    }
    tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
      getProfileData(mycgds,
        genename,
        geneticprofile,
        caselist)
    }
    EIF.tcga.RNAseq <- function(x, y) {
      EIF.tcga.RNAseq <- tcga.profiledata.RNAseq(x, geneticprofile.RNAseq(y),caselist.RNAseq(y))
      ## change the rownames into the first column
      setDT(EIF.tcga.RNAseq, keep.rownames = TRUE)[]
      return(EIF.tcga.RNAseq)
    }
    EIF.RNAseq.tcga.all <- function(x) {
      ## test[] to keep row.names by lapply function 
      test <- lapply(tcga.study.list, EIF.tcga.RNAseq, x = x)
      df2 <- melt(test)
      colnames(df2) <- c("SampleID", "EIFgene", "RNAseq","TCGAstudy")
      df2 <- data.frame(df2)
    }
    df2 <- EIF.RNAseq.tcga.all(EIF)
    df2$EIFgene <- as.factor(df2$EIFgene)
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
    df2 <- na.omit(df2)
    return(df2)
  }
  EIF.RNAseq.data <- EIF.RNAseq(EIF)
  Onco.Mut <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.Mut <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_sequenced",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
    geneticprofile.Mut <- function(x) {
      tcga.pan.geneticprofile[[x]][
        grep("tcga_pan_can_atlas_2018_mutations", 
          tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
    }
    tcga.profiledata.Mut <- function(genename, geneticprofile, caselist) {
      getProfileData(mycgds,
        genename,
        geneticprofile,
        caselist)
    }
    EIF.tcga.Mut <- function(x, y) {
      EIF.tcga.Mut <- tcga.profiledata.Mut(x, geneticprofile.Mut(y), caselist.Mut(y))
      ## change the rownames into the first column
      setDT(EIF.tcga.Mut, keep.rownames = TRUE)[]
      return(EIF.tcga.Mut)
    }
    EIF.tcga.Mut.all <- function(x) {
      test <- lapply(tcga.study.list, EIF.tcga.Mut, x = x)
      df2 <- melt(test)
      drops <- c("variable","value")
      df2 <- df2[ , !(names(df2) %in% drops)]
      colnames(df2) <- c("SampleID",x,"TCGAstudy")
      df2 <- data.frame(df2)
    }
    df2 <- EIF.tcga.Mut.all(EIF)
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
    return (df2)
  }
  Mut.data <- Onco.Mut(ONCO)
  Mut.data$Status <- NULL
  Mut.data$Status <- ifelse(is.na(Mut.data[[ONCO]]), "Wildtype", "Mutated")
  Mut.data <- as.data.frame(Mut.data)
  Mut.data$Status <- as.factor(Mut.data$Status)
  Mut.EIF.RNAseq <- merge(Mut.data, EIF.RNAseq.data, by = "SampleID", all = T)
  Mut.EIF.RNAseq <- Mut.EIF.RNAseq[complete.cases(Mut.EIF.RNAseq$Status), ]
  # levels(Mut.data$Status)
  # na.omit cannot eleminate NaN here!
  Wildtype.number <- 
    nrow(Mut.EIF.RNAseq[Mut.EIF.RNAseq$Status == "Wildtype", ])
  Mut.number <- 
    nrow(Mut.EIF.RNAseq[Mut.EIF.RNAseq$Status == "Mutated", ])
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  print(
    ggplot(Mut.EIF.RNAseq,
      aes(x     = Mut.EIF.RNAseq$Status,
          y     = log2(Mut.EIF.RNAseq$RNAseq),
          color = Mut.EIF.RNAseq$Status)) +    
      geom_violin(trim = FALSE) +
      geom_boxplot(alpha      = .01,
                   width      = .5) +
      labs(x = ONCO,
           y = paste("log2(", EIF, "RNA counts)")) +
      scale_x_discrete(limits = c("Wildtype", "Mutated"), # skip NaN data
                       labels = c("Wildtype" = paste(
                                                     "Wildtype \n n= ", 
                                                      Wildtype.number),
                                  "Mutated"  = paste(
                                                     "Mutated \n n= ", 
                                                      Mut.number)))+
      theme(axis.title      = element_text(face   = "bold",
                                           size   = 9,
                                           color  = "black"),
            axis.text       = element_text(size   = 9,
                                           hjust  = 1,
                                           face   = "bold",
                                           color  = "black"),
            axis.line.x     = element_line(color  = "black"),
            axis.line.y     = element_line(color  = "black"),
            panel.grid      = element_blank(),
            strip.text      = element_text(face   = "bold",
                                           size   = 9,
                                           colour = "black"),
                                           legend.position = "none"))
}
plot.Mut.RNAseq.all.tumor("EIF4E", "EIF4E")

########################################################################
##  Kaplan-Meier curve with clinic and mutation data from all tumors  ##
########################################################################
plot.km.mut.all.tumor <- function(ONCO) {
  #### retrieve mutation data from all tcga groups ####
  Onco.Mut <- function(EIF){
    tcga.pan.studies <- getCancerStudies(mycgds)[
      grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
    tcga.study.list <- tcga.pan.studies$cancer_study_id
    names(tcga.study.list) <- tcga.study.list
    caselist <- function(x) getCaseLists(mycgds, x)
    geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
    tcga.pan.caselist <- lapply(tcga.study.list, caselist)
    tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
    caselist.Mut <- function(x) {
      tcga.pan.caselist[[x]][
        grep("tcga_pan_can_atlas_2018_sequenced",
          tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
    }
    geneticprofile.Mut <- function(x) {
      tcga.pan.geneticprofile[[x]][
        grep("tcga_pan_can_atlas_2018_mutations", 
          tcga.pan.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
    }
    tcga.profiledata.Mut <- function(genename, geneticprofile, caselist) {
      getProfileData(mycgds,
        genename,
        geneticprofile,
        caselist)
    }
    EIF.tcga.Mut <- function(x, y) {
      EIF.tcga.Mut <- tcga.profiledata.Mut(x, geneticprofile.Mut(y), caselist.Mut(y))
      ## change the rownames into the first column
      setDT(EIF.tcga.Mut, keep.rownames = TRUE)[]
      return(EIF.tcga.Mut)
    }
    EIF.tcga.Mut.all <- function(x) {
      test <- lapply(tcga.study.list, EIF.tcga.Mut, x = x)
      df2 <- melt(test)
      drops <- c("variable","value")
      df2 <- df2[ , !(names(df2) %in% drops)]
      colnames(df2) <- c("case.id",x,"TCGAstudy")
      df2 <- data.frame(df2)
    }
    df2 <- EIF.tcga.Mut.all(EIF)
    df2$TCGAstudy <- as.factor(df2$TCGAstudy)
    return (df2)
  }
  Mut.data <- Onco.Mut(ONCO)
  Mut.data$Status <- NULL
  Mut.data$Status <- ifelse(is.na(Mut.data[[ONCO]]), "Wildtype", "Mutated")
  Mut.data <- as.data.frame(Mut.data)
  Mut.data$Status <- as.factor(Mut.data$Status)
  #### retrieve clinic data from all tcga groups ####
  tcga.clinic.data <- function(x) {
    print(x)
    url <- function(x){
      url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="
      url <- paste0(url, x, "_all")
      return(url)
    }
    # testurl <- url("acc_tcga")
    # tesereq <- GET(url("acc_tcga"))
    req <- function(x) {GET(url(x))}
    # req <- req("acc_tcga")
    clinical_data <- function(x) {content(req(x),
      type      = 'text/tab-separated-values',
      col_names = T,
      col_types = NULL)}
    data <- clinical_data(x)
    data <- data[c("OS_MONTHS",
      "OS_STATUS",
      "CASE_ID")]
  }
  pro.tcga.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "pro.tcga.study.list" contains all the tcga provisional cancer studies
  pro.tcga.study.list <- pro.tcga.studies$cancer_study_id
  names(pro.tcga.study.list) <- pro.tcga.study.list
  # three datasets donot have OS data and cause bugs remove them
  bug.data.set <- names(pro.tcga.study.list) %in% c("meso_tcga", 
                                                    "pcpg_tcga", 
                                                    "ucs_tcga")
  pro.tcga.study.list <- pro.tcga.study.list[!bug.data.set]
  all.tcga.clinic.data <- lapply(pro.tcga.study.list, tcga.clinic.data)
  all.tcga.clinic.data <- melt(all.tcga.clinic.data)
  all.tcga.clinic.data <- all.tcga.clinic.data[c("OS_STATUS",
                                                 "CASE_ID",
                                                 "value",
                                                 "L1")]
  colnames(all.tcga.clinic.data) <- c("OS_STATUS",
                                      "case.id",
                                      "OS_MONTHS",
                                      "TCGAstudy")
  all.tcga.clinic.data$case.id <- str_replace_all(all.tcga.clinic.data$case.id,
                                                  '-',
                                                  '.')
  message("clinical data retrieved")
  #### combine clinic data and mutation data from all tcga groups ####
  df <- join_all(list(all.tcga.clinic.data[c("OS_MONTHS",
                                             "OS_STATUS",
                                             "case.id")],
                      Mut.data[c("case.id",
                                 "TCGAstudy",
                                 "Status")]),
                      by   = "case.id")
  df <- merge(all.tcga.clinic.data[c("OS_MONTHS","OS_STATUS","case.id")],
              Mut.data[c("case.id", "TCGAstudy", "Status")], 
              by   = "case.id",  
              all  = T)
  df <- na.omit(df)
  message("clinical and Mutation data combined")
  #### print km plot ####
  Wildtype.number <- 
    nrow(df[df$Status == "Wildtype", ])
  Mut.number <- 
    nrow(df[df$Status == "Mutated", ])
  
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  km <- survfit(SurvObj ~ df$Status, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Status, data = df, rho = 0)
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  print(
    ggplot2::autoplot(
      km,
      xlab = "Months",
      ylab = "Survival Probability",
      main = paste("Kaplan-Meier plot", ONCO, "Mutations")) +
      theme(axis.title           = black.bold.12pt,
            axis.text            = black.bold.12pt,
            axis.line.x          = element_line(color  = "black"),
            axis.line.y          = element_line(color  = "black"),
            panel.grid           = element_blank(),
            strip.text           = black.bold.12pt,
            legend.text          = black.bold.12pt ,
            legend.title         = black.bold.12pt ,
            legend.justification = c(1,1),
            legend.position      = c(1,1))+
      guides(fill = FALSE) +
      scale_color_manual(values = c("red", "blue"),
                         name   = paste(ONCO, "Mutations"),
                         breaks = c("Wildtype", "Mutated"),
                         labels = c(paste("Wildtype n = ", Wildtype.number),
                                    paste("Mutated n = ", Mut.number))) +
      geom_point(size = 0.25) +
      annotate("text",
                x        = 400,
                y        = 0.70,
                label    = paste("log-rank test, p.val = ", p.val),
                size     = 4.5,
                hjust    = 1,
                fontface = "bold"))
  }

plot.km.mut.all.tumor ("EGFR")
lapply(ONCO, plot.km.mut.all.tumor)












#####################################################################
get.EIF.RNAseq.tcga <- function(x) {
  EIF.gene <- c("EIF4E","EIF4G1","EIF4EBP2","RPS6KB1")
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list,
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  # for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pro.caselist[[1]][
  # grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[1]]$case_list_id),
  # ][1,1]
  # b <- tcga.pro.geneticprofile[[1]][
  # grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  # tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  # how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.RNAseq <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      # double backslash \\ suppress the special meaning of ( )
      # in regular expression
      grep("mRNA Expression, RSEM", 
        tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  # test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  # caselist.RNAseq = caselist.RNAseq ('acc_tcga')
  # geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
  # Wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  # within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
      genename,
      geneticprofile,
      caselist)
  }
  EIF.RNAseq.tcga <- tcga.profiledata.RNAseq(EIF.gene,
                                             geneticprofile.RNAseq(x),
                                             caselist.RNAseq(x))
  EIF.RNAseq.tcga$rn <- rownames(EIF.RNAseq.tcga)
  return(EIF.RNAseq.tcga)
  }

get.EIF.score.tcga <- function(x){
  EIF.RNAseq.tcga <- get.EIF.RNAseq.tcga(x) 
  EIF.score.tcga <- EIF.RNAseq.tcga
  EIF.score.tcga$EIF4Escore <- EIF.RNAseq.tcga$EIF4E/EIF.RNAseq.tcga$EIF4E
  EIF.score.tcga$EIF4G1score <- EIF.RNAseq.tcga$EIF4G1/EIF.RNAseq.tcga$EIF4E
  EIF.score.tcga$EIF4EBP2score <- EIF.RNAseq.tcga$EIF4EBP2/EIF.RNAseq.tcga$EIF4E
#  EIF.score.tcga$RPS6KB1score <- EIF.RNAseq.tcga$RPS6KB1/EIF.RNAseq.tcga$EIF4E
  EIF.score.tcga <- EIF.score.tcga [, 5:8]
  return(EIF.score.tcga)
  }
### EIF.score.tcga$GEBPscore <- EIF.RNAseq.tcga$EIF4G1/EIF.RNAseq.tcga$EIF4EBP1

plot.EIF.RNAseq.score <- function (x) {
  EIF.RNAseq.tcga <- get.EIF.RNAseq.tcga(x)
  EIF.score.tcga <- get.EIF.score.tcga(x)
  par(mfrow=c(1,2))
  boxplot(log2(EIF.RNAseq.tcga[, 
                             c("EIF4E", "EIF4G1", 
                               "EIF4EBP2", "RPS6KB1")]),
          main= paste0("EIF RNAseq counts in ", x),
          las = 2)
  boxplot(log2(EIF.score.tcga[,
                        c("EIF4Escore","EIF4G1score",
                          "EIF4EBP2score")]),
          main= paste0("EIF scores in ", x),
          las = 2)
  }
lapply(tcga.study.list, plot.EIF.RNAseq.score)

### to be tested!
plot.EIF.score.all.tcga <- function(x) {
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  ### "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  EIF.score.tcga <- lapply(tcga.study.list, get.EIF.score.tcga)
  EIF.score.tcga.all.tumors <- melt(EIF.score.tcga)
  colnames(EIF.score.tcga.all.tumors) <- c("rn", "EIFgene", "Score", "TCGAstudy")
  EIF.score.tcga.all.tumors <- data.frame(EIF.score.tcga.all.tumors)
  EIF.score.tcga.all.tumors$EIFgene <- as.factor(EIF.score.tcga.all.tumors$EIFgene)
  EIF.score.tcga.all.tumors$TCGAstudy <- as.factor(EIF.score.tcga.all.tumors$TCGAstudy)
  ### to be tested!  
#  EIF.score.tcga.all.tumors <- EIF.score.tcga.all.tumors[EIF.score.tcga.all.tumors$EIFgene == x,]
#  mean <- within(EIF.score.tcga.all.tumors, TCGAstudy <- reorder(TCGAstudy, log2(x), median))
  y <- paste0(x, "score")
  median <- within(EIF.score.tcga.all.tumors[
    EIF.score.tcga.all.tumors$EIFgene == y,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(Score), median))
  print(
    ggplot(median,
           aes(x     = TCGAstudy,
               y     = log2(Score),
               color = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(EIF4G1/EIF4EBP1 ratio)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       angle  = 45,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
  
  }
plot.EIF.score.all.tcga("EIF4E")

##########################################################
##########################################################
plotEIF <-  function (x) {
#  name <- deparse(substitute(x))
  black_bold_tahoma_12 <- element_text(color  = "black", 
                                       face   = "bold",
                                       family = "Tahoma", 
                                       size   = 12)
  
  black_bold_tahoma_12_45 <- element_text(color  = "black",
                                          face   = "bold",
                                          family = "Tahoma", 
                                          size   = 12, 
                                          angle  = 45,
                                          hjust  = 1)
  ggplot(x,
         aes(x = variable,
             y = log2(value))) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .75,
                 position = position_dodge(width = .9)) +
    #    labs(title = paste0(name," n = 8555"),
    #         x     = "eIF4F complex components",
    #         y     = paste0("log2(value)")) +
    theme_bw() +
    theme(plot.title  = black_bold_tahoma_12,
          axis.title  = black_bold_tahoma_12,
          axis.text.x = black_bold_tahoma_12_45,
          axis.text.y = black_bold_tahoma_12,
          axis.line.x = element_line(color  = "black"),
          axis.line.y = element_line(color  = "black"),
          panel.grid  = element_blank(),
          strip.text  = black_bold_tahoma_12,
          legend.position = "none")
}

plot.EIFandScore.all.tumors <- function (){
  EIF.RNAseq.tcga <- lapply(tcga.study.list, get.EIF.RNAseq.tcga)
  EIF.score.tcga <- lapply(tcga.study.list, get.EIF.score.tcga)
  EIF.RNAseq.tcga.all.tumors <- melt(EIF.RNAseq.tcga)
#  x1  = factor(x, levels=c("B", "C", "A"))
#  levels(EIF.RNAseq.tcga.all.tumors$variable)
  EIF.RNAseq.tcga.all.tumors$variable <- ordered(EIF.RNAseq.tcga.all.tumors$variable, 
                                                 levels = c("EIF4E","EIF4G1",
                                                            "EIF4EBP2","RPS6KB1"))
  number <- nrow(EIF.RNAseq.tcga.all.tumors)/4
  EIF.score.tcga.all.tumors <- melt(EIF.score.tcga)
  my_comparison1 <- list( c("EIF4E", "EIF4G1"), 
                          c("EIF4G1", "EIF4EBP2"), 
                          c("EIF4E", "EIF4EBP2"),
                          c("EIF4E", "RPS6KB2"),
                          c("EIF4EBP1", "RPS6KB1"))
  my_comparison2 <- list( c("EIF4Escore", "EIF4G1score"), 
                          c("EIF4G1score", "EIF4EBP2score"), 
                          c("EIF4Escore", "EIF4EBP2score"),
                          c("EIF4Escore", "RPS6KB1score"),
                          c("EIF4EBP2score", "RPS6KB1score"))
  p1 <- plotEIF(EIF.RNAseq.tcga.all.tumors) +
    labs(title = paste0("All tumors n = ", number),
         x     = "eIF4F subunit RNAseq",
         y     = paste0("log2(value)")) +
    stat_compare_means(comparisons = my_comparison1, method = "t.test")
  p1$layers[[2]]$aes_params$textsize <- 5  
  p2 <- plotEIF(EIF.score.tcga.all.tumors) + 
    labs(title = paste0("All tumors n = ", number),
         x     = "eIF4E ratio score",
         y     = paste0("log2(value)")) +
    stat_compare_means(comparisons = my_comparison2, method = "t.test")
   p2$layers[[2]]$aes_params$textsize <- 5
  grid.arrange(p1, p2, ncol=2)
}
plot.EIFandScore.all.tumors()




