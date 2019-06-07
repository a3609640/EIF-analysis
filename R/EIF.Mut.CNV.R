if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("IsoformSwitchAnalyzeR")
BiocManager::install("maftools")
BiocManager::install("GenVisR")
BiocManager::install("graphics")
# devtools::install_github(repo = "PoisonAlien/TCGAmutations")
# devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")




library(GenVisR)
library(graphics)
library(IsoformSwitchAnalyzeR)
library(maftools)
library(TCGAmutations)
library(TCGAbiolinks)
library(dplyr)
library(DT)


vignette(package = "TCGAmutations", topic = "Introduction")
TCGA.all <- TCGAmutations::tcga_available()
TCGA.all <- TCGA.all$Study_Abbreviation
TCGA.all <- as.vector(TCGA.all)
# TCGAmutations::tcga_load(study = "SKCM")
# load MAF files for each tumor group
sapply(TCGA.all, TCGAmutations::tcga_load)
mafs.all <- c(tcga_acc_mc3, tcga_blca_mc3, tcga_brca_mc3,tcga_cesc_mc3,
              tcga_chol_mc3, tcga_coad_mc3, tcga_dlbc_mc3, tcga_esca_mc3,
              tcga_gbm_mc3,tcga_hnsc_mc3,tcga_kich_mc3,tcga_kirc_mc3,
              tcga_kirp_mc3,tcga_laml_mc3,tcga_lgg_mc3,tcga_lihc_mc3,
              tcga_luad_mc3, tcga_lusc_mc3,tcga_meso_mc3,tcga_ov_mc3,
              tcga_paad_mc3,tcga_pcpg_mc3,tcga_prad_mc3,tcga_read_mc3,
              tcga_sarc_mc3,tcga_skcm_mc3,tcga_stad_mc3,tcga_tgct_mc3,
              tcga_thca_mc3,tcga_thym_mc3,tcga_ucec_mc3,tcga_ucs_mc3,
              tcga_uvm_mc3)
all.cancer <- maftools::merge_mafs(mafs.all)
getSampleSummary(x = all.cancer)
# Shows gene summary
getGeneSummary(x = all.cancer)
# Clinical data; printing only first ten columns for display convenience
all.cancer.clin <- getClinicalData(x = all.cancer)[1:10, 1:10]
plotmafSummary(maf       = all.cancer, 
               rmOutlier = TRUE, 
               addStat   = 'median', 
               dashboard = TRUE, 
               titvRaw   = FALSE)

lollipopPlot(maf = all.cancer, 
             gene = 'EIF4A1',  
             labelPos = 'all',
             labPosAngle = 90)
lollipopPlot(maf = all.cancer, 
             gene = 'EIF4E',  
             labelPos = 'all',
             labPosAngle = 90)
lollipopPlot(maf = all.cancer, 
             gene = 'EIF4G1',  
             labelPos = 'all',
             labPosAngle = 90)
lollipopPlot(maf = all.cancer, 
             gene = 'EIF4EBP1',  
             labelPos = 'all',
             labPosAngle = 90)

oncoplot(maf = all.cancer, genes = c('EIF4E','EIF4G1','EIF4A1','EIF4EBP1','MYC'))

## drug interaction
eIF4F.dgi = drugInteractions(genes = "EIF4E", drugs = TRUE)
eIF4F.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
all.cancer.titv = titv(maf = all.cancer, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = all.cancer.titv)





laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
laml = read.maf(maf = laml.maf)
#Typing laml shows basic summary of MAF file.
laml
#Shows sample summry.
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#Shows all fields in MAF
getFields(laml)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')
plotmafSummary(maf = laml, rmOutlier = TRUE, 
  addStat = 'median', dashboard = TRUE, 
  titvRaw = FALSE)
#oncoplot for top ten mutated genes.
oncoplot(maf = skcm, genes = c('EIF4E','EIF4G1','EIF4A1','EIF4EBP1'))
oncostrip(maf = skcm, genes = c('EIF4E','EIF4G1','EIF4A1','EIF4EBP1'))
laml.titv = titv(maf = skcm, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(maf = skcm, gene = 'EIF4G1', showMutationRate = TRUE)

mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)

PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")

mafSurvival(maf    = laml, 
  genes  = 'DNMT3A', 
  time   = 'days_to_last_followup', 
  Status = 'vital_status', 
  isTCGA = TRUE)

laml.plus.gistic = read.maf(
  maf = tcga_laml_mc3,
  gisticAllLesionsFile = all.lesions,
  gisticAmpGenesFile = amp.genes,
  gisticDelGenesFile = del.genes,
  gisticScoresFile = scores.gis,
  isTCGA = TRUE,
  verbose = FALSE
)
oncoplot(maf = laml.plus.gistic, top = 10)



waterfall(brcaMAF, fileType="MAF")
# Load GenVisR and set seed
library(GenVisR)
set.seed(383)

# Plot only genes with mutations in 6% or more of samples
waterfall(brcaMAF, mainRecurCutoff = 0.06)

# Plot only the specified genes
waterfall(brcaMAF, plotGenes = c("PIK3CA", "TP53", "USH2A", "MLL3", "BRCA1"))

# Create clinical data
subtype <- c("lumA", "lumB", "her2", "basal", "normal")
subtype <- sample(subtype, 50, replace = TRUE)
age <- c("20-30", "31-50", "51-60", "61+")
age <- sample(age, 50, replace = TRUE)
sample <- as.character(unique(brcaMAF$Tumor_Sample_Barcode))
clinical <- as.data.frame(cbind(sample, subtype, age))

# Melt the clinical data into 'long' format.
library(reshape2)
clinical <- melt(clinical, id.vars = c("sample"))

# Run waterfall
waterfall(brcaMAF, 
          clinDat      = clinical, 
          clinVarCol   = c(lumA    = "blue4", 
                           lumB    = "deepskyblue", 
                           her2    = "hotpink2", 
                           basal   = "firebrick2", 
                           normal  = "green4", 
                           `20-30` = "#ddd1e7", 
                           `31-50` = "#bba3d0", 
                           `51-60` = "#9975b9", 
                           `61+`   = "#7647a2"), 
          plotGenes    = c("PIK3CA", "TP53", "USH2A", "MLL3", "BRCA1"), 
          clinLegCol   = 2, 
          clinVarOrder = c("lumA", "lumB", "her2", 
                           "basal", "normal", "20-30", 
                           "31-50", "51-60", "61+"))
# Create input data
data <- brcaMAF[brcaMAF$Hugo_Symbol == "TP53", c("Hugo_Symbol", "amino_acid_change_WU")]
data <- as.data.frame(cbind(data, "ENST00000269305"))
colnames(data) <- c("gene", "amino_acid_change", "transcript_name")

# Call lolliplot
lolliplot(data)






data("exampleSwitchListAnalyzed")
extractSwitchSummary(
  exampleSwitchListAnalyzed,
  filterForConsequences = TRUE
) 
subset(
  extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n=10,
    inEachComparison = TRUE
  )[,c('gene_name','condition_2','gene_switch_q_value','Rank')],
  gene_name == 'ZAK'
)
switchPlot(
  exampleSwitchListAnalyzed,
  gene='ZAK',
  condition1 = 'COAD_ctrl',
  condition2 = 'COAD_cancer'
)
extractSwitchOverlap(
  exampleSwitchListAnalyzed,
  filterForConsequences=TRUE
)

extractConsequenceSummary(
  exampleSwitchListAnalyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
