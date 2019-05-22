if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")
BiocManager::install("GenVisR")

library(maftools)
library(GenVisR)

#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
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
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10)
