BiocManager::install("RforProteomics", dependencies = TRUE)
library("RforProteomics")
pp <- proteomicsPackages()
display(pp)
