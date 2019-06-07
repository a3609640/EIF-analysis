BiocManager::install("RforProteomics")
library("ggplot2")  ## Convenient and nice plotting
library("RColorBrewer") ## Color palettes
library("RforProteomics")
library("reshape2") ## Flexibly reshape data



pp <- proteomicsPackages()
display(pp)
vignette(package = "RforProteomics")
