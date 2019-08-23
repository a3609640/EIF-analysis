library(dplyr)
library(DT)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v75)
library(EnvStats)
library(GenVisR)
library(ggplot2)
library(ggpubr)
library(graphics)
library(Gviz)
library(IsoformSwitchAnalyzeR)
library(maftools)
library(psichomics)
library(reshape2)
library(TCGAmutations)
library(TCGAbiolinks)


#############################
### TCGAmutations package ###
#############################
vignette(package = "TCGAmutations", topic = "Introduction")

EIF.gene <- c('EIF4E','EIF4G1','EIF4A1','EIF4EBP1','MYC')
names(EIF.gene) <- EIF.gene

TCGA.all <- TCGAmutations::tcga_available()
TCGA.all <- TCGA.all$Study_Abbreviation
TCGA.all <- as.vector(TCGA.all)
# TCGAmutations::tcga_load(study = "SKCM")
# load MAF files for all tumor groups
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

plot.mut.each.EIF <- function(x, y){
  getSampleSummary(x = x)
  getGeneSummary(x = x)
  plotmafSummary(maf       = x, 
                 rmOutlier = TRUE, 
                 addStat   = 'median', 
                 dashboard = TRUE, 
                 titvRaw   = FALSE)
  lollipopPlot(maf         = x, 
               gene        = y,  
               labelPos    = 'all',
               labPosAngle = 90)
  }
plot.mut.each.EIF(x = tcga_laml_mc3, y = 'MYC')
lapply(EIF.gene, plot.mut.each.EIF, x = all.cancer)

oncoplot(maf = all.cancer, genes = EIF.gene)


## drug interaction
eIF4F.dgi = drugInteractions(genes = "EIF4E", drugs = TRUE)
eIF4F.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
all.cancer.titv = titv(maf = all.cancer, plot = FALSE, useSyn = TRUE)
# plot titv summary
plotTiTv(res = all.cancer.titv)



##########################
### psichomics package ###
##########################
## Downloading and loading TCGA data
# Available tumour types
cohorts <- getFirebrowseCohorts()
# Available sample dates
date <- getFirebrowseDates()
# Available data types
dataTypes <- getFirebrowseDataTypes()



edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "8") 
gnm <- GRanges("8:128746765-128754082") # for MYC

edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "4") 
gnm <- GRanges("4:99788886-99854246") # for EIF4E


## Since we're using Ensembl chromosome names we have to set:
options(ucscChromosomeNames = FALSE)

## Define a genome axis track
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx, filter = GRangesFilter(gnm))
gtx <- GeneRegionTrack(gnm_gns, 
                       name       = "tx", 
                       geneSymbol = TRUE,
                       showId     = TRUE)

## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm)
## plot the region
plotTracks(list(ht))


#####################################################################
### map the translation start sites of MYC to genomic coordinates ###
#####################################################################
edb <- EnsDb.Hsapiens.v75
gns <- genes(edb, filter = ~ protein_domain_id == "PF01056" & seq_name == "8")
# gns

## Change chromosome naming style to UCSC
seqlevelsStyle(edb) <- "UCSC"
txs <- getGeneRegionTrackForGviz(
  edb, 
  filter = ~ genename == "MYC" & protein_domain_id == "PF01056")

pdoms <- proteins(
  edb, 
  filter  = ~ tx_id %in% txs$transcript & protein_domain_source == "pfam",
  columns = c("protein_domain_id", 
              "prot_dom_start",
              "prot_dom_end"))
# pdoms

pdoms_rng <- IRanges(start = pdoms$prot_dom_start, 
                     end   = pdoms$prot_dom_end,
                     names = pdoms$protein_id)
pdoms_gnm <- proteinToGenome(pdoms_rng, edb)

# pdoms_gnm

## Convert the list to a GRanges with grouping information
pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))
pdoms_gnm_grng$id <- rep(pdoms$protein_domain_id, lengths(pdoms_gnm))
pdoms_gnm_grng$grp <- rep(1:nrow(pdoms), lengths(pdoms_gnm))

# pdoms_gnm_grng

## Define the individual tracks:
## - Ideagram
ideo_track <- IdeogramTrack(genome = "hg19", chromosome = "chr8")
## - Genome axis
gaxis_track <- GenomeAxisTrack()
## - Transcripts
gene_track <- GeneRegionTrack(
  txs, 
  showId = TRUE, 
  just.group = "right",
  name = "Transcripts", 
  geneSymbol = TRUE, 
  # background.panel = "#FFFEDB",
  background.title = "brown",
  size = 0.5)
## - Protein domains
pdom_track <- AnnotationTrack(
  pdoms_gnm_grng, 
  group = pdoms_gnm_grng$grp,
  id = pdoms_gnm_grng$id, 
  groupAnnotation = "id",
  just.group = "right", 
  shape = "box",
  name = "Protein domains", 
  # background.panel = "#FFFEDB",
  background.title = "darkblue",
  size = 0.5)

# highlight
ht <- HighlightTrack(trackList  = list(gene_track),
                     start      = c(128747765,128748315,128750494), 
                     end        = c(128748082,128748869,128751265), 
                     chromosome = 8)
## Generate the plot  
plotTracks(list(gaxis_track, ht, pdom_track))



#####################################################################
### map the translation start sites of EIF4E to genomic coordinates ###
#####################################################################
edb <- EnsDb.Hsapiens.v75
gns <- genes(edb, filter = ~ protein_domain_id == "PF01652" & seq_name == "4")
# gns

## Change chromosome naming style to UCSC
seqlevelsStyle(edb) <- "UCSC"
txs <- getGeneRegionTrackForGviz(
  edb, 
  filter = ~ genename == "EIF4E" & protein_domain_id == "PF01652")

pdoms <- proteins(
  edb, 
  filter  = ~ tx_id %in% txs$transcript & protein_domain_source == "pfam",
  columns = c("protein_domain_id", 
    "prot_dom_start",
    "prot_dom_end"))
# pdoms

pdoms_rng <- IRanges(start = pdoms$prot_dom_start, 
  end   = pdoms$prot_dom_end,
  names = pdoms$protein_id)
pdoms_gnm <- proteinToGenome(pdoms_rng, edb)

# pdoms_gnm

## Convert the list to a GRanges with grouping information
pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))
pdoms_gnm_grng$id <- rep(pdoms$protein_domain_id, lengths(pdoms_gnm))
pdoms_gnm_grng$grp <- rep(1:nrow(pdoms), lengths(pdoms_gnm))

# pdoms_gnm_grng

## Define the individual tracks:
## - Ideagram
ideo_track <- IdeogramTrack(genome = "hg19", chromosome = "chr4")
## - Genome axis
gaxis_track <- GenomeAxisTrack()
## - Transcripts
gene_track <- GeneRegionTrack(
  txs, 
  showId = TRUE, 
  just.group = "right",
  name = "Transcripts", 
  geneSymbol = TRUE, 
  # background.panel = "#FFFEDB",
  background.title = "brown",
  size = 0.5)
## - Protein domains
pdom_track <- AnnotationTrack(
  pdoms_gnm_grng, 
  group = pdoms_gnm_grng$grp,
  id = pdoms_gnm_grng$id, 
  groupAnnotation = "id",
  just.group = "right", 
  shape = "box",
  name = "Protein domains", 
  # background.panel = "#FFFEDB",
  background.title = "darkblue",
  size = 0.5)

# highlight
ht <- HighlightTrack(trackList = list(gene_track),
  start = c(99850246,99850046,99823027,99812388,99809040,99807502,99806073,99799607,99792836), 
  end = c(99851786,99850243,99823133,99812483,99809103,99808343,99806212,99802293,99795886), 
  chromosome = 4)
## Generate the plot  
plotTracks(list(gaxis_track, ht, pdom_track), reverseStrand = TRUE)
plotTracks(list(gaxis_track, gene_track, pdom_track), reverseStrand = TRUE)



##################################
### Exon analysis on Xena data ###
##################################
plot.EIF.exon.all <- function() {
  Exon.data <- read.table(file = "~/Downloads/denseDataOnlyDownload.tsv", 
                          sep = '\t', 
                          header = TRUE)
  Exon.data.long <- melt(Exon.data)
  Exon.data.long <- subset(Exon.data.long, grepl("^chr8", Exon.data.long$variable))

  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
    #"Recurrent Tumor",
    "Solid Tissue Normal")
  Exon.data.long <- Exon.data.long[
    Exon.data.long$sample_type %in% tumor.type, ]
  Exon.data.long <- droplevels(Exon.data.long)
  Exon.data.long <- na.omit(Exon.data.long)
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  print(
    ggplot(data = Exon.data.long,
      #use[ ,genemutations] not $genemutations for variable in a function
      aes(x     = Exon.data.long$sample_type,
          y     = value,
          color = Exon.data.long$sample_type)) +
      facet_grid(~ variable,
                 scales = "free",
                 space  = "free") +
      facet_wrap(~ variable, 
                 ncol = 3) +
      geom_boxplot() +
      labs(x = "All TCGA Sample Type",
           y = paste("log2(RPKM + 1)")) +
      stat_n_text() + 
      theme_bw() +
      theme(
        plot.title      = black_bold_tahoma_12,
        axis.title      = black_bold_tahoma_12,
        axis.text.x     = black_bold_tahoma_12_45,
        axis.text.y     = black_bold_tahoma_12,
        axis.line.x     = element_line(color = "black"),
        axis.line.y     = element_line(color = "black"),
        panel.grid      = element_blank(),
        legend.position = "none",
        strip.text      = black_bold_tahoma_12) +
      ggpubr::stat_compare_means(
        comparisons = list(
          c("Metastatic", "Solid Tissue Normal"),
          c("Primary Tumor", "Solid Tissue Normal"),
          c("Metastatic", "Primary Tumor")), 
        method = "t.test", 
        label = "p.signif"))
  }
plot.EIF.exon.all ()

plot.EIF.exon <- function(x) {
  Exon.data <- read.table(file = "~/Downloads/denseDataOnlyDownload.tsv", 
                          sep = '\t', 
                          header = TRUE)
  Exon.data.long <- melt(Exon.data)
  Exon.data.long <- Exon.data.long[
    Exon.data.long$cancer.type.abbreviation == x, ]
  Exon.data.long <- subset(Exon.data.long, grepl("^chr4", Exon.data.long$variable))

  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
    #"Recurrent Tumor",
    "Solid Tissue Normal")
  Exon.data.long <- Exon.data.long[
    Exon.data.long$sample_type %in% tumor.type, ]
  Exon.data.long <- droplevels(Exon.data.long)
  Exon.data.long <- na.omit(Exon.data.long)
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  print(
    ggplot(data = Exon.data.long,
      #use[ ,genemutations] not $genemutations for variable in a function
      aes(x     = Exon.data.long$sample_type,
          y     = value,
          color = Exon.data.long$sample_type)) +
      stat_n_text() + 
      facet_grid(~ variable,
        scales = "free",
        space  = "free") +
      facet_wrap(~ variable, 
                 ncol = 3) +
      geom_boxplot() +
      labs(x = paste(x, " Sample Type"),
           y = paste("log2(RPKM + 1)")) +
      theme_bw() +
      theme(
        plot.title      = black_bold_tahoma_12,
        axis.title      = black_bold_tahoma_12,
        axis.text.x     = black_bold_tahoma_12_45,
        axis.text.y     = black_bold_tahoma_12,
        axis.line.x     = element_line(color = "black"),
        axis.line.y     = element_line(color = "black"),
        panel.grid      = element_blank(),
        legend.position = "none",
        strip.text      = black_bold_tahoma_12) +
      ggpubr::stat_compare_means(
        comparisons = list(
          c("Metastatic", "Solid Tissue Normal"),
          c("Primary Tumor", "Solid Tissue Normal"),
          c("Metastatic", "Primary Tumor")), 
        method = "t.test", 
        label = "p.signif"))
  }
plot.EIF.exon ("BRCA")








#####################################
### IsoformSwitchAnalyzeR package ###
#####################################
## acquire the genomic sequence to predict ORFs. 
library(BSgenome.Hsapiens.UCSC.hg19)
data("exampleSwitchList")
exampleSwitchList
exampleSwitchList <- isoformSwitchAnalysisPart1(
  input                   = exampleSwitchList, 
  genomeObject            = Hsapiens, 
  dIFcutoff               = 0.4,         # Set high for short runtime in example 
  dataoutputSequences     = FALSE, # keeps the function from outputting the fasta files from this 
  examplecalibratePvalues = FALSE
  )
switchPlot(exampleSwitchList, gene='HOXC13', condition1='Ctrl', condition2='KD1')
extractSwitchSummary(exampleSwitchList)

aSwitchList <- importGTF(pathToGTF = "Downloads/hg_ucsc.gtf")
aSwitchList
head(aSwitchList,2)
switchPlotTranscript(aSwitchList, gene = 'ENST00000488147.1')



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
  condition1 = 'LUAD_ctrl',
  condition2 = 'LUAD_cancer'
)
switchPlot(exampleSwitchList, gene='SRRM1')
