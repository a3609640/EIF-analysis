BiocManager::install("RforProteomics")
library("ggplot2")  ## Convenient and nice plotting
library("mzR")
library("RColorBrewer") ## Color palettes
library("RforProteomics")
library("reshape2") ## Flexibly reshape data
library("rpx")

###########################
## Importing experiments ##
###########################

# MSnbase is able to import raw MS data stored in XML-based formats, mzXML,
# mzData and mzML
file <- dir(system.file(package = "MSnbase", dir = "extdata"),
  full.names = TRUE, pattern = "mzXML$")
rawdata <- readMSData(file, msLevel = 2, verbose = FALSE)

library("MSnbase")
itraqdata
head(fData(itraqdata))

#####################
## Spectra objects ##
#####################
# The raw data is composed of the 55 MS spectra. The spectra are named individually (X1, X10, X11, X12, X13, X14, …) and stored in a environment
spectra(itraqdata)
sp <- itraqdata[["X1"]]
sp
peaksCount(sp)
head(peaksCount(itraqdata))
rtime(sp)
head(rtime(itraqdata))

###################
## Reporter ions ##
###################
# ReporterIons instances are required to quantify reporter peaks in MSnExp experiments 
iTRAQ4
TMT10

##########################
## Chromatogram objects ##
##########################










###################
## MS data space ##
###################

# a list of recent PX additions and updates
pxannounced()
# Pharmacoproteomic characterisation of human colon and rectal cancer - CPTAC Full Proteomes
# rpx package provids access to the ProteomeXchange (PX) central repository
px <- PXDataset("PXD005354")
px <- PXDataset("PXD000001")
px
pxtax(px)
pxurl(px)
pxref(px)
# All files available for the PX experiment 
pxfiles(px)

fn <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"
# download dataset with pxget function
mzf <- pxget(px, fn)
mzf
## reads the data
ms <- openMSfile(mzf)
ms
hd <- header(ms)
dim(hd)
names(hd)
hd[1000, ]
head(peaks(ms, 1000))
plot(peaks(ms, 1000), type = "h")

## a set of spectra of interest: MS1 spectra eluted
## between 30 and 35 minutes retention time
ms1 <- which(hd$msLevel == 1)
rtsel <- hd$retentionTime[ms1] / 60 > 30 &
  hd$retentionTime[ms1] / 60 < 35

## the heat map
M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd)
plot(M, aspect = 1, allTicks = FALSE)
plot3D(M)


i <- ms1[which(rtsel)][1]
j <- ms1[which(rtsel)][2]
M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
plot3D(M2)

plot(sp, reporters = iTRAQ4, full = TRUE)

#################
##  MS Spectra ##
#################
plot(sp, reporters = iTRAQ4, full = TRUE)
sel <- fData(itraqdata)$ProteinAccession == "BSA"
bsa <- itraqdata[sel]
bsa
as.character(fData(bsa)$ProteinAccession)
plot(bsa, reporters = iTRAQ4, full = FALSE) + theme_gray(8)

#####################
## MS Chromatogram ##
#####################





#########################
## Raw data processing ##
#########################
experiment <- removePeaks(itraqdata, t = 400, verbose = FALSE)
ionCount(itraqdata[["X55"]])

ionCount(experiment[["X55"]])

qnt <- quantify(experiment,
  method = "trap",
  reporters = iTRAQ4,
  strict = FALSE,
  verbose = FALSE)
qnt


