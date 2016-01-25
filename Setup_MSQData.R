########################################################
# Setup for creating figures for Manistique data paper #
########################################################
# Data output   
#   1. \samples\    Single data frame with all sample-level data: 
#                     - % lipid
#                     - % TOC
#                     - Total.PCB (sum of congeners)
#                     - Subset and normalized total.PCB values
#                     - Total PCB, lipid and TOC-normalized
#                     - Homolog totals
#   2. \means\      All data from \samples\, but averaged for each site (ie, each data point in \means\ is an average of ~4 data points in \samples\).
#   3. \means.w\    Total.PCB value for each sample in wide format (one column each for Ara, Tet, and Sed)
#   4. \rel.homologs\  Homolog proportions, calculated as a proportion of the mean Total.PCB for each site. [Note that the homolog totals used in this calculation are (like Total.PCB) site means, not measured values.]
#
# If standard deviation/standard error of Total.PCB [or other] values 
# are needed, run the commented-out code, adapt as necessary.

# Load packages
  library(plyr)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
#   library(gridExtra)
#   library(RColorBrewer)
#   library(scales)
#   library(vegan)
#   library(devEMF)   # for printing .emf format plots

# Define directory paths
  # R code
    DirCode <- "H:/DATA/Great Lakes Files/R Projects/MSQ_DataPaper_Anl/"
  # Data (input files and also write.table)
    DirData <- "H:/DATA/Great Lakes Files/R Projects/MSQ_DataPaper_Anl/Data/"
  # Output, temporary only, files may be overwitten on re-execution
    DirOut <- "H:/DATA/Great Lakes Files/R Projects/MSQ_DataPaper_Anl/Output/"

# Set WD
    setwd(DirData)


# Read in data
  congeners <- read.csv(paste(DirData, "CleanedData_Congeners.csv", sep=""))
  samples   <- read.csv(paste(DirData, "CleanedData_Samples.csv", sep=""), row.names=NULL)
  
# Calculate means by SITE 
  by.site <- group_by(samples, sample_year, category, site_number)
  sitemeans <- summarize(by.site,
                         smp.mass= mean(smp.mass), 
                         sum.pcb = mean(sum.pcb),
                         n.smp = n(), 
                         pct.toc = mean(toc, na.rm=TRUE),
                         pcb.tocnorm   = mean(pcb.tocnorm, na.rm=TRUE),
                         pct.lip = mean(pct.lip, na.rm=TRUE),
                         pcb.lipidnorm = mean(pcb.lipidnorm, na.rm=TRUE) )
  rm(by.site)
                         
# Create a wide-format data frame of PCB totals only
  means.w <- dcast(data = sitemeans,
                   formula = sample_year + site_number ~ category,
                   value.var = "sum.pcb")
  
