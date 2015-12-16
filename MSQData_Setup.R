##############################################################
# Setup for creating figures for Manistique 2015 report
##############################################################
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
  samples <- read.csv(paste(DirData, "AllSampleData.csv", sep=""), row.names=NULL)
  
# Calculate means by SITE (ie, create a new \means\ df)
  by.list <- list(sample_year = samples$sample_year,  #necessary for aggregate function
                  category = samples$category,
                  site_number = samples$site_number)
  
  means <- aggregate(x=samples[, 7:26],               #calculate mean for all data columns
                     by=by.list, 
                     FUN="mean", na.rm=TRUE)
  
  means <- means[with(means, order(sample_year, category, site_number)), ] #reorder
#   colnames(means)[11] <- "Total.PCB"                             
  means <- rename(means, Total.PCB = full)  #now, we are using [full] for the total.pcb values.
  row.names(means) <- NULL                  #get rid of row.names column
  rm(by.list)                               #unnecessary
  
# Create a wide-format data frame of PCB totals only
  means.w <- dcast(data = means,
                   formula = sample_year + site_number ~ category,
                   value.var = "Total.PCB")
  
# # Calculating st dev/st error of Tot.PCB if needed
#   totpcb.summ <- ddply(.data = samples, 
#                        .variables = c("sample_year", "category", "site_number"),
#                        .fun = summarise,
#                           mean = mean(full, na.rm=TRUE),
#                           sd = sd(full, na.rm=TRUE),
#                           n = sum(!is.na(full)),
#                           se = (sd(full, na.rm=TRUE)) / (sqrt(sum(!is.na(full)))) )
#   
#   means <- merge(means, totpcb.summ[, c(1:3, 7)], by=c("sample_year", "category", "site_number"), all=TRUE)
#   colnames(means)[24] <- "totpcb.se"                                       #better colname for new st error column
#   means <- means[with(means, order(sample_year, category, site_number)), ] #reorder
  
## Calculate relative homologs  
  
# Create a dataframe "rel.homologs" for new relative homolog values
  rel.homologs <- matrix(data=(rep(0, (nrow(means))*10)), ncol=10)  #create an empty data frame for writing new values into; there will be 10 values for every site (nrow(means)), and organize data into 10 columns
  rel.homologs <- as.data.frame(rel.homologs)
  colnames(rel.homologs) <- colnames(means)[14:23]   #[14:23] are the relevant, homolog columns of \means\
  
# calculate percent for each relevant column (i=column number)
  for (i in 14:23) {
    eachhom <- means[, i]/means$Total.PCB    
    rel.homologs[, (i-13)] <- eachhom 
  }
  
# turn proportions into percents for easier viewing
  rel.homologs <- rel.homologs*100
  
# Extract identifier columns from "means", and bind with the relative homolog values
  rel.homologs <- cbind((means[, 1:3]), rel.homologs)
  row.names(rel.homologs) <- NULL
  
# clean up unneeded extra variables
  rm(eachhom)
  rm(i)
