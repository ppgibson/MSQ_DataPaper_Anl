##################################################################
## MSQ Data paper, data analysis:                               ##
## Generate sum-PCB values for each sample according to several ##
## different approaches to handling censored data               ##
##################################################################
#
# Four different approaches to handling censored data:
# 1. substitute zero (standard approach)
# 2. substitute half the detection limit
# 3. estimate a sample sum-PCB value using the Kaplan-Meier approach,
#    as described by Helsel (2009)
# 4. estimate concentrations for censored data based on the concentration
#    of a set of common congeners 
#    (*for now, this method is applied to 2012-13 spider samples only)
#
# Input: \congeners\ df from setup  
# 
# Output- in the Output folder, a .csv table with a sum-PCB conc for
#  each sample according to the four methodsin the output folder:
#   - AllSampleData_CompareSums.csv


# Packages
  library(NADA)  #Helsel's package for estimating summary statistics including censored data

#### A. CALCULATE KAPLAN-MEIER SUMS ####
## As described by Helsel (2009), use the non-parametric Kaplan Meier approach
## to estimate a sum for each sample containing nondetects - 
## despite the issues, in our data, that
## (a) Many samples have < 0.5 rate of detection across all congeners.
## (b) Concentrations of congeners are not independent from one another- 
##     they occur according, loosely, to a characteristic congener profile.

# To calculate KM, need ss.MDL value for each nondetect result
  congeners <- mutate(congeners, resval.km=result_value) #For all detected obs, use the reported value.
    congeners$resval.km[congeners$censor.ind==1] <- congeners$mdl.ss[congeners$censor.ind==1]
  # Convert censor.ind to a TRUE/FALSE for use in the NADA function
  congeners$censor.ind <- as.logical(congeners$censor.ind)  

# Split the cong-level data into separate dfs, one for each sample
  smp.list <- split(x=congeners, f=congeners$sys_sample_code)  #n=307 congeners; shd be 119 or 162 rows (congeners) per df.

# Function to calculate the KM sum for a given data set 
  kmsum <- function(df){
    kmfit.cur <- cenfit(obs=df$resval.km, censored=df$censor.ind)
    mean.cur <- mean(kmfit.cur)[1]
    n.cur <- kmfit.cur@survfit$n
    est.sum <- mean.cur*n.cur
    return(est.sum)
  }

# Calculate KM-sum for each sample
  kmsums.l <- lapply(X=smp.list, FUN=kmsum)
  kmsums <- melt(unlist(kmsums.l))
  kmsums$sys_sample_code <- names(kmsums.l)
  kmsums <- rename(kmsums, sum.km=value)
  rownames(kmsums) <- NULL

# Clean up 
  rm(kmsums.l); rm(smp.list)

#### B. ESTIMATE VALUES FOR CENSORED DATA USING CONGENER PROFILES ####
# Calculate rate-of-detection by congener
# (i.e., which are the most frequently detected congeners, across 2012-13 spid smp?)
  bycon <- group_by(filter(congeners, sample_year!=2011 & category!="Sediment"), chemical_name)
  consum <- summarize(bycon,
                      n.det=sum(!is.na(result_value)),
                      n.smp=n())
  consum <- arrange(consum, desc(n.det))
  rm(bycon)

  # 6 congeners appear in all (n=119) or all but one (n=118) of the eligible samples 
  # (134 total samples minus 15 complete nondetects = 119 eligible samples). 
  # These 6 congners will be the reference congeners.
    refcons <- consum$chemical_name[consum$n.det>=118] 
  
# Calculate the sum total concentration of the refcons in each 2012-13 sample
  # A new column with result values only for the ref cons
  spid.cons <- filter(congeners, category!="Sediment" & sample_year!=2011)  #Separate data frame for 2012-13 spider data only.
  spid.cons$resval.ref <- NA  #most congeners should have NA in this column...
    spid.cons$resval.ref[spid.cons$chemical_name %in% refcons] <- spid.cons$result_value[spid.cons$chemical_name %in% refcons] #...but the ref cons should have a value (if anything was detected). 

  # Now sum the refcon result values for each sample
  bysmp <- group_by(spid.cons, sys_sample_code)
  refsums <- summarize(bysmp, 
                       sum.ref.cons = sum(resval.ref, na.rm=TRUE))
  
  # Clean up
  rm(bysmp); rm(consum); rm(spid.cons); rm(refcons)

#### C. CALCULATE PCB-SUMS ####
# New columns containing the desired type of result value for each sum-type
  # (Remove the current 'resval.km' (replaces res_val with ss.mdl for censored),
  #  as this is not for summing and it is just confusing)
    congeners <- select(congeners, -resval.km)
  # Substitute zero
    congeners <- mutate(congeners, resval.zero=result_value) #For all detected obs, use the reported value.
      congeners$resval.zero[congeners$censor.ind==1] <- 0     #For all censored obs (nondetect or below ssMDL), substitute zero.
  # Substitute half det lim
    congeners <- mutate(congeners, resval.halfdl=result_value) #For all detected obs, use the reported value.
      congeners$resval.halfdl[congeners$censor.ind==1] <- 0.5 * (congeners$mdl.ss[congeners$censor.ind==1])    #For all censored obs, substitute half the ssMDL.

# Calculate sums
  bysmp <- group_by(congeners, sample_year, category, site_number, stn_id, sys_sample_code)
  smp.sums <- summarize(bysmp, 
                        n.det = sum(!is.na(result_value)),
                        pct.det = (sum(!is.na(result_value)))/n(), 
                        avg.mdl = mean(mdl.ss),
                        smp.mass = min(smp.mass),       #all smp.mass values for a given sample should be same
                        sum.zero   = sum(resval.zero),   #sum of all observations above MDL
                        sum.halfdl = sum(resval.halfdl) ) #substitute half det lim for censored obs  #substitute zero for censored obs, but sum concentrations for good cons only.

# Add in KM-estimated sums and refcon sums
  smp.sums <- merge(smp.sums, kmsums,  by="sys_sample_code", all=TRUE)
  smp.sums <- merge(smp.sums, refsums, by="sys_sample_code", all=TRUE)  #only applies to 2012-13 spider smp, others will be NA.


#### D. WRITE STANDARD OUTPUT TABLE ####
  smp.sums <- arrange(ungroup(smp.sums), sample_year, category, site_number, stn_id)

  write.csv(smp.sums, paste(DirOut, "AllSampleData_CompareSums.csv", sep=""), row.names=FALSE)

##### END SCRIPT ####