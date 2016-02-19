################################################
## MSQ Data paper, data analysis:             ##
## Basic QAQC statistics reported in the text ##
################################################
# Inputs:
#  - \congeners\ df from setup
#  - <QAQC_USGS_Congener-Level_Results.csv>,
#     data file including QAQC smp and surrogate spike runs
#
# Output: basic statistics (in the console)

# Run setup first
# Prevent printing scientific notation
  options(scipen=8)

#### A. METHOD DETECTION LIMITS ####
# Use sample-adjusted method det. limits, already calculated for BDO during DataPrep.
# Use the reported MDLs for 2011 - which we have confirmed are, at least, sample specific.
  consbyyr <- group_by(congeners, sample_year)
  mdl.stats <- summarize(consbyyr, 
                        n.conruns = n(),
                        mean.mdl = mean(mdl.ss),
                        sd.mdl = sd(mdl.ss),
                        min.mdl = min(mdl.ss), 
                        max.mdl = max(mdl.ss))
  # Clean up
  rm(consbyyr)
  
# **Values to report:
  # What is the range in ssMDL values for 2011? (ALS, high-res anl)
  range(congeners$mdl.ss[congeners$sample_year==2011])
  
  # What is the range in ssMDL values for 2012-13?  (BDO, low-res anl)
  range(congeners$mdl.ss[congeners$sample_year!=2011])

  
#### B. SURROGATE SPIKES ####
# Read in large data file that includes (a) QAQC samples (the ones that apply
# to PCB data); and (b) surrogate congener data.  
  cons.qaqc <- read.csv("QAQC_USGS_Congener-Level_Results.csv")

  # Restrict to only the data we need for this analysis: 
  #  MSQ only (no ASH) and con_surr runs only.  
  cons.qaqc <- filter(cons.qaqc, river=="Manistique")
  sis <- filter(cons.qaqc, result_type_code=="SUR")
  
  # Are we working with the correct number of samples?
    # Numbers of samples with sis data
    sis.smp <- unique(select(sis, sample_year, category, sys_sample_code))
    table(sis.smp$sample_year, sis.smp$category)
    # Numbers of samples from 2012-13 in our standard data:
    table(samples$sample_year, samples$category)
      # perfect match! (except of course that sis includes ~20 QA smp from each yr)
    rm(sis.smp)
    
  # QA/QC check of the data: do the [qc_spike_recovery] values calculated by BDO
  # match the numbers reported in [qc_spike_added] and [qc_spike_recovered]?
    sis <- mutate(sis, spike.recov.man = 100*(round((qc_spike_measured/qc_spike_added), 2)))
    sis <- mutate(sis, recov.diff = qc_spike_recovery - spike.recov.man)
    range(sis$recov.diff)  # No differences, good.

# Calculate summary statistics by year/category groups   
  sis.bycatyr <- group_by(sis, sample_year, category)
  sis.stats <- summarize(sis.bycatyr,
                         n.samp = length(unique(sys_sample_code)),
                         n = n(), 
                         n.values = sum(!is.na(qc_spike_recovery)),  
                         min = min(qc_spike_recovery, na.rm=TRUE),
                         max = max(qc_spike_recovery, na.rm=TRUE),
                         mean = mean(qc_spike_recovery, na.rm=TRUE),
                         sd = sd(qc_spike_recovery, na.rm=TRUE))
  rm(sis.bycatyr)

# **Values to report:
  # What is the mean of sis recovery values acros ALL SAMPLES?
  mean(sis$qc_spike_recovery)
  # ...the standard deviation?
  sd(sis$qc_spike_recovery)
  # And how many values/how many samples are included in the mean? 
  length(unique(sis$sys_sample_code))
  length(!is.na(sis$qc_spike_recovery))
  
###################################################################