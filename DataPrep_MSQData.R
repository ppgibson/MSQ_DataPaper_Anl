## Purpose: generate cleaned/standardized data files for use in 
## analyses for the MSQ data paper
##
## Inputs: use data files exported directly from the database
##  (1) congener-level pcb data from MSQ, AllResults_Congeners_Man.csv
##  (2) sample-level lipid/TOC data, TotPCBLipTOC_bySample_Value.csv
##  (3) sample-level mass data for 2011 samples, MSQ_2011_SampleMasses.csv
## ...plus spreadsheet of the various lab-specific MDLs from Battelle:
##  (4) MDLs_BDOlab_AllYrs.csv
## 
## Outputs
##  (a) CleanedData_Congeners.csv: all concentrations in ng/g; 
##                                 all congeners have a sample-specific mdl value (manually calculated for BDO results), in ng/g;
##                                 all results associated with a sample mass (in g)
##  (b) CleanedData_Samples.csv: sum.pcb values are sums of all result values above MDL (ie, zero substituted for anything below MDL, including both NAs and reported values).
##                               each sample also has a value for smp.mass and pct.lip/toc and pcb conc normalized by lip/toc.


#### PREPARE CONGENER-LEVEL DATA ####
# Basic congener-level data
  # Read in AllResults_Congeners_Man query export file.  
    cons.raw <- read.csv("AllResults_Congeners_Man.csv")
  
  # Remove biphenyl values before doing anything else #?? Why are we removing biphenyl, again??
  # (There should be 40746 - 214 = 40532 rows remaining after)
    cons <- cons.raw[cons.raw$chemical_name!="Biphenyl", ]
  
  # Convert 2011 values (result, mdl, res_unit) from ng/kg to ng/g to match results from 2012-13
    cons$result_value[cons$result_unit=="ng/kg"] <- 0.001*(cons$result_value[cons$result_unit=="ng/kg"])
    cons$method_detection_limit[cons$result_unit=="ng/kg"] <- 0.001*(cons$method_detection_limit[cons$result_unit=="ng/kg"])
    cons$result_unit[cons$result_unit=="ng/kg"] <- "ng/g"

# Standardize MDLs: sample specific, manually calculated according to instructions 
# and lab-standard MDLs provided by Lisa Lefkovitz.  
  # Add a 'matrix' category for matching sed vs tissue MDLs
    cons <- mutate(cons, matrix="tissue")
    cons$matrix[cons$category=="Sediment"] <- "sediment"
  
  # Read in lab-MDL values
    mdls.lab <- read.csv("MDLs_BDOlab_AllYrs.csv")
    mdls.lab <- rename(mdls.lab, sediment=sed.MDL, tissue=tissue.MDL)
    mdls.lab.l <- melt(data=mdls.lab, id.vars=c("chemical_name", "mdl.year"), 
                       variable.name="matrix", value.name="mdl.lab.ul")  #mdl.lab.ul to indicate microliter units.
  # Relabel the values from the '2010' spreadsheet as '2012'; 
  # BDO data from 2012 used the '2010' lab standard MDL values. 
    mdls.lab.l$mdl.year[mdls.lab.l$mdl.year==2010] <- 2012
  # Add the standard mass values used for the lab-MDL tests
    mdls.lab.l$mass.mdl[mdls.lab.l$mdl.year==2012] <- 20.0  #both sed and tiss have std. mass of 20 g in the 2010 MDL list
    mdls.lab.l$mass.mdl[mdls.lab.l$mdl.year==2013 & mdls.lab.l$matrix=="sediment"] <- 16.64
    mdls.lab.l$mass.mdl[mdls.lab.l$mdl.year==2013 & mdls.lab.l$matrix=="tissue"]   <- 20.14
  # Convert the lab standard MDL values to units of ng/g (as explained by LL)
    mdls.lab.l <- mutate(mdls.lab.l, mdl.lab.ng = mdl.lab.ul * mass.mdl / 2)  #dilution factor for all is 2.0

  # Merge MDL data with cons data
    condata <- merge(cons, mdls.lab.l,
                  by.x=c("sample_year", "matrix", "chemical_name"),
                  by.y=c("mdl.year", "matrix", "chemical_name"),
                  all.x=TRUE, all.y=FALSE)
  # Remove no-longer-needed data
    rm(mdls.lab)
    rm(mdls.lab.l)

  # Calculate sample-specific MDLs for 2012/2013 (will be na for 2011)
    condata <- mutate(condata, mdl.ss = round(((mdl.lab.ng * dilution_factor) / subsample_amount), 2))
  # For 2011, use the reported MDLs (which we have confirmed are sample-specific, even if we don't know lab MDLs)
    condata$mdl.ss[condata$sample_year==2011] <- condata$method_detection_limit[condata$sample_year==2011]

# Sample mass
  # Need to read in 2011 mass data separately
    mass11 <- read.csv("MSQ_2011_SampleMasses.csv")  #Normal field samples only, no field splits or QC smp.
    colnames(mass11) <- c("sys_sample_code", "smp.mass")
  # Merge with condata - should be NAs for 2012-13 data
    condata <- merge(condata, mass11, by="sys_sample_code", all=TRUE)
  # Add 2012-13 information to [smp.mass] column so that all information is in one field
    condata$smp.mass[condata$sample_year!=2011] <- condata$subsample_amount[condata$sample_year!=2011]

# Indicate which results are censored
  condata <- mutate(condata, censor.ind=as.numeric(result_value<mdl.ss))  #any time reported result value is less than calculated mdl.ss
  condata$censor.ind[is.na(condata$result_value)] <- 1  #plus, any non-detected congener run counts as censored

# Remove unneeded columns
  condata <- condata[, c("sample_year", "matrix", "category", "site_number", "stn_id",
                         "sys_sample_code", "chemical_name", "result_value", 
                         "mdl.ss", "censor.ind", "smp.mass", "dilution_factor", 
                         "lab_testing", "CL_LEVEL", "DISPLAY_ORDER")]
  condata <- rename(condata, homolog=CL_LEVEL, con.order=DISPLAY_ORDER)
  condata <- arrange(condata, sample_year, category, site_number, stn_id, sys_sample_code, con.order)

# If desired, write an output file
  write.csv(condata, "CleanedData_Congeners.csv", row.names=FALSE)


#### PREPARE SAMPLE-LEVEL DATA ####
# Basic congener concentrations
  # Read in the cleaned congener data
  # (or use \condata\ from the previous section)
  # all units are ng/g (concentrations), g (sample masses)
    condata <- read.csv("CleanedData_Congeners.csv")
  
  # Add an extra [temporary] to enable easier calculation of sum-PCB
    condata <- mutate(condata, resval.cen=result_value) #For all detected obs, use the reported value.
      condata$resval.cen[condata$censor.ind==1] <- 0    #For all censored obs (nondetect or below ssMDL), substitute zero.
  
  # Calculate sum-PCB for each sample
    bysmp <- group_by(condata, sample_year, category, site_number, stn_id, sys_sample_code)
    smp.sums <- summarize(bysmp, 
                          sum.pcb  = sum(resval.cen),  #sum of all observations above MDL
                          smp.mass = min(smp.mass))    #all smp.mass values for a given sample should be same

# Lipid/TOC data
  # Read in data from database query
    liptoc <- read.csv("TotPCBLipTOC_bySample_Value.csv")
  # Manually correct some problems
    # Remove lip data for samples that have been disqualified due to QA problems (Man2012, spiders, 223 and 251)-
    # the lipid data itself should be fine, but we don't need it since we don't have reliable PCB data for these smp.
      liptoc <- liptoc[!(liptoc$sys_sample_code %in% c("223", "251")), ] 
    # Replace zeros for TOC
      # For four of the zero-samples, seems likely that the value was somewhere below MDL (0.02);
      # substitute half-MDL (i.e., 0.01) for these samples.
      liptoc$TOC[liptoc$TOC==0] <- 0.01
      # But, for the sample from site 5, it seems highly improbable that the true TOC was < 0.02.
      # Set this one to NA instead.
      liptoc$TOC[liptoc$sys_sample_code=="L1064557-18"] <- NA

  # Merge with sum-PCB totals 
    samples <- merge(smp.sums, liptoc[, c("sys_sample_code", "PCT_LIPID", "TOC")],
                     by="sys_sample_code", all=TRUE)
    samples <- rename(samples, pct.lip=PCT_LIPID, toc=TOC)

# Calculate lipid- and toc-normalized PCB concentrations
  samples <- mutate(samples, 
                    pcb.lipidnorm = sum.pcb/(pct.lip/100),
                    pcb.tocnorm   = sum.pcb/(toc/100))
#   # BUT! If we think that nondetects/zeros should not be used in normalized concentrations:
#   samples$pcb.lipidnorm[samples$sum.pcb==0] <- NA
#   samples$pcb.tocnorm[samples$sum.pcb==0] <- NA

# Write .csv, if desired
  # (Re-order first)
  samples <- arrange(samples, sample_year, category, site_number, stn_id, sys_sample_code)
  write.csv(samples, "CleanedData_Samples.csv", row.names=FALSE)
