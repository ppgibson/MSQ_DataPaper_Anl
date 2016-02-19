## Format sample data prepared in <Calculate_AlternateSumPCBs.R> for
## Table S4, sample data:
##  1. add replicate letters
##  2. covert pct.det to percentages
##  3. round pcb sums to one dec place
##  4. reorder columns/remove sys_smp_code and stn_id
## 
## Input: <AllSampleData_CompareSums.csv>
##
## Output: <TabVals_SampleData.csv>


# Read in base data, prepared with <Calculate_AlternateSumPCBs.R>
  smp.data.raw <- read.csv(paste(DirOut, "AllSampleData_CompareSums.csv", sep=""))
  smp.data <- smp.data.raw  #retain smp.data.raw as an unmodified df
  
# Generate 'replicate' column (to replace [stn_id])
  # Count number of samples in each 'site mean'
    sitelist <- split(x=smp.data, 
                      f=list(smp.data$site_number, smp.data$category, smp.data$sample_year), #note order of factors is important to get sequence in correct order
                      drop=TRUE)
    n.persite <- lapply(X=sitelist, FUN=nrow)
    n.vect    <- data.frame(do.call(rbind, n.persite))  #vector (sort of) of number of elements in each site mean level
  # Assign A-D letters to define replicates
    rep.letters <- c("A", "B", "C", "D")
  # Add replicate letters as a new column
    smp.data$rep <- rep.letters[sequence(n.vect[,1])]
  # Clean up
    rm(sitelist); rm(n.persite); rm(n.vect); rm(rep.letters)

# Reorder columns and remove unnec data
  smp.form <- select(smp.data, sample_year, category, site_number, rep,   #basic metadata
                     smp.mass, pct.det, avg.mdl,                          #sample data
                     sum.zero, sum.halfdl, sum.km, sum.ref.cons)          #four different sum-pcb values
  
# Convert pct column to percentages
  smp.form <- mutate(smp.form, pct.det=round(100*pct.det, 1))

# Round off PCB sum values
  smp.form <- mutate(smp.form, 
                     sum.zero = round(sum.zero, 1),
                     sum.halfdl = round(sum.halfdl, 1), 
                     sum.km = round(sum.halfdl, 1), 
                     sum.ref.cons = round(sum.ref.cons, 1))
  #! Wait and round the mdl values in excel
  
# Change K-M sums to NA for non-detect samples
  smp.form$sum.km[smp.form$sum.zero==0] <- NA

# Rearrange so that sediment is first
  smp.form$category <- factor(smp.form$category, levels=c("Sediment", "Araneid", "Tetragnathid"))
  smp.form <- arrange(smp.form, sample_year, category, site_number)

# Print output
  write.csv(smp.form, paste(DirOut, "TabVals_SampleData.csv", sep=""), row.names=FALSE)

#### END SCRIPT ####