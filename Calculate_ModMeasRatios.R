## Recalculate ratios between modeled sum-119 and measured sum-209
## PCB concentrations
## This version takes into account data censoring - removes these concentrations
## from calculation of the sums (*both full and subset)

## Calculate various sum-pcb totals 
  # Dataframe for manipulation
  cons <- congeners 
  
  # Extra column for excluding censored results
  cons$res.cen <- NA
    cons$res.cen[cons$censor.ind==0] <- cons$result_value[cons$censor.ind==0]
  
  # Extra column for conc of subset congeners only
  cons$res.sub <- NA
    cons$res.sub[cons$lab_testing=="both"] <- cons$result_value[cons$lab_testing=="both"]
    cons$res.sub[cons$censor.ind==1] <- NA  #!censored data shouldn't count toward subset totals, either. 
  
  # Now calculate sums
  by.smp <- group_by(cons, sample_year, matrix, sys_sample_code)
  smp.sums <- summarize(by.smp, 
                        sum.full = sum(res.cen, na.rm=TRUE),  #excluding results below det.lin
                        sum.sub  = sum(res.sub, na.rm=TRUE))
  
  # Add some extra columns for un-ambiguous ratio calculation
    # 2011
    smp.sums$sum.209 <- NA
      smp.sums$sum.209[smp.sums$sample_year==2011] <- smp.sums$sum.full[smp.sums$sample_year==2011]  #2011 sums should go in the sum.209 column
  
    # 2012-13
    smp.sums$sum.119 <- NA
      smp.sums$sum.119[smp.sums$sample_year>2011] <- smp.sums$sum.full[smp.sums$sample_year>2011]  #2012-13 sums should go in the sum.119 column
  
  
## Generate modeled values
# Sediment
  # Fit the lm sum.119 ~ sum.subset using 2012-13 data
  sed.bdo.rows <- which(smp.sums$sample_year>2011 & smp.sums$matrix=="sediment")
  sed.119.mod  <- lm(sum.full ~ 0 + sum.sub, data=smp.sums[sed.bdo.rows, ], na.action=na.omit)
  
  # Generate a predicted value for each sample - both 2011 and 2012-13
  smp.sums$mod.119[smp.sums$matrix=="sediment"] <- predict(sed.119.mod, newdata=smp.sums[smp.sums$matrix=="sediment", ])
  
# Spider
  # Fit the lm sum.119 ~ sum.subset using 2012-13 data
  spid.bdo.rows <- which(smp.sums$sample_year>2011 & smp.sums$matrix=="tissue")
  spid.119.mod  <- lm(sum.full ~ 0 + sum.sub, data=smp.sums[spid.bdo.rows, ], na.action=na.omit)
  
  # Generate a predicted value for each sample - both 2011 and 2012-13
  smp.sums$mod.119[smp.sums$matrix=="tissue"] <- predict(spid.119.mod, newdata=smp.sums[smp.sums$matrix=="tissue", ])


## Ratios
  # Mod.119 : Meas.209 [Full] - for all 2011 samples
  smp.sums <- mutate(smp.sums, ra.mo119.me209 = mod.119 / sum.209)
  
  # Summary stats: all 2011 smp
  mean(smp.sums$ra.mo119.me209, na.rm=TRUE)
  sd(smp.sums$ra.mo119.me209, na.rm=TRUE)
  
  # Summary stats: sed vs spider
  tapply(X=smp.sums$ra.mo119.me209, INDEX=smp.sums$matrix, 
         FUN=mean, na.rm=TRUE)
  
  tapply(X=smp.sums$ra.mo119.me209, INDEX=smp.sums$matrix, 
         FUN=sd, na.rm=TRUE)
  
#### END SCRIPT #### 