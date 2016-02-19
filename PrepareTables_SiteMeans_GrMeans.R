## Using standard sample-level data, calculate site means 
## and grand means for sum-PCB/lipid and TOC/normalized PCB
## (including SD and n), and format data for summary tables.  
## 
## Input: \samples\ df from setup
##
## Output: <TabVals_SiteMeans.csv>
##         <TabVals_GrMeans.csv>


#### SITE MEANS: tables S5-S6 ####
# Calculate site level mean/sd/n for Total PCB, % TOC, and % Lipid.
  bysiteyr <- group_by(samples, sample_year, category, site_number)
  site.totals  <- summarize(bysiteyr,  
                            mass.mean = mean(smp.mass),
                            pcb.mean = mean(sum.pcb, na.rm=TRUE),  #use pcb.composite - adjusted 2011 values - for PCB totals.
                            pcb.sd = sd(sum.pcb, na.rm=TRUE),
                            pcb.n = sum(!is.na(sum.pcb)),
                            toc.mean = mean(toc, na.rm=TRUE),
                            toc.sd = sd(toc, na.rm=TRUE),
                            toc.n = sum(!is.na(toc)),
                            lip.mean = mean(pct.lip, na.rm=TRUE),
                            lip.sd = sd(pct.lip, na.rm=TRUE),
                            lip.n = sum(!is.na(pct.lip)),
                            pcb.tocnorm.mean = mean(pcb.tocnorm, na.rm=TRUE),
                            pcb.tocnorm.sd = sd(pcb.tocnorm, na.rm=TRUE),
                            pcb.tocnorm.n = sum(!is.na(pcb.tocnorm)),
                            pcb.lipnorm.mean = mean(pcb.lipidnorm, na.rm=TRUE), 
                            pcb.lipnorm.sd = sd(pcb.lipidnorm, na.rm=TRUE), 
                            pcb.lipnorm.n = sum(!is.na(pcb.lipidnorm)) )
                           


# Format summary values 
  # Total.PCB
    pcb.mean <- format(round(site.totals$pcb.mean, 1), nsmall=1, trim=TRUE, big.mark=",")
    pcb.sd   <- paste(" (± ", 
                     format(round(site.totals$pcb.sd, 1), nsmall=1, trim=TRUE, big.mark=","),
                     ", ",
                     site.totals$pcb.n,
                     ")", sep="")

  # TOC
    toc.mean <- format(round(site.totals$toc.mean, 1), nsmall=1, trim=TRUE, big.mark=",")
    toc.sd   <- paste(" (± ",
                     format(round(site.totals$toc.sd, 1), nsmall=1, trim=TRUE, big.mark=","),
                     ", ",
                     site.totals$toc.n,
                     ")", sep="")

  # Total.PCB normalized to TOC (pcb.tocnorm)
    pcb.tocnorm.mean <- format(round(site.totals$pcb.tocnorm.mean, 0), nsmall=0, trim=TRUE, big.mark=",")                            
    pcb.tocnorm.sd   <- paste(" (± ", 
                          format(round(site.totals$pcb.tocnorm.sd, 0), nsmall=0, trim=TRUE, big.mark=","), 
                          ", ",
                          site.totals$pcb.tocnorm.n,
                          ")", sep="")

  # Lipids
    lip.mean <- format(round(site.totals$lip.mean, 1), nsmall=1, trim=TRUE, big.mark=",")  
    lip.sd   <- paste(" (± ",
                     format(round(site.totals$lip.sd, 1), nsmall=1, trim=TRUE, big.mark=","),
                     ", ",
                     site.totals$lip.n,
                     ")", sep="")

  # Total.PCB normalized to lipid content (pcb.lipidnorm)
    pcb.lipidnorm.mean <- format(round(site.totals$pcb.lipnorm.mean, 0), nsmall=0, trim=TRUE, big.mark=",")                            
    pcb.lipidnorm.sd   <- paste(" (± ", 
                              format(round(site.totals$pcb.lipnorm.sd, 0), nsmall=0, trim=TRUE, big.mark=","), 
                              ", ",
                              site.totals$pcb.lipnorm.n,
                              ")", sep="")

## Write output table
  # Combine all summary columns plus identifier columns in a single dataframe 
    table.sitemeans <- data.frame(cbind(site.totals[, 1:4], pcb.mean, pcb.sd, toc.mean, toc.sd, 
                                        pcb.tocnorm.mean, pcb.tocnorm.sd, lip.mean, lip.sd,
                                        pcb.lipidnorm.mean, pcb.lipidnorm.sd))

  # Export the table
    write.csv(x=table.sitemeans, file=paste(DirOut, "TabVals_SiteMeans.csv", sep=""), row.names=FALSE)

# Clean up
  rm(bysiteyr)
  rm(pcb.mean)
  rm(pcb.sd)
  rm(toc.mean)
  rm(toc.sd)
  rm(pcb.tocnorm.mean)
  rm(pcb.tocnorm.sd)
  rm(lip.mean)
  rm(lip.sd)
  rm(pcb.lipidnorm.mean)
  rm(pcb.lipidnorm.sd)


#### GRAND MEANS: Tables 1-2 ####

# Calculate the mean of means at each site, using values from site.totals
# n should =3 for most variables
  bysite <- group_by(site.totals, category, site_number)
  overall.means <- summarize(bysite, 
                             pcb.avg = mean(pcb.mean, na.rm=TRUE),
                             pcb.sd = sd(pcb.mean, na.rm=TRUE),
                             pcb.n = sum(!is.na(pcb.mean)),
                             toc.avg = mean(toc.mean, na.rm=TRUE),
                             toc.sd = sd(toc.mean, na.rm=TRUE),
                             toc.n = sum(!is.na(toc.mean)),
                             lip.avg = mean(lip.mean, na.rm=TRUE),
                             lip.sd = sd(lip.mean, na.rm=TRUE),
                             lip.n = sum(!is.na(lip.mean)),
                             pcb.toc.avg = mean(pcb.tocnorm.mean, na.rm=TRUE),
                             pcb.toc.sd = sd(pcb.tocnorm.mean, na.rm=TRUE),
                             pcb.toc.n = sum(!is.na(pcb.tocnorm.mean)),
                             pcb.lip.avg = mean(pcb.lipnorm.mean, na.rm=TRUE),
                             pcb.lip.sd = sd(pcb.lipnorm.mean, na.rm=TRUE),
                             pcb.lip.n = sum(!is.na(pcb.lipnorm.mean)) )


# Format summary values 
  # Total.PCB      #Note that this  calculation is based originally on [pcb.composite] - so 2011 values were BDO-adjusted before averageing.
    pcb.gm.mean <- format(round(overall.means$pcb.avg, 1), nsmall=1, trim=TRUE, big.mark=",")             
    pcb.gm.sd   <- paste(" (± ", 
                      format(round(overall.means$pcb.sd, 1), nsmall=1, trim=TRUE, big.mark=","),
                      ", ",
                      overall.means$pcb.n,
                      ")", sep="")
  
  # TOC
    toc.gm.mean <- format(round(overall.means$toc.avg, 1), nsmall=1, trim=TRUE, big.mark=",")
    toc.gm.sd   <- paste(" (± ",
                      format(round(overall.means$toc.sd, 1), nsmall=1, trim=TRUE, big.mark=","),
                      ", ",
                      overall.means$toc.n,
                      ")", sep="")
  
  # Lipids
    lip.gm.mean <- format(round(overall.means$lip.avg, 1), nsmall=1, trim=TRUE, big.mark=",")
    lip.gm.sd   <- paste(" (± ",
                      format(round(overall.means$lip.sd, 1), nsmall=1, trim=TRUE, big.mark=","),
                      ", ",
                      overall.means$lip.n,
                      ")", sep="")
  
  # Total.PCB normalized to TOC (pcb.tocnorm)
    pcb.toc.gm.mean <- format(round(overall.means$pcb.toc.avg, 0), nsmall=0, trim=TRUE, big.mark=",")
    pcb.toc.gm.sd   <- paste(" (± ", 
                             format(round(overall.means$pcb.toc.sd, 0), nsmall=0, trim=TRUE, big.mark=","), 
                             ", ",
                             overall.means$pcb.toc.n,
                             ")", sep="")

  # Total.PCB normalized to lipid content (pcb.lipidnorm)
    pcb.lip.gm.mean <- format(round(overall.means$pcb.lip.avg, 0), nsmall=0, trim=TRUE, big.mark=",")
    pcb.lip.gm.sd   <- paste(" (± ", 
                             format(round(overall.means$pcb.lip.sd, 0), nsmall=0, trim=TRUE, big.mark=","), 
                             ", ",
                             overall.means$pcb.lip.n,
                             ")", sep="")


# Combine all summary columns plus identifier columns in a single dataframe 
  table.grmeans <- data.frame(cbind(overall.means[, 1:2], pcb.gm.mean, pcb.gm.sd, toc.gm.mean, toc.gm.sd, 
                                    pcb.toc.gm.mean, pcb.toc.gm.sd, lip.gm.mean, lip.gm.sd,
                                    pcb.lip.gm.mean, pcb.lip.gm.sd))

# Export the table
  write.csv(x=table.grmeans, file=paste(DirOut, "TabVals_GrMeans.csv", sep=""), row.names=FALSE)


# Clean up
  rm(pcb.gm.mean)
  rm(pcb.gm.sd)
  rm(toc.gm.mean)
  rm(toc.gm.sd)
  rm(pcb.toc.gm.mean)
  rm(pcb.toc.gm.sd)
  rm(lip.gm.mean)
  rm(lip.gm.sd)
  rm(pcb.lip.gm.mean)
  rm(pcb.lip.gm.sd)
  

##### END SCRIPT ####