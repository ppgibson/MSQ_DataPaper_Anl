######################################
# Minor analyses for results section #
######################################

library(gridExtra)

#### 1. Site CVs ####
# Calculate cv among individual point samples/core samples within a site (ie, n<=4)
  byyrsite <- group_by(samples, sample_year, category, site_number)
  siteyrstats <- summarize(byyrsite, 
                         n.smp = n(), 
                         mean.conc = mean(sum.pcb, na.rm=TRUE), #yes we already have this data, but this is just a check.
                         sd.conc = sd(sum.pcb, na.rm=TRUE),
                         cv.conc = sd.conc/mean.conc,
                         range.conc = max(sum.pcb, na.rm=TRUE) - min(sum.pcb, na.rm=TRUE))

# Now calculate the mean cv across three sample years for each site-category combo
  bysite <- group_by(siteyrstats, category, site_number)
  sitestats <- summarize(bysite,
                         n.yrs = n(), 
                         grmean.conc = mean(mean.conc, na.rm=TRUE),
                         sd.sitemeans = sd(mean.conc),
                         mean.cv.conc = mean(cv.conc))

# What are the cv values in each year for sediment sites 3 and 5?
  arrange(filter(ungroup(siteyrstats), category=="Sediment" & site_number %in% c(6,7)), site_number, sample_year)

  #...and what are the mean CVs across all three years?
  filter(sitestats, category=="Sediment" & site_number %in% c(3,5))

  rm(bysite)
  rm(byyrsite)

# What is the CV among site means (all years combined - ?) for spiders and for sediment?
  # spiders
  sitemeans.spid <- filter(sitemeans, category!="Sediment")
  sd(sitemeans.spid$sum.pcb) / mean(sitemeans.spid$sum.pcb)  #CV = 0.883
  # sediment
  sitemeans.sed <- filter(sitemeans, category=="Sediment")
  sd(sitemeans.sed$sum.pcb) / mean(sitemeans.sed$sum.pcb)  #CV = 1.69

# Or, what if we calculate the CV among site means within a given year-category
  byyrcat <- group_by(sitemeans, sample_year, category)
  annual <- summarize(byyrcat,
                      yr.mean = mean(sum.pcb),
                      yr.sd = sd(sum.pcb),
                      yr.cv = yr.sd/yr.mean)

  mean(annual$yr.cv[annual$category=="Sediment"])  #mean CV across sites within a year is 1.76 for sed
  mean(annual$yr.cv[annual$category!="Sediment"])  #0.82 for spid (n=6, 3 tets and 3 ara)


#### 2. Normalization and site rankings ####
## How does toc-/lipid-normalization change the relative magnitude of 
## PCB concentration across sites?

# Grand means
  bysite  <- group_by(sitemeans, category, site_number)
  grmeans <- summarize(bysite, 
                       pcb.avg = mean(sum.pcb, na.rm=TRUE),
                       pcb.toc.avg = mean(pcb.tocnorm, na.rm=TRUE),
                       pcb.lip.avg = mean(pcb.lipidnorm, na.rm=TRUE) )

  # Make a single column for normalized pcb concentrations
  grmeans <- mutate(grmeans, pcb.norm=pcb.lip.avg)
  grmeans$pcb.norm[grmeans$category=="Sediment"] <- grmeans$pcb.toc.avg[grmeans$category=="Sediment"]

  pcb.bywt <- ggplot(data=grmeans) + 
    geom_bar(aes(x=factor(site_number), y=pcb.avg, fill=factor(site_number)), stat="identity") + 
    facet_wrap(~ category, scales="free_y") +
    ggtitle("Dry/Wet weight concentrations") +
    theme(legend.position="none")

  pcb.norm <- ggplot(data=grmeans) + 
    geom_bar(aes(x=factor(site_number), y=pcb.norm, fill=factor(site_number)), stat="identity") + 
    facet_wrap(~ category, scales="free_y") +
    ggtitle("Normalized concentrations") +
    theme(legend.position="none")

  grid.arrange(pcb.bywt, pcb.norm, ncol=1)
