#######################################
#### SIA DATA: EXPLORATORY FIGURES ####
#######################################

#### Read in and prepare data ####

# Read in DB export
  sia <- read.csv("Manistique_IsotopeData.csv")
  sia <- arrange(sia, sample_year, category, new_site_id)

# Condense life_hist into a single column
  sia$life_history[sia$life_history=="I"] <- sia$indeterminates[sia$life_history=="I"] 
  colnames(sia)
  sia <- sia[, -10] #get rid of [analysis] and [indeterminates]

# Create separate categories for aquatic and terr inverts
  sia$category <- as.character(sia$category)
  sia$category[sia$category=="Invertebrate" & sia$life_history=="T"] <- "Terrestrial"
  sia$category[sia$category=="Invertebrate" & sia$life_history=="A"] <- "Aquatic"
  sia$category[sia$category=="Invertebrate"] <- "Aquatic"  #deal with the 6 remaining aq insect comp records
#   sia$category[is.na(sia$category)]
#   sia$category <- factor(sia$category, 
#                          levels=c("Araneid", "Tetragnathid", "Seston", "Aquatic", "Terrestrial"))


#### Summary statistics ####

# Sample sizes
  by.site <- group_by(sia, sample_year, category, site_number)
  by.site <- filter(by.site, chemical_name=="d13C")
  smpsize <- summarize(by.site,
                     n.smp=sum(!is.na(result_value)) )
  # Convert to wide format
  samples <- dcast(data=smpsize, sample_year + site_number ~ category, value.var="n.smp")
  samples[is.na(samples)] <- 0
  # Export
  write.csv(samples, "Isotope_SampleSizes.csv", row.names=FALSE)

  # Totals by year
  by.yr <- group_by(smpsize, sample_year, category)
  totals <- summarize(by.yr, n.smp=sum(n.smp, na.rm=TRUE))

  rm(by.yr)
  rm(by.site)

# What is the range in values for the different isotopes?
  by.iso <- group_by(sia, category, chemical_name)
  ranges <- summarize(by.iso,
                     n.smp=sum(!is.na(result_value)),
                     min=min(result_value, na.rm=TRUE),
                     max=max(result_value, na.rm=TRUE),
                     med=median(result_value, na.rm=TRUE))


## Calculate summary stats for making figures

# By site and category, so n=4 or less per site (usually)
  by.iso.site <- group_by(sia, sample_year, category, chemical_name, site_number)
  sitemeans <- summarize(by.iso.site, 
                         n.smp=sum(!is.na(result_value)),
                         mean=mean(result_value, na.rm=TRUE),
                         sd=sd(result_value, na.rm=TRUE))
  rm(by.iso.site)
  
  # Add se
  sitemeans <- mutate(sitemeans, se=sd/sqrt(n.smp))

  # Convert to a mostly wide form for graphing
  sitemeans.l <- melt(sitemeans, measure.vars=c("mean", "se"))  #first have to fully melt
  sitemeans.w <- dcast(sitemeans.l, sample_year + category + site_number ~ chemical_name + variable,
                       value.var="value")


# Overall means by category, so n=40 or more for spiders or less per site (usually)
# !The other option would be to take the mean of site means...
  by.iso.cat <- group_by(sia, sample_year, category, chemical_name) #or, should this be mean of means?
  catmeans <- summarize(by.iso.cat, 
                        n.smp=sum(!is.na(result_value)),
                        mean=mean(result_value, na.rm=TRUE),
                        sd=sd(result_value, na.rm=TRUE))
  rm(by.iso.cat)

  # Add se
  catmeans <- mutate(catmeans, se=sd/sqrt(n.smp))

  # Convert to a mostly wide form for graphing
  catmeans.l <- melt(catmeans, measure.vars=c("mean", "se"))  #first have to fully melt
  catmeans.w <- dcast(catmeans.l, sample_year + category ~ chemical_name + variable,
                       value.var="value")


#### PLOTS ####

## Variability by site

# # Reorder factor levels, if desired
#   sia$category <- factor(sia$category, levels=c("Araneid", "Tetragnathid", "Invertebrate"))

# Make temp data frames
  allsmp <- sia
  spiders <- filter(sia, category=="Araneid" | category=="Tetragnathid")

  # Pick 1
  sia <- spiders
  sia <- allsmp

  ggplot(data=sia[sia$chemical_name=="d15N", ]) +
    geom_point(aes(x=as.factor(site_number), y=result_value, fill=category), 
               size=4, shape=21, color="black") +
    xlab("Site") + 
    ylab("delta 15N") + 
    theme(legend.position="none") + 
    ggtitle("Nitrogen") + 
    facet_grid(category ~ sample_year)  


## BIPLOTS 
# Plot

  # Summary - means by year
  ggplot(data=catmeans.w) +
    geom_point(aes(x=d13C_mean, y=d15N_mean, shape=category, color=factor(sample_year)), size=4) + 
    geom_errorbar (aes(x=d13C_mean, ymin=(d15N_mean - d15N_se), ymax=(d15N_mean + d15N_se)) ) + 
    geom_errorbarh(aes(x=d13C_mean, xmin=(d13C_mean - d13C_se), xmax=(d13C_mean + d13C_se), y=d15N_mean) ) +
    scale_shape_manual(name="Matrix", 
                         breaks=c("Araneid", "Tetragnathid", "Aquatic", "Terrestrial", "Basal"),
                         labels=c("Araneid", "Tetragnathid", "Aquatic", "Terrestrial", "Basal"),
                         values=c(1, 19, 11, 0, 15)) + 
    xlab("delta 13C") + 
    ylab("delta 15N") +
    ggtitle("Overall means by matrix")
  

  # Site-by-site
  ggplot(data=sitemeans.w) +
    geom_point(aes(x=d13C_mean, y=d15N_mean, shape=category, color=factor(sample_year)), size=4) + 
    geom_errorbar (aes(x=d13C_mean, ymin=(d15N_mean - d15N_se), ymax=(d15N_mean + d15N_se)) ) + 
    geom_errorbarh(aes(x=d13C_mean, xmin=(d13C_mean - d13C_se), xmax=(d13C_mean + d13C_se), y=d15N_mean) ) +
    scale_shape_manual(name="Matrix", 
                         breaks=c("Araneid", "Tetragnathid", "Aquatic", "Terrestrial", "Basal"),
                         labels=c("Araneid", "Tetragnathid", "Aquatic", "Terrestrial", "Basal"),
                         values=c(1, 19, 11, 0, 15)) + 
    xlab("delta 13C") + 
    ylab("delta 15N") +
    ggtitle("Site by site, all matrices") + 
    facet_wrap( ~ site_number, nrow=3)
 

  # Site-by-site again, this time for spiders only
  temp <- filter(sitemeans.w, (category=="Araneid" | category=="Tetragnathid"))
  ggplot(data=temp) +
    geom_point(aes(x=d13C_mean, y=d15N_mean, shape=category, color=factor(sample_year)), size=4) + 
    geom_errorbar (aes(x=d13C_mean, ymin=(d15N_mean - d15N_se), ymax=(d15N_mean + d15N_se)) ) + 
    geom_errorbarh(aes(x=d13C_mean, xmin=(d13C_mean - d13C_se), xmax=(d13C_mean + d13C_se), y=d15N_mean) ) +
    scale_shape_manual(name="Matrix", 
                         breaks=c("Araneid", "Tetragnathid", "Aquatic", "Terrestrial", "Basal"),
                         labels=c("Araneid", "Tetragnathid", "Aquatic", "Terrestrial", "Basal"),
                         values=c(19, 15)) + 
    xlab("delta 13C") + 
    ylab("delta 15N") +
    ggtitle("Site by site, spiders only") + 
    facet_wrap( ~ site_number, nrow=3)
 


## Correlation between spiders and other ecosystem compartments across sites
# Spiders vs seston
  noinvert <- filter(sitemeans.w, sample_year==2011 & 
                       (category=="Araneid" | category=="Tetragnathid" | category=="Basal") )
  noinvert$category <- factor(noinvert$category, levels=c("Araneid", "Tetragnathid", "Basal"))
  
  ggplot(data=noinvert[noinvert$category!="Araneid", ]) +
    geom_point(aes(x=factor(site_number), y=d13C_mean, color=category, shape=category), size=4) +
  #   coord_cartesian(ylim=c(3, 10)) +
    ggtitle("Carbon: seston vs Tetragnathids (2011)") + 
    geom_errorbar(aes(x=factor(site_number), ymin=(d13C_mean - d13C_se), ymax=(d13C_mean + d13C_se)))

# Spiders vs aquatic inverts
  aqonly <- filter(sitemeans.w, sample_year==2012 & 
                     category %in% c("Aquatic", "Araneid", "Tetragnathid"))
  aqonly$category <- factor(aqonly$category, levels=c("Araneid", "Tetragnathid", "Aquatic"))

  ggplot(data=aqonly) +
    geom_errorbar(aes(x=factor(site_number), 
                      ymin=(dD_mean - dD_se), 
                      ymax=(dD_mean + dD_se)), width=0.2) +
    geom_point(aes(x=factor(site_number), y=dD_mean, color=category, shape=category), size=4) +
    scale_color_manual(values=c("darkgreen", "lightgreen", "blue")) +
    ggtitle("DEUTERIUM: Aquatic inverts vs Spiders (2012)") 

## Regression of spider vs other compartments
# First need a new semi-wide df, with categories as columns
  catmeans.w <- dcast(sitemeans, sample_year + chemical_name + site_number ~ category,
                       value.var="mean")

  # Separate data frames for each spider type
  aras <- select(catmeans.w, -Tetragnathid)
  aras$spider.type <- "Araneid"
  aras <- rename(aras, spider.mean=Araneid)
  aras <- aras[, c(1:4, 6:7, 5, 8)]  #spid columns at the end, to match tets

  tets <- select(catmeans.w, -Araneid)
  tets$spider.type <- "Tetragnathid"
  tets <- rename(tets, spider.mean=Tetragnathid)

  catmeans.h <- rbind(aras, tets)
  rm(aras)
  rm(tets)

# Spiders vs Aquatic inverts
  ggplot(data=catmeans.h[catmeans.h$sample_year==2012, ]) +  #Aquatic invert data only available in 2012
    geom_point(aes(x=Aquatic, y=spider.mean, shape=chemical_name, size=chemical_name, fill=spider.type)) +
    scale_fill_manual(values=c("black", "white")) +
    scale_shape_manual(values=c(21, 22, 25)) + 
    scale_size_manual(values=c(3, 2, 2)) + 
    xlab("Site mean aquatic invert delta") +
    ylab("Site mean spider delta") + 
    facet_wrap(spider.type ~ chemical_name, scales="free") +
    theme(legend.position="none") +
    ggtitle("Correlations between spider and aquatic invert isotopic signatures")

  ggplot(data=catmeans.h[catmeans.h$sample_year==2011, ]) +  #Aquatic invert data only available in 2012
    geom_point(aes(x=Basal, y=spider.mean, shape=chemical_name, size=chemical_name, fill=spider.type)) +
    scale_fill_manual(values=c("black", "white")) +
    scale_shape_manual(values=c(21, 22, 25)) + 
    scale_size_manual(values=c(3, 2, 2)) + 
    xlab("Site mean basal (seston) delta") +
    ylab("Site mean spider delta") + 
    facet_wrap(spider.type ~ chemical_name, scales="free") +
    theme(legend.position="none") +
    ggtitle("Correlations between spider and seston isotopic signatures")




#### SCRATCH PAPER ##############################
# Calculate standard error
  ranges$se <- ranges$sd / sqrt(ranges$n.smp)

# Convert data to a format that will be usable for graphing
  # First have to be fully long...
    range.long <- melt(ranges, measure.vars=c("mean", "se"))
    range.long <- range.long[, c(1:4, 9, 10)]

  # Now back to wide
    means <- dcast(range.long, sample_year + category ~ chemical_name + variable, value.var="value")
