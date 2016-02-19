## Water Quality Data ##

# Read in summary w qual parameters for the MSQ sites
# Do a PCA to look for patterns by sites


#### READ IN AND PREPARE DATA ####
# Raw data files
  # Water quality
  raw.wq.data <- read.csv("WaterQualSummaryData.csv")
  
  # Horrible trace metals file for sodium
  metals.raw <- read.csv("11-30-11 Metals-undigested_MSQ.csv")
  
# Clean up the metals file in order to extract sodium data
  metals <- metals.raw[, 1:7]
  colnames(metals)[5:7] <- c("conc", "rsd", "sd.smp")
  # Filter to sodium rows only, don't want other metals for now.
    sodium <- filter(metals, substr(x=metals$Analyte.Name, start=1, stop=2)=="Na")  
  # Remove QA samples, only want standard field samples
    sodium <- filter(sodium, substr(x=sodium$Sample.ID, start=1, stop=2) %in% c("22", "23", "24"))  #should leave 45 rows
  # Split out separate fields for data and site number
    sodium$date.str <- substr(sodium$Sample.ID, start=1, stop=8)
    sodium$site.no  <- substr(sodium$Sample.ID, start=14, stop=15)  #take two digits, to include site 10
    sodium$site.no  <- gsub(pattern="-", replacement="", x=sodium$site.no, fixed=TRUE) #remove dashed from one-digit numbers
    sodium$site.no  <- as.numeric(sodium$site.no)

#### SUMMARIZE/CHECK DATA ####
# Water quality data - 
# Summary bar plots by parameter: am I getting same results as excel bar plots 
# shown in WaterQual file and in tech memo?
  param.barplot <- ggplot(data=raw.wq.data[raw.wq.data$parameter=="Temp", ], 
                          aes(x=factor(transect), y=mean)) + 
    geom_bar(stat="identity") + 
    geom_errorbar(aes(x=transect, ymax=(mean + stdev), ymin=(mean-stdev))) +
#     coord_cartesian(ylim=c(12.5, 16.5)) + 
    ggtitle("Temperature")

# Sodium data
  # How many data are there per site?
    table(sodium$site.no)  #n=5
  # Summary statistics for sodium
    na.sites <- summarize(group_by(sodium, site.no), 
                          n.dat= n(),
                          mean = mean(conc),
                          sd = sd(conc),
                          se = (sd/sqrt(n.dat)))
  # Summary sodium plot; am I getting same results as in the tech memo?
  # (Can't find any documentation of the data/procedure used to produce that plot.)
    na.means.plot <- ggplot(data=na.sites, aes(x=factor(site.no), y=mean)) +
      geom_point() +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))

# Reformat data
  # Convert main df to wide form (sites as rows, variables as columns) for PCA
  wq.w <- dcast(raw.wq.data, transect + location ~ parameter, value.var="mean")
  # Merge in sodium data
  na.sites <- rename(na.sites, Na=mean)
  wq.w <- merge(wq.w, na.sites[, c("site.no", "Na")],
                by.x="transect", by.y="site.no", 
                all=TRUE)
  
# Na vs SpCond:
# Is specific conductance, in fact, associated with Na conc as claimed in the tech memo text?
  ggplot(data=wq.w, aes(x=SpCond, y=Na)) + 
    geom_point(size=3) + 
    scale_y_log10()  #Note how linearity is improved with log.

  # Problem: missing Na data for site 8, this will be a problem for the PCA
  # Solution: calculate relationship between SpCond and Na, use rln to estimate an Na value for site 8.
    fit.na <- (lm(log10(Na) ~ SpCond, data=wq.w))
    site8.sp <- wq.w$SpCond[wq.w$transect==8]  #SpCond at site 8
    site8.na <- 10^(predict(object=fit.na, newdata=(data.frame(SpCond=0.2245))))  #...use this SpCond to predict Na.
    # Add this predicted value to the df
    wq.w[wq.w$transect==8, "Na"] <- site8.na


#### FORMAT AND PRINT SUMMARY TABLE ####
  wq.data <- raw.wq.data

# Rounding doesn't really work, number of dec places too inconsistent,
# I think we will need to just do this part in excel.
#   # Round to desired number of digits for mean and sd
#   wq.data$mean.sig <- signif(wq.data$mean, digits=3)
#   wq.data$sd.sig   <- paste("(", signif(wq.data$stdev, digits=2), ")", sep="")

  # Convert to wide format
    # But first have to melt to fully long form
    wq.long <- melt(wq.data, id.vars=c("transect", "location", "count", "parameter"),
                    measure.vars=c("mean", "stdev"), 
                    variable.name="stat", value.name="value")
    # Now cast
    wq.tab <- dcast(wq.long, transect+location+count ~ parameter+stat, value.var="value")

  # Merge in sodium data
#     # Round values
#     na.sites$na_mean.sig <- round(na.sites$mean, 2)
#     na.sites$na_sd.sig   <- round(na.sites$sd, 2)
    # Rename for consistent column names
    na.sites <- rename(na.sites, Na_mean=mean, Na_stdev=sd)
    # Add to wide-format table
    wq.tab <- merge(wq.tab, na.sites[, c("site.no", "Na_mean", "Na_stdev")],
                    by.x="transect", by.y="site.no", all=TRUE)
  
  # Rename columns for prettiness
    wq.tab <- rename(wq.tab, Site=transect, Location=location, n=count)

  # Write the table
  # (writing a .txt file in this case so that excel will not interpret parentheses as negatives)
    write.csv(wq.tab, paste(DirOut, "TabVals_WQuality.csv", sep=""), row.names=FALSE)
  
#### PCA ####
# Run the PCA
  # Without sodium (as in tech memo)
  env.pca.orig <- prcomp(x=wq.w[, c(4:6, 8:10)], scale=TRUE)  #this run (using [pH] and excluding [pHmV]) produces the same variance explained numbers as are reported in tech memo.
  # With sodium
  env.pca.na   <- prcomp(x=wq.w[, c(4:6, 8:11)], scale=TRUE, na.rm=TRUE)

# Variance explained by each component
  summary(env.pca.orig)  #essentiall no difference between the two
  summary(env.pca.na)

# Select which pca to use
  env.pca <- env.pca.orig
  #   env.pca <- env.pca.na

# Extract data in a form for plotting
  # Site scores (ie, coordinates on the component axes)
    env.scores <- as.data.frame(env.pca$x)
    env.scores$transect <- 1:10  #add site numbers for labeling
    env.scores <- cbind(env.scores, "location"=wq.w$location)  #add location factor for symbology
  # Variable vectors
    vects <- env.pca$rotation
    vects <- as.data.frame(vects)
    vects$parameter <- row.names(vects)


## Create a biplot 
  env.biplot <- ggplot(data=env.scores, aes(x=-PC1, y=PC2)) +
    geom_point(aes(shape=location, fill=location), size=4) +
    geom_text(aes(x=-PC1+0.1, y=PC2+0.1, label=transect)) + 
    geom_segment(data=vects,
                 aes(x=0, xend=-PC1*3,
                     y=0, yend= PC2*3),
                 arrow = arrow(length = unit(0.2, "cm")),
                 colour="grey") +   #, inherit_aes=FALSE) +
    xlab(paste("PC1 (", round(100*summary(env.pca)$importance[2,1], 0), "%)", sep="")) + 
    ylab(paste("PC2 (", round(100*summary(env.pca)$importance[2,2], 0), "%)", sep="")) + 
    geom_text(data=vects,
              aes(x=(-PC1*3), y=(PC2*3), label=parameter, vjust=0, hjust=1.1)) +
    scale_shape_manual(values=c(21, 0, 8, 24), name="Location") +
    scale_fill_manual(values=c("black", "grey", "black", "grey"), name="Location") +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(), 
          legend.key=element_rect(color=NA))

#### FORMAT AND PRINT A TABLE ####
# Table values

# pdf of biplot
  pdf(file=paste(DirOut, "FigSD_WQualBiplot.pdf", sep=""), width=8, height=6)
    env.biplot
  dev.off()

###########################################################
