#########################################################################
# Figure, boxplots showing distribution of Total.PCB for each site/year
#########################################################################
# 
# Creates a 3x3 panel of boxplots, showing the distribution of total PCB 
# concentration at each site, in each year, for each category.  Each 
# individual 'box' is composed of n=4 (or fewer) individual samples.
# This script uses total measured PCB values for all years (no normalization).
# 
# Input: \samples\ df from setup  
# 
# Output: two ggplot objects, 'sed' and 'spid', formatted boxplot figures
#   plus, some code for combining the two objects into a single figure
# 

# Rename data frame/main pcb column to avoid having to rename objects in code. 
  pcb <- samples    
  pcb <- rename(pcb, Total.PCB=full)  #We want to use the [full] column for EPCBs.

# # Add a column for location factor
#   pcb$location[pcb$site_number %in% c(1,2,6,7)] <- "River"
#   pcb$location[pcb$site_number %in% c(3,4,5)] <- "Backwater"
#   pcb$location[pcb$site_number %in% c(8,9,11)] <- "Harbor"
#   pcb$location[pcb$site_number %in% c(10)] <- "Lake Michigan"
#   # Put location levels in order
#     pcb$location <- factor(pcb$location, levels = c("River", "Backwater", "Harbor", "Lake Michigan"))

# Count sample size/quartile data for each boxplot (for displaying n= on plot)
  samplesize <- ddply(.data = samples, 
                       .variables = c("sample_year", "category", "site_number"),
                       .fun = summarise,
                        n = sum(!is.na(full)),
                        Q1 = quantile(full, 0.25, na.rm=TRUE),
                        Q2 = quantile(full, 0.50, na.rm=TRUE),
                        Q3 = quantile(full, 0.75, na.rm=TRUE),
                        IQR = IQR(full, na.rm=TRUE) )
  
  # Compute whisker ends separately, since it takes a special function
    whiskers <- ddply(.data = samples,
                      .variables = c("sample_year", "category", "site_number"),
                      function (x) max(fivenum(x$full, na.rm=TRUE)))
    colnames(whiskers)[4] <- "whisker"
  
  # Merge site-level information into one df  
    samplesize <- merge(samplesize, whiskers, by=c("sample_year", "category", "site_number"), all=TRUE)

  # Create a special column to adjust label position based on magnitude (to account for log-scaled axes) 
    samplesize$labelspot[samplesize$whisker >  10000] <- samplesize$whisker[samplesize$whisker >  10000] + log(samplesize$whisker[samplesize$whisker >  10000])*1000
    samplesize$labelspot[samplesize$whisker <= 10000] <- samplesize$whisker[samplesize$whisker <= 10000] + log(samplesize$whisker[samplesize$whisker <= 10000])*500
    samplesize$labelspot[samplesize$whisker <= 5000]  <- samplesize$whisker[samplesize$whisker <= 5000]  + log(samplesize$whisker[samplesize$whisker <= 5000]) *100
    samplesize$labelspot[samplesize$whisker <= 900]   <- samplesize$whisker[samplesize$whisker <= 900]   + log(samplesize$whisker[samplesize$whisker <= 900]) *50  #30 if it's log + 1
    samplesize$labelspot[samplesize$whisker <= 500]   <- samplesize$whisker[samplesize$whisker <= 500]   + log(samplesize$whisker[samplesize$whisker <= 500]) *30  #10 if it's log + 1
    samplesize$labelspot[samplesize$whisker <= 100]   <- samplesize$whisker[samplesize$whisker <= 100]   + log(samplesize$whisker[samplesize$whisker <= 100])  *4
    samplesize$labelspot[samplesize$whisker <= 10]    <- samplesize$whisker[samplesize$whisker <= 10]    + log(samplesize$whisker[samplesize$whisker <= 10])   *2
#     samplesize$labelspot[samplesize$whisker <= 10]    <- samplesize$whisker[samplesize$whisker <= 10]    + log(samplesize$whisker[samplesize$whisker <= 10])   *2
    # Manually set labelspot for -Inf 
    samplesize$labelspot[samplesize$whisker == 0]     <- 0.2   #1.2 if it's log + 1
    
# Create an index to subset overall datframe to category of interest 
# (have to do sediment and spiders separately because of different scales)
  # one for spiders
    spidindex <- which(pcb$category!="Sediment")
    spid.data <- pcb[spidindex, ]
    row.names(spid.data) <- NULL
  # and one for sediment
    sedindex <- which(pcb$category=="Sediment")


## Create and store separate plots for spiders and for sediment
## Plot Total.PCB + 1 to allow plotting 0 values (ie, non-detects)
# 1. Spiders
    spid <- ggplot(data=spid.data, 
                aes(x=as.factor(site_number), y=(Total.PCB + 0.1), fill=factor(site_number))) +
              geom_boxplot() + 
              
              xlab("Site") +
              ylab(expression(Sigma * "PCBs (ng g"^"-1"*" wet weight)")) +
              scale_x_discrete(limits = levels(factor(pcb$site_number)) ) +
              scale_y_log10(breaks = 10^(0:3), labels=comma) + 
              scale_fill_brewer(palette="Set3") +
              theme_bw() +
              theme(strip.text.x=element_text(size=13),
                    strip.text.y=element_text(size=13) ) +
#                     legend.position="top") +
#                     axis.title.x=element_text(size=12),
#                     axis.title.y=element_text(size=12) ) +
              guides(fill=FALSE) +                  #get rid of this to regain legend
              facet_grid(category ~ sample_year)

  # Hackishly add grey points to indicate 0 values
  spid <- 
    spid + geom_point(data=spid.data[c(68, 107, 138), ], shape=21, size=3, #grey circles for individual zero points
                      color="black", fill="gray33")  
#     geom_point(data=spid.data[c(54,87,92,187), ], shape="-", size=15, color="grey") #grey rectangles for when the whole distribution is zero.

  # Add annotation of sample size above each box (approx. at whisker end)
    spid <- 
      spid + geom_text(data=samplesize[(samplesize$category!="Sediment" & samplesize$n<4), ],
                       aes(x=as.factor(site_number), y=labelspot, label=n), size=4)


# 2. Sediment
    sed <- ggplot(data=pcb[sedindex, ], 
                 aes(x=as.factor(site_number), y=(Total.PCB + 0.1), fill=factor(site_number))) +
            geom_boxplot() +
            xlab("") +  #remove axis label, since it appears in the spider half, but retain the white space where the label would go.
            ylab(expression(Sigma * "PCBs (ng g"^"-1"*" dry weight)")) +
            scale_x_discrete(limits = levels(factor(pcb$site_number)) ) +
            scale_y_log10(breaks = 10^(0:4), labels=comma) +
            scale_fill_brewer(palette="Set3") +
            theme_bw() +
            theme(strip.text.x=element_text(size=13),
                  strip.text.y=element_text(size=13) ) +
#                   axis.title.x=element_text(size=12),
#                   axis.title.y=element_text(size=12)) +
            guides(fill=FALSE) +                  #get rid of this line in order to regain legend
            facet_grid(category ~ sample_year)    #years side by side (x-axis)

  # Add annotation of sample size above each box (approx. at whisker end)
    sed <- 
      sed + geom_text(data=samplesize[(samplesize$category=="Sediment" & samplesize$n<4), ],
                      aes(x=as.factor(site_number), y=labelspot, label=n), size=4)
  
# Adjust spacing to improve vertical lineup of the two parts of the figure
  sed  <- sed  + theme(plot.margin=unit(c(1,1,0,-0.2), "cm"), axis.title.y = element_text(vjust=0))
  spid <- spid + theme(plot.margin=unit(c(0,1,1,0), "cm"), axis.title.y = element_text(vjust=0.5))
    

## Create a combined figure
  # Option 1: use grid.arrange, lots of space between parts
    grid.arrange(sed, spid, ncol=1, heights=c(0.38, 0.62))   #for three approximately equal plot heights
    grid.arrange(sed, spid, ncol=1, heights=c(0.45, 0.55))   #for sediment taller than spider plots
    # save as image, 1000 x 800
  
  # Option 2: "print" the two parts separately and then manually arrange them later
    sed    # 900? X 308
    spid   # 900? x 503

# Print the file as a high quality .emf
  emf(file=paste(DirOut, "Fig1_PlusPoint1_v2.emf", sep=""), width=10, height=14)
#   emf(file=paste(DirOut, "Fig1_bysite_small.emf", sep=""), width=6.5, height=8.165, pointsize=10)
# 
#   tiff(filename="test_fig1.tiff",
#        width=10, height=13, units="in", pointsize=10, compression="lzw", res=300)

  grid.arrange(sed, spid, ncol=1, heights=c(0.38, 0.62))

  dev.off()

# Clean up
  rm(pcb)
  rm(sedindex)
  rm(spidindex)
  rm(sed)
  rm(spid)


tiff(filename=paste(DirOut, "tes_fig1.tiff", sep=""), 
     width=10, height=13, units="in", pointsize=10, compression="lzw", res=300)

################################
#   sed + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=4))
#   sed + scale_y_log10(breaks = 10^(0:4), labels=comma)
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=5),  #zoomed in: n=4
#                 labels = comma, limits=c(100, 1000000)) +               #zoomed in: c(1000, ...)
#     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=4),  #zoomed in: c(100, ...)
#                   labels = comma, limits=c(1, 10000)) +                   #zoomed in: n=3
#     

# #  Add points on top of boxplot, if desired
#     sed + geom_point(shape=21, size=4)  # add points on top of the boxplot