## Produce several standardized plots for figures in the MSQ
## data paper supplemental data.  Entire file can be sourced.
## 
## For each plot, a .pdf file and an .svg file are written to 
## the Output directory: 
## A. Fig SXa, FigSD_ZeroVsHalf.pdf 
##             Sum-half vs sum-zero for each year-category.
##             In two sections (three panels each) that will
##             have to be combined later.
## B. Fig SXb, FigSD_ConcVsMass.pdf
##             A single-panel plot containing all spider samples.
## C. Fig SXc, FigSD_SiteRanges_ConcMass.pdf
##             One six-panel plot, sample conc and sample mass for 
##             each site mean, grouped by year-cat for spider samples.
## D. Fig SXd, FigSD_SiteRankOrder.svg
##             One file combining two sets of four panels each.  
##             Site mean 2012-13 conc vs site mean 2011 conc, 
##             split by year-cat, and repeated for sum-zero totals
##             and sum-half totals.  
## E. Fig SXe, FigSD_FullVsRefConc.pdf
##             Two-panel plto, one for each spider cat, full vs ref
##             conc for each sample, combining data from 2012 and 2013. 


# For custom/standardized color sets 
  library(RColorBrewer)

  # Standard/custom color set for mapping site means
  # (designed for visual high contrast across all 11 sites)
    mycols <- c(brewer.pal(9, "Set1"), brewer.pal(4, "Dark2"))
    mycols <- mycols[c(1:10, 13)]

#### FIGURE SXa: CORRELATION, SUBS ZERO VS SUBS HALFDL ####
# Data including sums calculated with subs half mdl
  sums.data <- read.csv(paste(DirOut, "AllSampleData_CompareSums.csv", sep=""))
  sums.data$category <- factor(sums.data$category, levels=c("Sediment", "Araneid", "Tetragnathid"))

# Generate plots, separately for sed and spiders so axes can be different
# (will have to be combined in inkscape)
  baseplot <- ggplot(data=sums.data, aes(x=sum.zero, y=sum.halfdl)) +
    scale_fill_manual(values=mycols, name="Site") + 
    xlab(expression(Sigma * "PCB (ng g"^"-1"*"): zero substituted for non-detects")) + 
    ylab(expression(Sigma * "PCB (ng g"^"-1"*"): half detection limit substituted for non-detects")) + 
    facet_grid(category ~ sample_year, scales="free") + 
    theme_bw() + 
    theme(legend.position="right", 
          legend.key=element_rect(color=NA)) 

  # Sediment version
    plot.zero.half.sed <- baseplot + 
      geom_point(data=sums.data[sums.data$category=="Sediment", ], 
                 aes(fill=factor(site_number)), size=3, shape=21)  +
      geom_smooth(data=sums.data[sums.data$category=="Sediment", ], 
                  method=lm, se=FALSE, color="black", size=0.5, lty=5) +
      scale_x_continuous(breaks=c(10000, 20000, 30000), labels=comma) +
      scale_y_continuous(labels=comma)
  
  # Spider version
    plot.zero.half.spid <- baseplot + 
      geom_point(data=sums.data[sums.data$category!="Sediment", ], 
                 aes(fill=factor(site_number)), size=3, shape=21) +
      geom_smooth(data=sums.data[sums.data$category!="Sediment", ], 
                  method=lm, se=FALSE, color="black", size=0.5, lty=5)  

# Print output files
  # PDF (for easy viewing)
  pdf(file=paste(DirOut, "FigSD_ZeroVsHalf.pdf", sep=""), width=8, height=15)  #9 and 8 works well for single panel
    grid.arrange(plot.zero.half.sed, plot.zero.half.spid, ncol=1)
  dev.off()

  # SVG, for use in Inkscape
  svg(file=paste(DirOut, "FigSD_ZeroVsHalf.svg", sep=""), width=8, height=15)  #9 and 8 works well for single panel
    grid.arrange(plot.zero.half.sed, plot.zero.half.spid, ncol=1)
  dev.off()


#### FIGURE SXb: SAMPLE CONCENTRATION VS SAMPLE MASS ####
# The plot: sample conc as a function of sample mass, across all spider samples
  plot.conc.mass <-
    ggplot(data=samples[samples$category!="Sediment", ],
           aes(x=smp.mass, y=sum.pcb)) +
    geom_point(aes(fill=factor(sample_year)), shape=21, size=3) + 
    geom_smooth(method=lm, se=FALSE, color="black", size=0.5, lty=5) + 
    scale_fill_manual(values=c("gray", "white", "black"), name="Year") +
    coord_cartesian(xlim=c(0, 15)) + 
    xlab("Sample mass (g)") +
    ylab(expression("Sample "*Sigma * "PCB (ng g"^"-1"*" wet weight)")) +
    guides(shape=guide_legend(override.aes=list(size=3))) +   #control size of legend symbols
    theme_bw() + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.key=element_rect(color=NA))  #remove gray border around legend symbols


# Print output files
  # PDF (for easy viewing)
  pdf(file=paste(DirOut, "FigSD_ConcVsMass.pdf", sep=""), width=8, height=6)
    plot.conc.mass
  dev.off()

# SVG (for use in inkscape)
  svg(file=paste(DirOut, "FigSD_ConcVsMass.svg", sep=""), width=8, height=6)
    plot.conc.mass
  dev.off()

  
#### FIGURE SXc: RANGE IN SAMPLE MASS AND CONC PER SITE ####
# Generate a df with mean/min/max for smp.mass and sum.pcb
  spid.smp <- filter(samples, category!="Sediment")
  by.sitemean <- group_by(spid.smp, sample_year, category, site_number)
  spid.ranges <- summarize(by.sitemean, 
                         mean.conc = mean(sum.pcb, na.rm=TRUE),
                         min.conc  = min(sum.pcb, na.rm=TRUE),
                         max.conc  = max(sum.pcb, na.rm=TRUE), 
                         mean.mass = mean(smp.mass), 
                         sd.mass = sd(smp.mass),
                         min.mass = min(smp.mass),
                         max.mass = max(smp.mass)) 
  rm(by.sitemean); rm(spid.smp)
  

# Plot range in conc vs range in mass for each sitemean, by year/cat
  plot.siteranges <- ggplot(data=spid.ranges, aes(x=mean.mass, y=mean.conc)) + 
    geom_errorbar (aes(x=mean.mass, ymin=min.conc, ymax=max.conc)) +
    geom_errorbarh(aes(xmin=min.mass, xmax=max.mass)) + 
    geom_point(aes(fill=factor(site_number)), size=3, shape=21) + 
    #     geom_text(aes(x=mean.mass + 0.2, y=mean.conc-20, label=site_number), size=4) + 
    scale_fill_manual(values=mycols, name="Site") + 
    xlab("Sample mass (g)") + 
    ylab(expression("Sample "*Sigma * "PCB (ng g"^"-1"*" wet weight)")) +
    coord_cartesian(xlim=c(-0.1, 8), ylim=c(-20, 600)) + 
    facet_grid(category ~ sample_year) + 
    theme_bw() + 
    theme(legend.position="right", 
          legend.key=element_rect(color=NA))  

# Which max values are being cut off by the arbitrary plotting limits set on this graph?
  # Mass: max is set at 8g
  samples[(samples$smp.mass>8 & samples$category!="Sediment"), 1:7]
  # Conc: max is set at 600 ng/g
  samples[(samples$sum.pcb>600 & samples$category!="Sediment"), 1:7]

# Print output files
  # PDF (for easy viewing)
  pdf(file=paste(DirOut, "FigSD_SiteRanges_ConcMass.pdf", sep=""), width=10, height=6.5)
    plot.siteranges
  dev.off()

  # SVG for use in inkscape
  svg(file=paste(DirOut, "FigSD_SiteRanges_ConcMass.svg", sep=""), width=10, height=6.5)
    plot.siteranges
  dev.off()


#### FIGURE SXd: RANK ORDER CORRELATION IN SITE MEANS, 2011 vs other yrs ####
# Calculate site means using both subs.zero and subs.halfdl sum-PCB concentrations
# (!this section repeats code from <Analyses_CensoredData.R>)
  # Data including sums calculated with subs half mdl
    sums.data <- read.csv(paste(DirOut, "AllSampleData_CompareSums.csv", sep=""))
  # Use spider data only
    spidsums <- filter(sums.data, category!="Sediment")
  # Calculate site means
    by.site <- group_by(spidsums, sample_year, category, site_number)
    spid.means <- summarize(by.site,
                            sum.zero   = mean(sum.zero),
                            sum.halfdl = mean(sum.halfdl))
  
# Generate a hybrid wide/long df to match each 2012-13 sitemean with value from 2011
  # Full long form 
    spidmeans.l <- melt(spid.means, 
                        id.vars=c("sample_year", "category", "site_number"),
                        measure.vars=c("sum.zero", "sum.halfdl"), 
                        variable.name="sum.type", value.name="mean.conc") 
  
  # Split into separate dfs, in order to merge back into hybrid form
    spid11   <- filter(spidmeans.l, sample_year==2011)
    spid11   <- rename(spid11, conc.11 = mean.conc)
    spid1213 <- filter(spidmeans.l, sample_year!=2011)
    spidmeans.h <- merge(spid1213, spid11[, 2:5], 
                         by=c("category", "site_number", "sum.type"),
                         all.x=TRUE, all.y=FALSE)

  # Clean up 
    rm(spidsums); rm(by.site); rm(spid.means); rm(spidmeans.l); rm(spid11); rm(spid1213)

# Plots of the site means, 2011 site means vs equivalent from 2012/13
  baseplot <- ggplot(data=spidmeans.h, aes(x=conc.11, y=mean.conc)) +
    scale_fill_manual(values=mycols, name="Site") + 
    xlab(expression("Site mean 2011 "*Sigma * "PCB concentration (high resolution; ng g"^"-1"*" wet weight)")) + 
    ylab(expression("Site mean 2012-13 "*Sigma * "PCB concentration (low resolution; ng g"^"-1"*" wet weight)")) +
    facet_grid(sample_year ~ category) + 
    theme_bw() + 
    theme(legend.position="right", 
          legend.key=element_rect(color=NA))  

  # Subs-zero concentrations
    plot.rank.zero <- 
      baseplot +
      geom_point(data=spidmeans.h[spidmeans.h$sum.type=="sum.zero", ], 
                 aes(fill=factor(site_number)), size=3, shape=21) +
      theme(plot.margin= unit(c(1, 1, 2, 1), "lines"))

  # Subs-half-detlim concentrations
    plot.rank.halfdl <- 
      baseplot +
      geom_point(data=spidmeans.h[spidmeans.h$sum.type=="sum.halfdl", ], 
                 aes(fill=factor(site_number)), size=3, shape=21) +
      theme(plot.margin= unit(c(2, 1, 1, 1), "lines"))



# Print output files
  # PDF (for easy viewing; but ratios get messed up when importing to Inkscape)
    pdf(file=paste(DirOut, "FigSD_SiteRankOrder.pdf", sep=""), width=8, height=12)  #9 and 8 works well for single panel
      grid.arrange(plot.rank.zero, plot.rank.halfdl, ncol=1)
    dev.off()

  # SVG, for use in Inkscape
    svg(file=paste(DirOut, "FigSD_SiteRankeOrder.svg", sep=""), width=8, height=12)  #9 and 8 works well for single panel
      grid.arrange(plot.rank.zero, plot.rank.halfdl, ncol=1)
    dev.off()


#### FIGURE SXe: CORRELATION SUM.ZERO VS SUM.REF ####
# Data including sums of ref cons
  sums.data <- read.csv(paste(DirOut, "AllSampleData_CompareSums.csv", sep=""))
  # Use 2012-13 spider data only
  spidsums <- filter(sums.data, category!="Sediment" & sample_year!=2011)

# Two-panel plot, one for Ara and one for Tet
  plot.zerovsref <- ggplot(data=spidsums, aes(x=sum.ref.cons, y=sum.zero)) +
    geom_point(aes(fill=factor(sample_year)), shape=21, size=3) +
    geom_smooth(method=lm, se=FALSE, color="black", size=0.5, lty=5) + 
    scale_fill_manual(values=c("white", "black"), name="Year") +
    xlab(expression(Sigma * "PCB"[REF]*" (ng g"^-1*" wet weight)")) +
    ylab(expression(Sigma * "PCB"[FULL]*" (ng g"^-1*" wet weight)")) +
    theme_bw() + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.key=element_rect(color=NA),
          strip.text.x = element_text(size = 11)) + #remove gray border around legend symbols
    facet_wrap(~ category, scales="free")

# Print output files
  # PDF (for easy viewing; but ratios get messed up when importing to Inkscape)
  pdf(file=paste(DirOut, "FigSD_FullVsRefConc.pdf", sep=""), width=10, height=6)
    plot.zerovsref
  dev.off()

  # SVG, for use in Inkscape
  svg(file=paste(DirOut, "FigSD_FullVsRefConc.svg", sep=""), width=10, height=6)
    plot.zerovsref
  dev.off()



#### END SCRIPT ####