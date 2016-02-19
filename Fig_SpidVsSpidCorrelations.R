#########################################################################
# Fig "3", Regression of Araneid vs Tetragnathid pcb concentration
#########################################################################
# 
# Generates plot[s] of the correlation between Tet vs Ara total pcb concentration 
# for each year 2011-2013 and for all years combined, including a regression lines
# for each correlation.  
# - **plot shows raw data (with log-scaled axes), but the fitted coefficients (in the 
#     < Correlations_SpidVsSpid.R > script are for  [log10(Tet) ~ log10(Ara)]. 
# - Calculations use full measured values (everything below MDL is censored), no 
#   modeled values, zero substituted for nondetects.
#
# Note: a more elegant way to do regression/plotting by year is described at 
# <http://stackoverflow.com/questions/24983690/nls-and-log-scale-in-ggplot2>
#
# Inputs: \means.w\ df from setup
#
# Output: 
#  Two pdf files in Output folder:
#   - Fig_SpidVsSpid_LogScale.pdf
#   - Fig_SpidVsSpid_BaseScale.pdf


# For working with log-scale axes
  library(scales)

# Data set excluding zero-sitemeans
  means.nona <- filter(means.w, !is.na(Tetragnathid) & !is.na(Araneid)) #first, removes the 4 data points that have at least one NA obs.
  means.nozero <- filter(means.nona, Tetragnathid>0 & Araneid>0)        #then, remove the additional two points where one or both spiders is zero.

# Basic plot (data points plus reg line for each year)
  base.spidplot <- ggplot(data=means.nozero, aes(x=Araneid, y=Tetragnathid)) +
    geom_point(aes(shape=factor(sample_year), size=factor(sample_year))) + 
    geom_smooth(aes(linetype=factor(sample_year)), method=lm, se=FALSE, color="black", size=0.5) +
    scale_shape_manual(values=c(0,19,8), name="Year") +
    scale_linetype_manual(values=c(2,4,5)) +
    scale_size_manual(values=c(2.5,3.50,3), guide="none") + 
    xlab(expression("Araneid "*Sigma * "PCB (ng g"^"-1"*" wet weight)")) +
    ylab(expression("Tetragnathid "*Sigma * "PCB (ng g"^"-1"*" wet weight)")) +
    guides(shape=guide_legend(override.aes=list(size=3))) +   #control size of legend symbols
    theme_bw() + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.key=element_rect(color=NA),  #remove gray border around legend symbols
          legend.key.width=unit(2, "cm"))    #make linetypes in legend longer so pattern is visible.
  #         panel.border=element_blank(),
  #         axis.line=element_line(color="black"), )
    
  # Add overall regression line
    spidplot <- 
      base.spidplot + 
      geom_smooth(method=lm, se=FALSE, size=0.5, lty=1, color="black")
  
  # Add one:one line [?]
   spidplot <- 
    spidplot + 
    geom_segment(aes(x=7, y=7, xend=500, yend=500), color="black", lty=3)

  # Add log scales
    log.spidplot <- 
      spidplot + 
      scale_x_log10(limits=c(6,600), breaks=c(10,50,100,500), labels=c(10, "", 100, ""), oob=squish) +
      scale_y_log10(limits=c(6,600), breaks=c(10,50,100,500), labels=c(10, "", 100, ""), oob=squish) 

  # Add excluded 'zero' points
    # Log plot
    log.spidplot <- 
      log.spidplot +  
      geom_point(data=means.nona[means.nona$Tetragnathid==0 | means.nona$Araneid==0, ], 
                 shape=21, fill="grey", size=4) 
    # Untransformed-axes plot (formatting of pt is slightly different)
    spidplot <- 
      spidplot +  
      geom_point(data=means.nona[means.nona$Tetragnathid==0 | means.nona$Araneid==0, ], 
                 shape=19, color="grey", size=3.5) +
      coord_cartesian(xlim=c(20, 500), ylim=c(20, 550))  #to avoid lots of extra space past 0
 

## Print pdf files
## In inkscape, these files will need to be shrunk down to an apprpriate size
## (width ~ 4") and annotated with regression coefficients (etc) from 
## < Correlations_SpidVsSpid.R >.

# Log plot (the main-text figure)
  # PDF for easy viewing
  pdf(file=paste(DirOut, "Fig_SpidVsSpid_LogScale.pdf", sep=""), width=8, height=6)
    log.spidplot
  dev.off()

  # SVG for inskscape
  svg(file=paste(DirOut, "Fig_SpidVsSpid_LogScale.svg", sep=""), width=8, height=6)
    log.spidplot
  dev.off()

# Untransformed axes (possible appendix figure)
  # PDF for easy viewing
  pdf(file=paste(DirOut, "Fig_SpidVsSpid_BaseScale.pdf", sep=""), width=8, height=6)
    spidplot
  dev.off()

  # SVG for Inkscape
  svg(file=paste(DirOut, "Fig_SpidVsSpid_BaseScale.svg", sep=""), width=8, height=6)
    spidplot
  dev.off()

#####################################################################################