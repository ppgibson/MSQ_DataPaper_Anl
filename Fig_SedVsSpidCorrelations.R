#########################################################################
# Fig "4", Regressions of Spider vs Sediment pcb concentrations
#########################################################################
# 
# Generates plots of the correlations between spider and sediment for different
# combinations of normalized vs -wt concentrations.  Plots have log scale
# and regression lines.
# Note: plots show raw data (with log-scaled axes), but the fitted coefficients  
#       (in the < Correlations_SedVsSpid.R > script) are for log10 values. 
# 
# Input: \sitemeans\ df from setup  
# 
# 
# Output: 
#  Three pdf files in Output folder:
#   - Fig_SedVsSpid_A_WetTOC.pdf
#   - Fig_SedVsSpid_B_LipTOC.pdf
#   - Fig_SedVsSpid_C_LipDry.pdf



# For working with log-scale axes
  library(scales)

## A. Data prep ##
## First, have to create a hybrid wide/long format data frame, 
## with sediment data and in column and matching spider data in other columns

# Split 'means' into two data frames, one for sediment and one for spiders
  means.sed  <- subset(sitemeans, category=="Sediment", 
                       select=c("sample_year", "site_number", 
                                "sum.pcb", "pcb.tocnorm"))
  means.sed <- rename(means.sed, sed.pcb=sum.pcb)  

  means.spid <- subset(sitemeans, category!="Sediment", 
                       select=c("sample_year", "category", "site_number", 
                                "sum.pcb", "pcb.lipidnorm"))
  means.spid <- rename(means.spid, spid.pcb=sum.pcb)

# Merge the spider and sediment data frames (each sediment data point will appear twice, once for Tet and once for Ara)
  means.h <- merge(means.spid, means.sed, by=c("sample_year", "site_number"), all=TRUE)   #"h" for hybrid, vs "l" for long

# Add a factor to control shape aesthetic by location type
  means.h$location[means.h$site_number %in% c(1,2,6,7)] <- "River"
  means.h$location[means.h$site_number %in% c(3,4,5)] <- "Backwater"
  means.h$location[means.h$site_number %in% c(8,9,11)] <- "Harbor"
  means.h$location[means.h$site_number %in% c(10)] <- "Lake"

# Sitemeans of zero are being excluded from the regressions
  means.nozero <- means.h
  means.nozero[means.nozero==0] <- NA  #note this doesn't delete the whole row, just replaces zeros with NAs - so that the other non-zero values can still be used.
  # But, we also want a df solely of zeros for adding the excluded points to the plot
  # (Sediment values are never zero, only spiders)
  means.zero <- filter(means.h, spid.pcb==0 | pcb.lipidnorm==0)

# Clean up
  rm(means.sed); rm(means.spid)

## B. Plotting ##
# Formatted labels for axes
  dry.label <- expression("Sediment " * Sigma * "PCB (ng g"^-1*" dry weight)")   #label for x-axis in plot
  toc.label <- expression("Sediment " * Sigma * "PCB"[TOC]*" (ng g"^-1*" organic carbon)")   #label for x-axis in plot
  wet.label <- expression("Spider " * Sigma * "PCB (ng g"^-1*" wet weight)")   #label for x-axis in plot
  lip.label <- expression("Spider " * Sigma * "PCB" [lipid] * " (ng g"^-1*" lipids)")   #label for y-axis in plot

# Baseplot to establish standard formatting parameters
  baseplot <- ggplot(data=means.nozero) +  
        scale_fill_manual(values=c("black", "white"), labels=c("Araneid", "Tetragnathid")) +
        scale_shape_manual(values=c(21, 22, 24, 25)) +
        guides(fill=guide_legend(override.aes=list(shape=21))) + 
        theme_bw() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour="black"),
              legend.title= element_blank(), 
              legend.key=element_rect(color=NA),  #remove gray border around legend symbols
              legend.position = "right")   

# Panel A: wet weight spiders vs TOC-normalized sediment
  wetplot <- baseplot +       
        geom_point(aes(x=pcb.tocnorm, y=spid.pcb, fill=category, shape=location), color="black", size=3.5) +
        geom_smooth(aes(x=pcb.tocnorm, y=spid.pcb), 
                    method=lm, se=FALSE, color="black", size=0.5, lty=5) + 
        geom_point(data=means.zero, 
                     aes(x=pcb.tocnorm, y=spid.pcb), shape=24, fill="grey", size=4) + #the zero points; they show up as tiny pts on x-axis, will need to fix this in inkscape.
        xlab(toc.label) + ylab(wet.label)
  
  (wetlog <- wetplot +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=4),  
                      labels = comma, limits=c(100, 1000000)) + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=4),  
                      labels = comma, limits=c(1, 1000)))  

# Panel B: lipid-norm spiders vs TOC-normalized sediment
  lipplot <- baseplot +       
        geom_point(aes(x=pcb.tocnorm, y=pcb.lipidnorm, fill=category, shape=location), color="black", size=3.5) +
        geom_smooth(aes(x=pcb.tocnorm, y=pcb.lipidnorm), 
                    method=lm, se=FALSE, color="black", size=0.5, lty=5) + 
        geom_point(data=means.zero, 
                     aes(x=pcb.tocnorm, y=pcb.lipidnorm), shape=24, fill="grey", size=4) + #the zero points; they show up as tiny pts on x-axis, will need to fix this in inkscape.
        xlab(toc.label) + ylab(lip.label)
  
  (liplog <- lipplot +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=4),  
                      labels = comma, limits=c(100, 1000000)) + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=3),  
                      labels = comma, limits=c(10, 12000)))  

# Panel C: alternate to panel A, lipid-norm spiders vs dry-wt sediment
  dryplot <- baseplot +       
        geom_point(aes(x=sed.pcb, y=pcb.lipidnorm, fill=category, shape=location), color="black", size=3.5) +
        geom_smooth(aes(x=sed.pcb, y=pcb.lipidnorm), 
                    method=lm, se=FALSE, color="black", size=0.5, lty=5) + 
        geom_point(data=means.zero, 
                     aes(x=sed.pcb, y=pcb.lipidnorm), shape=24, fill="grey", size=4) + #the zero points; they show up as tiny pts on x-axis, will need to fix this in inkscape.
        xlab(dry.label) + ylab(lip.label)
  
  (drylog <- dryplot +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=6),  
                      labels = comma, limits=c(0.8, 12000)) + 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=3),  
                      labels = comma, limits=c(10, 12000)))  
 
## C. Printing ##
## In inkscape, these files will need to be shrunk down to an apprpriate size
## (width ~ 4") and annotated with regression coefficients (etc) from 
## < Correlations_SedVsSpid.R >.

# Use obscure grid functions to ensure all plots take same plot area 
# despite different units/n.digits on y-axis 
# (code copied directly from http://stackoverflow.com/questions/24709307/keep-all-plot-components-same-size-in-ggplot2-between-two-plots)
  gl <- lapply(list(wetlog, liplog, drylog), ggplotGrob)
  library(grid)
  widths <- do.call(unit.pmax, lapply(gl, "[[", "widths"))
  heights <- do.call(unit.pmax, lapply(gl, "[[", "heights"))
  lg <- lapply(gl, function(g) {g$widths <- widths; g$heights <- heights; g})
  # extract elements 1, 2, and 3 of lg to draw plots (using grid.draw)  

# Panel A: wet wt vs TOC-norm
  pdf(file=paste(DirOut, "Fig_SedVsSpid_A_WetTOC.pdf", sep=""), width=8, height=5)
    grid.draw(lg[[1]])
  dev.off()

  svg(file=paste(DirOut, "Fig_SedVsSpid_A_WetTOC.svg", sep=""), width=8, height=5)
    grid.draw(lg[[1]])
  dev.off()

# Panel B: lip-norm wt vs TOC-norm
  pdf(file=paste(DirOut, "Fig_SedVsSpid_B_LipTOC.pdf", sep=""), width=8, height=5)
    grid.draw(lg[[2]])
  dev.off()

  svg(file=paste(DirOut, "Fig_SedVsSpid_B_LipTOC.svg", sep=""), width=8, height=5)
    grid.draw(lg[[2]])
  dev.off()

# Panel 'C': lip-norm vs dry wt
  pdf(file=paste(DirOut, "Fig_SedVsSpid_C_LipDry.pdf", sep=""), width=8, height=5)
    grid.draw(lg[[3]])
  dev.off()

  svg(file=paste(DirOut, "Fig_SedVsSpid_C_LipDry.svg", sep=""), width=8, height=5)
    grid.draw(lg[[3]])
  dev.off()

## D. Extra stuff ##
# # 'Standard' parameters for the various plots  
#   # Establish some variables based on function inputs
#     if (sed=="dry") {
#       xcol <- 6  #column to use for x (ie, sediment) data
#       xcap <- expression("Sediment " * Sigma * "PCBs (ng g"^-1*" dry weight)")   #label for x-axis in plot
#       xbreaks <- 6
#       xlimits <- c(0.8, 12000)
#       xpos <- 1000
#      } 
#   
#     if (sed=="toc") {
#       xcol <- 7  #column to use for x (ie, sediment) data
#       xcap <- expression("Sediment " * Sigma * "PCBs"[TOC]*" (ng g"^-1*" organic carbon)")   #label for x-axis in plot
#       xbreaks <- 4
#       xlimits <- c(100, 1000000)
#       xpos <- 100000
#     } 
#   
#     if (spid=="wet") {
#       ycol <- 4  #column to use for y (ie, spider) data
#       ycap <- expression("Spider " * Sigma * "PCBs (ng g"^-1*" wet weight)")   #label for x-axis in plot
#       ybreaks <- 4
#       ylimits <- c(1, 1000)
#       ypos <- 5
#     } 
#       
#     if (spid=="lip") {
#       ycol <- 5  #column to use for y (ie, spider) data
#       ycap <- expression("Spider " * Sigma * "PCBs" [lipid] * " (ng g"^-1*" lipids)")   #label for y-axis in plot
#       ybreaks <- 3
#       ylimits <- c(100, 12000)
#       ypos <- 200
#     } 

###############################################################################################