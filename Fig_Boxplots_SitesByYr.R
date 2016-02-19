#########################################################################
# Figure '1', points (or boxplots) showing the distribution of sample 
# sum.PCB values for each site/year.
#########################################################################
# 
# Creates a 3x3 panel of boxplot-like dotplots, showing the distribution of 
# sum-PCB concentration at each site, in each year, for each category. Each 
# individual site/year/category combination is composed of n=4 or fewer individual 
# samples/data points.  Data points are colored according to location category. 
# 
# Input: \samples\ df from setup  
# 
# Output is one pdf file in Output directory:
#   - Fig1_SampleConc_SitesByYr.pdf

## A. Data prep ##
# Rename data frame to avoid having to rename objects in code. 
  pcb <- samples    

# Add a column for location factor
  pcb$location[pcb$site_number %in% c(1,2,6,7)] <- "River"
  pcb$location[pcb$site_number %in% c(3,4,5)] <- "Backwater"
  pcb$location[pcb$site_number %in% c(8,9,11)] <- "Harbor"
  pcb$location[pcb$site_number %in% c(10)] <- "Lake"
  # Put location levels in order
    pcb$location <- factor(pcb$location, levels = c("River", "Backwater", "Harbor", "Lake"))

# Reorder category levels for appropriate display
  pcb$category <- factor(pcb$category, levels=c("Sediment", "Araneid", "Tetragnathid"))


## B. Create plots ##
## This approach uses facet_grid function to ensure equal spacing for all panels; however,
## spider vs. sediment panels require slightly different y-axis formatting, so first all 
## three rows are plotted using spider formatting parameters, then all three rows are 
## plotted using sediment formatting parameters, then the sediment row in the spider plot 
## is replaced, piece-by-piece, with the equivalent row from the sediment plot. 
## Reference: http://stackoverflow.com/questions/27929549/different-y-axis-labels-facet-grid-and-sizes?answertab=active#tab-top

# Baseplot to establish standard formatting parameters
  baseplot.samples <- ggplot(data=pcb[pcb$sum.pcb!=0, ], #we'll plot zeros separately
                aes(x=as.factor(site_number), y=(sum.pcb), fill=location)) +
              geom_point(shape=21, size=4, alpha=0.6) +    
              geom_point(data=pcb[pcb$sum.pcb==0, ],     #add the nondetect (zero) samples with their own formatting.
                         shape=21, size=4, 
                         alpha=0.6, position=position_jitter(width=0.3, height=0)) + 
              xlab("Site number") +
              ylab(expression(Sigma * "PCB (ng g"^"-1"*")")) +
              scale_fill_brewer(palette="Set3", name="Location") +
              annotation_logticks(sides="l") +
              theme_bw() +
              theme(strip.text.x=element_text(size=11),
                    strip.text.y=element_text(size=11), 
                    legend.position="right",
                    legend.key=element_rect(color=NA),
                    axis.title.x=element_text(size=10.5),
                    axis.title.y=element_text(size=10.5) ) +
              facet_grid(category ~ sample_year)

# All three rows with axis formatting for spiders
  plot.spidformat <- 
    baseplot.samples + 
    scale_x_discrete(limits = levels(factor(pcb$site_number)) ) +
    scale_y_log10(breaks = 10^(0:3), labels=comma) + 
    coord_cartesian(ylim=c(0.7, 1000))  

# All three rows with axis formatting for sediment
  plot.sedformat <- 
    baseplot.samples + 
    scale_x_discrete(limits = levels(factor(pcb$site_number)) ) +
    scale_y_log10(breaks = 10^(0:4), labels=comma) +
    coord_cartesian(ylim=c(0.2, 50000)) 
  
# Hackishly superimpose sediment features onto top row of spider plot
  gspid <- ggplotGrob(plot.spidformat)
  gsed  <- ggplotGrob(plot.sedformat)
  gboth <- gspid

  gboth[["grobs"]][[c(2)]]  <- gsed[["grobs"]][[2]]
  gboth[["grobs"]][[c(5)]]  <- gsed[["grobs"]][[5]]
  gboth[["grobs"]][[c(8)]]  <- gsed[["grobs"]][[8]]
  gboth[["grobs"]][[c(11)]] <- gsed[["grobs"]][[11]]
  gboth[["grobs"]][[c(14)]] <- gsed[["grobs"]][[14]]
  gboth[["grobs"]][[c(17)]] <- gsed[["grobs"]][[17]]
  gboth[["grobs"]][[c(20)]] <- gsed[["grobs"]][[20]]

# Print output files of the Frankenstein figure
  # PDF for easy viewing
  pdf(paste(DirOut, "Fig1_SampleConc_SitesByYr.pdf", sep=""), width=10, height=10)
    grid.newpage()
    grid.draw(gboth)
  dev.off()

  # SVG for use in Inkscape
  svg(paste(DirOut, "Fig1_SampleConc_SitesByYr.svg", sep=""), width=10, height=10)
    grid.newpage()
    grid.draw(gboth)
  dev.off()

#############################################################################