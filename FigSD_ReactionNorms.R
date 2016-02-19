#############################
# Fig SD, "Reaction Norms" ##
#############################
# 
# Generates a 3 panel plot of site mean sum-PCB concentration at each site
# over time (x-axis = time); each site has a different shape, with color of 
# the points used to indicate location type.   
# Input: \sitemeans\ df from setup  
# 
# Output - one pdf/svg in Output folder: 
#   FigSD_ReactionNorms.pdf 
#   FigSD_ReactionNorms.svg 


# For custom color scheme
  library(RColorBrewer)

## A. Data prep ##
# Rename data frame to avoid having to rename objects in code. 
  pcb <- sitemeans   

# Add a column for location factor
  pcb$location[pcb$site_number %in% c(1,2,6,7)] <- "River"
  pcb$location[pcb$site_number %in% c(3,4,5)] <- "Backwater"
  pcb$location[pcb$site_number %in% c(8,9,11)] <- "Harbor"
  pcb$location[pcb$site_number %in% c(10)] <- "Lake"
  # Put location and category levels in order for display
    pcb$location <- factor(pcb$location, levels = c("River", "Backwater", "Harbor", "Lake"))
    pcb$category <- factor(pcb$category, levels=c("Sediment", "Araneid", "Tetragnathid"))

# Add a new site code (non-numeric) for aes mapping
  pcb$site <- paste("s", pcb$site_number, sep="")
  pcb$site <- factor(pcb$site, levels=c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11"))

# Create custom color scheme (should match colors in Fig. 1)
  set3.cols <- brewer.pal(4, "Set3")
  set3.cols[2] <- "gold"  #replace existing yellow (very pale) with a stronger color
  loc.cols <- set3.cols[c(1,1,4,4,4,1,1,3,3,2,3)]  #switched red and yellow so that backwater sites will be red and stand out more. 
  
  
## B. Plot
# Manually establish symbol shapes, colors, and labels...
  site.shapes <- c("s1"= 0, "s2"= 1, 
                   "s3"=21, "s4"=22, "s5"=23, 
                   "s6"= 2, "s7"= 6, 
                   "s8"=24, "s9"=25, 
                   "s10"=8, 
                   "s11"=5)

  site.colors <- c("s1"=loc.cols[1], "s2"=loc.cols[2], 
                   "s3"=loc.cols[3], "s4"=loc.cols[4], "s5"=loc.cols[5], 
                   "s6"=loc.cols[6], "s7"=loc.cols[7], 
                   "s8"=loc.cols[8], "s9"=loc.cols[9], 
                   "s10"=loc.cols[10], 
                   "s11"=loc.cols[11])

  site.labels <- c("s1"=1, "s2"=2, 
                   "s3"=3, "s4"=4, "s5"=5, 
                   "s6"=6, "s7"=7, 
                   "s8"=8, "s9"=9, 
                   "s10"=10, 
                   "s11"=11)

# Baseplot with standard formatting
  baseplot <- ggplot(data=pcb, aes(x=factor(sample_year), y=sum.pcb)) +
    geom_line (aes(group=site, lty=site, color=site)) + 
    scale_shape_manual(name="Site", values=site.shapes, labels=site.labels) +
    scale_color_manual(name="Site", values=site.colors, labels=site.labels) +
    scale_fill_manual (name="Site", values=site.colors, labels=site.labels) +
    scale_x_discrete(expand=c(0, 0.25)) + #moves data closer to y-axis. Not clear on what first value does.
    xlab("Year") + 
    guides(linetype=FALSE) + 
    theme_bw() +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          strip.text.x=element_text(size=11),
          strip.text.y=element_text(size=11), 
          legend.position="right",
          legend.key=element_rect(color=NA),
          axis.title.x=element_text(size=10.5),
          axis.title.y=element_text(size=10.5) ) +
      facet_wrap(~category)

  # Legend version (will still need modifying in Inkscape)
    leg.plot <- baseplot + 
      geom_point(aes(shape=site, 
                     fill=site, 
                     color=site), size=3) 
  
  # Sediment-scale version (for left panel)
    sed.plot <- baseplot + 
      geom_point(data=pcb[pcb$site_number %in% c(3,4,5,8,9), ], 
                 aes(shape=site, fill=site), 
                 color="black", size=3) + 
      geom_point(data=pcb[pcb$site_number %in% c(1,2,6,7,10,11), ], 
                 aes(shape=site, color=site), 
                 fill="white", size=3) +
      ylab(expression(Sigma * "PCB (ng g"^"-1"*" dry weight)")) + 
      theme(legend.position="none")
    
  # Spider-scale version (for right two panels)
    spid.plot <- sed.plot +
      coord_cartesian(ylim=c(0,540)) +
      ylab(expression(Sigma * "PCB (ng g"^"-1"*" wet weight)"))  
  

## C. Print output
  # PDF for easy viewing
    pdf(file=paste(DirOut, "FigSD_ReactionNorms.pdf", sep=""), width=8, height=12)
      grid.arrange(sed.plot, spid.plot, leg.plot, ncol=1)
    dev.off()
    
  # SVG for use in inkscape
    svg(file=paste(DirOut, "FigSD_ReactionNorms.svg", sep=""), width=8, height=12)
      grid.arrange(sed.plot, spid.plot, leg.plot, ncol=1)
    dev.off()

##### END SCRIPT ####