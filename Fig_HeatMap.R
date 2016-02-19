##############################################################
# Fig "2", Grand mean concentrations overlaid on spatial map #
##############################################################
# 
# - Brings in the latest google earth satellite imagery of Manistique,
# coverts it into a ggplot layer; 
# - calculates grand mean concentrations for each site/category; 
# - splits grand mean concentrations (dry/wet weight and toc-normalized)
#   into semi-arbitrary categories - currently orders of magnitude, except
#   for final spider category;
# - plot symbol indicating concentration (by color and size) in the location
#   of each site, for sediment; toc-normalized sediment; araneid; tetragnathid.  
#
# Input: \sitemeans\ df from setup  
#   plus <Man_FW_SiteCoords.csv">, lat/long values for each MSQ FW site
#    and satellite imagery from google earth
# 
# Output: 
#  Four large .tiff files in Output folder:
#   - Fig2_Heatmap_A_Sed_13.tif
#   - Fig2_Heatmap_B_TOC_13.tif
#   - Fig3_Heatmap_C_Ara_13.tif
#   - Fig4_Heatmap_D_Tet_13.tif


# For working with google earth imagery
  library(ggmap)

## A. Data prep ##
# Read in coordinate data
  coords <- read.csv("Man_FW_SiteCoords.csv")
  coords <- rename(coords, long=alt_x_coord, lat=alt_y_coord)
  coords[1, 4:5] <- c(-86.254, 45.959)   #replace site 1 with fake coords so that the symbol for site 1 will plot within the range of the main map.  
  
# Get GoogleEarth imagery and save as a ggplot object
  msq.zoom14.bw  <- get_googlemap(center=c(long=-86.252769, lat=45.958286), 
                       zoom=14, maptype="satellite", color="bw")   #or 'color' instead of 'bw'
  msq.map.bw  <- ggmap(msq.zoom14.bw)

# Calculate grand means
  bysite  <- group_by(sitemeans, category, site_number)
  grmeans <- summarize(bysite, 
            mean.sum.pcb = mean(sum.pcb), 
            mean.pct.toc = mean(pct.toc, na.rm=TRUE), 
            mean.pcb.tocnorm = mean(pcb.tocnorm, na.rm=TRUE),
            mean.pct.lip = mean(pct.lip, na.rm=TRUE), 
            mean.pcb.lipnorm = mean(pcb.lipidnorm, na.rm=TRUE) )
  rm(bysite)      
  
  # Add site coordinates to each grand mean
    grmeans <- merge(grmeans, coords, by="site_number", all.x=TRUE, all.y=FALSE)

# Convert concentrations to factor levels
  # Split out separate data frames for spiders and sediment for simpler manipulation
  # and plotting (since each facet has different scale)
    grmeans.spid <- filter(grmeans, category!="Sediment")
    grmeans.sed  <- filter(grmeans, category=="Sediment")

  # Establish categorical levels for concentration symbology, 
  # depending on concentration type
    sedbreaks.dry  <- c(0, 10, 100, 1000, 10000)
    sedbreaks.toc  <- c(0,     100, 1000, 10000, 100000, 150000)
    spidbreaks.wet <- c(0, 10, 100, 200, 300)
  
  # Add new columns for factorized concentration values
    grmeans.spid <- mutate(grmeans.spid, 
                           sum.pcb.cat=cut(mean.sum.pcb, breaks=spidbreaks.wet),
                           norm.pcb.cat=NA)  #leave normalized column as NA for spiders; we don't need actual data here, just a column so it can merge with the sediment data.
    grmeans.sed  <- mutate(grmeans.sed,  
                           sum.pcb.cat=cut(mean.sum.pcb, breaks=sedbreaks.dry),
                           norm.pcb.cat=cut(mean.pcb.tocnorm, breaks=sedbreaks.toc))
  # Clean up
    rm(sedbreaks.dry); rm(sedbreaks.toc); rm(spidbreaks.wet); 

       
## B. Plotting ##
# Annotation/scales standard values
  # Titles for legends
  sedtitle <- bquote(atop("Sediment " * Sigma * "PCB", 
                          "(ng g"^-1*" dry weight)"))
  toctitle <- bquote(atop("Sediment " * Sigma * "PCB", 
                          "(" * mu * "g g"^-1*" organic carbon)"))
  aratitle <- bquote(atop("Araneid " * Sigma * "PCB", 
                          "(ng g"^-1*" wet weight)"))
  tettitle <- bquote(atop("Tetragnathid " * Sigma * "PCB", 
                          "(ng g"^-1*" wet weight)"))

  # Names for concentration categories
  sedcats <- c("< 10", "10 - 100", "100 - 1000", "> 1000")
  toccats <- c("< 1", "1 - 10", "10 - 100", "> 100")
  spidcats<- c("< 10", "10 - 100", "100 - 200", "> 200")
  
  # Consistent color scale
  colvalues  <- c("green", "yellow", "orange", "red")
  sizevalues <- c(1.5,2,3,4)
  
# Standard base layer/formatting, on which to plot points for each category
  basemap <- msq.map.bw + 
    geom_text (data=grmeans.sed,
               aes(x=long-0.0003, y=lat+0.0005, label=site_number), 
               color="white", size=3) + 
    theme_bw() + 
    theme(legend.text.align=1,
          legend.text=element_text(size=8),
          legend.title=element_text(size=9, face="bold"), 
          legend.key=element_rect(color=NA) )
  
# Now specific plots for the four category/normalization combos, one at a time.
  # Sediment dry weight
  sed.dry.map <- 
    basemap +
    geom_point(data=grmeans.sed, 
               aes(x=long, y=lat, fill=sum.pcb.cat, size=sum.pcb.cat), 
               shape=21) +
    scale_fill_manual(name=sedtitle, values=colvalues,  labels=sedcats) +
    scale_size_manual(name=sedtitle, values=sizevalues, labels=sedcats)

  # Sediment TOC-norm
  sed.toc.map <- 
    basemap +
    geom_point(data=grmeans.sed, 
               aes(x=long, y=lat, fill=norm.pcb.cat, size=norm.pcb.cat), 
               shape=21) +
    scale_fill_manual(name=toctitle, values=colvalues,  labels=toccats) +
    scale_size_manual(name=toctitle, values=sizevalues, labels=toccats)

  # Araneid wet weight
  ara.wet.map <- 
    basemap +
    geom_point(data=grmeans.spid[grmeans.spid$category=="Araneid", ], 
               aes(x=long, y=lat, fill=sum.pcb.cat, size=sum.pcb.cat), 
               shape=21) +
    scale_fill_manual(name=aratitle, values=colvalues,  labels=spidcats) +
    scale_size_manual(name=aratitle, values=sizevalues, labels=spidcats)

  # Tetragnathid wet weight
  tet.wet.map <- 
    basemap +
    geom_point(data=grmeans.spid[grmeans.spid$category=="Tetragnathid", ], 
               aes(x=long, y=lat, fill=sum.pcb.cat, size=sum.pcb.cat), 
               shape=21) +
    scale_fill_manual(name=tettitle, values=colvalues,  labels=spidcats) +
    scale_size_manual(name=tettitle, values=sizevalues, labels=spidcats)
    
## C. Write output files. (large, high res) ##
  # Sediment dry weight
  ggsave(filename=paste(DirOut, "Fig2_Heatmap_A_Sed_13.tif", sep=""),
         plot=sed.dry.map,
         width=13, height=13, units="in", 
         device="tiff")
  
  # Sediment TOC-norm
  ggsave(filename=paste(DirOut, "Fig2_Heatmap_B_TOC_13.tif", sep=""),
         plot=sed.toc.map,
         width=13, height=13, units="in", 
         device="tiff")
  
  # Araneid
  ggsave(filename=paste(DirOut, "Fig2_Heatmap_C_Ara_13.tif", sep=""),
         plot=ara.wet.map,
         width=13, height=13, units="in", 
         device="tiff")

  # Tetragnathid
  ggsave(filename=paste(DirOut, "Fig2_Heatmap_D_Ara_13.tif", sep=""),
         plot=tet.wet.map,
         width=13, height=13, units="in", 
         device="tiff")
  
##############################################