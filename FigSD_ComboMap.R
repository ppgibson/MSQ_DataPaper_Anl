########################################################################
# Fig S1, Base map for combined figure of Manistique area and features #
########################################################################
# 
# - Brings in the latest google earth satellite imagery of Manistique,
#   coverts it into a ggplot layer; 
# - overlays shapefile of AOC boundary and points for each sample site;  
# - adds a scale bar.  
# - Also, creates a simple outline map of the US/Canada to show location.
#
# Inputs: 
#   - <Man_FW_SiteCoords.csv">, lat/long values for each MSQ FW site
#   - satellite imagery from google earth (pulled by ggmap package)
#   - shapefile of MSQ AOC boundary: Output/aoc_mi_manistique/  
# 
# Output: 
#   - Raster of the annotated MSQ site image (incluing sites, AOC outline, 
#     scale bar): FigSD_ComboMap.tif
#   - Simple map of North America w/ Great Lakes: FigSD_NAm_Map.pdf / .svg


# Special spatial packages
  library(rgeos)
  library(maptools)
  library(ggmap)
# For the USA map
  library(maps)
  library(mapproj)


## A. Data prep ##
# For base map: get GoogleEarth imagery and save as a ggplot object
  msq.zoom14.bw  <- get_googlemap(center=c(long=-86.252769, lat=45.958286), 
                                  zoom=14, maptype="satellite", color="bw")   #or 'color' instead of 'bw'
  msq.map.bw  <- ggmap(msq.zoom14.bw)
  bb <- attr(msq.zoom14.bw, "bb")  #extract bounding box [bb] coords; from http://stackoverflow.com/questions/18136468/is-there-a-way-to-add-a-scale-bar-for-linear-distances-to-ggmap

# Site coordinate data
  coords <- read.csv("Man_FW_SiteCoords.csv")
  coords <- rename(coords, long=alt_x_coord, lat=alt_y_coord)

# Shapefile for AOC boundary polygon
  aoc.bd <- readShapeSpatial(fn=paste(DirData, "aoc_mi_manistique\\AOC_MI_MANISTIQUE", sep=""))
  aoc.bd.df <- fortify(aoc.bd) # ?somehow converts object to usable df. Took this http://stackoverflow.com/questions/18084609/im-having-trouble-adding-a-shapefile-to-my-ggmap-due-to-differing-geographic-un


## B. Plotting ##
# Basic map with background, sample sites, and AOC bdry
  anno.map <- msq.map.bw +
    geom_polygon(data=aoc.bd.df[aoc.bd.df$group==0.1, ], aes(x=long, y=lat),   #group==0.1 gets it to use AOC bdry only, excluding MSQ watershed bdry.
                 color="green", size=0.5, fill=NA) + 
    geom_point(data=coords, aes(x=long, y=lat), shape=21, fill="orange", size=3) +
    geom_text (data=coords,
               aes(x=long-0.0004, y=lat+0.0006, label=site_number), 
               color="white", size=5) 

# Laboriously add scale bar: see http://stackoverflow.com/questions/18136468/is-there-a-way-to-add-a-scale-bar-for-linear-distances-to-ggmap
  # Simple function to estimate distance between two pts on the map (possibly inaccurate? But good enough for my needs.)
    distHaversine <- function(long, lat){
      
      long <- long*pi/180
      lat <- lat*pi/180  
      dlong = (long[2] - long[1])
      dlat  = (lat[2] - lat[1])
      
      # Haversine formula:
      R = 6371;
      a = sin(dlat/2)*sin(dlat/2) + cos(lat[1])*cos(lat[2])*sin(dlong/2)*sin(dlong/2)
      c = 2 * atan2( sqrt(a), sqrt(1-a) )
      d = R * c
      return(d) # in km
    }

  # Distance conversion factor
    ptspermm <- 2.83464567  # need this because geom_text uses mm, and themes use pts. Urgh.

  # Extract data from the bb, as needed for a scale bar (multiples edited to plot in the right location on my map)
    sbar <- data.frame(lon.start = c(bb$ll.lon + 0.4*(bb$ur.lon - bb$ll.lon)),
                       lon.end = c(bb$ll.lon + 0.65*(bb$ur.lon - bb$ll.lon)),
                       lat.start = c(bb$ll.lat + .15*(bb$ur.lat - bb$ll.lat)),
                       lat.end = c(bb$ll.lat + .15*(bb$ur.lat - bb$ll.lat)))

  # Use the function to calculate the distance between the start and end coords 
  # of the future scale bar (should be about 1km)
    sbar$distance = distHaversine(long = c(sbar$lon.start,sbar$lon.end),  #units=km
                                  lat = c(sbar$lat.start,sbar$lat.end))

  # Set desired length (distance) to be shown by scale bar
    scalebar.length <- 0.5 #in km
  
  # Calculate new end coord for scale bar based on above distance
  # (note I edited this slightly from the internet, making it a new column,  
  #  so that the scale bar length can be changed.)
    sbar$new.lon.end <- sbar$lon.start +
      ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length

  # Now add the scale bar to the map as manual annotation
    scale.map <- anno.map + 
      geom_segment(data = sbar,
                   aes(x = lon.start,
                       xend = new.lon.end,
                       y = lat.start,
                       yend = lat.end),
                   arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                               ends = "both", type = "open"), color="white") +
      geom_text(data = sbar,
                aes(x = (lon.start + new.lon.end)/2,
                    y = lat.start + 0.005*(bb$ur.lat - bb$ll.lat),
                    label = paste(format(scalebar.length), 'km')),
                hjust = 0.5,
                vjust = 0,
                size = 8/ptspermm, color="white")

## C. Print output ##
  ggsave(filename=paste(DirOut, "FigSD_ComboMap.tif", sep=""),
         plot=scale.map,
         width=13, height=13, units="in", 
         device="tiff")


## D. Area map of USA/Canada
# Pull data
  # USA state outlines
  states.dat <- map_data("state")
    # To get one state only...
    mi.only <- subset(states.dat, region %in% "michigan")  #to do only some states; for multiple states, c("michigan", "illinois", etc.
  # Canada outline
  canada.dat <- map_data(map="world", region=c("canada"))

# Create a map with a point showing Manistique
  noram.map <- ggplot(states.dat, aes(long, lat, group = group)) +
    geom_polygon(fill = "white", colour = "black") +  #outlines of all usa states (could be limited to a subset only)
    geom_polygon(data=canada.dat, fill="white", colour="black") +   #add canada outline
    geom_point(aes(x=-86.25, y=45.95), color="red") +   #point showing location of MSQ
    coord_map("gilbert") +                              #map projection
    theme(panel.background=element_rect(fill="lightsteelblue1"), #make the background blue so that the lakes stand out.
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())  

# Print output
  pdf(file=paste(DirOut, "FigSD_NAm_Map.pdf", sep=""), width=10, height=10)
    noram.map
  dev.off()

  svg(file=paste(DirOut, "FigSD_NAm_Map.svg", sep=""), width=10, height=10)
    noram.map
  dev.off()

#### END SCRIPT ####