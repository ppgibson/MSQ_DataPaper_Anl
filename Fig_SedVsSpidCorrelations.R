#########################################################################
# Fig "5", relationship between sediment and spider pcb concentrations
#########################################################################
# 
# Plot lipid-normalized spider PCB concentrations as a function of  
# TOC-normalized sediment concentration. 
# Three parts to the main script:
# A. Create a hybrid-format data frame with sed (toc-norm) vs both spiders (lipid-norm).
# B. Complicated function to calculate lm/remove outliers/plot figure/return list
# C. Code to run the function (including exclude all of sites 1 and 10).
# 
# Main script calculates regression for all site/year data points combined.
# Code to calculate the same thing with grand means (years combined) is still mostly here -
# part 2 / Fig5b.
# 
# Input: \means\ df from setup  
# 
# Output: function returns a list with ggplot object + a variety of data re outliers, etc.
#   byyr.aoc -  list object, after running function with sites 1 and 10 excluded
#   Fig5 - ggplot object



## A. ##
## First, have to create a hybrid wide/long format data frame, 
## with sediment data and in column and matching spider data in other columns

# Split 'means' into two data frames, one for sediment and one for spiders
  means.sed  <- subset(means, category=="Sediment", 
                       select=c("sample_year", "site_number", "Total.PCB"))
  means.sed <- rename(means.sed, sed.pcb=Total.PCB)  

  means.spid <- subset(means, category!="Sediment", 
                       select=c("sample_year", "category", "site_number", "Total.PCB"))
  means.spid <- rename(means.spid, spid.pcb=Total.PCB)

# Merge the spider and sediment data frames (each sediment data point will appear twice, once for Tet and once for Ara)
  means.h <- merge(means.spid, means.sed, by=c("sample_year", "site_number"), all=TRUE)   #"h" for hybrid, vs "l" for long

# Merge normalized data with the hybrid 'means.h' dataframe
  # toc-normalized
    means.h <- merge(means.h, means[means$category=="Sediment", c(1, 3, 13)],  #id columns plus [pcb.tocnorm] from \means\
                     by=c("sample_year", "site_number"), all=TRUE)  
  # lipid-normalized
    means.h <- merge(means.h, means[means$category!="Sediment", c(1:3, 12)], 
                     by=c("sample_year", "site_number", "category"), all=TRUE)

# Add a new column to get a single aes category for ggplot mapping (Ara vs Tet vs OutOfAOC)
  means.h$plotshape[means.h$category=="Araneid"] <- "Araneid"
  means.h$plotshape[means.h$category=="Tetragnathid"] <- "Tetragnathid"  #fill in spider values everywhere FIRST,
  means.h$plotshape[(means.h$site_number==1 | means.h$site_number==10)] <- "outofaoc"
  # Reorder factor levels so that they will plot in correct order
  means.h$plotshape <- factor(means.h$plotshape, levels=c("Araneid", "Tetragnathid", "outofaoc"))


# # Convert spider 0 values to 1 for plotting purposes
#   means.h$pcb.lipidnorm[means.h$pcb.lipidnorm==0] <- 1  #currently 4 data pts
#   means.h$spid.pcb[means.h$spid.pcb==0] <- 1            #currently 4 data pts
#   # There are no sed.pcb = 0 (at least for now).

# Clean up
  rm(means.sed)
  rm(means.spid)



# ## Calculate grand means for use in Part 2  #currently not being used! 
# 
# # Spider grand means
#   # With spider taxa separate
#     spid.sep.gm <- ddply(.data=means.h, .variables=c("site_number", "category"), 
#                          .fun=summarise, gms.lipnorm = mean(pcb.lipidnorm, na.rm=TRUE))
#   # With spider taxa combined
#     spid.comb.gm <- ddply(.data=means.h, .variables=c("site_number"), 
#                           .fun=summarise, gmc.lipnorm = mean(pcb.lipidnorm, na.rm=TRUE))
# 
# # Sediment grand means
#   sed.gm <- ddply(.data=means.h, .variables=c("site_number"), 
#                   .fun=summarise, gm.toc = mean(pcb.tocnorm, na.rm=TRUE)) 
# 
# # Merge sediment + spiders into a single gm data frame
#   all.comb.gm <- merge(sed.gm, spid.comb.gm, all=TRUE)
#   all.sep.gm  <- merge(sed.gm, spid.sep.gm, by="site_number", all=TRUE)
# #   # Fix columns/column names for all.sep.gm
# #     all.sep.gm <- all.sep.gm[, -2]           #why would I want to do this?
# #     colnames(all.sep.gm)[3] <- "category"    #not needed
##########################################################################


## B. ##
# Function to calculate regression, remove outliers in preparation for plotting
# Inputs are a) the dataframe: typically means.h, all.sep.gm, or all.comb.gm
#            b) column indexes for the x and y columns in the dataframe: 
#                 means.h = 6, 7
#                 all.sep = 2, 4
#                 all.com = 3, 4
#            c) remove.1.10 = TRUE : should the data for sites 1 and 10 be removed?
# 
# Output is a list with numerous elements:
#           1) fit - the basic lm output
#           2) summary - summary of the lm output
#           3) outliers - a data frame with one row for each outlier that was removed
#           4) n.data - the number of data points used to calculate the regression (so it excludes rows with nas; outliers; and reference site points)
#           5) n.out - the number of outlier data points removed
#           6) plot - the complete plot, with regression line and r2 annotation

#sed = "dry" or "toc"; spid="wet" or "lip"
calcnew <- function(df, sed, spid, remove.1.10=TRUE){

  # Establish local environment, otherwise ggplot won't recognize variables defined within function
    localenv <- environment()
  
  # Establish some variables based on function inputs
    if (sed=="dry") {
      xcol <- 5  #column to use for x (ie, sediment) data
      xcap <- expression("Sediment " * Sigma * "PCBs (ng g"^-1*" dry weight)")   #label for x-axis in plot
      xbreaks <- 6
      xlimits <- c(0.8, 12000)
    } 
  
    if (sed=="toc") {
      xcol <- 6  #column to use for x (ie, sediment) data
      xcap <- expression("Sediment " * Sigma * "PCBs"[TOC]*" (ng g"^-1*" organic carbon)")   #label for x-axis in plot
      xbreaks <- 4
      xlimits <- c(100, 1000000)
    } 
  
    if (spid=="wet") {
      ycol <- 4  #column to use for y (ie, spider) data
      ycap <- expression("Spider " * Sigma * "PCBs (ng g"^-1*" wet weight)")   #label for x-axis in plot
      ybreaks <- 4
      ylimits <- c(1, 1000)
    } 
      
    if (spid=="lip") {
      ycol <- 7  #column to use for y (ie, spider) data
      ycap <- expression("Spider " * Sigma * "PCBs" [lipid] * " (ng g"^-1*" lipids)")   #label for y-axis in plot
      ybreaks <- 3
      ylimits <- c(100, 12000)
    } 
  
  # Create a data frame limited to only the complete data points used in the regression
    means.data <- df[complete.cases(df[, c(xcol, ycol)]), ]  
  
  # Remove reference site and site 10 (ie, outside AOC) if desired (specified in function)
    if (remove.1.10==TRUE) {
      means.data <- means.data[(means.data$site_number!=1 & means.data$site_number!=10), ]  #if desired, remove site #1!  (not formally an outlier, but looks like one)
    }
  
  # Calcualte an initial lm regression using all data
    fit <- lm(log10(means.data[, ycol] + 1) ~ log10(means.data[, xcol] + 1), na.action=na.omit)
  
  # Calculate residuals
    fit.res    <- resid(fit)     # plain residuals
    fit.stres  <- rstandard(fit) # standardized residuals
    fit.stures <- rstudent(fit)  # studentized residuals
  

#   # Remove outliers (stresid > 2) from the data frame
#     means.out <- means.data[-(as.numeric(names(which(abs(fit.stres) > 2)))), ]  #subset means.h df to remove outlier rows (rows where stand resid > 2)

    means.out <- means.data  #instead of eliminating outliers, I will just rename complete df as means.out

  # Sample sizes and outliers
    # Data frame of outliers only, in case data on specific outliers is desired 
      outliers <- means.data[(as.numeric(names(which(abs(fit.stres) > 2)))), ]  #subset means.h df to remove outlier rows (rows where stand resid > 2)
    # Sample sizes for number of data points used in the regression and number excluded as outliers
      n.data <- nrow(means.out[complete.cases(means.out[, c(xcol, ycol)]), ])
      n.data.simp <- nrow(means.out)
      n.out  <- nrow(outliers)

  
  # Run the lm again
    fit.out <- lm(log10(means.out[, ycol] + 1) ~ log10(means.out[, xcol] + 1), na.action=na.omit)
    summary(fit.out)  
  
  # Print sample size
    print(paste("n=", n.data, "; n outliers = ", n.out))
  
  # Extract components from the lm output
    r.squared.out <- round(summary(fit.out)$r.squared, 3)
    p.value.out   <- round(summary(fit.out)$coefficients[2,4], 5)
  
  
  # Create the plot! (using 'category' as the classification variable, for now) 
    plot.cat <- ggplot(data=means.h, environment=localenv) + 
      geom_point(aes(x=means.h[, xcol], y=means.h[, ycol], fill=plotshape, shape=plotshape), 
                 color="black", size=3.5) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=xbreaks),  #zoomed in: n=4
                    labels = comma, limits=xlimits) +               #zoomed in: c(1000, ...)
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=ybreaks),  #zoomed in: c(100, ...)
                    labels = comma, limits=ylimits) +                   #zoomed in: n=3
      scale_fill_manual(values=c("black", "white", "grey"), labels=c("Araneid", "Tetragnathid", "samples from\noutside the AOC")) +
      scale_shape_manual(values=c(21, 21, 24), labels=c("Araneid", "Tetragnathid", "samples from\noutside the AOC")) +
      xlab(xcap) +
      ylab(ycap) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line = element_line(colour="black"),
            legend.title= element_blank(), 
            legend.position = "none") #c(0.85, 0.2)) #c(0.9, 0.2) or "none"
  
  # Add regression line to the base plot
    plot <- plot.cat +
      geom_smooth(data=means.out, 
                  aes(means.out[, xcol], y=means.out[, ycol]), method=lm, se=FALSE)
#       annotate(geom="text", x=1000, y=5, 
#                label=paste("r2 = ", r.squared.out, ", p = ", (format(p.value.out, scientific=FALSE)), sep=""))
    
  # Create a list with elements to return
  l <- list(fit = fit.out,
            summary = summary(fit.out),
            outliers = outliers,
            n.data = n.data,
            n.out = n.out,
            means.out = means.out,
            plot = plot,
            oldfit = fit)
  
  return(l)
}

############################################################

## C. ##

## Run the functions/create the plots

# Fig 5a: Sediment dry weight vs spider wet weight (no normalization)
  drywet.l <- calcnew(df=means.h, sed="dry", spid="wet", remove.1.10=TRUE)
  drywet.p <- drywet.l$plot

# Fig 5b: Sediment toc-normed vs spider wet weight 
  tocwet.l <- calcnew(df=means.h, sed="toc", spid="wet", remove.1.10=TRUE)
  tocwet.p <- tocwet.l$plot + theme(plot.margin=unit(c(0.5,0.5,0.5,0.6), "cm"))

# Fig 5c: Sediment dry weight vs spider lipid-normed
  drylip.l <- calcnew(df=means.h, sed="dry", spid="lip", remove.1.10=TRUE)
  drylip.p <- drylip.l$plot
  # note that two spid=0 values from site 10 are hidden

# Fig 5d: Sediment toc-normed vs spider lipid-normed
  toclip.l <- calcnew(df=means.h, sed="toc", spid="lip", remove.1.10=TRUE)
  toclip.p <- toclip.l$plot
  # note that two spid=0 values from site 10 are (again) hidden


## Export an .emf of both figures, lined up
#   emf(file=paste(DirOut, "Fig5_wetdry_toclip.emf", sep=""), width=8, height=10)
#   emf(file=paste(DirOut, "Fig5_tocwet_drylip.emf", sep=""), width=8, height=10)
  emf(file=paste(DirOut, "Fig5_tocwet_toclip.emf", sep=""), width=8, height=10)
 
  grid.arrange(tocwet.p, toclip.p, ncol=1)
 
  dev.off()

# Output is a list with numerous elements:
#           1) fit - the basic lm output
#           2) summary - summary of the lm output
#           3) outliers - a data frame with one row for each outlier that was removed
#           4) n.data - the number of data points used to calculate the regression (so it excludes rows with nas; outliers; and reference site points)
#           5) n.out - the number of outlier data points removed
#           6) plot - the complete plot, with regression line and r2 annotation

    
# Clean up
  rm(means.aoc)
  rm(means.h)
  rm(Fig5)
  rm(byyr.aoc)
  rm(calclm)


# END MAIN SCRIPT
########################