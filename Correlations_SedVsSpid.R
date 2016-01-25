##################################################################
## MSQ Data paper, data analysis:                               ##
## Correlations between Sediment and Spider concentrations      ##
##################################################################

# Run setup first

# Packages
  library(gridExtra)   #for arranging plots
  library(scales)      #for making formatted plots (with 'calcnew' function)

#### RESULTS NUMBERS ####
# Analysis decisions:
# - Present log-log correlations.
# - Although it is problematic, for now we decided to include zeros
#   when calculating site means for normalized concentrations - this low-bias
#   seems preferable to the potentially strong high-bias of eliminating zero 
#   values.  However, we are still eliminating zero-sitemeans from the 
#   regressions - they become serious/problematic outliers when using a log-
#   log scale - but this decision requires justification - ?.  
# - Currently (2016-01-22) the text and figures use the toc-wet and toc-lip 
#   correlations, although I think it would be better to report the dry-lip 
#   toc-lip - at least, to make the argument later for bioavailability.

## Prepare the clean dataset
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
  # Clean up
  rm(means.sed); rm(means.spid)

# Sitemeans of zero are being excluded from the regressions
  means.nozero <- means.h
  means.nozero[means.nozero==0] <- NA  #note this doesn't delete the whole row, just replaces zeros with NAs - so that the other non-zero values can still be used.

  # Set of the different possible combinations of concentration types
    conc.combos <- list(c("sed.pcb", "spid.pcb"),
         c("pcb.tocnorm", "spid.pcb"),
         c("sed.pcb", "pcb.lipidnorm"), 
         c("pcb.tocnorm", "pcb.lipidnorm"))


# Function to calculate desired regression parameters for a given pair of columns
  fit.log <- function(df, sed.conc, spid.conc, add.const=0){
    m <- lm(log10(df[, spid.conc] +add.const) ~ log10(df[, sed.conc]+add.const))
    cf <- coef(m)
    tinfo <- summary(m)$coefficients[2, c(2, 4)]
    r2 <- round(summary(m)$r.squared, 3)
    dat.cur <- data.frame(intercept = cf[1], slope = cf[2], n = length(resid(m)),
               slope.se = tinfo[1], pval = tinfo[2], Rsq = r2)
    dat.cur <- mutate(dat.cur, sed.var=sed.conc, spid.var=spid.conc)
    return(dat.cur)
  }

# Now, loop through the 4 desired pairs (lapply/mapply would be better,  
# but the pairs of columns makes this too difficult for now).
  fitlist=list()
  for (i in 1:length(conc.combos)){
    cols.cur <- conc.combos[[i]]
    sed.conc <- cols.cur[1]
    spid.conc <- cols.cur[2]
    print(cols.cur)
    param.cur <- fit.log(df=means.nozero, sed.conc=sed.conc, spid.conc=spid.conc, add.const=0)
    fitlist[[i]] <- param.cur
  }
  # Convert output list to one consistent data frame
  (sed.spid.coefs <- ldply(fitlist))

  # Clean up
  rm(param.cur); rm(cols.cur); rm(fitlist); rm(i); rm(sed.conc); rm(spid.conc)


#### DATA EXPLORATION ####
# Add a factor to control shape aesthetic by location type
  means.h$location[means.h$site_number %in% c(1,2,6,7)] <- "River"
  means.h$location[means.h$site_number %in% c(3,4,5)] <- "Backwater"
  means.h$location[means.h$site_number %in% c(8,9,11)] <- "Harbor"
  means.h$location[means.h$site_number %in% c(10)] <- "Lake"


# So how many NAs/zeros are there?
  for (i in colnames(means.h)[4:7]){
    print(i)
    print(paste("n.na=", sum(is.na(means.h[, i])),
                " n.zero=", sum(means.h[, i]==0, na.rm=TRUE), sep=""))
  }

# Sample sizes for each comparison
  # Set of the different possible combinations
  conc.combos <- list(c("sed.pcb", "spid.pcb"),
       c("pcb.tocnorm", "spid.pcb"),
       c("sed.pcb", "pcb.lipidnorm"), 
       c("pcb.tocnorm", "pcb.lipidnorm"))

  # Print the sample size with and without zeros included, for each comparison
  # (dry-wet, toc-wet, dry-lip, toc-lip)
  for(i in 1:(length(conc.combos))) {
    cols.cur <- conc.combos[[i]]
    print(cols.cur)
    means.nona.cur <- means.h[complete.cases(means.h[, cols.cur]), ]
    n.zero.inc <- nrow(means.nona.cur)  #df with no NAs in specified columns.
    n.zero.exc <- sum((means.nona.cur[, cols.cur[1]]!=0 )&(means.nona.cur[, cols.cur[2]]!=0))
    print(paste("n with zero =", n.zero.inc, 
                "   n without zero =", n.zero.exc))
  }

  # Clean up
  rm(cols.cur); rm(i); rm(means.nona.cur); rm(n.zero.inc); rm(n.zero.exc)

# Datasets to use: either all (non-na) rows, or with all zeros excluded
  means.all <- means.h
  means.nozero <- means.h
  means.nozero[means.nozero==0] <- NA  #note this doesn't delete the whole row, just replaces zeros with NAs - so that the other non-zero values can still be used.

# Quick functions for extracting the R2 of the correlation between
# spider concentrations and sediment concentrations (by weight or 
# normalized) from a data frame in the format of /means.h/
  cor.log <- function(df, sed.conc, spid.conc, add.const=0){
    fit <- lm(log10(df[, spid.conc] +add.const) ~ log10(df[, sed.conc]+add.const))
    rsq <- summary(fit)$r.squared
    return(rsq)
  }

  cor.nolog <- function(df, sed.conc, spid.conc, add.const=0){
    fit <- lm((df[, spid.conc]) ~ (df[, sed.conc]))
    rsq <- summary(fit)$r.squared
    return(rsq)
  }


# Loop to calculate the r2 for each combination, given a df/add.const
  data.cur <- means.nozero
  add.cur <- 0   

  for (i in 1:(length(conc.combos))){
    cols.cur <- conc.combos[[i]]
    cor.cur  <- cor.nolog(df=data.cur, sed.conc=cols.cur[1], spid.conc=cols.cur[2], add.const=add.cur)
    cor.cur <- round(cor.cur, 3)
    print(paste(cols.cur[1], cols.cur[2], cor.cur))
  #   print(cor.cur)
  }
  
  # Clean up
  rm(data.cur); rm(add.cur); rm(cols.cur); rm(cor.cur); rm(i)


## Exploratory Plots ##
# Function to do a basic plot
  plotcor <- function(data.cur, add.cur, sed.conc, spid.conc, logscale=TRUE) {

    # Establish local environment, otherwise ggplot won't recognize variables defined within function
    localenv <- environment()

    baseplot <- ggplot(data=data.cur, environment=localenv, 
                       aes(x=(data.cur[, sed.conc] +add.cur), 
                           y=(data.cur[, spid.conc]+add.cur)) ) +
      geom_point(aes(shape=category, color=location), size=3) +
      xlab(sed.conc) +
      ylab(spid.conc) +
      geom_smooth(method=lm, se=FALSE, color="black")
    
    logplot <- 
      baseplot + 
      scale_x_log10() +
      scale_y_log10()
    
    if(logscale==TRUE){
      return(logplot)
    }
    if(logscale==FALSE){
      return(baseplot)
    }
  }

# Create plots
  # Set parameters
  data.cur <- means.all 
  add.cur  <- 0   #use 1 or 0.1 to show the plots at a 'real' scale; use 0 to 'squish' the zeros in with whatever the scale of the non-zero data points is.
  sed.conc <- "pcb.tocnorm" 
  spid.conc<- "pcb.lipidnorm" 
  logscale <- FALSE

  # Produce a single plot using given parameters
#   plotcor(data.cur, add.cur, sed.conc, spid.conc, logscale) 

  # Loop to do one plot for each conc-combo, using the same parameters
  plotlist=list()
  for (i in 1:length(conc.combos)){
    cols.cur <- conc.combos[[i]]
    sed.conc <- cols.cur[1]
    spid.conc <- cols.cur[2]
    print(cols.cur)
    plot.cur <- plotcor(data.cur, add.cur, sed.conc, spid.conc, logscale=logscale)
    plotlist[[i]] <- plot.cur
  }
    # Arrange the plots in a 2X2 grid
    grid.arrange(plotlist[[1]], plotlist[[2]], plotlist[[3]], plotlist[[4]], ncol=2)

  # Clean up
  rm(data.cur); rm(add.cur); rm(sed.conc); rm(spid.conc); rm(logscale); rm(plotlist); rm(i); rm(cols.cur); rm(plot.cur)

##### EXTRA PLOTTING FUNCTION ####

# Function to calculate the correlation, make a standard plot that indicates
# outlier points and includes annotation of function and regression parameters.  
# Inputs are a) the dataframe: typically means.all or means.nozero
#            b) sediment data type: {dry, toc}
#            c) spider data type:   {wet, lip}
#            d) if retaining zeros, have to add some value; default is 0.1
#            c) remove.zeros=TRUE; should zero-sitemeans be included in plots and regressions?
# 
# Output is a list with numerous elements:
#           1) fit - the basic lm output
#           2) summary - summary of the lm output
#           3) outliers - a data frame with one row for each 'outlier' data point (these points are not excluded however)
#           4) n.data - the number of data points used to calculate the regression (so it excludes rows with nas)
#           5) n.out - the number of outlier data points 
#           6) plot - the complete plot, with regression line and r2 annotation

test <- calcnew(df=means.nozero, sed="toc", spid="lip", add.const=0, remove.zeros=TRUE)
  test$plot


calcnew <- function(df, sed, spid, add.const=0.1, remove.zeros=FALSE){

  # Establish local environment, otherwise ggplot won't recognize variables defined within function
    localenv <- environment()
  
  # Establish some variables based on function inputs
    if (sed=="dry") {
      xcol <- 6  #column to use for x (ie, sediment) data
      xcap <- expression("Sediment " * Sigma * "PCBs (ng g"^-1*" dry weight)")   #label for x-axis in plot
      xbreaks <- 6
      xlimits <- c(0.8, 12000)
      xpos <- 1000
     } 
  
    if (sed=="toc") {
      xcol <- 7  #column to use for x (ie, sediment) data
      xcap <- expression("Sediment " * Sigma * "PCBs"[TOC]*" (ng g"^-1*" organic carbon)")   #label for x-axis in plot
      xbreaks <- 4
      xlimits <- c(100, 1000000)
      xpos <- 100000
    } 
  
    if (spid=="wet") {
      ycol <- 4  #column to use for y (ie, spider) data
      ycap <- expression("Spider " * Sigma * "PCBs (ng g"^-1*" wet weight)")   #label for x-axis in plot
      ybreaks <- 4
      ylimits <- c(1, 1000)
      ypos <- 5
    } 
      
    if (spid=="lip") {
      ycol <- 5  #column to use for y (ie, spider) data
      ycap <- expression("Spider " * Sigma * "PCBs" [lipid] * " (ng g"^-1*" lipids)")   #label for y-axis in plot
      ybreaks <- 3
      ylimits <- c(100, 12000)
      ypos <- 200
    } 
  
  # Create a data frame limited to only the complete data points used in the regression (no NAs)
    means.data <- df[complete.cases(df[, c(xcol, ycol)]), ]  
  
  # Remove data points where a concentration is zero - 
  #  necessary for normalized values, optional for dry/wet wt concentrations
    if (remove.zeros==TRUE) {
      means.data <- means.data[(means.data[, xcol]>0) & (means.data[, ycol]>0), ]  #if desired, remove site #1!  (not formally an outlier, but looks like one)
    }
  
  # Calcualte an initial lm regression using all data
    fit <- lm(log10(means.data[, ycol] + add.const) ~ log10(means.data[, xcol] + add.const), na.action=na.omit)
  
  # Calculate studentized residuals, for indicating which pts might function as outliers.
    fit.stres  <- rstandard(fit) # standardized residuals

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
    fit.out <- lm(log10(means.out[, ycol] + add.const) ~ log10(means.out[, xcol] + add.const), na.action=na.omit)
  
  # Print sample size
    print(paste("n=", n.data, "; n.outliers=", n.out, sep=""))
  
  # Extract components from the lm output
    r.squared.out <- round(summary(fit.out)$r.squared, 3)
    p.value.out   <- round(summary(fit.out)$coefficients[2,4], 5)
    # Create a custom annotation for the plot
    zero.note <- "zeros inc."
    if(remove.zeros==TRUE){
      zero.note <- "zeros exc."
    }
    fit.notes <- paste("r2 = ", r.squared.out, ", p = ", (format(p.value.out, scientific=FALSE)), "
  add ", add.const, ", ", zero.note, sep="")
  
  # Create the plot! (using 'category' as the classification variable, for now) 
    plot.cat <- ggplot(data=means.h, environment=localenv) + 
      geom_point(aes(x=means.h[, xcol], y=means.h[, ycol], fill=category, shape=location), 
                 color="black", size=3.5) +
#   coord_cartesian(xlim=xlimits, ylim=ylimits) +   #temporary solution, but it screws up the axes.
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=xbreaks),  #zoomed in: n=4
                    labels = comma) + 
#                     , limits=xlimits) +               #zoomed in: c(1000, ...)
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=ybreaks),  #zoomed in: c(100, ...)
                    labels = comma) + 
#                     , limits=ylimits) +                   #zoomed in: n=3
      scale_fill_manual(values=c("black", "white"), labels=c("Araneid", "Tetragnathid")) +
      scale_shape_manual(values=c(21, 22, 23, 24)) +
      xlab(xcap) +
      ylab(ycap) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line = element_line(colour="black"),
            legend.title= element_blank(), 
            legend.position = "right") + #c(0.85, 0.2)) #c(0.9, 0.2) or "none" +
    guides(fill=guide_legend(override.aes=list(shape=21)))
  
  # Add regression line and annotations re the fit to the base plot
    plot <- plot.cat +
      geom_smooth(data=means.out, 
                  aes(x=means.out[, xcol] + 0.1, y=means.out[, ycol] + 0.1), method=lm, se=FALSE) +
      annotate(geom="text", x=xpos, y=ypos, label=fit.notes)

  # Circle the 'outlier' plots
    plot <- plot + 
      geom_point(data=outliers, aes(x=(outliers[, xcol] + add.const), y=(outliers[,ycol]+add.const)), 
                 shape=1, size=6, color="green") 
#       geom_text(data=outliers,    #currently this is giving a weird error message when add.const<1; not sure why, can't fix it, so forget it for now.
#                      aes(x=(outliers[, xcol]+add.const), 
#                          y=(outliers[, ycol]+add.const), 
#                          label=site_number) , size=8) 
                                   

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

#### END SCRIPT ####