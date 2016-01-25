##################################################################
## MSQ Data paper, data analysis:                               ##
## Correlations between Araneid and Tetragnathid concentrations ##
##################################################################

# Run setup first

# For arranging multiple plots:
  library(gridExtra)

#### DATA EXPLORATION ####
# Quick functions for extracting the R2 of the correlation between
# Tet and Ara from a data frame in the format of /means.w/
  r2.nolog <- function(df){
    fit <- lm(Tetragnathid ~ Araneid, data=df)
    rsq <- summary(fit)$r.squared
    return(rsq)
  }
  
  r2.log <- function(df, add.const){
    fit <- lm(log10(df$Tetragnathid+add.const) ~ log10(df$Araneid+add.const))
    rsq <- summary(fit)$r.squared
    return(rsq)
  }

# Different possible datasets for use with the r2. functions
  means.all <- means.w  #all data points, inc zeros
  means.aoc <- means.w[means.w$site!=1 & means.w$site!=10, ]  #exclude sites 1 and 10...no good reason to do this.
  means.nozero <- means.w[means.w$Tetragnathid!=0 & means.w$Araneid!=0, ] #exclude all zero values (only two data points, in fact)

# Correlation among all spider data points, data points from all years combined
  r2.nolog(means.all)
  r2.nolog(means.nozero)
  r2.nolog(means.aoc)
  
  r2.log(means.all, add.const=0.01)
  r2.log(means.nozero, add.const=0)
  r2.log(means.aoc, add.const=0)

# Correlation coefficients among spider data points within each year
  # no log
  yr.list <- split(x=means.aoc, f=means.aoc$sample_year)
  sapply(X=yr.list, FUN=r2.nolog)

  # log
  yr.list <- split(x=means.nozero, f=means.nozero$sample_year)
  sapply(X=yr.list, FUN=r2.log, add.const=0)

# Are the slopes significantly different between years?
  # Utility function for use in plyr's ldply() function.
  extractfun <- function(m) {
    cf <- coef(m)
    tinfo <- summary(m)$coefficients[2, c(2, 4)]
    r2 <- round(summary(m)$r.squared, 3)
    data.frame(intercept = cf[1], slope = cf[2], n = length(resid(m)),
               slope.se = tinfo[1], pval = tinfo[2], Rsq = r2)
  }

  # LOG-LOG SCALE
    logfits <- dlply(means.nozero, 
                     "sample_year", 
                     function(d) lm(log10(Tetragnathid) ~ log10(Araneid), data=d) ) #no need to add 0.1 because zeros have alreayd been removed.
    # Compare the 95% CIs for the slopes
    ldply(logfits, extractfun)
    sapply(X=logfits, FUN=confint, parm="log10(Araneid)", level=0.95)  #suggests poss sig difference between 2011, 2013?
    # Fit a full model with interaction; is the interaction term significant?
    anova(lm(log10(Tetragnathid)~log10(Araneid)*sample_year, data=means.nozero))
    
  # UNTRANSFORMED SCALE
    basefits <- dlply(means.nozero, 
                     "sample_year", 
                     function(d) lm(Tetragnathid ~ Araneid, data=d) ) #no need to add 0.1 because zeros have alreayd been removed.
    # Compare the 95% CIs for the slopes
    ldply(basefits, extractfun)
    sapply(X=basefits, FUN=confint, parm="Araneid", level=0.95)  #suggests poss sig difference between 2011, 2013?
    # Fit a full model with interaction; is the interaction term significant?
    anova(lm(Tetragnathid~Araneid*sample_year, data=means.nozero))
 

  
## Exploratory Plots
baseplot <- ggplot(data=means.all, 
                   aes(x=Araneid, y=Tetragnathid)) +
  geom_point(aes(color=factor(sample_year)), size=3) 

  # just the data points
    baseplot
  
  # Add regression lines for each year
    baseplot <- 
      baseplot + 
        geom_smooth(aes(color=factor(sample_year)), method=lm, se=FALSE)

  # Add the overall regression line
    baseplot <- 
      baseplot + 
        geom_smooth(method=lm, se=FALSE, color="black") 

  # Add a 1:1 line
    baseplot <- 
      baseplot + 
        geom_segment(aes(x=0, xend=400, y=0, yend=400), lty=3)  #1:1 line
  
  # Log-scale the axes
#     baseplot <- 
      baseplot + 
        scale_x_log10(limits=c(10, 1000)) +  #limits=c(10, 1000)             
        scale_y_log10(limits=c( 5, 1000))    #limits=c(5, 1000)

# Logging of axes vs logging of pts - do they produce the same figure/
# regression lines?  
  # Logging of axes (reg line is on un-transformed data points)
  logaxes <- ggplot(data=means.nozero, aes(x=(Araneid), y=(Tetragnathid)) ) +
    geom_point(aes(shape=factor(sample_year)), size=3) + 
    scale_shape_manual(values=c(0, 19, 8)) +
    geom_smooth(aes(linetype=factor(sample_year)), method=lm, se=FALSE, color="black") +
    theme(legend.position="none") +
        scale_x_log10() +  #limits=c(10, 1000)             
        scale_y_log10()    #limits=c(5, 1000)

  # Logging of points (regression is on logged points)
  logpts <- ggplot(data=means.nozero, aes(x=log10(Araneid), y=log10(Tetragnathid)) ) +
    geom_point(aes(shape=factor(sample_year)), size=3) + 
    scale_shape_manual(values=c(0, 19, 8)) +
    geom_smooth(aes(linetype=factor(sample_year)), method=lm, se=FALSE, color="black") +
    theme(legend.position="none")

  # Now plot them together to compare/confirm they are same shapes
  grid.arrange(logaxes, logpts, ncol=2)


#### RESULTS NUMBERS ####
# Analysis decisions:
# - Present log-log correlation (log-log does a better job of 
#   handling the large values, accounts for the increase in 
#   spread/variance when looking at all data points together).
# - But, also prepare a figure showing the unscaled plot and 
#   regression coefficients - for the appendix.
# - Exclude zero-sitemeans from the regression (as is done 
#   for the sed-spid figure), and for the same justification.

# Data set excluding zero-sitemeans
  means.nozero <- filter(means.w, !is.na(Tetragnathid) & !is.na(Araneid)) #first, removes the 4 data points that have at least one NA obs.
  means.nozero <- filter(means.nozero, Tetragnathid>0 & Araneid>0)        #then, remove the additional two points where one or both spiders is zero.
  # The no-zero policy excludes 3 rows/4 sitemeans: 
  #     2012 site 1 is zero for both Ara and Tet; 
  #     2012 site 10 is zero for Ara, 38 for Tet
  #     2013 site 10 is NA for Ara, zero for Tet - so this point is exluced in any case.

# Correlations for individual years (2011, 2012, 2013)
# (Code taken from https://stat.ethz.ch/pipermail/r-help/2011-September/291070.html)
  # Calculate the lm-fit object for each year (return a list of these fit-objects).
  year.l <- dlply(means.nozero, 
                  "sample_year", 
                  function(d) lm(log10(Tetragnathid) ~ log10(Araneid), data=d) ) #no need to add 0.1 because zeros have alreayd been removed.
 
  # Utility function for use in plyr's ldply() function.
  extractfun <- function(m) {
    cf <- coef(m)
    tinfo <- summary(m)$coefficients[2, c(2, 4)]
    r2 <- round(summary(m)$r.squared, 3)
    data.frame(intercept = cf[1], slope = cf[2], n = length(resid(m)),
               slope.se = tinfo[1], pval = tinfo[2], Rsq = r2)
  }

  # Data frame giving the regression coeffecients for each year
  (coefs <- ldply(year.l, extractfun))

# Correlation of all spider points combined
  # Regression coefficients
  fit.all <- lm(log10(Tetragnathid) ~ log10(Araneid), data=means.nozero) 
  extractfun(fit.all)

  # 95% CI for the slope of the all-points regression
  confint(object=fit.all, parm="log10(Araneid)", level=0.95)

## To generate figure(s), run <Fig_SpidVsSpidCorrelations.R > 
## Use the results generated here to annotate the figures.  
##############################################################