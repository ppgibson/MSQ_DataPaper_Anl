#########################################################################
# Fig "4", Regression of Araneid vs Tetragnathid pcb concentration
#########################################################################
# 
# Generates a plot of the correlation between Tet vs Ara total pcb concentration 
# for each year 2011-2013 and for all years combined; calculates lm regressions
# for each year and annotates the plot with these data.
# - **plot shows raw data (with scaled axes), but coefficients in the current 
#     script are for [log10(Tet) ~ log10(Ara)]; to change this, edit lm calculations
#     in 'Calculations' section, then adjust annotation locations in the ggplot script.
# - Calculations use full measured values, no normalization.
# - Note that the r2 annotation should have a superscript.
# - For final versions, remove plot annotations and manually add later.
#
# Note: a more elegant way to do regression/plotting by year is described at 
# <http://stackoverflow.com/questions/24983690/nls-and-log-scale-in-ggplot2>
#
# Inputs: \means.w\ df from setup
#
# Output: 
#  Version A: 'Fig4' (ggplot object)
#     plots the untransformed data (with log-scaled axes) and regression coefficients
#  (Version B plots log Tet ~ log Ara; this part of the script hasn't been updated)



###################
## Version A: Un-transformed data

## Calculations
  # Adjust means.w for use in the functions: remove sites 1 and 10
    # Note that it's no longer necessary to adjust for zeros - no zeroes without sites 1 and 10.
    means.temp <- means.w[means.w$site!=1 & means.w$site!=10, ]

  # Make a dummy data table for plotting 1:1 line (geom_segment doesn't work properly with scale_log10())
    onetooneline <- data.frame(cbind(c(10, 200, 600), c(10, 200, 600)))
  
  # Calculate the regression for all years, to extract r2 for annotation
    allyears.lm <- lm(log10(Tetragnathid) ~ log10(Araneid), data=means.temp)
    rsq.all <- summary(allyears.lm)$r.squared 
  
  # Calculate mid- and end-values for each year of data, in order to add annotation later
    max <- unlist(dlply(means.temp, .(sample_year), 
                        function(dat) max(dat$Araneid, na.rm=TRUE) ))  #get the maximum 'x' (ie, Araneid) value for each year
    mid <- unlist(dlply(means.temp, .(sample_year), 
                        function(dat) ( (max(dat$Araneid, na.rm=TRUE) + min(dat$Araneid, na.rm=TRUE))/2) ) )  #get the middle(mean) 'x' (ie, Araneid) value for each year
  
  # Calculate lm-fits for each year
  # code taken from https://stat.ethz.ch/pipermail/r-help/2011-September/291070.html
    # Input a data frame, output a list of lm objects
      year.l <- dlply(means.temp, .(sample_year), function(d) lm(log10(Tetragnathid) ~ log10(Araneid), data=d) )
    
    # Utility function for use in plyr's ldply() function.
      extractfun <- function(m) {
        cf <- coef(m)
        tinfo <- summary(m)$coefficients[2, c(2, 4)]
        r2 <- summary(m)$r.squared
        data.frame(intercept = cf[1], slope = cf[2], n = length(resid(m)),
                   slope.se = tinfo[1], pval = tinfo[2], Rsq = r2)
      }
    
  # Take a list (of models) as input and output a data frame:
    coefs <- ldply(year.l, extractfun)
    coefs <- cbind(coefs, max, mid)   #add the previously calculated max and mid to the coefs df
    
## Plotting, in a series of steps
  # The basic plot, Tet vs Ara different colors for each year
    Fig4 <- ggplot(data=(means.temp), aes(x=Araneid, y=Tetragnathid)) +   #remove the color=factor(sample_year) to get all years to plot as one color
      geom_point(shape=19, size=3, aes(color=factor(sample_year))) + 
      geom_smooth(aes(color=factor(sample_year)), method=lm, se=FALSE) +
#       geom_text(aes(label=site_number, vjust=-1)) +                     #label each point with site number 
      xlab(expression("Araneid "*Sigma * "PCBs (ng g"^"-1"*" wet weight)")) +
      ylab(expression("Tetragnathid "*Sigma * "PCBs (ng g"^"-1"*" wet weight)")) +
#       scale_x_continuous(limits=c(0, 500), breaks=c(0, 200, 400)) +
#       scale_y_continuous(limits=c(0, 600), breaks=c(0, 200, 400, 600)) +
      theme_bw() + 
      theme(panel.grid.major=element_blank(),          
            panel.grid.minor=element_blank(), 
            axis.line=element_line(colour="black"), 
            panel.border = element_blank(),
            legend.position="none")  

  # Add the 1: 1 line and anotation
    Fig4 <- Fig4 + 
            geom_smooth(data=onetooneline, aes(x=X1, y=X2), 
                        method=lm, se=FALSE, color="black", linetype="dashed") + 
            annotate(geom="text", x=650, y=650, label="1:1 line")
    
  # Add the regression line for all years combined, plus annotation
    Fig4 <- Fig4 + 
            geom_smooth(method=lm, se=FALSE, color="black", size=1) +
                         annotate(geom="text", x=475, y=350, label="all years") 
#                          annotate(geom="text", x=450, y=300, label=paste("r2 = ", round(rsq.all, 2)))
    
#   # Make x and y axes equal        #applying this does strange things to the plot, better to do without.
#     Fig4 <- Fig4 + coord_fixed()   
    
  # Add plot annotations
    # 1. Direct-label each regression line with the year
        Fig4 <- Fig4 + geom_text(data=coefs, 
#                        aes(x=max+50, y=(slope*(max+20) + intercept),   #values for non-log regression
                       aes(x=(max+c(0, -90, 0)), y=c(200, 80, 600), 
                           color=factor(sample_year), label=sample_year) )  
                                         
#     # 2. Direct-label reg lines with r2 values
#         Fig4 <- Fig4 + geom_text(data=coefs, 
# #                          aes(x=c(400, 180, 130), y=(slope*(mid) + intercept + c(50, -10, 25)),  #values for non-log regression (or, do x=mid)
#                          aes(x=c(400, 100, 130), y=c(170, 40, 300),  
#                                  label=paste("r2 = ", round(Rsq, 2)), color=factor(sample_year) ) )

  # Add log-log scaling if desired
    Fig4 <-
      Fig4 + 
      scale_x_log10(limits=c(10, 1000)) +  
      scale_y_log10(limits=c(5, 1000))    #be careful about setting limits any higher than this - lowest y val is ~7.5
      
  # Plot it!
    Fig4


# Export -> save as image -> metafile  800x800

# Print the file as a high quality .emf
  emf(file=paste(DirOut, "Fig4_plain.emf", sep=""), width=8.5, height=8.5)

  Fig4

  dev.off()



 
###############################################################################
# SCRATCH PAPER
## Version B: log-transformed data 

  means.w -> means.w.orig
  means.w[(means.w$site_number!=1 & means.w$site_number!=10), ] -> means.w
  
# The basic plot, Tet vs Ara different colors for each year
  Fig4.log10 <- ggplot(data=means.w, aes(x=log10(Araneid + 1), y=log10(Tetragnathid + 1))) +   #remove the color=factor(sample_year) to get all years to plot as one color
    geom_point(shape=19, size=3, aes(color=factor(sample_year))) + 
    geom_smooth(aes(color=factor(sample_year)), method=lm, se=FALSE) +
    xlab(expression("log Araneid "*Sigma * "PCBs (ng g"^"-1"*"wet weight)")) +
    ylab(expression("log Tetragnathid "*Sigma * "PCBs (ng g"^"-1"*"wet weight)")) +
    coord_fixed() + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(),          #remove all the background junk
          panel.grid.minor=element_blank(), 
          axis.line=element_line(colour="black"), 
          panel.border = element_blank(),
          legend.position="none") +
    ggtitle("log-transformed data")

# Add the 1: 1 line and anotation
  Fig4.log10 <- Fig4.log10 + 
    geom_segment(x=0, xend=3.0, y=0, yend=3.0, color="black", linetype="dashed") +  #log-e values: x=6.5, y=6.5
    annotate(geom="text", x=3.1, y=3.1, label="1:1 line")  #log-e values: 6.5, 6.5


# Add the regression line for all years combined, plus annotation
  # Calculate the regression for all years, to extract r2 for annotation
    all.log.lm  <- lm(log10(Tetragnathid + 1) ~ log10(Araneid + 1), data=means.w)
    all.log.rsq <- summary(all.log.lm)$r.squared    
  # Now add line and annotation to the plot
    Fig4.log10 <- Fig4.log10 + 
                geom_smooth(method=lm, se=FALSE, color="black", size=1) + 
                annotate(geom="text", x=2.75, y=2.3, label="all years") +  #log-e values: x=6.3, y=5.3
                annotate(geom="text", x=1.4, y=1.5, label=paste("r2 = ", round(all.log.rsq, 2)))  #log-e values: x=1.1, y=2.5

# Calculate mid- and end-values for each year of data, in order to add annotation later
  max.log <- unlist(dlply(means.w, .(sample_year), 
                      function(dat) max(log10(dat$Araneid + 1), na.rm=TRUE) ))  #get the maximum 'x' (ie, Araneid) value for each year
  mid.log <- unlist(dlply(means.w, .(sample_year), 
                      function(dat) ( (max(log10(dat$Araneid+1), na.rm=TRUE) + min(log10(dat$Araneid + 1), na.rm=TRUE))/2) ) )  #get the middle(mean) 'x' (ie, Araneid) value for each year


## Calculate lm-fits for each year
# code taken from https://stat.ethz.ch/pipermail/r-help/2011-September/291070.html
# Input a data frame, output a list of lm objects
  year.l.log <- dlply(means.w, .(sample_year), function(d) lm(log10(Tetragnathid + 1) ~ log10(Araneid + 1), data=d) )

  # Utility function for use in plyr's ldply() function. #Same as version in part a
    extractfun <- function(m) {
      cf <- coef(m)
      tinfo <- summary(m)$coefficients[2, c(2, 4)]
      r2 <- summary(m)$r.squared
      data.frame(intercept = cf[1], slope = cf[2], n = length(resid(m)),
                 slope.se = tinfo[1], pval = tinfo[2], Rsq = r2)
    }

  # Take a list (of models) as input, and output a data frame of relevant model coefficients:
    coefs.log <- ldply(year.l.log, extractfun)
    coefs.log <- cbind(coefs.log, max.log, mid.log)   #add the previously calculated max and mid to the coefs df

# Add plot annotations
  # 1. Direct-label each regression line with the year
    Fig4.log10 <- Fig4.log10 + geom_text(data=coefs.log, 
                             aes(x=max.log+0.12, y=(slope*(max.log+0.12) + intercept), # log-e values: 0.2? 2?
                                 color=factor(sample_year), label=sample_year) )
  
  # 2. Direct-label reg lines with r2 values
  # sadly, doesn't really work - r2 should be superscript
    Fig4.log10 <- Fig4.log10 + geom_text(data=coefs.log, 
                             aes(x=(mid.log + c(0.3, 0.3, -0.3)), y=(slope*(mid.log) + intercept),  #log-e values: (-0.8, 0.5, -0.7)
                                 label=paste("r2 = ", round(Rsq, 2)), color=factor(sample_year) ) )
# Plot it!
  Fig4.log10 + ggtitle("log10") + coord_fixed() + scale_y_continuous(limits=c(1, 3))

  Fig4.log10 + geom_text(aes(label=site_number)) + scale_x_continuous(limits=c(1, .5))

#########################################
## Scratch paper

# Add log-scaled axes to plot of untransformed data
Fig4 <- Fig4 + scale_x_log10(breaks=c(5, 10, 50, 100, 500), limits=c(4, 600), expand=c(0,0)) +
   scale_y_log10(breaks=c(5, 10, 50, 100, 500), limits=c(4, 700), expand=c(0,0))


# Another way to calculate model fits, and add them to the plot:
# See <http://stackoverflow.com/questions/24983690/nls-and-log-scale-in-ggplot2>
lm(log(Tetragnathid+1) ~ log(Araneid + 1), data=means.w) -> logfit
range(means.w[, c(3,5)], na.rm=TRUE) -> xr
seq(xr[1], xr[2], length.out=100)->xx
data.frame(Araneid=xx)->temp
predict(logfit, newdata=temp)->yy
data.frame(Araneid=xx,Tetragnathid=yy)->all.log.fit

Fig4 + geom_line(data=all.log.fit, color="black")

#########################################
## Scratch paper

rsquaredlab <- paste("r^2", round(coefs$Rsq, 2))
  expression(r^{2}==(round(coefs$Rsq, 2)))
label = "r^2 == 0.585~~~p < 0.001"), parse = TRUE) 
paste("r^2", round(coefs$Rsq, 2))
expression("r"^"2"*" = "*round(coefs$Rsq, 2)) -> rsquaredlab


lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

charformat = function(data) {
  l <- list(a = format(data$slope, digits = 2),
           b = format(data$intercept, digits = 2), 
           r2 = format(data$Rsq, digits = 3));
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  
  as.character(as.expression(eq));
  
}

format(summary(year.l$r.squared, digits = 3))

Fig4 + geom_text(aes(x = 25, y = 300, label = lm_eqn(lm(y ~ x, df))), parse = TRUE)



paste(italic('r'), " = 0.01" )))
ylab(expression(Sigma * "PCBs (ng g"^"-1"*"wet weight)")) 
expression(italic(r)^"2"*" = ") -> expr
           round(Rsq, 2)) -> test
expression(A,italic(A),bolditalic(A),Delta*italic(D)))
expression(r^{2}~ = ~round(Rsq, 2))

                       paste(~~italic(r)^2~"="~r2"R2 =", round(Rsq, 2)) )
Fig4 + annotate(geom="text", aes(x=mid, y=(slope*(mid) + intercept + c(20, -20, 25)), 
                                label=test)) 
                 # (another option)
                 Fig4 + annotate(geom="text", aes(x=max+15, y=(coefs$slope*(max+15) + coefs$intercept), color=coefs$sample_year), 
                    label=coefs$sample_year ), 
                    color=c("red", "green", "blue"))

######################
# scratch paper

    mid <- dlply(means.w, .(sample_year), fun=function(x) (max(x) + min(x))/2, na.rm=TRUE)
Fig4 + stat_function(fun=function(x) 3*(x) + 2, color="black")

  # Annotate each regression line with series name: - currently failed, will have to add them by hand.
  # Red=2011 (middle), green=2012(bottom), blue=2013 (top).
  # geom_dl() is adding annotation ot the last point, rather than the end of the geom_smooth line;
  # being able to pass "method = lm" to stat=smooth might fix the problem.
  # See <http://stackoverflow.com/questions/10065196/how-to-show-directlabels-after-geom-smooth-and-not-after-geom-line>
    # p + geom_dl(aes(label=sample_year),method="last.points", stat="smooth", span="span")

  # Annotate each regression line with the associated r-squared value
    # r-squared for each year
      lms <- lapply(years.l, function(x) lm(Tetragnathid ~ Araneid, data=x))  # output is a list of the three lm objects, one for each year
      rsquareds <- sapply(lms, function(x) {
                              summary(x)$r.squared } )
      intercepts <- sapply(lms, function (x) {
                              summary(x)$
      })
      lapply(lms, fun(x) {myFun()}) -> tester

for (i in length(lms)) {
  myFun(lms[i])
}
    # add the values to the plot
      Fig4 <- Fig4 + geom_dl(aes(label=rsquareds), method="smart.grid")

fit <- lm(Tetragnathid~Araneid+sample_year, data=means.w)
years.l <- split(x=means.w, f=means.w$sample_year)

lapply(X=years.l$site_number, FUN=mean)
sapply(years.l, function(x) mean(x$site_number))
lapply(years.l, function(x) lm(Tetragnathid ~ Araneid, data=x)) -> test

sapply(test, function(x) {
  summary(x)$r.squared
})

test[1] -> fit2011
r.squared <- round(summary(fit)$r.squared, 2)
p.value   <- round(summary(fit)$coefficients[2,4], 4)

myFun <-
  function(lm)
  {
    out <- c(lm$coefficients[1],
             lm$coefficients[2],
             length(lm$model$y),
             summary(lm)$coefficients[2,2],
             pf(summary(lm)$fstatistic[1], summary(lm)$fstatistic[2],
                summary(lm)$fstatistic[3], lower.tail = FALSE),
             summary(lm)$r.squared)
    names(out) <- c("intercept","slope","n","slope.SE","p.value","r.squared")
    return(out)}


# Find the end points of each line in order to annotate them
  # Get the maximum araneid value for each year  
    max.ara <- tapply(X=means.w$Araneid, INDEX=means.w$sample_year, FUN=max, na.rm=TRUE)
  # Calculate a row index for these data points
    label.index <- which(means.w$Araneid %in% test)
  # Add the data series annotations to the plot
    p3 + geom_text(data=means.w[label.index, ], aes(label=sample_year), hjust=0, vjust=0.5)
    p3 + annotate(geom="text", x=420, y=420, label="1:1 line")

p <- ggplot(data=means.w, aes(x=Araneid, y=Tetragnathid, color=factor(sample_year))) +   #remove the color=factor(sample_year) to get all years to plot as one color
  geom_point(shape=19) + 
  #   geom_smooth(method=lm, se=FALSE) +
  stat_smooth(method=loess, se=FALSE, ) +
  geom_dl(aes(label=sample_year), method="last.points", (stat="smooth", method=lm))


p + geom_dl(aes(label=sample_year), method="last.points", args)
p + geom_dl(aes(label=sample_year), list(method="last.points", stat="smooth", lm))

test <- tapply(X=means.w$Araneid, INDEX=means.w$sample_year, FUN=max, na.rm=TRUE)
which(tapply(X=means.w$Araneid, INDEX=means.w$sample_year, FUN=max, na.rm=TRUE))
means.w$Araneid %in% test
which(means.w$Araneid %in% test)

model <- lm(means.w$Tetragnathid ~ means.w$Araneid)
