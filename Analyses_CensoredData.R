## Formally alculate the statistics re censored data 
## that are reported in the main text or in the supp 
## data section.
##
## A. Correlations between sum-zero and sum-halfdl pcb totals.
## B. Correlations between sample mass and sample sum-pcb (sum-zero).
## C. Rank order correlations between site means - 2011 vs other yrs.
## D. Correlations between ref-con sums vs full sums.

#### A. SUBS ZERO VS SUBS HALFDL ####
# Data including sums calculated with subs half mdl
  sums.data <- read.csv(paste(DirOut, "AllSampleData_CompareSums.csv", sep=""))

# Function to fit a lm and return parameters for a given data set (df)
  fitlm <- function(df) {
    fit.cur <- lm(sum.halfdl ~ sum.zero, data=df)
    r.sq  <- summary(fit.cur)$r.squared
    p.val <- summary(fit.cur)$coefficients[2, 4]
    n.dat <- length(fit.cur$residuals)
    dat.cur <- data.frame(r.sq=r.sq, p.val=p.val, n.dat=n.dat)
    return(dat.cur)
  }

# Apply to each year/category grp of samples
  smpgrp.l <- split(sums.data, f=list(sums.data$category, sums.data$sample_year))
  zerohalf.mods <- do.call("rbind", lapply(X=smpgrp.l, FUN=fitlm))

  rm(smpgrp.l)


#### B. SAMPLE MASS VS SAMPLE CONCENTRATION ####
# Correlation between mass and concentration for all spider samples
  fit.spid.all <- lm(sum.pcb ~ smp.mass, 
                     data=samples[samples$category!="Sediment", ])
  summary(fit.spid.all)
  
  # Numbers for text
  print(paste("Relationship for all spider samples:", 
              " r.sq=", round(summary(fit.spid.all)$r.squared, 3),
              "; p=", round(summary(fit.spid.all)$coefficients[2, 4], 3),
              "; n=", length(fit.spid.all$residuals), sep=""))
      
# Correlation between mass and concentration for spider samples from 2012-13 only
  fit.spid.1213 <- lm(sum.pcb ~ smp.mass, 
                      data=samples[samples$category!="Sediment" & samples$sample_year!=2011, ])
  summary(fit.spid.1213)
  
  # Numbers for text
  print(paste("Relationship for 2012-13 spider samples:", 
              " r.sq=", round(summary(fit.spid.1213)$r.squared, 3),
              "; p=", round(summary(fit.spid.1213)$coefficients[2, 4], 3),
              "; n=", length(fit.spid.1213$residuals), sep=""))


#### C. RANK ORDER CORRELATION IN SITE MEANS, 2011 vs other yrs ####
# Calculate site means using both subs.zero and subs.halfdl sum-PCB concentrations
  # Data including sums calculated with subs half mdl
    sums.data <- read.csv(paste(DirOut, "AllSampleData_CompareSums.csv", sep=""))
  # Use spider data only
    spidsums <- filter(sums.data, category!="Sediment")
  # Calculate site means
    by.site <- group_by(spidsums, sample_year, category, site_number)
    spid.means <- summarize(by.site,
                            sum.zero   = mean(sum.zero),
                            sum.halfdl = mean(sum.halfdl))

# Basic stats for reporting in text
  # 2012-13 spider smp:
  spidsums.1213 <- filter(spidsums, sample_year!=2011)
    # n. samples
    paste("n =", nrow(spidsums.1213))  
    # Sample mass (mean/sd)
    paste("mean mass=", round(mean(spidsums.1213$smp.mass), 2),
          ", StDev=", round(sd(spidsums.1213$smp.mass), 3), sep="")#avg sample mass 
    # ssMDL (mean/sd)
    spidcons.1213 <- filter(congeners, category!="Sediment" & sample_year!=2011)
    paste("mean MDL=", round(mean(spidcons.1213$mdl.ss), 2),
          ", StDev=", round(sd(spidcons.1213$mdl.ss), 3), sep="")#avg sample mass 

  # 2011 spider smp:
  spidsums.11 <- filter(spidsums, sample_year==2011)
    # n. samples
    paste("n =", nrow(spidsums.11))  
    # Sample mass (mean/sd)
    paste("mean mass=", round(mean(spidsums.11$smp.mass), 2),
          ", StDev=", round(sd(spidsums.11$smp.mass), 3), sep="")#avg sample mass 
    # ssMDL (mean/sd)
    spidcons.11 <- filter(congeners, category!="Sediment" & sample_year==2011)
    paste("mean MDL=", round(mean(spidcons.11$mdl.ss), 2),
          ", StDev=", round(sd(spidcons.11$mdl.ss), 3), sep="")#avg sample mass 
  # Clean up
    rm(spidsums.1213); rm(spidcons.1213); rm(spidsums.11); rm(spidcons.11)

# Generate a hybrid wide/long df to match each 2012-13 sitemean with value from 2011
  # Full long form 
    spidmeans.l <- melt(spid.means, 
                       id.vars=c("sample_year", "category", "site_number"),
                       measure.vars=c("sum.zero", "sum.halfdl"), 
                       variable.name="sum.type", value.name="mean.conc") 

  # Split into separate dfs, in order to merge back into hybrid form
    spid11   <- filter(spidmeans.l, sample_year==2011)
    spid11   <- rename(spid11, conc.11 = mean.conc)
    spid1213 <- filter(spidmeans.l, sample_year!=2011)
    spidmeans.h <- merge(spid1213, spid11[, 2:5], 
                         by=c("category", "site_number", "sum.type"),
                         all.x=TRUE, all.y=FALSE)

  # Clean up 
    rm(spidsums); rm(by.site); rm(spid.means); rm(spidmeans.l); rm(spid11); rm(spid1213)

# Calculate rank-order correlations
  # Data frame with all rho/pval/n
    by.grp <- group_by(spidmeans.h, category, sample_year, sum.type)
    spear.coefs <- summarize(by.grp,   #ignore warnings about cannot compute exact pvals with ties
                             rho   = round(cor.test(mean.conc, conc.11, method="spearman")[[4]], 3),
                             p.val = round(cor.test(mean.conc, conc.11, method="spearman")[[3]], 3),
                             n.dat = sum(!is.na(mean.conc) & !is.na(conc.11)) )
  # Data frame to compare rho values only
    (rho.vals <- dcast(data=spear.coefs, category + sample_year ~ sum.type, value.var="rho"))

  # Average rho value across the four groups
    mean(rho.vals$sum.zero)   #for the subs zero grp
    mean(rho.vals$sum.halfdl) #for the subs half detlim grp

  # Clean up
    rm(by.grp)


#### D. COMPARE SUMS OF VERY COMMON CONGENERS ONLY ####
# Use \sums.data\, read in during part B

# Magnitude of sum-ref as a percentage of sum-full
  spidsums <- filter(sums.data, category!="Sediment" & sample_year!=2011 & sum.zero>0) 
  spidsums$ref.pct <- spidsums$sum.ref.cons / spidsums$sum.zero
  # Min and max pct
  range(spidsums$ref.pct, na.rm=TRUE)
  # Mean and SD
  paste("mean=", round(mean(spidsums$ref.pct, na.rm=TRUE), 3), 
        "; sd=", round(sd(spidsums$ref.pct, na.rm=TRUE), 3) ) 
        
# Linear correlation between sum.ref and sum.zero
  # Araneids
  summary(lm(sum.zero ~ sum.ref.cons, data=spidsums[spidsums$category=="Araneid", ]))
  sum(complete.cases(spidsums[spidsums$category=="Araneid", c("sum.zero", "sum.ref.cons")]))  #n = df + 2
  # Tetragnathids
  summary(lm(sum.zero ~ sum.ref.cons, data=spidsums[spidsums$category=="Tetragnathid", ]))


#### END SCRIPT ####