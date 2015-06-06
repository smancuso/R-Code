# -------------------------------------------------------------------------
# Project: Bootstrap AN(C)OVA     
# Coded by: Sam Mancuso
# Mail: sammancuso.wordpress.com
# Date: 18 May 2015
# Version: 1.0
# Revision Date: NA
# -------------------------------------------------------------------------

# Load Libraries ----------------------------------------------------------

if (!require("car")) {install.packages("car"); library("car")}
if (!require("boot")) {install.packages("boot"); library("boot")}

# Bootstrapped AN(C)OVA Options -------------------------------------------
#
# formula:      Formula for AN(C)OVA model
#
# conf.int:     Confidence intervals (default = 0.95)
#
# dec:          Number of decimal places (default = 3)
#
# reps:         Number of bootstrap replications (default = 1000)

# Bootstrap statistic -----------------------------------------------------

anboot <- function(formula, data, i){
 
  data.resamp <- data[i,] # Resample rows
  
  fit <- lm(formula, data.resamp)
  
  f.coef <- Anova(fit, type = 3)[["F value"]] # Extract F-values
  
  return <- c(f.coef)
}

# Function to run bootstrap and provide output ----------------------------
boot.anova <- function(data, formula, conf.int = 0.95, dec = 2, reps = 1000) {
  
  # Fit non-bootstrapped ANOVA to obtain Sum Sq and df
  lm.fit <- lm(formula, data = data)
  # Get DV name
  dv.name <- colnames(model.frame(lm.fit))[1]                    
  # Get IV names
  iv.names <- attr(lm.fit$terms , "term.labels")
  
  # Fit ANOVA model with Type III Sum of Squares
  anova.fit <- Anova(lm.fit, type = 3)
  
  
  # Call anboot function
  anboot.run <- boot(statistic = anboot, formula = formula, data = data, R = reps)
  
  f.ssq <- anova.fit[["Sum Sq"]]           # Sum of Squares (Type III)    
  f.values <- summary(anboot.run)$original # F-values
  f.se <- summary(anboot.run)$bootSE       # Bootstrap SE
  f.z <- f.values/f.se                     # z-value
  f.z.p <- 2*pnorm(-abs(f.z))              # Pr(>|z|)
  f.df <- anova.fit[["Df"]]                # df
  
  f.crit <- NULL # To prevent errors
  
  # Critical F-Value (excludes the 'Residuals' term)
  for(i in 1:(length(f.df)-1)){
    f.crit[i] <- qf(conf.int, f.df[i], f.df[length(f.df)])
  }
    
  # Get Bootstrapped confidence intervals
  ci.lower <- NULL
  ci.upper <- NULL
  
  for(i in 1:length(f.values)){
    boot.bca <- boot.ci(anboot.run, index = i, type = "bca", conf = conf.int)$bca
    ci.lower[i] <- boot.bca[4]
    ci.upper[i] <- boot.bca[5]     
  }
  
  # P-value stars
  fstar <- NULL
 
  for(i in 1:length(na.omit(f.z.p))) {
    if(f.z.p[i] < .001) {
      fstar[i] <- "***"
    } else if(f.z.p[i] < .01) {
        fstar[i] <- "**"
    } else if(f.z.p[i] < .05) {
        fstar[i] <- "*"
    } else if(f.z.p[i] < .10){
        fstar[i] <- "."
    } else {
        fstar[i] <- ""
    }
  }
  
  # Column names
  cnames <- c("Sum Sq", "Df", "F value", "SE", "LB", "UB", "Crit F", "z", "Pr(>|z|)", "")
  
  # Row names
  rnames <- c("(Intercept)", iv.names, "Residuals")
  
  # Turn warnings off as the variables are not equal length when creating output.df
  options(warn = -1)
  
  # Create data frame to use as output
  output.df <- as.data.frame(cbind(round(f.ssq, dec), 
                                   f.df, 
                                   format(round(f.values, dec), dec), 
                                   format(round(f.se, dec), dec), 
                                   format(round(ci.lower, dec), dec), 
                                   format(round(ci.upper, dec), dec), 
                                   format(round(f.crit, dec), dec),
                                   format(round(f.z, dec), dec), 
                                   signif(f.z.p, 3),
                                   fstar))
  
  
  # Remove redundant values for the Residuals row
  output.df[nrow(output.df), 3:10] <- NA
  
  # Set column and row names
  colnames(output.df) <- cnames
  rownames(output.df) <- rnames
  
  
  # Print results
  cat("BCa Bootstrap Anova Table (Type III tests)", 
      "\n",
      "\n", "Response: ", dv.name,
      "\n",
      sep ="")
  
  print(as.matrix(output.df), justify = "right", na.print = "" , quote = FALSE )
   
  cat("---",
      "\n", "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1)",
      sep = "")
  
  # Turn warnings on
  options(warn = 0)
}