# =================================================================
#
# Analysis of repeated single-case data
#
# Mueser, K. T., Yarnold, P. R., & Foy, D. W. (1991).
# Statistical Analysis for Single-Case Designs Evaluating Outcome 
# of Imaginal Exposure Treatment of Chronic PTSD. 
# Behavior Modification, 15(2), 134-155.
# 
# Coded by: Sam Mancuso
# Blog: sammancuso.wordpress.com
# 
# Version: 2.00
#
# Date: 9 April 2014
# Updated: 28 January 2015
#
# =================================================================
#
# Options
# -------
# 
# baseline: select comparison assessment point (default = 1)
#
# level: significance level (default = .05)
#
# two.sided: two-sided test (default = TRUE)
#
# =================================================================

single.case <- function(data, level = .05, two.sided = T) {
  
  #----1-lag Autocorrelation----
  
  # Get length of array
  j <- length(data)
  
  if(j >= 4) {    # Check minimum of 4 data points
    
    # Create arrays
    x <- data[(2):j]
    y <- data[1:(j - 1)]
    
    # Calculate ACF(1)
    acf <- abs(cor(x, y))
    
    # Get critical z value
    sides <- ifelse(two.sided == T, 2, 1)
    
    crit.zval <- abs(qnorm(level/sides))
    
    # Calculate critical difference (CD)
    cd <- crit.zval*(sqrt(j*(1-acf)))
    
    #----Ipsative Data----
    
    # Calculate ipsative z-scores
    zi <- scale(data)[, 1] # Ipsative Z-scores
    
    # Calculate Difference Scores
    x.1 <- NULL
    x.2 <- NULL
    z.x.1 <- NULL
    z.x.2 <- NULL
    df.z <- NULL
    df.raw <- NULL
    comp.label <- NULL
    k <- 0  # Initialise counter
    
    for(m in 1:(j - 1)) {
      
      for(n in (m + 1):j) {
          
        k <- k + 1
        
        # Raw scores
        x.1[k] <- data[m]
        x.2[k] <- data[n]
        
        # Raw score difference
        df.raw[k] <- data[n] - data[m]
        
        # Ipsative z-scores
        z.x.1[k] <- zi[m]
        z.x.2[k] <- zi[n]
        
        # Ipsative z-score difference
        df.z[k] <- zi[n] - zi[m]
        
        comp.label[k] <- paste(m, " vs ", n, sep = "")
      }
    }
    
    # Identify significant differences
    df.sig <- abs(df.z) > cd
    df.sig[df.sig == T] <- "*"
    df.sig[df.sig == F] <- "ns"
      
    # Round zi and df.z scores to 2 decimal places
    x.1 <- format(round(x.1, 2), 2)
    x.2 <- format(round(x.2, 2), 2)
    df.raw <- format(round(df.raw, 2), 2)
    z.x.1 <- format(round(z.x.1, 2), 2)
    z.x.2 <- format(round(z.x.2, 2), 2)
    df.z <- format(round(df.z, 2), 2)
    
    #----Create table----
    table.out <- cbind(x.1, x.2, df.raw, z.x.1, z.x.2, df.z, df.sig)
    
    # Column names and row names
    colnames(table.out) <- c("x1", "x2", "d.x", "Z.x1", "Z.x2", "d.Z", "Sig")
    rownames(table.out) <- comp.label
    
    # Change table to data frame
    table.out <- as.data.frame(table.out)
    
    # Round acf and cd to 2 decimal places for output
    acf <- format(round(acf, 2), 2)
    cd <- format(round(cd, 2), 2)
    
    # Output
    cat("\n",
        "Single-case data analysis", "\n",
        "\n", sep="")
    
    print(table.out)
    
    cat("\n", "ACF(1) = ", acf, "\n",
        "Critical Difference = ", cd, sep ="")
  } else {
      
    cat("\n",
        "Single-case data analysis", "\n", "\n",
        "Error: ", j, " data points entered. Minimum of 4 required",
        "\n", sep = "")
  }
}
