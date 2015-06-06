# -------------------------------------------------------------------------
# Sample Size Estimation for Comparisons between Two Groups
# in Experimental Designs
# Coded by: Sam Mancuso
# Date: 03 March 2015
# Version: 2.0
# Revision Date: 05 March 2015
# Blog: sammancuso.wordpress.com
# -------------------------------------------------------------------------


expsample <- function(num = 1,       # Number of outcome measures (default = 1)
                      alpha = 0.05,  # Study-wide two-tailed alpha (default = 0.05)
                      points = 1,    # Number of follow-up assessments (default = 1)
                      power = 0.8,   # Power (default = 0.8)
                      avg = F,       # Compare at a single-time point (F) or average over the total follow up period (T)
                      d,             # Cohen's d effect size
                      rho = 0.5,     # Correlation of repeated measures (default = 0.5)
                      attr = 0) {    # Attrition rate (default = 0), can be decimal (e.g.., 0.25) or percentage (e.g., 25)

                      # Applies Bonferroni correction for multiple outcome measures
                      alphaadjust = (alpha/2)/num
                     
                      # Standardised score cutting off alpha/2 proportion of each tail (standard normal distribution)
                      za = qnorm(1 - alphaadjust)
                     
                      # Standardised score cutting off the upper Beta proportion
                      zb = qnorm(power)
                     
                      # Average over total follow up period (and more than one follow-up point)
                      if(avg == T & points > 1) {   
                        numerator = 2*((za + zb)^2)*(1 + ((points - 1)*rho))
                        denominator = points*(d^2)
                      } else {
                          # Difference between single time point
                          numerator = 2*((za + zb)^2)   
                          denominator = d^2
                      }
                      
                      # Number per group
                      npergroup = ceiling(numerator/denominator)
                      Ntotal = 2*npergroup
                      
                      # Inflated sample size for attrition
                    
                      # Change attrition from percentage to decimal if required
                      if(attr > 1){
                        attr = attr/100
                      }
                     
                      # Inflate sample size for attrition
                      nattr = ceiling(npergroup * (1/(1 - attr)))
                      Nattr = 2*nattr
                     
                      # Output
                      cat("\n",
                          "Parameters Specified",
                          "\n", "--------------------",
                          "\n", "Alpha: ", alpha, 
                          "\n", "Power: ", power, 
                          "\n", "Effect size (Cohen's d): ", d,
                          "\n", "Number of outcome measures: ", num, sep = "")
                      
                      # If averaging over the total follow-up period
                      if(avg == T){
                        cat("\n", "Average over the total follow-up period: ", avg,
                            "\n", "Numer of follow-up measurements: ", points, 
                            "\n", "Correlation between repeated measures: ", rho, sep = "")
                      }
                          
                      cat("\n",
                          "\n", "Required Sample Size",
                          "\n", "--------------------",
                          "\n", "Number of participants per group: ", npergroup,
                          "\n", "Total participants required: ", Ntotal, sep ="")
                     
                      # Following output is displayed if attrition rate specified
                      if(attr > 0){
                        cat("\n",
                            "\n", "Sample size adjusted for attrition rate",
                            "\n", "---------------------------------------",
                            "\n", "Expected attrition rate: ", attr,
                            "\n", "Number of participants per group: ", nattr,
                            "\n", "Total participants required: ", Nattr,
                            sep ="")
                      }
}