require(boot)
require(car)
require(MASS)
require(QuantPsyc) # lm.beta() function

#---- Bootstrap Mediation Function ----

boot.med <- function(data, formula.m, formula.y, iv, mv, i) {
  data.resamp <- data[i,]
  reg.m <- lm(formula.m, data=data.resamp)
  reg.y <- lm(formula.y, data=data.resamp)
  # Make sure that you pull the correct coefficients from each regression equation
  # With interactions, quadratics, etc, R may reorder the equation so you need to check it.
  a <- coefficients(reg.m)[iv] # the a path coefficient
  b <- coefficients(reg.y)[mv] # the b path coefficient
  ind <- a*b # the indirect effect (product of a*b)
  dir <- coefficients(reg.y)[iv] # the direct effect (c')
  total <- ind + dir # teturns total effect
  boot.med <- return(c(ind, dir, total))
}

#---- Mediation Analyses ----
#
# data = data frame
# iv = independent variable (X)
# dv = dependent variable (Y)
# mv = mediating variable (M)
# cv = covariate(s) if any. Use: cv = c("var.1", "var.2", "var.3", ..., "var.n")
# reps = number of bootstrap replications (default = 1000)
#

med <- function(data, iv, dv, mv, cv, dec = 4, reps = 1000) {
  # Formulas for Regressions
  if(missing(cv)){
    iv.formula <- as.formula(paste(mv, paste("~"), paste(iv)))
    mv.formula <- as.formula(paste(dv, paste("~"), paste(mv)))
    dv.formula <- as.formula(paste(dv, paste("~"), paste(iv)))
    dvmv.formula <- as.formula(paste(dv, paste("~"), paste(iv), paste("+"), paste(mv)))
  } else {
    iv.formula <- as.formula(paste(mv, paste("~"), paste(iv), paste("+"), paste(cv, collapse="+")))
    mv.formula <- as.formula(paste(dv, paste("~"), paste(mv), paste("+"), paste(cv, collapse="+")))
    dv.formula <- as.formula(paste(dv, paste("~"), paste(iv), paste("+"), paste(cv, collapse="+")))
    dvmv.formula <- as.formula(paste(dv, paste("~"), paste(iv), paste("+"), paste(mv), paste("+"), paste(cv, collapse="+")))
  }
  ## Regressions

  # IV on MV
  iv.out <- lm(iv.formula, data = data)
  iv.coef <- coef(iv.out) # Coefficients
  iv.se <- summary(iv.out)$coefficients[,2] # SEs
  iv.pr <- 2 * pt(abs(iv.coef/iv.se), df=nrow(iv.out$model)-2, lower.tail = FALSE) # p-values
  iv.ll <- iv.coef - 1.96 * iv.se # lower 95% CIs
  iv.ul <- iv.coef + 1.96 * iv.se # upper 95% CIs
  # Standardised Betas
  iv.beta <- lm.beta(iv.out)
  iv.beta <- cbind(NA, iv.beta)
  iv.beta <- rbind(NA, iv.beta)[,2]
  iv.rsquared <- summary(iv.out)$r.squared #R2
  iv.rsquared.pr <- pf(summary(iv.out)$fstatistic[1], summary(iv.out)$fstatistic[2],
                       summary(iv.out)$fstatistic[3], lower.tail = F) # p-value for R2
  # Output Table
  iv.table <- cbind("Estimate" = round(iv.coef, dec),
                    "SE" = round(iv.se, dec),
                    "Pr(>|z|)" = round(2 * pt(abs(iv.coef/iv.se), df=nrow(iv.out$model)-2,
                                              lower.tail = FALSE), dec),
                    "LL" = round(iv.coef - 1.96 * iv.se, dec),
                    "UL" = round(iv.coef + 1.96 * iv.se, dec),
                    "Beta" = round(iv.beta, dec))  
  
  # MV on DV
  mv.out <- lm(mv.formula, data = data)
  mv.coef <- coef(mv.out) # Coefficients
  mv.se <- summary(mv.out)$coefficients[,2] # SEs
  mv.pr <- 2 * pt(abs(mv.coef/mv.se), df=nrow(mv.out$model)-2, lower.tail = FALSE) # p-values
  mv.ll <- mv.coef - 1.96 * mv.se # lower 95% CIs
  mv.ul <- mv.coef + 1.96 * mv.se # upper 95% CIs
  # Standardised Betas
  mv.beta <- lm.beta(mv.out)
  mv.beta <- cbind(NA, mv.beta)
  mv.beta <- rbind(NA, mv.beta)[,2]
  mv.rsquared <- summary(mv.out)$r.squared #R2
  mv.rsquared.pr <- pf(summary(mv.out)$fstatistic[1], summary(mv.out)$fstatistic[2],
                       summary(mv.out)$fstatistic[3], lower.tail = F) # p-value for R2
  # Output Table
  mv.table <- cbind("Estimate" = round(mv.coef, dec),
                    "SE" = round(mv.se, dec),
                    "Pr(>|z|)" = round(2 * pt(abs(mv.coef/mv.se), df=nrow(mv.out$model)-2,
                                              lower.tail = FALSE), dec),
                    "LL" = round(mv.coef - 1.96 * mv.se, dec),
                    "UL" = round(mv.coef + 1.96 * mv.se, dec),
                    "Beta" = round(mv.beta, dec))  
  
  
  # IV on DV
  dv.out <- lm(dv.formula, data = data)
  dv.coef <- coef(dv.out) # Coefficients
  dv.se <- summary(dv.out)$coefficients[,2] # SEs
  dv.pr <- 2 * pt(abs(dv.coef/dv.se), df=nrow(dv.out$model)-2, lower.tail = FALSE) # p-values
  dv.ll <- dv.coef - 1.96 * dv.se # lower 95% CIs
  dv.ul <- dv.coef + 1.96 * dv.se # upper 95% CIs
  # Standardised Betas
  dv.beta <- lm.beta(dv.out)
  dv.beta <- cbind(NA, dv.beta)
  dv.beta <- rbind(NA, dv.beta)[,2]
  dv.rsquared <- summary(dv.out)$r.squared #R2
  dv.rsquared.pr <- pf(summary(dv.out)$fstatistic[1], summary(dv.out)$fstatistic[2],
                       summary(dv.out)$fstatistic[3], lower.tail = F) # p-value for R2
  # Output Table
  dv.table <- cbind("Estimate" = round(dv.coef, dec),
                    "SE" = round(dv.se, dec),
                    "Pr(>|z|)" = round(2 * pt(abs(dv.coef/dv.se), df=nrow(dv.out$model)-2,
                                              lower.tail = FALSE), dec),
                    "LL" = round(dv.coef - 1.96 * dv.se, dec),
                    "UL" = round(dv.coef + 1.96 * dv.se, dec),
                    "Beta" = round(dv.beta, dec))
  
  # IV and MV on DV
  dvmv.out <- lm(dvmv.formula, data = data)
  dvmv.coef <- coef(dvmv.out)
  dvmv.se <- summary(dvmv.out)$coefficients[,2]
  dvmv.pr <- 2 * pt(abs(dvmv.coef/dvmv.se), df=nrow(dvmv.out$model)-2, lower.tail = FALSE)
  dvmv.ll <- dvmv.coef - 1.96 * dvmv.se
  dvmv.ul <- dvmv.coef + 1.96 * dvmv.se
  # Standardised Betas
  dvmv.beta <- lm.beta(dvmv.out)
  dvmv.beta <- cbind(NA, dvmv.beta)
  dvmv.beta <- rbind(NA, dvmv.beta)[,2]
  dvmv.rsquared <- summary(dvmv.out)$r.squared
  dvmv.rsquared.pr <- pf(summary(dvmv.out)$fstatistic[1], summary(dvmv.out)$fstatistic[2],
                         summary(dvmv.out)$fstatistic[3], lower.tail = F)
  # Output Table
  dvmv.table <- cbind("Estimate" = round(dvmv.coef, dec),
                      "SE" = round(dvmv.se, dec),
                      "Pr(>|z|)" = round(2 * pt(abs(dvmv.coef/dvmv.se), df=nrow(data)-2,
                                          lower.tail = FALSE), dec),
                      "LL" = round(dvmv.coef - 1.96 * dvmv.se, dec),
                      "UL" = round(dvmv.coef + 1.96 * dvmv.se, dec),
                      "Beta" = round(dvmv.beta, dec))
  ## Mediation Bootstrap
  medboot.out <- boot(statistic = boot.med,
                      data = data,
                      formula.m = iv.formula,
                      formula.y = dvmv.formula,
                      iv = iv,
                      mv = mv,
                      R = reps,
                      parallel="multicore")
  indirect.b <- round(summary(medboot.out)$original[1], dec) # Indirect B
  indirect.se <- round(summary(medboot.out)$bootSE[1], dec) # Indirect SE
  indirect.bca <- boot.ci(medboot.out, type="bca", index = 1) # Indirect BCA
  indirect.bca.l <- round(indirect.bca$bca[4], dec) # Extract lower BCa CI
  indirect.bca.u <- round(indirect.bca$bca[5], dec) # Extract upper BCa CI
  direct.b <- round(summary(medboot.out)$original[2], dec) # Direct B
  direct.se <- round(summary(medboot.out)$bootSE[2], dec) # Direct SE
  direct.bca <- boot.ci(medboot.out, type="bca", index = 2) # indirect
  direct.bca.l <- round(direct.bca$bca[4], dec) # Extract lower BCa CI
  direct.bca.u <- round(direct.bca$bca[5], dec) # Extract upper BCa CI
  total.b <- round(summary(medboot.out)$original[3], dec) # Direct B
  total.se <- round(summary(medboot.out)$bootSE[3], dec) # Direct SE
  total.bca <- boot.ci(medboot.out, type="bca", index = 3) # Total Effect
  total.bca.l <- round(total.bca$bca[4], dec) # Extract lower BCa CI
  total.bca.u <- round(total.bca$bca[5], dec) # Extract upper BCa CI
  ## Output

  cat(paste("\n",
            "Independent Variable: ", iv,
            "\n", "Dependent Variable: ", dv,
            "\n", "Mediator Variable: ", mv,
            "\n",
            "\n",
            "IV (and CVs) on MV: Path a",
            "\n",
            sep = ""))
    
  print(iv.table)

  cat(paste("\n",
            "R-Squared: ", round(iv.rsquared, dec), ", p-value: ", round(iv.rsquared.pr, dec),
            "\n",
            sep=""))
  
  cat(paste("\n",
            "\n",
            "MV (and CVs) on DV: Path b",
            "\n",
            sep = ""))
  print(mv.table)
  
  cat(paste("\n",
            "R-Squared: ", round(mv.rsquared, dec), ", p-value: ", round(mv.rsquared.pr, dec),
            "\n",
            sep=""))
  
  cat(paste("\n",
            "\n",
            "IV (and CVs) on DV: Path c",
            "\n",
            sep = ""))
  print(dv.table)
  cat(paste("\n",
            "R-Squared: ", round(dv.rsquared, dec), ", p-value: ", round(dv.rsquared.pr, dec),
            sep=""))
  
  cat(paste("\n",
            "\n",
            "IV and MV (and CVs) on DV",
            "\n",
            sep = ""))
  print(dvmv.table)
  cat(paste("\n",
            "R-Squared: ", round(dvmv.rsquared, dec), ", p-value: ", round(dvmv.rsquared.pr, dec),
            "\n",
            sep=""))
  cat("\n",
      "============================================", "\n",
      "Bootstrap Mediation Effects with 95% BCa CIs", "\n",
      "============================================", "\n",
      "\n",
      "Indirect effect", "\n",
      "---------------", "\n",
      paste("B = ", indirect.b, ", SE = ", indirect.se, ", 95% CI [", indirect.bca.l, ", ", indirect.bca.u, "]", sep=""),
      "\n", "\n",
      "Direct effect", "\n",
      "-------------", "\n",
      paste("B = ", direct.b, ", SE = ", direct.se, ", 95% CI [", direct.bca.l, ", ", direct.bca.u, "]", sep=""),
      "\n", "\n",
      "Total effect", "\n",
      "-------------", "\n",
      paste("B = ", total.b, ", SE = ", total.se, ", 95% CI [", total.bca.l, ", ", total.bca.u, "]", sep=""),
      "\n",
      sep ="")
} 