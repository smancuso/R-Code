# Load required packages --------------------------------------------------
if (!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(car, boot)

# Function ----------------------------------------------------------------
nullboot.Anova <- function(null.model, full.model,
                           B = 1000, scaled = TRUE, seed = 1234,
                           ci = TRUE, cent = .95, dec = 5,
                           ci.type = c("bca", "perc")) {

  # Returns error if models are not 'lm' class
  stopifnot(class(null.model) == "lm" | class(full.model) == "lm")

  # Set seed
  set.seed(seed)

  # Extract data from the two models
  data_full <- model.frame(full.model)
  data_null <- model.frame(null.model)

  # Get residuals from null model and recenter and rescale based on
  # Davidson and MacKinnon (2004) if requested.
  if (isTRUE(scaled)) {

    # Number of observations
    N <- nrow(data_null)

    # Number of predictors (from full model)
    K <- length(labels(terms(full.model)))

    # Extract and scale residuals (if requested) from null model
    Er <- (residuals(null.model) - mean(residuals(null.model))) * sqrt(N / (N - K))
  } else {
    Er <- residuals(null.model)
  }

  # Predicted values from null model
  Yhat <- fitted(null.model)

  # Get ANOVA statistics for full model
  mod_ANOVA <- Anova(full.model, type = 3)

  # Get observed F-values (excluding Intercept and NA [Residuals])
  Fobs <- as.vector(na.omit(mod_ANOVA[-1, 3]))

  # Updated formula for bootstapping with Ystar as DV
  formula_resample <- update.formula(formula(full.model), Ystar ~ .)

  # Bootstrap function
  bs.Anova <- function(data, i) {
    data_full$Ystar <- Yhat + Er[i]

    # Fit linear model
    model_resample <- lm(formula_resample, data = data_full)

    # Get F-values from ANOVA
    Fstar <- as.vector(na.omit(Anova(model_resample, type = 3)[-1, 3]))

    return(Fstar)
  }

  # Run bootstrapping
  bootAnova <- boot(
    data_full, statistic = bs.Anova, R = B,
    parallel = "snow"
  )

  # Bootstrapped F-values
  Fstar <- bootAnova$t

  # The proportion of resampled F-values greater than or equal to
  # the observed F-values

  pBoot <- vector()

  for (i in 1:length(Fobs)) {
    pBoot[i] <- (sum(Fstar[, i] >= Fobs[i]) + 1) / (length(Fstar[, i]) + 1)
  }

  # Create data.frame for outputting
  varNames <- labels(terms(full.model))
  pRaw <- as.vector(na.omit(mod_ANOVA[-1, 4]))

  # Calculate confidence intervals if requested

  if (isTRUE(ci)) {
    ci.type <- match.arg(ci.type)
    FstarCI <- list()
    ciLB <- vector()
    ciUB <- vector()
    bcaLB <- vector()
    bcaUB <- vector()

    for (i in 1:length(Fobs)) {
      suppressWarnings(FstarCI[[i]] <- boot.ci(
        bootAnova,
        type = ci.type,
        index = i,
        conf = cent
      ))

      # Percentile confidence intervals
      if (ci.type == "perc") {
        ci_desc <- "Percentile"
        ciLB[i] <- FstarCI[[i]]$perc[, 4]
        ciUB[i] <- FstarCI[[i]]$perc[, 5]
      }

      # BCa confidence intervals
      if (ci.type == "bca") {
        ci_desc <- "BCa"
        ciLB[i] <- FstarCI[[i]]$bca[, 4]
        ciUB[i] <- FstarCI[[i]]$bca[, 5]
      }
    }

    bootOut <- data.frame(
      "F.value" = Fobs,
      "p.value" = pRaw,
      "p.boot" = pBoot,
      "CI.LB" = ciLB,
      "CI.UB" = ciUB,
      row.names = varNames
    )
  } else {
    ci_desc <- "Not requested"
    bootOut <- data.frame(
      "F.value" = Fobs,
      "p.value" = pRaw,
      "p.boot" = pBoot,
      row.names = varNames
    )
  }

  cat(
    "\n", "Bootstrapped ANOVA with Type III tests",
    "\n",
    "\n", "Resampling type: ",
    ifelse(isTRUE(scaled),
      "rescaled residuals",
      "residuals"
    ),
    "\n", "Number of bootstrap resamples: ", B,
    "\n", "Bootstrapped confidence interval type: ", ci_desc,
    "\n", "Confidence interval: ", 100 * cent, "%",
    "\n",
    "\n", sep = ""
  )

  # Turn off scientific
  options(scipen = 999)

  print(round(bootOut, dec))

  # Turn on scientific
  options(scipen = 0)
}
