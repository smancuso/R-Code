SDadjP.lm <- function(model, B = 10000, seed = 1, rob = FALSE) {
  # Stop if model not 'lm' object
  stopifnot(class(model)[1] == "lm")
  # Set seed for reproducibility
  set.seed(seed)
  # Set robust estimator
  est <- ifelse(rob == TRUE, "HC3", "const")
  # Set covariance matrix
  covMat <- vcovHC(model, est)
  # Get observed t-values
  tObs <- as.vector(abs(coeftest(model, vcov. = covMat)[-1, 3]))
  # Sort observed t-values
  tObsSort <- as.vector(sort(tObs))
  # Save order of the observed t-values
  tOrder <- order(tObs)
  # Get observed p-values
  pObs <- as.vector(coeftest(model, vcov. = covMat)[-1, 4])
  # Sort observed p-values in descending order
  pObsSort <- as.vector(sort(pObs, decreasing = TRUE))
  # Save order of the observed p-values
  pOrder <- order(pObs)
  # Number of terms
  nPred <- length(pObs)
  # Create blank matrices
  tResample <- matrix(
    nrow = B,
    ncol = nPred
  )
  pResample <- matrix(
    nrow = B,
    ncol = nPred
  )
  # Number of observations
  nObs <- nrow(model$model)
  # Resampling loop
  for (i in 1:B) {
    # Resample residuals
    Ystar <- sample(
      model$residuals, 
      size = nObs, 
      replace = TRUE
    )
    # Run model with resampled residuals as DV
    ModResample <- lm(
      as.formula(paste0("Ystar~", paste0(labels(terms(model)), collapse = "+"))),
      data = model$model[, -1]
    )
    # Covariance matrix of resampled data
    covResample <- vcovHC(ModResample, est)
    # Get resampled t-values based on order of the observed
    # t-values and store in matrix
    tResample[i, ] <- as.vector(
      abs(coeftest(ModResample, vcov. = covResample)[-1, 3])[tOrder]
    )
    # Get resampled p-values based on order of the observed
    # p-values and store in matrix
    pResample[i, ] <- as.vector(
      coeftest(ModResample, vcov. = covResample)[-1, 4][pOrder])
  }
  # Get the cumulative maxima
  Qtmat <- t(apply(tResample, 1, cummax))
  # Get the cumulative minima
  QPmat <- t(apply(pResample, 1, cummin))
  # Compute adjusted p-values:
  # maxT
  p.maxT <- apply(t(matrix(rep(tObsSort, B), nPred)) < Qtmat, 2, mean)
  # minP
  p.minP <- apply(t(matrix(rep(pObsSort, B), nPred)) > QPmat, 2, mean)
  # Turn off scientific notation
  options(scipen = 999)
  # Sort the p-values
  p.order.maxT <- match(tObs, tObsSort)
  p.order.minP <- match(pObs, pObsSort)
  # Get predictor names
  varnames <- labels(terms(model))
  # Outputs the raw and adjusted p-values
  # Create data.frame to output
  adjP <- data.frame(
    p.raw = pObs,
    p.maxT = p.maxT[p.order.maxT],
    p.minP = p.minP[p.order.minP],
    row.names = varnames
  )
  # Print data.frame
  cat("\n", "Step-down residuals-based resampling p-value adjustment (Westfall & Young, 1993)",
    "\n",
    "\n", "Number of bootstrapped resamples: ", B,
    "\n",
    "\n",
    sep = ""
  )
  print(adjP)
  # Turn back on scientific notation
  options(scipen = 0)
}
