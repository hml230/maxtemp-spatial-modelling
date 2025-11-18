# Fit OLS using basis functions
fit_ols_single <- function(lyr_idx, train_data, S) 
{
  cat("Fitting OLS for layer ", lyr_idx, "\n")
  train_layer <- train_data[[lyr_idx]]
  vals <- values(train_layer)[1:nrow(S)]
  valid_idx <- which(!is.na(vals))

  y <- vals[valid_idx]
  X <- S[valid_idx, , drop = FALSE] # basis for observed locations

  # Create data frame
  X_df <- as.data.frame(X)
  ols_model <- lm(y ~ ., data = X_df)

  # Predict for all locations
  S_df <- as.data.frame(S)
  names(S_df) <- names(X_df)
  preds <- predict(ols_model, newdata = S_df)

  result_raster <- train_layer
  values(result_raster) <- preds

  return(list(model = ols_model, raster = result_raster))
}

# OLS evaluation
calculate_ols_fit_metrics <- function(lyr_idx, ols_results, temp_data, holdout_masks) 
{
  mask <- holdout_masks[[lyr_idx]]
  pred_vals <- values(ols_results[[lyr_idx]])[mask]
  obs_vals <- values(temp_data[[lyr_idx]])[mask]

  complete_idx <- complete.cases(pred_vals, obs_vals)
  pred_vals <- pred_vals[complete_idx]
  obs_vals <- obs_vals[complete_idx]

  return(
    data.frame(
    observed = obs_vals,
    predicted = pred_vals,
    pred_mean = mean(pred_vals),
    residual = obs_vals - pred_vals,
    month = lyr_idx,
    mse = mean((pred_vals - obs_vals)^2),
    rmse = sqrt(mean((pred_vals - obs_vals)^2)),
    mae = mean(abs(pred_vals - obs_vals))
    )
  )
}

# Calculate residuals
calculate_ols_lyr_resid <- function(lyr_idx, ols_results, temp_data) {
  ols_vals <- values(ols_results[[lyr_idx]])
  orig_vals <- values(temp_data[[lyr_idx]])
  residuals <-  orig_vals - ols_vals
  
  df <- as.data.frame(temp_data[[lyr_idx]], xy = TRUE, na.rm = FALSE)
  df$residual <- residuals
  df$layer <- lyr_idx
  df$method <- "TSR"
  
  return(df)
}

calculate_crps <- function(lyr_idx, results_df) 
{
  month_data  <- results_df %>% dplyr::filter(month == lyr_idx)
  
  # Extract numeric vectors
  obs <- month_data$observed
  pred_mean <- month_data$pred_mean
  
  # Calculate CRPS for Gaussian distribution
  # Mean = pred_mean, SD = standard deviation of residuals
  sd_resid <- sd(month_data$residual)
  
  crps_scores <- verification::crps(obs, cbind(pred_mean, sd_resid))
  
  return(as.data.frame(
    tibble(
    month = lyr_idx,
    mean_crps = mean(crps_scores$CRPS),
    median_crps = median(crps_scores$CRPS)
    )
  )
  )
}
