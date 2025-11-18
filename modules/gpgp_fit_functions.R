fit_gpgp_single <- function(month_col, train_df, temp_data) {
  cat("Fitting GpGp for", month_col, "...\n")
  
  # PREPROCESSING TRAIN DATA
  train_complete <- train_df %>%
    filter(!is.na(.data[[month_col]])) %>%
    select(x, y, all_of(month_col))
  
  y_train <- train_complete[[month_col]]
  locs_train <- as.matrix(train_complete[, c("x", "y")])
  X_train <- matrix(1, nrow = length(y_train), ncol = 1)
  
  # GET PREDICTION (MASKED) DATA
  pred_complete <- train_df %>%
    filter(is.na(.data[[month_col]])) %>%
    select(x, y)
  
  locs_pred <- as.matrix(pred_complete[, c("x", "y")])
  X_pred <- matrix(1, nrow = nrow(pred_complete), ncol = 1)
  
  # Fit model
  fit <- gpgp::fit_model(
    y = y_train,
    locs = locs_train,
    X = X_train,
    covfun_name = "matern_isotropic",
    reorder = TRUE,
    m = 30,
    silent = TRUE
  )
  
  # Predict at the masked locations
  preds <- predictions(
    fit = fit,
    locs_pred = locs_pred,
    X_pred = X_pred,
    m = 30
  )
  
  # Create full grid predictions dataframe
  pred_df <- data.frame(
    x = pred_complete$x,
    y = pred_complete$y,
    mu = preds$means,
    sd = sqrt(preds$vars)
  )
  
  # Add training data predictions (use observations)
  train_pred_df <- data.frame(
    x = train_complete$x,
    y = train_complete$y,
    mu = y_train,
    sd = 0  # Zero variance at observed locations
  )
  
  # Combine for full grid
  full_pred_df <- bind_rows(pred_df, train_pred_df)
  
  return(list(
    fit = fit,
    predictions = full_pred_df,
    pred_only = pred_df
  ))
}

# Apply the fit to first 5 months (Jan-May)
cat("Performing GpGp fitting for first 5 months...\n")
gpgp_fits <- map(month_names, ~ fit_gpgp_single(.x, train_df, temp_data))
names(gpgp_fits) <- month_names

# Extract predictions for evaluation
preds_list <- map(gpgp_fits, ~ .x$predictions)

# Extract true values at masked locations
true_values_at_mask <- map(1:5, function(i) {
  month_col <- month_names[i]
  return(temp_data[[i]][holdout_masks[[month_col]]])
})

# Extract predictions at masked locations
predictions_at_mask <- map(1:5, function(i) {
  month_col <- month_names[i]
  
  # Get prediction dataframe
  pred_df <- preds_list[[i]]
  
  # Convert to temp_data format
  pred_raster <- rast(
    xmin = ext(temp_data)[1], xmax = ext(temp_data)[2],
    ymin = ext(temp_data)[3], ymax = ext(temp_data)[4],
    resolution = res(temp_data)
  )
  
  # Rasterize
  pred_raster <- rasterize(
    vect(pred_df, geom = c("x", "y")),
    pred_raster,
    field = "mu"
  )
  
  # Extract using the boolean mask
  return(pred_raster[holdout_masks[[month_col]]])
})