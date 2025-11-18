rm(list = ls())
library(tidyr)
library(dplyr)
library(terra)
library(purrr)
library(GpGp)
library(ggplot2)
library(verification)

source("./modules/preprocessing_functions.R")
set.seed(1111)

# ============================================================================
# DATA LOADING AND PREPROCESSING
# ============================================================================

cat("Loading and aggregating temperature data...\n")
data_path <- "./data/agcd_v1_tmax_mean_r005_monthly_2017.nc"
temp_data <- rast(data_path)
temp_data <- aggregate(temp_data, fact = 3)

# Create train/test splits for each month (80/20 split)
cat("Creating train/test splits...\n")
masked <- map(1:nlyr(temp_data), ~ mask_random_cells(temp_data[[.x]], 0.2))

# Extract training layers and holdout masks
train_layers <- map(masked, "train")
train_raster <- rast(train_layers)
names(train_raster) <- names(temp_data)

holdout_masks <- map(masked, "mask")
names(holdout_masks) <- names(temp_data)

# Convert to dataframe for GpGp
train_df <- as.data.frame(train_raster, xy = TRUE)

# Focus on first 5 months (Jan-May)
month_names <- names(holdout_masks)[1:5]

# ============================================================================
# GAUSSIAN PROCESS MODELING
# ============================================================================

# Fit GpGp model for a single month
fit_gpgp_single <- function(month_col, train_df, temp_data) {
  cat("  Fitting", month_col, "...\n")
  
  # Prepare training data
  train_complete <- train_df %>%
    filter(!is.na(.data[[month_col]])) %>%
    dplyr::select(x, y, all_of(month_col))
  
  y_train <- train_complete[[month_col]]
  locs_train <- as.matrix(train_complete[, c("x", "y")])
  X_train <- matrix(1, nrow = length(y_train), ncol = 1)
  
  # Prepare prediction locations (masked cells)
  pred_complete <- train_df %>%
    filter(is.na(.data[[month_col]])) %>%
    dplyr::select(x, y)
  
  locs_pred <- as.matrix(pred_complete[, c("x", "y")])
  X_pred <- matrix(1, nrow = nrow(pred_complete), ncol = 1)
  
  # Fit Matern Gaussian process model
  fit <- GpGp::fit_model(
    y = y_train,
    locs = locs_train,
    X = X_train,
    covfun_name = "matern_isotropic"
  )
  
  # Generate predictions using the predictions() function
  pred_mean <- GpGp::predictions(
    fit = fit,
    locs_pred = locs_pred,
    X_pred = X_pred,
    m = 30
  )
  
  # MC Sim to get the predictive
  sims <- GpGp::cond_sim(
    fit = fit,
    locs_pred = locs_pred,
    X_pred = X_pred,
    m = 30,
    nsims = 500   # adjust (100–500 for most cases)
  )

  pred_sd <- apply(sims, 1, sd)
  
  # Create prediction dataframe
  pred_df <- data.frame(
    x = pred_complete$x,
    y = pred_complete$y,
    mu = as.vector(pred_mean),
    sd = as.vector(pred_sd)
  )
  
  # Add training observations (zero uncertainty)
  train_pred_df <- data.frame(
    x = train_complete$x,
    y = train_complete$y,
    mu = y_train,
    sd = 0
  )
  
  # Combine for full spatial grid
  full_pred_df <- bind_rows(pred_df, train_pred_df)
  
    return(list(
    fit = fit,
    predictions = full_pred_df,
    pred_only = pred_df
  ))
}

# Fit models for first 5 months
cat("\nFitting GpGp models for Jan-May...\n")
gpgp_fits <- map(month_names, ~ fit_gpgp_single(.x, train_df, temp_data))
names(gpgp_fits) <- month_names

# ============================================================================
# MODEL EVALUATION
# ============================================================================

cat("\nEvaluating model performance...\n")

# Extract predictions
preds_list <- map(gpgp_fits, ~ .x$predictions)

# Get true values at holdout locations
true_values_at_mask <- map(1:5, function(i) {
  temp_data[[i]][holdout_masks[[month_names[i]]]]
})

# Get predictions at holdout locations
predictions_at_mask <- map(1:5, function(i) {
  pred_raster <- predictions_to_raster(preds_list[[i]], temp_data)
  pred_raster[holdout_masks[[month_names[i]]]]
})

# Calculate RMSE and MAE
rmse_results <- map2_df(
  predictions_at_mask, 
  true_values_at_mask, 
  calculate_metrics,
  .id = "month"
)
rmse_results$month <- month_names

cat("\nPrediction Errors at Holdout Locations:\n")
print(rmse_results)

# Calculate CRPS (Continuous Ranked Probability Score)
crps_results <- map_df(1:5, function(i) {
  month_col <- month_names[i]
  pred_df <- preds_list[[i]]
  
  # Convert predictions to raster
  mu_raster <- predictions_to_raster(pred_df, temp_data)
  
  sd_raster <- rast(
    xmin = ext(temp_data)[1], xmax = ext(temp_data)[2],
    ymin = ext(temp_data)[3], ymax = ext(temp_data)[4],
    resolution = res(temp_data)
  )
  sd_raster <- rasterize(
    vect(pred_df, geom = c("x", "y")),
    sd_raster,
    field = "sd"
  )
  
  # Extract at holdout locations
  mu_masked <- mu_raster[holdout_masks[[month_col]]]
  sd_masked <- sd_raster[holdout_masks[[month_col]]]
  true_masked <- temp_data[[i]][holdout_masks[[month_col]]]
  
  # Remove NAs
  valid_idx <- !is.na(mu_masked) & !is.na(sd_masked) & !is.na(true_masked)
  mu_masked <- mu_masked[valid_idx]
  sd_masked <- sd_masked[valid_idx]
  true_masked <- true_masked[valid_idx]
  
  # Calculate CRPS assuming Gaussian predictive distribution
  crps_vals <- vapply(seq_along(mu_masked), function(idx) {
    mu <- mu_masked[idx]
    sd <- max(sd_masked[idx], 1e-10)  # Avoid division by zero
    obs <- true_masked[idx]
    verification::crps(obs, c(mu, sd))[[1]]
  }, numeric(1))
  
  tibble(
    month = month_col,
    CRPS = mean(crps_vals, na.rm = TRUE),
    n_holdout = length(crps_vals)
  )
})

cat("\nCRPS at Holdout Locations:\n")
print(crps_results)

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("\nGenerating visualizations...\n")

# Prepare GpGp predictions for plotting
gpgp_plot_data <- map_dfr(1:5, function(i) {
  df <- preds_list[[i]]
  df$month <- paste0("Month ", i)
  df
})
gpgp_plot_data$method <- "GpGp"

# Prepare actual data
actual_data <- as.data.frame(temp_data[[1:5]], xy = TRUE) %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "month_raw",
    values_to = "mu"
  ) %>%
  mutate(
    month = paste0("Month ", as.integer(gsub("tmax_", "", month_raw))),
    method = "Actual"
  ) %>%
  dplyr::select(x, y, mu, month, method)

# Combine data
all_plot_data <- bind_rows(gpgp_plot_data, actual_data)
all_plot_data$method <- factor(all_plot_data$method, levels = c("Actual", "GpGp"))

# Plot 1: Predictions vs Actual
p1 <- ggplot(all_plot_data, aes(x = x, y = y, fill = mu)) +
  geom_tile(width = 0.15, height = 0.15) +
  scale_fill_distiller(
    palette = "Spectral",
    name = "Max Temp (°C)"
  ) +
  facet_grid(method ~ month) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

# Prepare residual data
error_plot_data <- map_dfr(1:5, function(i) {
  month_col <- month_names[i]
  
  holdout_coords <- xyFromCell(temp_data[[i]], which(holdout_masks[[month_col]]))
  pred <- predictions_at_mask[[i]]
  true <- true_values_at_mask[[i]]
  
  valid_idx <- !is.na(pred) & !is.na(true)
  errors <- pred[valid_idx] - true[valid_idx]
  
  data.frame(
    x = holdout_coords[valid_idx, 1],
    y = holdout_coords[valid_idx, 2],
    error = errors,
    month = paste0("Month ", i),
    method = "GpGp"
  )
})

# Plot 2: Spatial residuals at holdout locations
p2 <- ggplot(error_plot_data, aes(x = x, y = y, fill = error)) +
  geom_raster() +
  facet_grid(method ~ month) +
  scale_fill_distiller(
    palette = "Spectral",
    name = "Residual (°C)",
    limits = c(-max(abs(error_plot_data$error)), max(abs(error_plot_data$error)))
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

print(p1)
print(p2)

# ============================================================================
# MODEL SUMMARIES
# ============================================================================

cat("\n=== Model Summaries ===\n")
for (i in 1:5) {
  cat("\n", month_names[i], ":\n")
  print(summary(gpgp_fits[[i]]$fit))
}

cat("\nWorkflow complete!\n")