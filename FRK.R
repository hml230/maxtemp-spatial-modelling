rm(list = ls())
library(tibble)
library(sf)
library(tidyr)
library(dplyr)
library(raster)
library(terra)
library(purrr)
library(FRK)
library(ggplot2)

source("./modules/preprocessing_functions.R")
source("./modules/frk_fit_functions.R")
set.seed(1111)

# Load and aggregate data
data_path <- "./data/agcd_v1_tmax_mean_r005_monthly_2017.nc"
temp_data <- rast(data_path)
temp_data <- aggregate(temp_data, fact = 3)

# Create a test set for each month/layer (20% random sampling)
masked <- map(1:nlyr(temp_data), ~ mask_random_cells(temp_data[[.x]], 0.2))

# Apply extract the masked data and convert to raster
train_layers <- map(masked, "train")
train_raster <- rast(train_layers)
names(train_raster) <- names(temp_data)

# Convert to a DF then a SPDF
train_df <- as.data.frame(train_raster, xy = TRUE)
train_spdf <- train_df
coordinates(train_spdf) <- ~ x + y

# Extract holdout masks (for RMSE calculation)
holdout_masks <- map(masked, "mask")
names(holdout_masks) <- names(temp_data)

# Create BAUs and spatial basis functions 
BAUs <- auto_BAUs(
  manifold = plane(),
  data = train_spdf,
  type = "grid",
  cellsize = c(0.15, 0.15),
  nonconvex_hull = 0
)
BAUs$fs <- 1

G_spatial <- auto_basis(
  manifold = plane(),
  data = train_spdf,
  nres = 2,
  type = "bisquare",
  regular = 0
)

month_names <- names(holdout_masks)
month_names <- month_names[1:5]

# Apply the fit to 1 to 5th layers (Jan-May)
cat("Performing FRK Interpolation for first 5 months...\n")
grid_BAUs <- map(month_names, ~ fit_frk_single(
  .x, train_spdf,
  BAUs, G_spatial
))
preds_list <- purrr::map(grid_BAUs, as.data.frame)

# Calculate RMSE for masked locations
# Extract true values at masked locations
true_values_at_mask <- map(1:5, function(i) {
  month_col <- month_names[i]
  return(temp_data[[i]][holdout_masks[[month_col]]])
})

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
  
  # Rasterise
  pred_raster <- rasterize(
    vect(pred_df, geom = c("x", "y")),
    pred_raster,
    field = "mu"
  )
  
  # Extract using the boolean mask
  return(pred_raster[holdout_masks[[month_col]]])
})

# Calculate RMSE and MAE for masked locations only
rmse_results <- map2_df(predictions_at_mask, true_values_at_mask, function(pred, true) {
  # Ensure same length
  valid_idx <- !is.na(pred) & !is.na(true)
  pred <- pred[valid_idx]
  true <- true[valid_idx]
  
  return(tibble(
    RMSE = sqrt(mean((pred - true)^2)),
    MAE = mean(abs(pred - true)),
    n_holdout = length(true)
    )
  )
}, .id = "month")

rmse_results$month <- month_names
cat("\nFRK Prediction Errors on unsampled locations:\n")
print(rmse_results)

# Calculate CRPS for masked locations
crps_results <- map_df(1:5, ~calculate_frk_crps(.x, month_names, preds_list, temp_data, holdout_masks))
  
# Convert to a combined dataframe for plots (predictions on full grid)
frk_plot_data <- purrr::map_dfr(1:5, ~ {
  df <- preds_list[[.x]]
  df$month <- paste0("tmax_", .x)
  return(df)
})
frk_plot_data$method <- "FRK"

actual_data <- as.data.frame(temp_data[[1:5]], xy = TRUE) %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "month",
    values_to = "mu"
  ) %>%
  mutate(method = "Actual")

# Relabel months for plotting
frk_plot_data <- frk_plot_data %>%
  mutate(month = factor(month, 
                        levels = paste0("tmax_", 1:5),
                        labels = paste0("Month ", 1:5)))

actual_data <- actual_data %>%
  mutate(month = factor(month,
                        levels = paste0("tmax_", 1:5),
                        labels = paste0("Month ", 1:5)))

# Assemble the full data
all_plot_data <- bind_rows(frk_plot_data, actual_data)
all_plot_data$method <- factor(all_plot_data$method, levels = c("Actual", "FRK"))

# Plot
p1 <- ggplot(all_plot_data) +
  geom_tile(aes(x = x, y = y, fill = mu), width = 0.15, height = 0.15) +
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

# Extract errors at holdout locations
error_plot_data <- map_dfr(1:5, function(i) {
  month_col <- month_names[i]
  
  # Get coordinates of holdout cells
  holdout_coords <- raster::xyFromCell(temp_data[[i]], which(holdout_masks[[month_col]]))
  
  # Get prediction errors
  pred <- predictions_at_mask[[i]]
  true <- true_values_at_mask[[i]]
  
  # Calculate errors (remove NAs)
  valid_idx <- !is.na(pred) & !is.na(true)
  errors <- pred[valid_idx] - true[valid_idx]
  
  # Combine into dataframe
  data.frame(
    x = holdout_coords[valid_idx, 1],
    y = holdout_coords[valid_idx, 2],
    error = errors,
    month = paste0("Month ", i),
    method = "FRK"
  )
})

# Plot spatial distribution of errors
p2 <- ggplot(error_plot_data,aes(x = x, y = y, fill = error)) +
  geom_raster() +
  facet_grid(method ~ month) +
  scale_fill_distiller(
    palette = "Spectral",
    name = "Residual (°C)"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

# print(p1)
# print(p2)