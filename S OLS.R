rm(list = ls())
library(sf)
library(tidyr)
library(dplyr)
library(raster)
library(terra)
library(purrr)
library(FRK)
library(ggplot2)
library(verification)

source("./modules/preprocessing_functions.R")
source("./modules/ols_fit_functions.R")

# SELECTED MONTHS FOR MODELLING
SELECTED_MONTHS <- 1:5
set.seed(1111)

# Load and aggregate data
data_path <- "./data/agcd_v1_tmax_mean_r005_monthly_2017.nc"
temp_data <- rast(data_path)
temp_data <- aggregate(temp_data, fact = 3)

# Prepare coordinates and basis functions
coords_df <- as.data.frame(crds(temp_data), xy = TRUE)
G <- auto_basis(data = SpatialPoints(coords_df), nres = 2, type = "Gaussian")
S <- eval_basis(basis = G, s = as.matrix(coords_df)) %>% as.matrix()
colnames(S) <- paste0("B", 1:ncol(S))

# Create a test set for each month/layer (20% random sampling)
masked <- map(1:nlyr(temp_data), ~ mask_random_cells(temp_data[[.x]], 0.2))
train_data <- map(masked, "train") %>% rast()
holdout_masks <- map(masked, "mask")
names(train_data) <- names(temp_data)

# ============ TSR FIT ============
cat("Performing TSR Interpolation for first 5 months...\n")

# Apply OLS to all layers using the fit function
all_layers_results_list <- map(1:5, ~ fit_ols_single(.x, train_data, S))
ols_results <- map(all_layers_results_list, "raster") %>% rast()
names(ols_results) <- names(temp_data[[1:5]])

# EVALUATION
ols_metrics <- map_dfr(SELECTED_MONTHS, ~ calculate_ols_fit_metrics(.x, ols_results, temp_data, holdout_masks))
cat("TSR Mean RMSE:", round(mean(ols_metrics$rmse, na.rm = T), 3), "\n")

# SR plotting data
ols_plot_data <- map_dfr(SELECTED_MONTHS, ~ raster_to_df(ols_results, .x))
ols_plot_data$method <- "TSR"

# Original data for comparison
orig_plot_data <- map_dfr(SELECTED_MONTHS, ~ raster_to_df(temp_data, .x))
orig_plot_data$method <- "Actual"

# Combine all data
all_plot_data <- bind_rows(ols_plot_data, orig_plot_data)
all_plot_data$method <- factor(all_plot_data$method, levels = c("Actual", "TSR"))

# Create prediction plots
p1 <- ggplot(all_plot_data, aes(x = x, y = y, fill = tmax)) +
  geom_raster() +
  facet_grid(method ~ paste("Month", layer)) +
  scale_fill_distiller(palette = "Spectral", name = "Max Temp (°C)") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(
    x = "Longitude", y = "Latitude"
  )

# Calculate residuals for error mapping
ols_residuals <- map_dfr(SELECTED_MONTHS, ~ calculate_ols_lyr_resid(.x, ols_results, temp_data))

# CRPS for TSR, assuming Gaussian noises with respective parameters of predictions residuals and means
crps_gaussian <- map(SELECTED_MONTHS, ~ calculate_crps(.x, ols_metrics))

# Plot residuals
p2 <- ggplot(ols_residuals, aes(x = x, y = y, fill = residual)) +
  geom_raster() +
  facet_grid(method ~ paste("Month", layer)) +
  scale_fill_distiller(palette = "Spectral", name = "Residual (°C)") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(
    x = "Longitude", y = "Latitude"
  )

print(p1)
print(p2)

# Write to files
combined_df <- do.call(rbind, crps_gaussian)
write.csv(combined_df, "crps_gaussian_combined.csv", row.names = FALSE)

sink(file = "fitted_models_summary.txt")
print(
  map(SELECTED_MONTHS, ~ summary(all_layers_results_list[[.x]]$model))
)
sink(file = NULL)

