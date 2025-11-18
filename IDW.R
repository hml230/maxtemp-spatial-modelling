rm(list = ls())
library(sf)
library(tidyr)
library(dplyr)
library(terra)
library(purrr)
library(gstat)
library(ggplot2)

source("./modules/preprocessing_functions.R")
source("./modules/idw_fit_functions.R")

# SELECTED MONTHS FOR MODELLING
SELECTED_MONTHS <- 1:5
set.seed(1111)

# Load and aggregate data
data_path <- "./data/agcd_v1_tmax_mean_r005_monthly_2017.nc"
temp_data <- rast(data_path)
temp_data <- aggregate(temp_data, fact = 3)

# Get coordinates
coords_df <- as.data.frame(crds(temp_data), xy = TRUE)

# Create a test set for each month/layer (20% random sampling)
masked <- map(1:nlyr(temp_data), ~ mask_random_cells(temp_data[[.x]], 0.2))
train_data <- map(masked, "train") %>% rast()
holdout_masks <- map(masked, "mask")
names(train_data) <- names(temp_data)

# ============ IDW FIT ============
cat("Performing IDW Interpolation...\n")

# Apply IDW to all layers
idw_results <- map(1:5, ~fit_idw_single(.x, train_data))
names(idw_results) <- names(train_data[[1:5]])

# IDW evaluation
idw_metrics <- map_dfr(1:5, ~calculate_idw_fit_metrics(.x, idw_results, temp_data, holdout_masks))

cat("IDW Mean RMSE:", mean(idw_metrics$rmse, na.rm = T), "\n")

# ============ IDW PLOT ============

# IDW plotting data
idw_plot_data <- map_dfr(SELECTED_MONTHS, ~ raster_to_df(idw_results, .x))
idw_plot_data$method <- "IDW"

# Original data for comparison
orig_plot_data <- map_dfr(SELECTED_MONTHS, ~ raster_to_df(temp_data, .x))
orig_plot_data$method <- "Actual"

# Combine all data
all_plot_data <- bind_rows(idw_plot_data, orig_plot_data)
all_plot_data$method <- factor(all_plot_data$method, levels = c("Actual", "IDW"))

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
idw_residuals <- map_dfr(SELECTED_MONTHS, ~calculate_idw_lyr_resid(.x, idw_results, temp_data))

# Plot residuals (proxy for standard errors)
p2 <- ggplot(idw_residuals, aes(x = x, y = y, fill = residual)) +
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
