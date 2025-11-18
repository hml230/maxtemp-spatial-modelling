rm(list = ls())
library(sf)
library(tidyr)
library(dplyr)
library(raster)
library(terra)
library(purrr)
library(ggplot2)

source("./modules/preprocessing_functions.R")
source("./modules/ols_fit_functions.R")

# SELECTED MONTHS FOR MODELLING
SELECTED_MONTHS <- 1:5
set.seed(1111)

# Load and aggregate data
data_path <- "./data/agcd_v1_tmax_mean_r005_monthly_2017.nc"
temp_data <- rast(data_path)
temp_data <- aggregate(temp_data, fact = 3)

# Create a test set for each month/layer (20% random sampling)
masked <- map(1:nlyr(temp_data), ~ mask_random_cells(temp_data[[.x]], 0.2))
train_data <- map(masked, "train") %>% rast()
holdout_masks <- map(masked, "mask")
names(train_data) <- names(temp_data)

# Fill NA values with extreme value (-999) in the SpatRaster object
train_data_filled <- train_data
train_data_filled[is.na(train_data_filled)] <- 999

# Select a single month for visualization (e.g., January)
month_to_plot <- 1

# Convert SpatRaster layers to data.frames for plotting
original_df <- as.data.frame(temp_data[[month_to_plot]], xy = TRUE)
colnames(original_df)[3] <- "value"
original_df$type <- "Original"

masked_df <- as.data.frame(train_data_filled[[month_to_plot]], xy = TRUE)
colnames(masked_df)[3] <- "value"
masked_df$type <- "Simulated MAR (filled with 999)"

# Combine for faceted plot
combined_df <- rbind(original_df, masked_df)
combined_df$type <- factor(combined_df$type,
                           levels = c("Original", "Simulated MAR (filled with 999)"))

# Get the actual temperature range (excluding fill values)
temp_range <- range(values(temp_data[[month_to_plot]]), na.rm = TRUE)
# Create side-by-side plot
p <- ggplot(combined_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  facet_wrap(~type, ncol = 2) +
  scale_fill_distiller(palette = "Spectral", name = "Max Temp (°C)",
                       limits = temp_range,
                       oob = scales::squish) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

print(p)