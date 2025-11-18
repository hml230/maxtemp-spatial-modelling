rm(list = ls())
rm(list = ls())
library(tibble)
library(sf)
library(raster)
library(terra)
library(purrr)
library(FRK)


source("./modules/preprocessing_functions.R")
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


G_spatial <- auto_basis(
  manifold = plane(),
  data = train_spdf,
  nres = 2,
  type = "Gaussian",
  regular = 0
)


FRK::show_basis(G_spatial)