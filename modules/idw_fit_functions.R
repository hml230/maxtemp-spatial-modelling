# ============ Fit Functions ============
fit_idw_single <- function(lyr_idx, train_data)
{
  cat("Fitting IDW for layer ", lyr_idx, "\n")
  train_layer <- train_data[[lyr_idx]]
  
  # Convert raster to points
  train_points <- as.points(train_layer, values = TRUE, na.rm = TRUE)
  train_sf <- st_as_sf(train_points)
  
  # Create prediction grid from all cells
  pred_sf <- st_as_sf(as.points(train_layer, na.rm = FALSE))

  idw_result <- gstat::idw(formula = as.formula(paste(names(train_layer), "~ 1")),
                           locations = train_sf,
                           newdata = pred_sf,
                           idp = 2)

  result_raster <- rasterize(vect(idw_result), train_layer, field = "var1.pred")
  return(result_raster)
}

# ============ Metrics Calculation Functions ============

calculate_idw_fit_metrics <- function(lyr_idx, idw_results, temp_data, holdout_masks)
{
  mask <- holdout_masks[[lyr_idx]]
  pred_vals <- values(idw_results[[lyr_idx]])[mask]
  obs_vals <- values(temp_data[[lyr_idx]])[mask]

  complete_idx <- complete.cases(pred_vals, obs_vals)
  pred_vals <- pred_vals[complete_idx]
  obs_vals <- obs_vals[complete_idx]

  return(
    data.frame(
    layer = lyr_idx,
    rmse = sqrt(mean((pred_vals - obs_vals)^2)),
    mae = mean(abs(pred_vals - obs_vals))
    )
  )
}

calculate_idw_lyr_resid <- function(lyr_idx, idw_results, temp_data)
{
  idw_vals <- values(idw_results[[lyr_idx]])
  orig_vals <- values(temp_data[[lyr_idx]])
  residuals <- idw_vals - orig_vals

  # Include ALL cells (na.rm = FALSE)
  df <- as.data.frame(temp_data[[lyr_idx]], xy = TRUE, na.rm = FALSE)
  df$residual <- residuals
  df$layer <- lyr_idx
  df$method <- "IDW"

  return(df)
}
