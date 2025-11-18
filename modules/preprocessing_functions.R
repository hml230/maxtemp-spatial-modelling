# Function to randomly set proportion of cells to NA in one layer
mask_random_cells <- function(layer, proportion = 0.2) {
  vals <- values(layer) # get Tmax values
  valid_idx <- which(!is.na(vals)) # filter NA values

  sample_size <- round(length(valid_idx) * proportion)
  drop_idx <- sample(valid_idx, sample_size) # sample() returns index vector

  test_mask <- rep(FALSE, ncell(layer))
  test_mask[drop_idx] <- TRUE

  vals[drop_idx] <- NA
  values(layer) <- vals

  return(list(train = layer, mask = test_mask))
}


# Convert raster layers to data frames for ggplot
raster_to_df <- function(raster_stack, layer_num) {
  r <- raster_stack[[layer_num]]
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "tmax"
  df$layer <- layer_num

  return(df)
}

# Convert predictions dataframe to raster
predictions_to_raster <- function(pred_df, template_raster) {
  pred_raster <- rast(
    xmin = ext(template_raster)[1], 
    xmax = ext(template_raster)[2],
    ymin = ext(template_raster)[3], 
    ymax = ext(template_raster)[4],
    resolution = res(template_raster)
  )
  
  rasterize(
    vect(pred_df, geom = c("x", "y")),
    pred_raster,
    field = "mu"
  )
}

# Calculate prediction metrics
calculate_metrics <- function(predictions, true_values) {
  valid_idx <- !is.na(predictions) & !is.na(true_values)
  pred <- predictions[valid_idx]
  true <- true_values[valid_idx]
  
  tibble(
    RMSE = sqrt(mean((pred - true)^2)),
    MAE = mean(abs(pred - true)),
    n_holdout = length(true)
  )
}