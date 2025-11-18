# Define function to fit and predict single month
fit_frk_single <- function(month_var, full_df,
                           created_BAUs, basis_functions) {

  single_month_data <- as.data.frame(full_df)[c("x", "y", month_var)]
  non_na_data <- single_month_data[!is.na(single_month_data[[month_var]]), ]
  coordinates(non_na_data) <- ~ x + y

  f <- as.formula(paste(month_var, "~ x + y + 1"))

  S <- FRK(
    f = f,
    data = list(non_na_data),
    basis = basis_functions,
    BAUs = created_BAUs,
    n_EM = 35,
    tol = 1e-12
  )

  return(predict(S))
}

calculate_frk_crps <- function(i, month_names, preds_list, temp_data, holdout_masks)
{
  month_col <- month_names[i]
  cat("Processing month:", month_col, "\n")
  
  # Prediction DF (x, y, mu, sd)
  pred_df <- preds_list[[i]]
  
  # Holdout mask and observed data
  holdout_mask <- holdout_masks[[month_col]]
  obs_raster <- temp_data[[i]]
  
  # Reference raster with same extent/resolution as observed raster
  ref_rast <- rast(
    xmin = ext(obs_raster)[1], xmax = ext(obs_raster)[2],
    ymin = ext(obs_raster)[3], ymax = ext(obs_raster)[4],
    resolution = res(obs_raster)
  )
  
  # Rasterise predicted mean and standard error
  pred_mu_rast <- rasterize(vect(pred_df, geom = c("x", "y")), ref_rast, field = "mu")
  pred_se_rast <- rasterize(vect(pred_df, geom = c("x", "y")), ref_rast, field = "sd")
  
  # Extract predictions at holdout (masked) cells
  pred_mu_vals <- pred_mu_rast[holdout_mask]
  pred_se_vals <- pred_se_rast[holdout_mask]
  
  # Extract true observations at those cells
  true_vals <- obs_raster[holdout_mask]
  
  # Remove any NAs (e.g. due to incomplete raster coverage)
  valid_idx <- !is.na(pred_mu_vals) & !is.na(pred_se_vals) & !is.na(true_vals)
  obs <- as.numeric(true_vals[valid_idx])
  pred_mean <- as.numeric(pred_mu_vals[valid_idx])
  pred_se <- as.numeric(pred_se_vals[valid_idx])
  
  # Compute CRPS under Gaussian assumption
  crps_scores <- verification::crps(obs, cbind(pred_mean, pred_se))
  
  return(tibble(
    month = month_col,
    mean_crps = mean(crps_scores$CRPS, na.rm = TRUE),
    median_crps = median(crps_scores$CRPS, na.rm = TRUE),
    n_holdout = length(obs)
    )
  )
}