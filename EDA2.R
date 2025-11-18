library(raster)
library(ggplot2)
library(tidyr)
library(dplyr)
library(maps)

# Read NetCDF file as a raster stack
data_path <- "./data/agcd_v1_tmax_mean_r005_monthly_2017.nc"
nc_data <- stack(data_path)

# Aggregate the raster stack with fact=3
nc_data_agg <- aggregate(nc_data, fact=3, fun=mean, na.rm=TRUE)

# Extract layer names (they should be tmax_1, tmax_2, ..., tmax_12)
# If names are different, we'll rename them
layer_names <- names(nc_data_agg)
# Rename layers to proper month names for facet plot
month_names <- c("January", "February", "March", "April", "May", "June",
                 "July", "August", "September", "October", "November", "December")
names(nc_data_agg) <- month_names

# Convert aggregated raster stack to data frame for ggplot2
tmax_df <- as.data.frame(nc_data_agg, xy=TRUE)

# Reshape to long format for easier manipulation
tmax_long <- tmax_df %>%
  pivot_longer(cols = -c(x, y), 
               names_to = "time", 
               values_to = "Tmax") %>%
  drop_na()

# Convert time to factor with proper ordering
tmax_long$time <- factor(tmax_long$time, levels = month_names)

# 1. TIME SERIES PLOT
# Average over all locations for each time period
tmax_time_av <- tmax_long %>%
  group_by(time) %>%
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE),
            .groups = "drop")

time_series_plot <- ggplot(data = tmax_time_av,
                           aes(x = time, y = mean_tmax, group=1)) +
  geom_line() +
  geom_point(color = "darkred") +
  labs(x = "Time",
       y = "Average Max Temp (°C)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 2. HOVMOLLER PLOTS
# Hovmoller plot: time vs latitude
hovmoller_lat <- tmax_long %>%
  group_by(time, y) %>%
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE), .groups = 'drop')

plot_hovmoller_lat <- ggplot(hovmoller_lat, aes(x = time, y = y, fill = mean_tmax)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", name = "Mean Max Temp (°C)") + 
  labs(x = "Time (Month)",
       y = "Latitude (°)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Hovmoller plot: time vs longitude
hovmoller_lon <- tmax_long %>%
  group_by(time, x) %>%
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE), .groups = 'drop')

plot_hovmoller_lon <- ggplot(hovmoller_lon, aes(x = time, y = x, fill = mean_tmax)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", name = "Mean Max Temp (°C)") +
  labs(x = "Time (Month)",
       y = "Longitude (°)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


# 3. SPATIAL AVERAGED PLOT
# Average over all time periods for each location
tmax_spat_av <- tmax_long %>% 
  group_by(y, x) %>%
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE),
            .groups = "drop")

# Get Australia map outline
aus_map <- map_data("world", region = "Australia")

spatial_mean_heatmap <- ggplot() +
  geom_tile(data = tmax_spat_av, aes(x = x, y = y, fill = mean_tmax)) +
  geom_polygon(data = aus_map, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  scale_fill_distiller(palette = "Spectral", name = "Mean Max Temp (°C)") +
  labs(x = "Longitude (°)",
       y = "Latitude (°)") +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0), limits = range(tmax_spat_av$x)) +
  scale_y_continuous(expand = c(0, 0), limits = range(tmax_spat_av$y))


# 4. FACET GRID PLOT
# All months in one faceted plot
facet_plot <- ggplot(tmax_long) + 
  geom_tile(aes(x = x, y = y, fill = Tmax)) +
  facet_wrap(vars(time), ncol=4) +
  scale_fill_distiller(palette = "Spectral", name = "Max Temperature (°C)") +
  labs(x = "Longitude (°)", 
       y = "Latitude (°)") +
  theme_bw() +
  theme(panel.grid = element_blank())


# Display all plots
print(time_series_plot)
print(plot_hovmoller_lat)
print(plot_hovmoller_lon)
print(spatial_mean_heatmap)
print(facet_plot)