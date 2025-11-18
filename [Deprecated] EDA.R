library(ncdf4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(maps)
library(sf)

data_path <- "./data/agcd_v1_tmax_mean_r005_monthly_2017.nc"
nc_temp_data <- nc_open(data_path)
print(nc_temp_data)


lat <- ncvar_get(nc_temp_data, "lat")
lon <- ncvar_get(nc_temp_data, "lon")
time <- ncvar_get(nc_temp_data, "time")
time <- as.Date(time, origin = "1850-01-01")

# Tmax is a 3-D array, ordering: Tmax[lon, lat, time]
Tmax <- ncvar_get(nc_temp_data,
                      "tmax",
                      count = c(-1, -1, -1))
nc_close(nc_temp_data)


# Convert the Tmax 3D array to a long data frame
tmax_df <- cbind(Tmax = c(Tmax),
            expand.grid(lon = lon,
                        lat = lat,
                        time = time))
tmax_df <- drop_na(tmax_df)


# Average over each time period for all locations, for time series analysis
tmax_time_av <- tmax_df %>%
  group_by(time) %>%
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE),
            .groups = "drop")

# Average over all observation in each location in the 4 periods
tmax_1 <- drop_na(subset(tmax_df, time %in% c("2017-01-16", "2017-04-16",
                                      "2017-07-16", "2017-10-16")))
tmax_spat_av <- tmax_1 %>% 
  group_by(lat, lon) %>% # group by lat-lon pair
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE),
            .groups = "drop")
rm(tmax_1)


# Create aggregate df for Hovmoller Plots
hovmoller_lat <- tmax_df %>%
  group_by(time, lat) %>%
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE), .groups = 'drop')

hovmoller_lon <- tmax_df %>%
  group_by(time, lon) %>%
  summarise(mean_tmax = mean(Tmax, na.rm = TRUE), .groups = 'drop')


# Create Hovmoller plots, both averaged by space dim
plot_hovmoller_lat <- ggplot(hovmoller_lat, aes(x = time, y = lat, fill = mean_tmax)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", name = "Mean Max Temp (°C)") + 
  labs(
    x = "Time (Month)",
    y = "Latitude (°)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
  )

plot_hovmoller_lon <- ggplot(hovmoller_lon, aes(x = time, y = lon, fill = mean_tmax)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", name = "Mean Max Temp (°C)") +
  labs(
    x = "Time (Month)",
    y = "Longitude (°)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
  )

# Create a facet grid with time as the grid variable
facet_plot <- ggplot(tmax_df) + 
  geom_tile(aes(x = lon, y = lat, fill = Tmax)) +
  facet_wrap(vars(time)) +
  scale_fill_distiller(palette = "Spectral") +
  labs(x = "Longitude", 
       y = "Latitude", 
       fill = "Max Temperature (°C)") +
  theme_bw()


# Create a line plot to show a time series
time_series_plot <- ggplot(data = tmax_time_av,
       aes(x = time, y = mean_tmax)) +
  geom_line() +
  geom_point(color = "darkred") +
  labs(x = "Time",
       y = "Average Maximum Temperature (°C)",
       title = "Averaged Maximum Temperature over Time") +
  theme_bw()


# Create an spatial average plot
aus_map <- map_data("world", region = "Australia")
spatial_mean_heatmap <- ggplot() +
  geom_tile(data = tmax_spat_av, aes(x = lon, y = lat, fill = mean_tmax)) +
  geom_polygon(data = aus_map, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  scale_fill_distiller(palette = "Spectral", name = "Mean Max Temp (°C)") +
  labs(
    x = "Longitude (°)",
    y = "Latitude (°)"
  ) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0), limits = range(tmax_spat_av$lon)) +
  scale_y_continuous(expand = c(0, 0), limits = range(tmax_spat_av$lat))

print(plot_hovmoller_lat)
print(plot_hovmoller_lon)
print(spatial_mean_heatmap)
print(time_series_plot)
print(facet_plot)
