# Filtering ROS to arctic land points only and exporting plots
# 14-01-24
# Maya Lemaire

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(raster)
library(sp)
library(ggspatial)
library(terra)
library(geosphere)
library(vegan)

# Read the shapefile for the Arctic
arctic_shape <- st_read("Arctic_Climate_Project/data/evidence-map-scope/evidence-map-scope.shp")

# Clean the geometry of the Arctic shapefile
arctic_shape <- st_make_valid(arctic_shape)
arctic_shape <- st_union(arctic_shape)

# Create a bounding box and subtract the arctic_shape to get the complement
bounding_box <- st_bbox(arctic_shape) %>%
  st_as_sfc() %>%
  st_sf(geometry = .)

arctic_shape <- st_difference(bounding_box, arctic_shape)

# Importing ocean shapefile 
ocean_polygons <- read_sf("Arctic_Climate_Project/data/ne_10m_ocean/ne_10m_ocean.shp")

coastline <- read_sf("Arctic_Climate_Project/data/ne_10m_coastline/ne_10m_coastline.shp")

# Define bounding box
bounding_box <- st_bbox(c(xmin = -180, xmax = 180, ymin = 50, ymax = 90), crs = st_crs(coastline))
# Crop coastline data using bounding box
coastline <- st_crop(coastline, bounding_box)
coastline <- st_union(coastline)

# Function to calculate average ROS events per year
yearly_gridded_count <- function(ROS_file) {
  
  # Extracting year
  filename_without_path <- basename(ROS_file)
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  year <- as.integer(filename_parts[length(filename_parts)])
  
  print(paste("processing", year, sep = " "))
    
  # Read the CSV file with lon-lat points
  ROS_data <- read.csv(ROS_file)
    
  # Convert lon values to -180 to 180 range
  ROS_data$lon <- ifelse(ROS_data$lon > 180, ROS_data$lon - 360, ROS_data$lon)
  
  ROS_data <- ROS_data[ROS_data$lat > 50, ]
    
  total_rf <- ROS_data %>%
    group_by(lat, lon) %>%
    summarise(total_rf = sum(rf))
    
    # Add year information to the result
  total_rf$year <- year
    
  return(total_rf)
  
}

# Specify the path to the folder containing CSV files
ROS_folder <- "Arctic_Climate_Project/output/ROS_3mm_1t2m_sd"

# Get a list of all CSV files in the folder
ROS_files <- list.files(path = ROS_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
ROS_data_list <- lapply(ROS_files, yearly_gridded_count)

# Combine the results into a single data frame
ROS_data <- bind_rows(ROS_data_list)

ROS_data$total_rf <- ROS_data$total_rf * 1000

## Filtering out land and out of arctic data for baseline
filter_fun <- function(global_data, ocean_polygons, arctic_shape) {
  
  coords_df <- global_data[, c("lon", "lat")]
  coords_df <- distinct(coords_df, lat, lon)
  
  coords_df <- st_as_sf(coords_df, coords=c("lon","lat"))
  
  st_crs(coords_df) <- st_crs(ocean_polygons)
  sf_use_s2(FALSE)
  
  ##find where out points intersect with the ocean
  tmp <- sapply(st_intersects(coords_df, ocean_polygons), function(z) if (length(z)==0) NA_integer_ else z[1])
  
  if (sum(!is.na(tmp))>0) {
    coords_df<-data.frame(st_coordinates(coords_df[is.na(tmp),]))} else {
      coords_df<-data.frame(st_coordinates(coords_df))}
  
  colnames(coords_df) <- c("lon","lat")
  
  coords_df <- st_as_sf(coords_df, coords=c("lon","lat"))
  
  st_crs(coords_df) <- st_crs(arctic_shape)
  sf_use_s2(FALSE)
  
  ##find where out points intersect with the arctic area
  tmp <- sapply(st_intersects(coords_df, arctic_shape), function(z) if (length(z)==0) NA_integer_ else z[1])
  
  if (sum(!is.na(tmp))>0) {
    coords_df <- data.frame(st_coordinates(coords_df[is.na(tmp),]))} else {
      coords_df <- data.frame(st_coordinates(coords_df))}
  
  colnames(coords_df) <- c("lon","lat")
  
  filtered_data <- global_data %>%
    inner_join(coords_df, by = c("lat", "lon"))  
  
  return(filtered_data)

}

filtered_ROS <- filter_fun(ROS_data, ocean_polygons, arctic_shape)
rm(ROS_data)

filtered_baseline <- filtered_ROS %>%
  subset(year > 1950 & year <= 1980) %>%
  group_by(lon, lat) %>%
  summarise(ROS_baseline = sum(total_rf)/30)

filtered_current <- filtered_ROS %>%
  subset(year > 1990 & year <= 2020) %>%
  group_by(lon, lat) %>%
  summarise(ROS_current = sum(total_rf)/30)

ROS_absolute_change <- full_join(filtered_baseline, filtered_current, by = c("lon", "lat"))
rm(filtered_baseline, filtered_current)

ROS_absolute_change[is.na(ROS_absolute_change)] <- 0

ROS_absolute_change <- ROS_absolute_change %>%
  mutate(Absolute_diff = ROS_current - ROS_baseline)

bounding_box_alaska <- st_bbox(c(xmin = -168, xmax = -141, ymin = 66, ymax = 72), crs = st_crs(coastline))
alaska_coastline <- st_crop(coastline, bounding_box_alaska)

alaska_change <- subset(ROS_absolute_change, 
                        lat >= 66 & lat <= 72 & 
                          lon >= -168 & lon <= -141)

bounding_box_northern_fennoscandia <- st_bbox(c(xmin = 16, xmax = 32, ymin = 67, ymax = 71), crs = st_crs(coastline))
northern_fennoscandia_coastline <- st_crop(coastline, bounding_box_northern_fennoscandia)

northern_fennoscandia_change <- subset(ROS_absolute_change, 
                         lat >= 67 & lat <= 71 & 
                           lon >= 16 & lon <= 32)

bounding_box_svalbard <- st_bbox(c(xmin = 10, xmax = 35, ymin = 76, ymax = 82), crs = st_crs(coastline))
svalbard_coastline <- st_crop(coastline, bounding_box_svalbard)

svalbard_change <- subset(ROS_absolute_change, 
                          lat >= 76 & lat <= 82 & 
                            lon >= 10 & lon <= 35)

bounding_box_yamal <- st_bbox(c(xmin = 64, xmax = 82, ymin = 64, ymax = 74), crs = st_crs(coastline))
yamal_coastline <- st_crop(coastline, bounding_box_yamal)

yamal_change <- subset(ROS_absolute_change, 
                       lat >= 64 & lat <= 73 & 
                         lon >= 64 & lon <= 82)

df_sf <- st_as_sf(ROS_absolute_change, coords = c("lon", "lat"), crs = 4326)
raster_template <- rast(ext(df_sf), resolution = 1)

baseline_raster <- rasterize(df_sf, raster_template, field = "ROS_baseline", fun = mean)
current_raster <- rasterize(df_sf, raster_template, field = "ROS_current", fun = mean)
change_raster <- rasterize(df_sf, raster_template, field = "Absolute_diff", fun = mean)

# Stack the rasters
multi_layer_raster <- c(baseline_raster, current_raster, change_raster)

# Set layer names
names(multi_layer_raster) <- c("baseline", "current", "change")

# Define the polar stereographic projection
polar_proj <- "+proj=stere +lat_0=90 +lon_0=-45 +lat_ts=70 +datum=WGS84"

# Project the raster
multi_layer_raster_polar <- project(multi_layer_raster, polar_proj)

# Extract each layer
baseline_layer <- multi_layer_raster_polar[[1]]
current_layer <- multi_layer_raster_polar[[2]]
change_layer <- multi_layer_raster_polar[[3]]

coastline_mapping <- read_sf("Arctic_Climate_Project/data/ne_10m_coastline/ne_10m_coastline.shp")

graticule <- read_sf("Arctic_Climate_Project/data/ne_10d_graticule/ne_10m_graticules_10.shp")
graticule <- graticule[, "degrees"]

extent <- ext(baseline_layer)
extent <- as.vector(extent)

# Create a bounding box polygon covering the entire area of interest
bounding_polygon <- st_polygon(list(rbind(c(extent[1], extent[3]), 
                                          c(extent[1], extent[4]), 
                                          c(extent[2], extent[4]), 
                                          c(extent[2], extent[3]), 
                                          c(extent[1], extent[3]))))

arctic_shape <- st_transform(arctic_shape, crs = st_crs(baseline_layer))

# Convert the bounding polygon to the same CRS as arctic_shape
bounding_polygon <- st_sfc(bounding_polygon, crs = st_crs(baseline_layer))

# Plot inverse of arctic_mapping in grey
arctic_inverse <- st_difference(bounding_polygon, arctic_shape)

coastline_mapping <- st_as_sfc(coastline_mapping, crs = st_crs(baseline_layer))
coastline_mapping <- st_transform(coastline_mapping, st_crs(bounding_polygon))
coastline_mapping <- st_intersection(coastline_mapping, bounding_polygon)

graticule <- st_as_sfc(graticule, crs = st_crs(baseline_layer))
graticule <- st_transform(graticule, st_crs(bounding_polygon))
graticule <- st_intersection(graticule, bounding_polygon)

# Define color palettes and breaks for baseline and current layers
breaks <- c(0, 10, 20, 30, 40, 50, Inf)
color_palette <- brewer.pal(9, "PuBu")
continuous_palette <- colorRampPalette(color_palette)

# Plot baseline_layer with graticule and other layers
plot(baseline_layer, col = continuous_palette(100), breaks = breaks, 
     axis.args = list(at = breaks), legend = TRUE)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)

# Plot current_layer with graticule and other layers
plot(current_layer, col = continuous_palette(100), breaks = breaks, 
     axis.args = list(at = breaks), legend = TRUE)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)

# Define color palette and breaks
color_palette <- brewer.pal(9, "RdBu")
continuous_palette <- colorRampPalette(color_palette)

# Reverse the color palette
continuous_palette_reversed <- function(n) rev(continuous_palette(n))

custom_palette <- c("#5698c5", "#a9cde3", "#e0f2fe", "#feddcc", "#e6795a", "#b2172c")
breaks <- c(-Inf, -10, -5, 0, 5, 10, Inf)

# Project baseline_zero_points to match the projection of change_layer
baseline_zero_points <- ROS_absolute_change[ROS_absolute_change$ROS_baseline == 0, ]
baseline_zero_points_sf <- st_as_sf(baseline_zero_points, coords = c("lon", "lat"), crs = 4326)
baseline_zero_points_projected <- st_transform(baseline_zero_points_sf, crs = crs(change_layer))

current_zero_points <- ROS_absolute_change[ROS_absolute_change$ROS_current == 0, ]
current_zero_points_sf <- st_as_sf(current_zero_points, coords = c("lon", "lat"), crs = 4326)
current_zero_points_projected <- st_transform(current_zero_points_sf, crs = crs(change_layer))

# Define a bounding box polygon for Alaska
alaska_bbox <- st_bbox(c(xmin = -168, xmax = -141, ymin = 66, ymax = 72), crs = 4326)
northern_fennoscandia_bbox <- st_bbox(c(xmin = 16, xmax = 32, ymin = 67, ymax = 71), crs = 4326)
svalbard_bbox <- st_bbox(c(xmin = 10, xmax = 35, ymin = 76, ymax = 82), crs = 4326)
yamal_bbox <- st_bbox(c(xmin = 64, xmax = 82, ymin = 64, ymax = 74), crs = 4326)

# Transform the bounding box to match the projection of change_layer
alaska_bbox <- st_as_sfc(alaska_bbox)
alaska_bbox <- st_transform(alaska_bbox, st_crs(change_layer))

northern_fennoscandia_bbox <- st_as_sfc(northern_fennoscandia_bbox)
northern_fennoscandia_bbox <- st_transform(northern_fennoscandia_bbox, st_crs(change_layer))

svalbard_bbox <- st_as_sfc(svalbard_bbox)
svalbard_bbox <- st_transform(svalbard_bbox, st_crs(change_layer))

yamal_bbox <- st_as_sfc(yamal_bbox)
yamal_bbox <- st_transform(yamal_bbox, st_crs(change_layer))

# Plot change_layer with reversed color scale and customized legend labels
plot(change_layer, col = custom_palette, breaks = breaks)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)
points(baseline_zero_points_projected, col = "darkred", pch = "*", cex = 1.5)
points(current_zero_points_projected, col = "darkblue", pch = "*", cex = 1.5)
plot(alaska_bbox, col = alpha("grey", 0.3), border = "black", lwd = 2, lty = "dashed", add = TRUE)
plot(northern_fennoscandia_bbox, col = alpha("grey", 0.3), border = "black", lwd = 2, lty = "dashed", add = TRUE)
plot(svalbard_bbox, col = alpha("grey", 0.3), border = "black", lwd = 2, lty = "dashed", add = TRUE)
plot(yamal_bbox, col = alpha("grey", 0.3), border = "black", lwd = 2, lty = "dashed", add = TRUE)

# ggplot with reversed color palette
gg_change_alaska <- ggplot() +
  geom_tile(data = subset(alaska_change, Absolute_diff <= -10), 
            aes(x = lon, y = lat), 
            fill = "#5698c5", width = 1, height = 1) +
  geom_tile(data = subset(alaska_change, Absolute_diff > -10 & Absolute_diff <= -5), 
            aes(x = lon, y = lat), 
            fill = "#a9cde3", width = 1, height = 1) +
  geom_tile(data = subset(alaska_change, Absolute_diff > -5 & Absolute_diff <= -0), 
            aes(x = lon, y = lat), 
            fill = "#e0f2fe", width = 1, height = 1) +
  geom_tile(data = subset(alaska_change, Absolute_diff > 0 & Absolute_diff <= 5), 
            aes(x = lon, y = lat), 
            fill = "#feddcc", width = 1, height = 1) +
  geom_tile(data = subset(alaska_change, Absolute_diff > 5 & Absolute_diff <= 10), 
            aes(x = lon, y = lat), 
            fill = "#e6795a", width = 1, height = 1) +
  geom_tile(data = subset(alaska_change, Absolute_diff > 10), 
            aes(x = lon, y = lat), 
            fill = "#b2172c", width = 1, height = 1) +
  geom_point(data = subset(alaska_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  geom_point(data = subset(alaska_change, ROS_current == 0), 
             aes(x = lon, y = lat), color = "darkblue", size = 15, shape = "*") +
  geom_sf(data = alaska_coastline, fill = "transparent", color = "black") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )

plot(gg_change_alaska)

gg_change_northern_fennoscandia <- ggplot() +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff <= -10), 
            aes(x = lon, y = lat), 
            fill = "#5698c5", width = 1, height = 1) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff > -10 & Absolute_diff <= -5), 
            aes(x = lon, y = lat), 
            fill = "#a9cde3", width = 1, height = 1) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff > -5 & Absolute_diff <= -0), 
            aes(x = lon, y = lat), 
            fill = "#e0f2fe", width = 1, height = 1) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff > 0 & Absolute_diff <= 5), 
            aes(x = lon, y = lat), 
            fill = "#feddcc", width = 1, height = 1) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff > 5 & Absolute_diff <= 10), 
            aes(x = lon, y = lat), 
            fill = "#e6795a", width = 1, height = 1) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff > 10), 
            aes(x = lon, y = lat), 
            fill = "#b2172c", width = 1, height = 1) +
  geom_point(data = subset(northern_fennoscandia_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  geom_point(data = subset(northern_fennoscandia_change, ROS_current == 0), 
             aes(x = lon, y = lat), color = "darkblue", size = 15, shape = "*") +
  geom_sf(data = northern_fennoscandia_coastline, fill = "transparent", color = "black") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Add dashed grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_blank(),
    axis.text.x = element_blank(),  # Rotate x-axis labels
    axis.text.y = element_blank(),  # Adjust y-axis labels alignment
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "none"  # Remove the legend
  ) 

plot(gg_change_northern_fennoscandia)

gg_change_svalbard <- ggplot() +
  geom_tile(data = subset(svalbard_change, Absolute_diff <= -10), 
            aes(x = lon, y = lat), 
            fill = "#5698c5", width = 1, height = 1) +
  geom_tile(data = subset(svalbard_change, Absolute_diff > -10 & Absolute_diff <= -5), 
            aes(x = lon, y = lat), 
            fill = "#a9cde3", width = 1, height = 1) +
  geom_tile(data = subset(svalbard_change, Absolute_diff > -5 & Absolute_diff <= -0), 
            aes(x = lon, y = lat), 
            fill = "#e0f2fe", width = 1, height = 1) +
  geom_tile(data = subset(svalbard_change, Absolute_diff > 0 & Absolute_diff <= 5), 
            aes(x = lon, y = lat), 
            fill = "#feddcc", width = 1, height = 1) +
  geom_tile(data = subset(svalbard_change, Absolute_diff > 5 & Absolute_diff <= 10), 
            aes(x = lon, y = lat), 
            fill = "#e6795a", width = 1, height = 1) +
  geom_tile(data = subset(svalbard_change, Absolute_diff > 10), 
            aes(x = lon, y = lat), 
            fill = "#b2172c", width = 1, height = 1) +
  geom_point(data = subset(svalbard_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  geom_point(data = subset(svalbard_change, ROS_current == 0), 
             aes(x = lon, y = lat), color = "darkblue", size = 15, shape = "*") +
  geom_sf(data = svalbard_coastline, fill = "transparent", color = "black") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Add dashed grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_blank(),
    axis.text.x = element_blank(),  # Rotate x-axis labels
    axis.text.y = element_blank(),  # Adjust y-axis labels alignment
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "none"  # Remove the legend
  )

plot(gg_change_svalbard)

gg_change_yamal <- ggplot() +
  geom_tile(data = subset(yamal_change, Absolute_diff <= -10), 
            aes(x = lon, y = lat), 
            fill = "#5698c5", width = 1, height = 1) +
  geom_tile(data = subset(yamal_change, Absolute_diff > -10 & Absolute_diff <= -5), 
            aes(x = lon, y = lat), 
            fill = "#a9cde3", width = 1, height = 1) +
  geom_tile(data = subset(yamal_change, Absolute_diff > -5 & Absolute_diff <= -0), 
            aes(x = lon, y = lat), 
            fill = "#e0f2fe", width = 1, height = 1) +
  geom_tile(data = subset(yamal_change, Absolute_diff > 0 & Absolute_diff <= 5), 
            aes(x = lon, y = lat), 
            fill = "#feddcc", width = 1, height = 1) +
  geom_tile(data = subset(yamal_change, Absolute_diff > 5 & Absolute_diff <= 10), 
            aes(x = lon, y = lat), 
            fill = "#e6795a", width = 1, height = 1) +
  geom_tile(data = subset(yamal_change, Absolute_diff > 10), 
            aes(x = lon, y = lat), 
            fill = "#b2172c", width = 1, height = 1) +
  geom_point(data = subset(yamal_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  geom_point(data = subset(yamal_change, ROS_current == 0), 
             aes(x = lon, y = lat), color = "darkblue", size = 15, shape = "*") +
  geom_sf(data = yamal_coastline, fill = "transparent", color = "black") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Add dashed grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_blank(),
    axis.text.x = element_blank(),  # Rotate x-axis labels
    axis.text.y = element_blank(),  # Adjust y-axis labels alignment
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "none"  # Remove the legend
  )

gg_change_yamal

## MANTEL TEST

coords <- ROS_absolute_change[c("lon", "lat")]
spatial_dist <- distm(coords)

# Compute distance matrix for ROS and time
ROS_absolute_dist <- dist(ROS_absolute_change$Absolute_diff)

mantel_test_result <- mantel(ROS_absolute_dist, spatial_dist, method = "pearson", permutations = 99)

print(mantel_test_result)

densities <- densCols(space_dist_vector, ROS_absolute_dist, colramp = colorRampPalette(c("black", "white")))

# Extract densities
dens_values <- col2rgb(densities)[1,] + 1L

# Create a dataframe with spatial distance, exceedance distance, and densities
df <- data.frame(space_dist = space_dist_vector,
                 ROS_absolute_dist = ROS_absolute_dist,
                 dens = dens_values)

# Plot it
mantel_plot <- ggplot(df, aes(x = space_dist, y = ROS_absolute_dist, color = dens)) +
  geom_point(size = 2) +
  scale_color_viridis(option = "inferno") +  # Using viridis color palette
  labs(
    x = "Spatial Distance",
    y = "Rain-On-Snow Absolute Change Distance"
    ) +
  theme_minimal()

plot(mantel_plot)

## PAIRED T-TEST

result <- t.test(ROS_absolute_change$ROS_current, ROS_absolute_change$ROS_baseline, paired = TRUE)
result_alaska <- t.test(alaska_change$ROS_current, alaska_change$ROS_baseline, paired = TRUE)
result_northern_fennoscandia <- t.test(northern_fennoscandia_change$ROS_current, northern_fennoscandia_change$ROS_baseline, paired = TRUE)
result_svalbard <- t.test(svalbard_change$ROS_current, svalbard_change$ROS_baseline, paired = TRUE)
result_yamal <- t.test(yamal_change$ROS_current, yamal_change$ROS_baseline, paired = TRUE)

print(result)
print(result_alaska)
print(result_northern_fennoscandia)
print(result_svalbard)
print(result_yamal)
