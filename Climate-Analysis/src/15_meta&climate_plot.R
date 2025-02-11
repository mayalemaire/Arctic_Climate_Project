# Filtering ROS to arctic land points only and exporting plots
# 22-01-25
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

## Meta analysis results ##
# Only including overall effects and phyla effects with statistical significance
ROS_effects <- -1.476
WW_effects <- -0.655
chordata_effects <- -1.400
angiosperm_effects <- -0.987
study_locations <- read.csv("Arctic_Climate_Project/output/study_locations.csv")
# Filter WW and ROS studies
WW_studies <- study_locations %>% filter(Extreme_event == "WW")
ROS_studies <- study_locations %>% filter(Extreme_event == "ROS")

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
yearly_gridded_count_ROS <- function(ROS_file) {
  
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
ROS_folder <- "Arctic_Climate_Project/output/ROS_3mm_1t2m_sd_m3d"

# Get a list of all CSV files in the folder
ROS_files <- list.files(path = ROS_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
ROS_data_list <- lapply(ROS_files, yearly_gridded_count_ROS)

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

filtered_baseline_ROS <- filtered_ROS %>%
  subset(year > 1950 & year <= 1980) %>%
  group_by(lon, lat) %>%
  summarise(ROS_baseline = sum(total_rf)/30)

filtered_current_ROS <- filtered_ROS %>%
  subset(year > 1990 & year <= 2020) %>%
  group_by(lon, lat) %>%
  summarise(ROS_current = sum(total_rf)/30)

ROS_absolute_change <- full_join(filtered_baseline_ROS, filtered_current_ROS, by = c("lon", "lat"))
rm(filtered_baseline_ROS, filtered_current_ROS)

ROS_absolute_change[is.na(ROS_absolute_change)] <- 0

ROS_absolute_change <- ROS_absolute_change %>%
  mutate(Absolute_diff = ROS_current - ROS_baseline)

ROS_absolute_change <- ROS_absolute_change %>%
  mutate(
    ROS_diff = Absolute_diff * ROS_effects,
    Chordata_diff = Absolute_diff * chordata_effects,
    Angiosperm_diff = Absolute_diff * angiosperm_effects
  ) %>%
  dplyr::select(-ROS_baseline, -ROS_current, -Absolute_diff)

# Function to calculate average WW events per year
yearly_gridded_count_WW <- function(WW_file) {
  
  # Extracting year
  filename_without_path <- basename(WW_file)
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  year <- as.integer(filename_parts[length(filename_parts)])
  
  print(paste("processing", year, sep = " "))
  
  # Read the CSV file with lon-lat points
  WW_data <- read.csv(WW_file)
  
  # Convert lon values to -180 to 180 range
  WW_data$lon <- ifelse(WW_data$lon > 180, WW_data$lon - 360, WW_data$lon)
  
  WW_data <- WW_data[WW_data$lat > 50, ]
  
  total_WW <- WW_data %>%
    group_by(lat, lon) %>%
    summarise(total_exceedance = sum(exceedance))
  
  # Add year information to the result
  total_WW$year <- year
  
  return(total_WW)
  
}

# Specify the path to the folder containing CSV files
WW_folder <- "Arctic_Climate_Project/output/WW_gridded_t2m"

# Get a list of all CSV files in the folder
WW_files <- list.files(path = WW_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
WW_data_list <- lapply(WW_files, yearly_gridded_count_WW)

# Combine the results into a single data frame
WW_data <- bind_rows(WW_data_list)

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

filtered_WW <- filter_fun(WW_data, ocean_polygons, arctic_shape)
rm(WW_data)

filtered_baseline_WW <- filtered_WW %>%
  subset(year > 1950 & year <= 1980) %>%
  group_by(lon, lat) %>%
  summarise(WW_baseline = sum(total_exceedance)/30)

filtered_current_WW <- filtered_WW %>%
  subset(year > 1990 & year <= 2020) %>%
  group_by(lon, lat) %>%
  summarise(WW_current = sum(total_exceedance)/30)

WW_absolute_change <- full_join(filtered_baseline_WW, filtered_current_WW, by = c("lon", "lat"))
rm(filtered_baseline_WW, filtered_current_WW)

WW_absolute_change[is.na(WW_absolute_change)] <- 0

WW_absolute_change <- WW_absolute_change %>%
  mutate(Absolute_diff = WW_current - WW_baseline)

WW_absolute_change <- WW_absolute_change %>%
  mutate(
    WW_diff = Absolute_diff * WW_effects,
    Chordata_diff = Absolute_diff * chordata_effects,
    Angiosperm_diff = Absolute_diff * angiosperm_effects
  ) %>%
  dplyr::select(-WW_baseline, -WW_current, -Absolute_diff)

df_sf <- st_as_sf(WW_absolute_change, coords = c("lon", "lat"), crs = 4326)
raster_template <- rast(ext(df_sf), resolution = 1)

biota_raster_WW <- rasterize(df_sf, raster_template, field = "WW_diff", fun = mean)
chordata_raster_WW <- rasterize(df_sf, raster_template, field = "Chordata_diff", fun = mean)
angiosperm_raster_WW <- rasterize(df_sf, raster_template, field = "Angiosperm_diff", fun = mean)

# Stack the rasters
multi_layer_raster <- c(biota_raster_WW, chordata_raster_WW, angiosperm_raster_WW)

# Set layer names
names(multi_layer_raster) <- c("biota", "chordata", "angiosperm")

# Define the polar stereographic projection
polar_proj <- "+proj=stere +lat_0=90 +lon_0=-45 +lat_ts=70 +datum=WGS84"

# Project the raster
multi_layer_raster_polar <- project(multi_layer_raster, polar_proj)

# Extract each layer
biota_layer <- multi_layer_raster_polar[[1]]
chordata_layer <- multi_layer_raster_polar[[2]]
angiosperm_layer <- multi_layer_raster_polar[[3]]

coastline_mapping <- read_sf("Arctic_Climate_Project/data/ne_10m_coastline/ne_10m_coastline.shp")

graticule <- read_sf("Arctic_Climate_Project/data/ne_10d_graticule/ne_10m_graticules_10.shp")
graticule <- graticule[, "degrees"]

extent <- ext(biota_layer)
extent <- as.vector(extent)

# Create a bounding box polygon covering the entire area of interest
bounding_polygon <- st_polygon(list(rbind(c(extent[1], extent[3]), 
                                          c(extent[1], extent[4]), 
                                          c(extent[2], extent[4]), 
                                          c(extent[2], extent[3]), 
                                          c(extent[1], extent[3]))))

arctic_shape <- st_transform(arctic_shape, crs = st_crs(biota_layer))

# Convert the bounding polygon to the same CRS as arctic_shape
bounding_polygon <- st_sfc(bounding_polygon, crs = st_crs(biota_layer))

# Plot inverse of arctic_mapping in grey
arctic_inverse <- st_difference(bounding_polygon, arctic_shape)

coastline_mapping <- st_as_sfc(coastline_mapping, crs = st_crs(biota_layer))
coastline_mapping <- st_transform(coastline_mapping, st_crs(bounding_polygon))
coastline_mapping <- st_intersection(coastline_mapping, bounding_polygon)

graticule <- st_as_sfc(graticule, crs = st_crs(biota_layer))
graticule <- st_transform(graticule, st_crs(bounding_polygon))
graticule <- st_intersection(graticule, bounding_polygon)

# Convert study_locations to an sf object
study_locations_sf <- st_as_sf(study_locations, coords = c("Lon", "Lat"), crs = 4326)

# Reproject to polar projection
study_locations_polar <- st_transform(study_locations_sf, crs = polar_proj)

# Separate WW and ROS studies
WW_studies <- study_locations_polar %>% filter(Extreme_event == "WW")
ROS_studies <- study_locations_polar %>% filter(Extreme_event == "ROS")

# Define a custom color palette with specific colors
custom_palette <- rev(c("#3c546a", "#6B96B9", "#6596c1", "#b1cce1", "#e2f4fd", "#ec9172", "#de543e", "#ad2d24", "#5e1412"))

# Create a continuous palette function using your custom colors
continuous_palette <- colorRampPalette(custom_palette)
breaks <- c(-Inf, -9, -6, -3, 0, 3, 6, 9, 12, Inf)

# Generate a color map that matches the breaks
num_colors <- length(breaks) - 1
colors <- continuous_palette(num_colors)

# Plot biota_layer with reversed color scale and customized legend labels
plot(biota_layer, col = continuous_palette(100), breaks = breaks)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)
# Add WW_studies as thick black crosses
plot(
  WW_studies["Num_Studies"], 
  pch = 24,  # Cross symbol
  col = "black",  # Black color
  bg = "yellow",  # Fill color to contrast with background
  cex = 2,  # Increase size of crosses
  lwd = 2,  # Increase line thickness
  add = TRUE
)

# Plot chordata_layer with reversed color scale and customized legend labels
plot(chordata_layer, col = continuous_palette(100), breaks = breaks)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)

# Plot change_layer with reversed color scale and customized legend labels
plot(angiosperm_layer, col = continuous_palette(100), breaks = breaks)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)

df_sf <- st_as_sf(ROS_absolute_change, coords = c("lon", "lat"), crs = 4326)
raster_template <- rast(ext(df_sf), resolution = 1)

biota_raster_ROS <- rasterize(df_sf, raster_template, field = "ROS_diff", fun = mean)
chordata_raster_ROS <- rasterize(df_sf, raster_template, field = "Chordata_diff", fun = mean)
angiosperm_raster_ROS <- rasterize(df_sf, raster_template, field = "Angiosperm_diff", fun = mean)

# Stack the rasters
multi_layer_raster_ROS <- c(biota_raster_ROS, chordata_raster_ROS, angiosperm_raster_ROS)

# Set layer names
names(multi_layer_raster_ROS) <- c("biota", "chordata", "angiosperm")

# Define the polar stereographic projection
polar_proj <- "+proj=stere +lat_0=90 +lon_0=-45 +lat_ts=70 +datum=WGS84"

# Project the raster
multi_layer_raster_polar_ROS <- project(multi_layer_raster_ROS, polar_proj)

# Extract each layer
biota_layer_ROS <- multi_layer_raster_polar_ROS[[1]]
chordata_layer_ROS <- multi_layer_raster_polar_ROS[[2]]
angiosperm_layer_ROS <- multi_layer_raster_polar_ROS[[3]]

# Plot biota_layer with reversed color scale and customized legend labels
plot(biota_layer_ROS, col = continuous_palette(100), breaks = breaks)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)
# Add WW_studies as thick black crosses
plot(
  ROS_studies["Num_Studies"], 
  pch = 24,  # Cross symbol
  col = "black",  # Black color
  bg = "black",  # Fill color to contrast with background
  cex = 2,  # Increase size of crosses
  lwd = 2,  # Increase line thickness
  add = TRUE
)

# Plot chordata_layer with reversed color scale and customized legend labels
plot(chordata_layer_ROS, col = continuous_palette(100), breaks = breaks)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)

# Plot change_layer with reversed color scale and customized legend labels
plot(angiosperm_layer_ROS, col = continuous_palette(100), breaks = breaks)
plot(arctic_inverse, col = "gray", border = NA, add = TRUE)
plot(coastline_mapping, add = TRUE, col = "black", lwd = 0.5)
plot(graticule, add = TRUE, col = "black", lwd = 0.5)

