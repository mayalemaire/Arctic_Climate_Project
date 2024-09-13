# Adapting plot from Martin et al., 2020
# 13-09-24
# Maya Lemaire

# Load necessary libraries
library(sf)          # For handling spatial data
library(dplyr)       # For data manipulation
library(RColorBrewer) # For color palettes

# Read the shapefile for the Arctic
arctic_shape <- st_read("Arctic_Climate_Project/data/evidence-map-scope/evidence-map-scope.shp")

# Filter the shapes with OBJECTID 1 to 4
selected_shapes <- arctic_shape %>% filter(OBJECTID_1 %in% c(1, 2, 3, 4))

# Set the CRS to polar stereographic (EPSG:3413 is commonly used for Arctic regions)
selected_shapes <- st_transform(selected_shapes, crs = 3413)

# Read the coastline and graticule layers
coastline_mapping <- read_sf("Arctic_Climate_Project/data/ne_10m_coastline/ne_10m_coastline.shp")
graticule <- read_sf("Arctic_Climate_Project/data/ne_10d_graticule/ne_10m_graticules_10.shp")
graticule <- graticule[, "degrees"]

# Ensure all layers have the same CRS
coastline_mapping <- st_transform(coastline_mapping, st_crs(selected_shapes))
graticule <- st_transform(graticule, st_crs(selected_shapes))

# Plot each Arctic shape with a different color
plot(st_geometry(selected_shapes[selected_shapes$Zone == "Sub arctic",]), col = "#42A5F5", main = "Arctic Shapes with Coastline and Graticule")
plot(st_geometry(selected_shapes[selected_shapes$Zone == "CAFF",]), col = "#FF8F00", add = TRUE)
plot(st_geometry(selected_shapes[selected_shapes$Zone == "High arctic",]), col = "#4db368", add = TRUE)
plot(st_geometry(selected_shapes[selected_shapes$Zone == "Low arctic",]), col = "#D32F2F", add = TRUE)

# Overlay the coastline
plot(st_geometry(coastline_mapping), col = "black", lwd = 0.5, add = TRUE)

# Overlay the graticule
plot(st_geometry(graticule), col = "black", lwd = 0.5, add = TRUE)

# Add a simplified legend for Arctic Zones
legend("topright", legend = c("CAFF", "High Arctic", "Low Arctic", "Sub Arctic"), 
       fill = c("#FF8F00", "#4db368", "#D32F2F", "#42A5F5"), 
       cex = 0.7, bty = "o", bg = "white")

