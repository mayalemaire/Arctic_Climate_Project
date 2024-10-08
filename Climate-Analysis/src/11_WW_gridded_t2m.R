# Producing a data frame with the WW events for every year
# Maya Lemaire
# 14-03-24

library(dplyr)

baseline_percentiles <- read.csv("Arctic_Climate_Project/output/1951_1980_gridded_0.99_percentile.csv")

ww_count <- function(t2m_file) {
  
  year_data <- read.csv(t2m_file)
  
  # Extracting year
  filename_without_path <- basename(t2m_file)
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  year <- filename_parts[length(filename_parts)]
  
  print(paste("processing", year, sep = " "))
  
  # Iterate through each unique grid point in testing_data
  unique_grid_points <- unique(year_data[, c("lon", "lat")])
  
  result_list <- list()
  
  for (i in 1:nrow(unique_grid_points)) {
    current_lon <- unique_grid_points$lon[i]
    current_lat <- unique_grid_points$lat[i]
    
    # Find the corresponding baseline percentile for the grid point
    baseline_percentile <- baseline_percentiles %>%
      filter(lon == current_lon, lat == current_lat) %>%
      pull(p99_temperature)
    
    # Filter year_data for the current grid point
    current_grid_data <- year_data %>%
      filter(lon == current_lon, lat == current_lat)
    
    # Filter to just have WW events
    current_grid_data <- filter(current_grid_data, t2m > baseline_percentile)
    current_grid_data$exceedance <- current_grid_data$t2m - baseline_percentile
    
    # Append the filtered dataframe to result_list
    result_list[[i]] <- current_grid_data
  }
  
  # Combine all subsets into a single data frame
  result_data <- do.call(rbind, result_list)
  
  csv_fname <- paste("Arctic_Climate_Project/output/WW_gridded_t2m/WW_gridded_t2m_", year, ".csv", sep="")
  
  print("writing")
  
  # Write the filtered dataframe to a CSV file
  write.table(result_data, file = csv_fname, row.names = FALSE, sep = ",")
}

# Specify the path to the folder containing NetCDF files
t2m_folder <- "Arctic_Climate_Project/data/2m_temperature"

# Get a list of all NetCDF files in the folder
t2m_files <- list.files(path = t2m_folder, pattern = "\\.csv$", full.names = TRUE)

mapply(FUN = ww_count, t2m_files)
