# Producing a data frame with the ROS events for criteria:
# t2m > 274.15, rf > 3mm, sd > 3, and mean max temp over next 3 days < 273.24
# Scale - 1 by 1 degree grid boxes for entire arctic region
# 24-01-25
# Maya Lemaire

library(dplyr)

ROS_count <- function(ROS_file) {

  result_data <- data.frame(lon = integer(), lat = integer(), date = character(), rf = numeric()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 )
  
  # Importing data
  year_data <- read.csv(ROS_file) %>%
    na.omit()  # Remove NAs
  
  # Extracting year
  filename_without_path <- basename(ROS_file)
  
  # Remove the ".nc" extension
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  
  # Split the filename using "_" as a separator
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  
  # Extract the year (assuming it's always in the same format
  year <- filename_parts[length(filename_parts)]
  
  print(paste("processing", year, sep = " "))
  
  # Filter data based on conditions for ROS event
  result_data <- year_data %>%
    filter(rf > 0.003, t2m > 274.15, sd > 0.003, mean_t2m_3days < 273.15) %>%
    dplyr::select(lon, lat, date, rf)
  
  if (nrow(result_data) > 0) {
    csv_fname <- paste("Arctic_Climate_Project/output/ROS_3mm_1t2m_sd_m3d_2/ROS_gridded_3mm_1t2m_sd", year, ".csv", sep="")
    
    print("writing")
    
    # Write the filtered dataframe to a CSV file
    write.table(result_data, file = csv_fname, row.names = FALSE, sep = ",")
    
  }

}

# Specify the path to the folder containing NetCDF files
ROS_folder <- "Arctic_Climate_Project/data/t2m_rf_sd_m3d"

# Get a list of all NetCDF files in the folder
ROS_files <- list.files(path = ROS_folder, pattern = "\\.csv$", full.names = TRUE)

mapply(FUN = ROS_count, ROS_files)

