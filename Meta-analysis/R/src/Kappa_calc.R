# Calculating Kappa value to determine agreement between author and independent
# Maya Lemaire
# 04-09-24

# Load necessary package
library(irr)

# Replace 'your_file.csv' with the path to your CSV file
file_path <- 'data/Kappa_calc.csv'

# Load the data from the CSV file
df <- read.csv(file_path, sep = ";")

# Use the 'First_author' and 'Independent' columns as the two raters
rater1 <- df$First_author
rater2 <- df$Independent

# Create a data frame of the ratings
ratings <- data.frame(rater1, rater2)

# Calculate Cohen's Kappa
kappa_result <- kappa2(ratings)

# Print the Kappa result
print(kappa_result)
