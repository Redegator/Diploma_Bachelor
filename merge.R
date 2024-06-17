setwd("D:/Диплом")
getwd()

library(purrr)
library(dplyr)

file_paths = "D:/Диплом/data"

data <- read.table(file_paths, header = FALSE, col.names = c("V1", gsm_id))

# library(GEOquery)
# library(tidyverse)


# Define a vector of GSM IDs
gsm_ids <- c("GSM2102531",	"GSM2102532",	"GSM2102533",	"GSM2102534",	"GSM2102535",	"GSM2102536",	"GSM2102537",	"GSM2102538",	"GSM2102539",	"GSM2102540",	"GSM2102541",	"GSM2102542",	"GSM2102543",	"GSM2102544",	"GSM2102545",	"GSM2102546",	"GSM2102547",	"GSM2102548",	"GSM2102549",	"GSM2102550",	"GSM2102551",	"GSM2102552",	"GSM2102553",	"GSM2102554",	"GSM2102555",	"GSM2102556",	"GSM2102557",	"GSM2102558",	"GSM2102559")


# Function to merge data from a file based on gsm_id
merge_data <- function(gsm_id) {
  # Find all files matching gsm_id pattern in the "data" directory
  file_paths <- dir("data", pattern = paste0(gsm_id, "_.*\\.txt$"), full.names = TRUE)
  # Read data from each file and assign column names
  data <- read.table(file_paths, header = FALSE, col.names = c("V1", gsm_id))
  return(data)
}

# Apply the function to each gsm_id
list_of_data <- map(gsm_ids, merge_data)

# Merge dataframes based on the "V1" column
df_merged <- reduce(list_of_data, left_join, by = "V1")


head(df_merged)

# Save the results
output_file <- "df_merged.txt"
write.table(df_merged, file = output_file, sep = "\t", row.names = FALSE)

