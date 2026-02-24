#---------------------------------------
# Set GBIF credentials as env variables
#---------------------------------------
Sys.setenv(GBIF_USER = "your_username")
Sys.setenv(GBIF_PWD = "your_password")
Sys.setenv(GBIF_EMAIL = "your_email")


#--------------------------------------
# Load libraries
#--------------------------------------

library(rgbif)

#--------------------------------------
# Download Florabank2 dataset
#--------------------------------------

dataset_key <- "1e9b6eff-af44-4e48-90f0-35ca8d2cdb7b"

download_key <- occ_download(
  pred("datasetKey", dataset_key),
  format = "SIMPLE_CSV",
  user = Sys.getenv("GBIF_USER"),
  pwd = Sys.getenv("GBIF_PWD"),
  email = Sys.getenv("GBIF_EMAIL")
)

# Wait for download to finish
occ_download_wait(download_key)

# Download to local folder
occ_download_get(download_key, path = "data/")

# Import into R
data <- occ_download_import(download_key)