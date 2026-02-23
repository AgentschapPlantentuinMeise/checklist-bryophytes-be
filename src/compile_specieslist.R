###########################
# Load necessary packages #
###########################

library(readr)
library(dplyr)
library(rgbif)
library(stringr)
library(purrr)
library(httr)
library(jsonlite)
library(tibble)
library(tidyr)

#################
# Read in files #
#################

# Florabank2 occurrence dataset
mosses_be_occ <- read_tsv('data/raw/sources/florabank2/0009938-251120083545085/occurrence.txt')
df <- mosses_be_occ
mosses_be_specieslist <- mosses_be_occ[mosses_be_occ$taxonRank=="SPECIES",] %>%
  select(
    taxonKey,
    acceptedTaxonKey,
    scientificName,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    taxonRank
  ) %>%
  distinct(taxonKey, .keep_all = TRUE) %>%
  mutate(
    bibliographicReference = paste0("Van Landuyt W, Brosens D, De Beer D",
                                    "(2025). Florabank2: a grid-based database on distribution of bryophytes",
                                    " in the northern part of Belgium (Flanders and the Brussels Capital",
                                    " region). Version 1.18. Research Institute for Nature and Forest (INBO).",
                                    " Occurrence dataset https://doi.org/10.15468/385t22")
  )

# Quentin's list of bryophytes from publications

bryo_pub <- read_csv('data/raw/sources/bryophytes_published.csv')

######################################
# Cleaning and merging species lists #
######################################

# Define match_species function to replace rgbif::name_backbone_specieslist()
match_species <- function(sp_name, dataset_key, pause = 0.2) {
  Sys.sleep(pause)  # avoid API throttling
  
  res <- GET(
    "https://api.gbif.org/v1/species/match",
    query = list(name = sp_name, datasetKey = dataset_key)
  )
  
  if (res$status_code != 200) {
    warning("API request failed for ", sp_name)
    return(tibble(submittedName = sp_name, key = NA_character_))
  }
  
  data <- fromJSON(content(res, "text", encoding = "UTF-8"), simplifyVector =
                     TRUE)
  
  tibble(
    submittedName     = sp_name,
    usageKey          = data$usageKey,
    scientificName    = data$scientificName,
    status            = data$status,
    rank              = data$rank,
    matchType         = data$matchType,
    confidence        = data$confidence,
    acceptedUsageKey  = data$acceptedUsageKey,
    kingdomKey        = data$kingdomKey,
    phylumKey         = data$phylumKey,
    classKey          = data$classKey,
    orderKey          = data$orderKey,
    familyKey         = data$familyKey
  )
}

#====================================================
# Prep Quentin's bryophyte list from 3 publications
#====================================================

# Match scientific names in Quentin's spreadsheet with the GBIF backbone

## Using rgbif (doesn't return acceptedUsageKey at the moment)

# bryo_pub_matched <- name_backbone_checklist(bryo_pub$taxon)

## Using match_species()

species_list <- bryo_pub$taxon
dataset_key <- "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c" #gbif backbone

# Loop over all species
bryo_pub_matched <- bind_rows(lapply(species_list, match_species,
                                     dataset_key = dataset_key))

# Inspect result
bryo_pub_matched
table(bryo_pub_matched$matchType)
table(bryo_pub_matched$status)

# Keep only EXACT matches
bryo_pub_matched_exact <- bryo_pub_matched %>% filter(matchType == "EXACT")

# Make new column which contains the acceptedUsageKey for synonyms and the
# usageKey for accepted species

bryo_pub_matched_exact <- bryo_pub_matched_exact %>%
  mutate(
    taxonKey = case_when(
      status == "ACCEPTED" ~ usageKey,
      status == "DOUBTFUL" ~ usageKey,
      status == "SYNONYM"  ~ acceptedUsageKey,
      TRUE ~ NA_integer_
    )
  )

# Link the GBIF usageKey and taxonKey (=accepted taxonKey) to the names in the
# list
bryo_pub_wkeys <-
  bryo_pub %>%
  left_join(
    bryo_pub_matched_exact %>% select(usageKey, taxonKey, submittedName),
    by = c("taxon" = "submittedName")
  )

bryo_pub_wkeys

# Look up extra columns to fill out template

extra_col <- map_df(bryo_pub_wkeys$usageKey, ~ {
  dat <- name_usage(key = .x)$data
  
  dat %>%
    select(
      key,
      kingdom,
      phylum,
      class,
      order,
      family,
      genus,
      rank
    )
})

bryo_pub_wkeys <-
  bryo_pub_wkeys %>%
  left_join(extra_col, by= c("usageKey"="key"))

# Re-order columns to fit template
bryo_pub_wkeys_curated <-
  bryo_pub_wkeys %>%
  select(
    usageKey,
    taxonKey,
    taxon,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    rank,
    everything()
  ) %>%
  select(-Division) %>%
  rename(acceptedTaxonKey = taxonKey) %>%
  rename(taxonKey = usageKey) %>%
  rename(scientificName = taxon) %>%
  rename(taxonRank = rank)

# write.csv(bryo_pub_wkeys,"./output/bryo_pub_wkeys.csv")

# Combine the two datasets into one checklist

## Check which species of the florabank2 dataset are not yet in
## Quentin's dataset

taxonKeys_to_add <- setdiff(mosses_be_specieslist$taxonKey,
                            bryo_pub_wkeys$usageKey)
species_to_add <-
  mosses_be_specieslist[mosses_be_specieslist$taxonKey %in% taxonKeys_to_add,]

# Write.csv(species_to_add,"./output/species_to_add.csv")

## Add extra species data to bryo_pub_wkeys

bryophyte_list_belgium <- bind_rows(bryo_pub_wkeys_curated, species_to_add)

## Retrieve first and last dates of observation in Belgium

# Generate cube with min and max eventDate
accepted_keys <- unique(bryophyte_list_belgium$acceptedTaxonKey)

sql_query <- sprintf("
  SELECT
    acceptedTaxonKey,
    MIN(eventDate) AS first_observation_date,
    MAX(eventDate) AS last_observation_date
  FROM occurrence
  WHERE countryCode = 'BE'
    AND acceptedTaxonKey IN (%s)
    AND eventDate IS NOT NULL
  GROUP BY acceptedTaxonKey
  ORDER BY acceptedTaxonKey
", paste(accepted_keys, collapse = ", "))

download_key <- occ_download_sql(
  q = sql_query,
  user = "username",
  pwd = "passw",
  email = "email"
)

firstlastdates <- occ_download_get(download_key) %>%
  occ_download_import()

# read in cube with first and last observation dates
firstlastdates <- read.csv("data/raw/sources/cube/0046111-251120083545085/0046111-251120083545085.csv",
                           header=TRUE, sep="\t")

firstlastdates <- firstlastdates %>%
  rename("acceptedTaxonKey"="acceptedtaxonkey")

# Join back to species list
bryophyte_list_belgium <- bryophyte_list_belgium %>%
  left_join(firstlastdates, by = "acceptedTaxonKey")

#===============================================================================
# The following code is abandoned as it was decided to use exact eventDates from
# GBIF occurrences for all taxa and not compare with date from literature.
#===============================================================================

# Cleaning up dates

## Convert year and yearmonth and iso datetime to iso date

# bryophyte_list_belgium$date_fixed_first <- dplyr::case_when(
#   # year only
#   grepl("^\\d{4}$", bryophyte_list_belgium$first_observation_date) ~
#     paste0(bryophyte_list_belgium$first_observation_date, "-01-01"),
#
#   # year-month (YYYY-MM)
#   grepl("^\\d{4}-\\d{2}$", bryophyte_list_belgium$first_observation_date) ~
#     paste0(bryophyte_list_belgium$first_observation_date, "-01"),
#
#   # ISO date or datetime (YYYY-MM-DD or YYYY-MM-DDTHH:MM...)
#   grepl("^\\d{4}-\\d{2}-\\d{2}", bryophyte_list_belgium$first_observation_date) ~
#     substr(bryophyte_list_belgium$first_observation_date, 1, 10),
#
#   # otherwise leave as-is
#   TRUE ~ NA_character_
# )

# bryophyte_list_belgium$date_fixed_last <- dplyr::case_when(
#   # year only
#   grepl("^\\d{4}$", bryophyte_list_belgium$last_observation_date) ~
#     paste0(bryophyte_list_belgium$last_observation_date, "-01-01"),
#
#   # year-month (YYYY-MM)
#   grepl("^\\d{4}-\\d{2}$", bryophyte_list_belgium$last_observation_date) ~
#     paste0(bryophyte_list_belgium$last_observation_date, "-01"),
#
#   # ISO date or datetime (YYYY-MM-DD or YYYY-MM-DDTHH:MM...)
#   grepl("^\\d{4}-\\d{2}-\\d{2}", bryophyte_list_belgium$last_observation_date) ~
#     substr(bryophyte_list_belgium$last_observation_date, 1, 10),
#
#   # otherwise leave as-is
#   TRUE ~ NA_character_
# )

# bryophyte_list_belgium$first_obs_date_gbif <-
#   as.Date(bryophyte_list_belgium$date_fixed_first)
# bryophyte_list_belgium$last_obs_date_gbif <-
#   as.Date(bryophyte_list_belgium$date_fixed_last)
# bryophyte_list_belgium <-
#   bryophyte_list_belgium %>%
#   select(-first_observation_date,-last_observation_date,
#          -date_fixed_first, -date_fixed_last)

## Rename literature date columns
bryophyte_list_belgium <-
  bryophyte_list_belgium %>%
  rename(date_from_lit='date from', date_to_lit='date to')

# ## Clean up dates from publication
#
# bryophyte_list_belgium$clean_from_date <- dplyr::case_when(
#
#   # true NA
#   is.na(bryophyte_list_belgium$date_from_lit) ~ NA_character_,
#
#   # text "NA"
#   bryophyte_list_belgium$date_from_lit == "NA" ~ NA_character_,
#
#   # c.year → NA
#   grepl("^c\\.?\\d{4}$", bryophyte_list_belgium$date_from_lit) ~ NA_character_,
#
#   # <year → extract the year and convert
#   grepl("^<\\d{4}$", bryophyte_list_belgium$date_from_lit) ~
#     paste0(sub("^<", "", bryophyte_list_belgium$date_from_lit), "-01-01"),
#
#   # pure 4-digit year
#   grepl("^\\d{4}$", bryophyte_list_belgium$date_from_lit) ~
#     paste0(bryophyte_list_belgium$date_from_lit, "-01-01"),
#
#   TRUE ~ NA_character_
# )

# bryophyte_list_belgium$clean_from_date <-
#   as.Date(bryophyte_list_belgium$clean_from_date)
#
# bryophyte_list_belgium$clean_to_date <- dplyr::case_when(
#
#   # true NA
#   is.na(bryophyte_list_belgium$date_to_lit) ~ NA_character_,
#
#   # text "NA"
#   bryophyte_list_belgium$date_to_lit == "NA" ~ NA_character_,
#
#   # c.year → NA
#   grepl("^c\\.?\\d{4}$", bryophyte_list_belgium$date_to_lit) ~ NA_character_,
#
#   # <year → extract the year and convert
#   grepl("^<\\d{4}$", bryophyte_list_belgium$date_to_lit) ~
#     paste0(sub("^<", "", bryophyte_list_belgium$date_to_lit), "-01-01"),
#
#   # pure 4-digit year
#   grepl("^\\d{4}$", bryophyte_list_belgium$date_to_lit) ~
#     paste0(bryophyte_list_belgium$date_to_lit, "-01-01"),
#
#   TRUE ~ NA_character_
# )

# bryophyte_list_belgium$clean_to_date <-
#   as.Date(bryophyte_list_belgium$clean_to_date)
# bryophyte_list_belgium$date_from_lit <-
#   bryophyte_list_belgium$clean_from_date
# bryophyte_list_belgium$date_to_lit <-
#   bryophyte_list_belgium$clean_to_date
# bryophyte_list_belgium <-
#   bryophyte_list_belgium %>% select(-clean_from_date, -clean_to_date)

## Create final date column with earliest and latest date

# bryophyte_list_belgium <- bryophyte_list_belgium %>%
#   mutate(
#     first_date = case_when(
#       # first_obs_date_gbif is NA or empty → use date_from_lit
#       is.na(first_obs_date_gbif) | first_obs_date_gbif == "" ~ date_from_lit,
#
#       # first_obs_date_gbif is later than date_from_lit → use date_from_lit
#       first_obs_date_gbif > date_from_lit ~ date_from_lit,
#
#       # Otherwise use first_obs_date_gbif
#       TRUE ~ first_obs_date_gbif
#     ),
#     last_date = case_when(
#       # last_obs_date_gbif is NA or empty → use date_from_lit
#       is.na(last_obs_date_gbif) | last_obs_date_gbif == "" ~ date_to_lit,
#
#       # last_obs_date_gbif is earlier than date_to_lit → use date_to_lit
#       last_obs_date_gbif < date_to_lit ~ date_to_lit,
#
#       # Otherwise use last_obs_date_gbif
#       TRUE ~ last_obs_date_gbif
#     )
#
#   )

## Add distribution data for additional mosses from florabank2 (Quentin)

bryophyte_list_belgium <- bryophyte_list_belgium %>%
  mutate(
    Flanders = as.integer(Flanders),  # convert from character to integer
    Brussels = na_if(Brussels, "?"),   # convert "?" to NA
    Brussels = as.integer(Brussels),  # convert from character to integer
    Wallonia = na_if(Wallonia, "?"),   # convert "?" to NA
    Wallonia = as.integer(Wallonia))  # convert from character to integer


add_moss_distr <-
  read.csv("data/raw/sources/florabank2/additional_mosses_distr.csv", header=TRUE)
add_moss_distr <-
  add_moss_distr %>% select(key, Flanders, Brussels, Wallonia) %>%
  rename(taxonKey=key)
bryophyte_list_belgium <- bryophyte_list_belgium  %>%
  rows_update(add_moss_distr, by = "taxonKey")

## Create distribution column Belgium

bryophyte_list_belgium <- bryophyte_list_belgium %>% mutate(Belgium = 1) %>%
  relocate(Belgium, .after = Brussels)

# Final table
head(bryophyte_list_belgium)
# write.csv(bryophyte_list_belgium,"./output/bryophyte_list_belgium.csv")

## vertical version

# Pivot longer and filter presence = 1
bryophyte_list_be_longv <- bryophyte_list_belgium %>%
  pivot_longer(
    cols = c(Belgium, Flanders, Wallonia, Brussels),
    names_to = "location",
    values_to = "presence"
  ) %>%
  filter(presence == 1) %>%
  select(-presence)

# write.csv(bryophyte_list_be_longv,"./output/bryophyte_list_be_longv.csv")

# mosses_be_occ <- read_csv('')

# If date is a range, take first date

keep_first_date <- function(x) {
  sapply(strsplit(x, "/"), `[`, 1)
}

cl_belgium$first_observation_date <-
  keep_first_date(cl_belgium$date.first.observation)
cl_belgium$last_observation_date <-
  keep_first_date(cl_belgium$date.last.observation)

