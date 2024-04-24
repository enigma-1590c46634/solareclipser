## Description:
##   Whether the subject was ascertained using Siemens (1) or
##   non-Siemens (0) scanner. Add a column with these binary values.



source("./R/install_packages.R")

a0_wave_01   <- read_csv("./data/csv/input/wave_01.csv") #bl
a0_wave_02   <- read_csv("./data/csv/input/wave_02.csv") #yr2
a0_wave_03   <- read_csv("./data/csv/input/wave_03.csv") #yr4
a0_site_info <- read_csv("./data/csv/input/site_info.csv")
a0_mri_info  <- read_csv("./data/csv/input/mri_info.csv")
a0_depression_score <- read_csv("./data/csv/input/depression_score.csv")
a0_ple_score        <- read_csv("./data/csv/input/PLE_score.csv")
a0_ple <- read_csv("./data/csv/input/PLE.csv")

## New Data 2024-01-03
subjs <- read_csv("./temp/subjects.csv")

## --- Tidy ---
## Rename columns
a0_depression_score <-
  a0_depression_score %>% rename(subject_id = 1)
a0_mri_info <-
  a0_mri_info %>% rename(subject_id = 1)
a0_ple <-
  a0_ple %>% rename(subject_id = 1)
a0_ple_score <-
  a0_ple_score %>% rename(subject_id = 1)
a0_site_info <-
  a0_site_info %>% rename(subject_id = 1)
subjs <-
  subjs %>% rename(subject_id = 1)

## --- Transform ---
## Return only rows with event name "baseline_year_1_arm_1"
a1_site_info_bl <-
  a0_site_info[a0_site_info$eventname == "baseline_year_1_arm_1", ]

a1_mri_info_bl <-
  a0_mri_info[a0_mri_info$eventname == "baseline_year_1_arm_1", ]

## Get the common subject_id from subjs and a0_wave_01
a1_com_subjs <-
  intersect(subjs$subject_id, a0_wave_01$subject_id)

## Create a data frame matching subject_id from
## a1_com_subjs and a1_mri_info_bl. On match, return the
## mri_info_manufacturer and mri_info_manufacturersmn from a1_mri_info_bl.
a2_mri_info_bl <-
  a1_mri_info_bl[a1_mri_info_bl$subject_id %in% a1_com_subjs,
                 c("subject_id",
                   "mri_info_manufacturer",
                   "mri_info_manufacturersmn")]

## Add a column to a2_mri_info_bl under the following conditions:
##  - If mri_info_manufacturer is "SIEMENS", then "1"
##  - else then "0"
a2_mri_info_bl$siemens <-
  ifelse(a2_mri_info_bl$mri_info_manufacturer == "SIEMENS", 1, 0)

## Rename the column to is_siemens
a3_result <-
  a2_mri_info_bl %>% rename(is_siemens = siemens)

## Drop the mri_info_manufacturer and mri_info_manufacturersmn columns
a3_result <-
  a3_result %>% select(-mri_info_manufacturer, -mri_info_manufacturersmn)

## Write to csv
write_csv(a3_result, "./data/csv/output/subj_scanner_is_siemens.csv")
