##Program name: 9_summary_df_ibsi.r
##Purpose: To summarize from the image headers the information needed
##        for ibsi guideline reporting
# Path to your data
PATH_TO_P_DRIVE = "/DataLibrary"
library("tidyverse")
# Used to read the RDS files from P-drive in parallel.
library("future")
library("future.apply")
plan(multisession)

all_header_rds <- Sys.glob(file.path(PATH_TO_P_DRIVE, "/GRADS_CT/DataProcessed/1_ct_nifti/*.rds"))
all_header_rds_df <- tibble(path=all_header_rds, sid= str_to_upper(str_extract(all_header_rds, ".*/([sS0-9]{6})_header.rds",group=1)))
sprintf("There are %d files in total", dim(all_header_rds_df)[1])

# Read in list of GRADS subjects, should be 321.
grads_subjs <- read_csv(paste0(PATH_TO_P_DRIVE,"GRADS_CT/Subjects_to_Include.csv"))
sprintf("There are %d GRADS subjects.", length(grads_subjs$gradsid))

use_paths <- (all_header_rds_df %>% filter(sid %in% grads_subjs$gradsid))$path
sprintf("Reading headers from %d paths.", length(use_paths))

# Split up DICOM fields into string, numeric, and numeric list types. 
# We need to parse them separately for everything to work.
string_fields = c("PatientName",
                  "SeriesDescription",
                  "SoftwareVersions",
                  "Manufacturer",
                  "ManufacturerModelName",
                  "ConvolutionKernel")
numeric_fields = c("KVP",
                   "XRayTubeCurrent",
                   "Exposure",
                   "ExposureTime",
                   "SliceThickness",
                   "SpacingBetweenSlices",
                   "InstanceNumber")
numeric_list_fields = c("ImagePositionPatient",
                        "PixelSpacing")
fields_to_extract <- c(string_fields, numeric_fields, numeric_list_fields)

# Function to read a RDS file, filter for any of the columns we want (not all files have all columns), and then add a path variable)
read_and_extract_RDS <- function(path) {
    readRDS(path) %>% filter(`name` %in% fields_to_extract)  %>% mutate(path = path)
}

# Read files in parallel and conatenate into a dataframe
all_header_futures_df <- do.call(rbind, future_lapply(use_paths, FUN = read_and_extract_RDS))

#write_rds(x=all_header_futures_df, file="all_header_futures_df.RDS")

object.size(all_header_futures_df)

parse_single_str <- function(x) {
    str_remove_all(x, "[\\[\\]]")
}
parse_single_num <- function(x) {
    as.numeric(parse_single_str(x))
}

out_df <- all_header_futures_df %>% 
            pivot_wider(id_cols = all_of(c("path","file")), names_from = `name`, values_from = `value`) %>% 
            mutate(
                across(all_of(c(string_fields, numeric_list_fields)), parse_single_str),
                across(all_of(numeric_fields), parse_single_num)        
            ) %>%
            separate_wider_delim(
                all_of(numeric_list_fields), delim="\\", names_sep=""
            ) %>%
            mutate(
                across(any_of(c(starts_with("PixelSpacing"),starts_with("ImagePositionPatient"))), as.numeric),
            ) %>%
            group_by(`PatientName`) %>%
            arrange(`InstanceNumber`, by_group=T) %>%
            mutate(                
                PixelSpacing3 = lag(ImagePositionPatient3) - ImagePositionPatient3
            ) %>%
            fill(PixelSpacing3, .direction = "up") %>% ungroup

get_stats <- function(x) {
    tibble(`mean` = mean(x, na.rm = T),
           `std` = sd(x, na.rm = T),
           `min` = min(x, na.rm = T),
           `max` = max(x, na.rm = T),
           `num.non.na` = sum(!is.na(x)),
           `pct.non.na` = sum(!is.na(x))/length(x))
}
get_unq_vals <- function(x) {
    tibble(`unq.vals` = paste0(unique(x),collapse=", "))
}

global_summary <- out_df %>% reframe(
    across(all_of(c(numeric_fields, c("PixelSpacing1","PixelSpacing2","PixelSpacing3"))), get_stats, .unpack=T),
    across(all_of(setdiff(string_fields, c("PatientName", "SeriesDescription"))), get_unq_vals, .unpack=T)
) %>% pivot_longer(everything(), names_to = c("keyword",".value"), names_sep="_")

per_subj_summary <- out_df %>% reframe(
    across(all_of(c(numeric_fields, c("PixelSpacing1","PixelSpacing2","PixelSpacing3"))), get_stats, .unpack=T),
    across(all_of(setdiff(string_fields, c("PatientName", "SeriesDescription"))), get_unq_vals, .unpack=T),
    .by="PatientName"
) %>% pivot_longer(-any_of(c("PatientName")), names_to = c("keyword",".value"), names_sep="_")

write_csv(x = global_summary, file = file.path(PATH_TO_P_DRIVE,"/all_files_summary.csv"))

write_csv(x = per_subj_summary, file = file.path(PATH_TO_P_DRIVE,"/per_subj_summary.csv"))


