#-------- Load libraries --------
library(sf)
library(spdep)
library(spatstat)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(fixest)
library(did)
library(ggplot2)

#-------- Load shapefiles --------

blckgrps_2020 <- st_read("nhgis0014_shape\\nhgis0014_shape\\nhgis0014_shapefile_tl2020_us_blck_grp_2020\\US_blck_grp_2020.shp")
study_area <- st_read("study_area.shp")

denver_bgs_2020 <- st_filter(
  blckgrps_2020,
  study_area,
  .predicate = st_within
)

RTD <- st_read("Local_transit\\Local_transit.shp")
RTD_stations <- st_read("LightrailStations\\LightrailStations.shp")

#-------- Filter commuter rail stations --------

RTD_stations_commuter <- RTD_stations %>%
  filter(str_detect(RAIL_LINE,"A|B|G|N")) %>%
  mutate(
    year = case_when(
      str_detect(RAIL_LINE, "A") ~ 2016,
      str_detect(RAIL_LINE, "N") & !str_detect(RAIL_LINE, "A") ~ 2020,
      CITY == "Westminster" ~ 2016,
      TRUE ~ 2019
    )
  )

#-------- Load OD Cost Matrix Network Analysis --------

matrices <- list.files(
  "OD_Cost_Matrices",
  pattern = "\\.csv$",
  full.names = TRUE
)

matrix_list <- setNames(
  lapply(matrices, read.csv, stringsAsFactors = FALSE),
  tools::file_path_sans_ext(basename(matrices))
)

matrix_list <- lapply(matrix_list, function(df) {
  df %>%
    mutate(
      OriginID = str_extract(Name, "^[^ ]+(?= - )"),
      DestinationID = str_extract(Name, "(?<= - ).*$")
    )
})

#-------- Load LODES LEHD Data --------

lodes <- list.files(
  "OnTheMap",
  pattern = "^polygon_.*\\.csv$",
  full.names = TRUE
)

lodes_list <- setNames(
  lapply(lodes, read.csv, stringsAsFactors = FALSE),
  str_sub(tools::file_path_sans_ext(basename(lodes)),9,12)
)

lodes_list <- lapply(lodes_list, function(df) {
  df %>%
    mutate(
      GISJOIN = str_replace(id,"8","G080"),
      GISJOIN = str_c(
        str_sub(GISJOIN, 1, 7),
        "0",
        str_sub(GISJOIN, 8)
      ),
      county_fips = str_c(
        str_sub(GISJOIN, 3, 3),   # 3rd character
        str_sub(GISJOIN, 5, 7)    # 5th–7th characters
      )
    )
})

#-------- Load County-level QCEW Data --------

county_fips <- c("08001","08005","08014","08031","08035","08059")
NAICS <- c("11","21","22","23","31-33","42","44-45",
           "48-49","51","52","53","54","55","56",
           "61","62","71","72","81","92")

qcew_files <- list.files(
  "QCEW",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

fips_pattern <- paste0(county_fips, collapse = "|")

qcew_files <- qcew_files[str_detect(qcew_files, fips_pattern)]

qcew_list <- setNames(
  lapply(qcew_files, read.csv, stringsAsFactors = FALSE),
  tools::file_path_sans_ext(basename(qcew_files))
)

#-------- Prep QCEW Data for Wage Estimate using NAICS Sectors --------

qcew_wage <- lapply(qcew_list, function(df){
  df %>% 
    filter(
      qtr == 2,
      industry_code %in% NAICS
    ) %>%
    mutate(area_fips = as.character(area_fips)) %>%
    group_by(area_fips,industry_code) %>%
    summarize(
      avg_empl = sum(month1_emplvl,month2_emplvl,month3_emplvl)/3,
      total_wages = sum(total_qtrly_wages),
      avg_weekly_wage = (total_wages/avg_empl)/13
    )
})

#-------- Load and clean BDS Firm Entry Data --------

bds <- read.csv("bds2023_st_cty_sec.csv")

bds.clean <- bds %>%
  filter(
    year >= 2013,
    year <= 2023,
    st == 8,
    cty %in% c(1, 5, 14, 31, 35, 59)
  ) %>%
  mutate(
    estabs_entry = as.character(estabs_entry),
    estabs_entry = ifelse(
      estabs_entry == "D",
      0,
      as.numeric(estabs_entry)
    ),
    area_fips = paste0(
      sprintf("%01d", st),
      sprintf("%03d", cty)
    )
  )

bds_list <- bds.clean %>%
  mutate(year = as.character(year)) %>%   # important for list indexing
  group_split(year) %>%
  setNames(unique(bds.clean$year))

#-------- Define helper functions for LODES, BDS & QCEW joining --------

pivot_lodes <- function(df) {
  # rename cns columns to NAICS codes
  cns_cols <- paste0("cns", str_pad(1:20, 2, pad = "0"))
  colnames(df)[match(cns_cols, colnames(df))] <- NAICS
  
  df %>%
    pivot_longer(
      cols = NAICS,
      names_to = "NAICS2",
      values_to = "jobs"
    ) %>%
    filter(jobs > 0)
}

compute_bg_employment_weights <- function(lodes_long) {
  
  lodes_long %>%
    filter(jobs > 0) %>%
    group_by(county_fips, NAICS2) %>%
    mutate(
      county_jobs = sum(jobs, na.rm = TRUE),
      emp_weight = jobs / county_jobs
    ) %>%
    ungroup()
}

allocate_bds_births_bg <- function(bg_weights, bds_cty) {
  
  bg_weights %>%
    inner_join(
      bds_cty %>%
        select(year, area_fips, sector, estabs_entry) %>%
        rename(NAICS2 = sector),
      by = c("county_fips" = "area_fips", "NAICS2")
    ) %>%
    mutate(
      estabs_entry_bg = estabs_entry * emp_weight
    )
}

estimate_bg_births <- function(lodes_long, bds_cty) {
  
  lodes_long %>%
    compute_bg_employment_weights() %>%
    allocate_bds_births_bg(bds_cty)
}

join_lodes_qcew <- function(lodes_df, qcew_df) {
  lodes_df %>%
    left_join(
      qcew_df %>%
        select(
          county_fips = area_fips,
          NAICS2 = industry_code,
          avg_weekly_wage
        ),
      by = c("county_fips", "NAICS2")
    )
}

compute_shift_share <- function(df) {
  df %>%
    group_by(GISJOIN) %>%
    mutate(
      total_jobs = sum(jobs, na.rm = TRUE),
      weight = jobs / total_jobs,
      weighted_wage = weight * avg_weekly_wage
    ) %>%
    summarise(
      bg_avg_weekly_wage = sum(weighted_wage, na.rm = TRUE),
      .groups = "drop"
    )
}

#-------- Join LODES to bgs and estimate firm entry --------

lodes_st_list <- lapply(lodes_list,function(df){
  df %>% inner_join(denver_bgs_2020, by = "GISJOIN")
})

lodes_st_list_long <- lapply((lodes_st_list),pivot_lodes)

bg_births <- lapply(names(lodes_st_list_long), function(y) {
  estimate_bg_births(
    lodes_long = lodes_st_list_long[[y]],
    bds_cty    = bds_list[[y]]
  )
})

names(bg_births) <- names(lodes_st_list_long)

#-------- Join QCEW to LODES and estimate wages --------

wages_destination <- lapply(names(lodes_list), function(l_name) {
  y <- str_extract(l_name, "\\d{4}")
  
  lodes_list[[l_name]] %>%
    pivot_lodes() %>%
    join_lodes_qcew(qcew_yearly[[y]]) %>%
    compute_shift_share()
})

names(wages_destination) <- names(qcew_yearly)

#-------- Join wages to OD Cost Matrix & calculate RCMA --------

compute_RCMA <- function(od_df, wages_df, epsilon = 1e-3) {
  od_df %>%
    left_join(
      wages_df %>% rename(DestinationID = GISJOIN),
      by = "DestinationID"
    ) %>%
    mutate(RCMA_contrib = bg_avg_weekly_wage / if_else(Total_TravelTime == 0, epsilon, Total_TravelTime))
}

RCMA_list <- lapply(names(wages_destination), function(y) {
  
  od_mat <- matrix_list[[paste0(y,"CostMatrix")]]
  wages <- wages_destination[[y]]
  
  compute_RCMA(od_mat, wages) %>%
    group_by(OriginID) %>%
    summarise(RCMA = (sum(RCMA_contrib, na.rm = TRUE))^(-3.4)) %>%
    mutate(year = y)
})

names(RCMA_list) <- names(wages_destination)

#-------- Point-Pattern Analysis using Centroids --------

control.treatment <- function(blckgrp, rail_stations, stdevs=2) {
  # ---- CRS setup ----
  target_crs <- 26913  # UTM Zone 13N
  blckgrp <- st_transform(blckgrp, target_crs)
  rail_stations <- st_transform(rail_stations, target_crs)
  study_area <- st_transform(study_area, target_crs)
  
  # ---- Geometry prep ----
  center <- st_centroid(blckgrp)
  study_area_valid <- st_make_valid(st_union(study_area))
  rail_stations_valid <- st_make_valid(rail_stations)
  rail_buffer <- st_buffer(rail_stations_valid, 1180) # buffer 1 standard deviation
  
  study_buffer <- st_buffer(study_area_valid, 200)
  win <- as.owin(st_union(study_buffer))
  
  # ---- Convert to planar points ----
  center.xy <- st_coordinates(center)
  rail.xy <- st_coordinates(rail_stations_valid)
  
  center.ppp <- ppp(center.xy[,1], center.xy[,2], window = win, check = FALSE)
  rail.ppp <- ppp(rail.xy[,1], rail.xy[,2], window = win, check = FALSE)
  
  # ---- Compute nearest station index and distance ----
  nn.df <- nncross(center.ppp, rail.ppp)  # returns $which and $dist
  
  # ---- Attach route, distance, and ref date ----
  rail_attrs <- st_drop_geometry(rail_stations_valid) %>%
    mutate(station_id = row_number())
  
  blckgrp.dist <- blckgrp %>%
    mutate(
      center_dist_to_rail = nn.df$dist,
      nearest_station = rail_attrs$NAME[nn.df$which],
      nearest_route = rail_attrs$RAIL_LINE[nn.df$which],
      treat = case_when(
        center_dist_to_rail <= 805 + 375 * stdevs | lengths(st_intersects(blckgrp, rail_buffer)) > 0 ~ rail_attrs$year[nn.df$which],
        TRUE ~ 0
      )
    )
}

denver_controlled <- control.treatment(denver_bgs_2020,RTD_stations_commuter)

#-------- Join spatially-assigned data to regression var --------

RCMA_spatial <- lapply(RCMA_list,function(df){
  inner_join(df, denver_controlled, by = c("OriginID" = "GISJOIN"))
})

RCMA_regressionset <- bind_rows(RCMA_spatial) %>%
  mutate(
    lRCMA = log(RCMA),
    year = as.numeric(year),
    GEOID = as.numeric(GEOID)
  )
