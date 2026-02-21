#-------- Load libraries --------
library(sf)
library(spdep)
library(spatstat)
library(dplyr)
library(data.table)
library(did)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)
library(units)
library(fixest)
library(ggplot2)
library(viridis)

#-------- Load shapefiles --------

blckgrps_2010 <- st_read("nhgis0014_shape\\nhgis0014_shape\\nhgis0014_shapefile_tl2010_us_blck_grp_2010\\US_blck_grp_2010.shp")
blckgrps_2020 <- st_read("nhgis0014_shape\\nhgis0014_shape\\nhgis0014_shapefile_tl2020_us_blck_grp_2020\\US_blck_grp_2020.shp")

denver_bgs_2010 <- blckgrps_2010 %>%
  filter(
    STATEFP10 == "08",
    COUNTYFP10 %in% c("001","005","013","014","031","035","059")
  )

denver_bgs_2020 <- blckgrps_2020 %>%
  filter(
    STATEFP == "08",
    COUNTYFP %in% c("001","005","013","014","031","035","059")
  )

RTD <- st_read("Local_transit\\Local_transit.shp")
RTD_boundaries <- st_read("RTDBoundary.shp")
RTD_boundaries <- RTD_boundaries %>%
  st_transform(crs = st_crs(counties))
RTD_stations <- st_read("LightrailStations\\LightrailStations.shp")
highways <- st_read("Highways\\Highways.shp")

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

RTD_commuter <- RTD %>%
  filter(agency_id == "RTD",route_shor %in% c("A","B","G","N")) %>%
  mutate(year = case_when(
    route_shor == "A" ~ 2016,
    route_shor == "B" ~ 2016,
    route_shor == "G" ~ 2019,
    route_shor == "N" ~ 2020
  ))

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

county_fips <- c("08001","08005","08013","08014","08031","08035","08059")
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

qcew_yearly <- qcew_wage %>%
  bind_rows(.id = "source") %>%
  mutate(year = str_extract(source, "\\d{4}")) %>%
  split(.$year) %>%
  lapply(\(df) select(df, -source, -year))

#-------- Load and clean BDS Firm Entry Data --------

bds_sec <- read.csv("bds2023_st_cty_sec.csv")
bds_size <- read.csv("bds2023_st_cty_fzc.csv")

clean_bds <- function(bds) {
  bds %>%
   filter(
      year >= 2013,
      year <= 2023,
      st == 8,
      cty %in% c(1, 5, 13, 14, 31, 35, 59)
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
}

bds_sec_list <- bds_sec %>%
  clean_bds() %>%
  mutate(year = as.character(year)) %>%   # important for list indexing
  group_split(year)

bds_size <- bds_size %>%
  clean_bds() %>%
  mutate(year = as.character(year))

names(bds_sec_list) <- sapply(bds_sec_list, function(df) unique(df$year))

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
    bds_cty    = bds_sec_list[[y]]
  )
})

names(bg_births) <- names(lodes_st_list_long)

births_by_naics <- bg_births %>%
  bind_rows() %>%
  select(year, GISJOIN, NAICS2, estabs_entry_bg) %>%
  pivot_wider(
    names_from  = NAICS2,
    values_from = estabs_entry_bg,
    names_prefix = "estabs_entry_"
  )

births <- bg_births %>%
  bind_rows() %>%
  group_by(year,GISJOIN) %>%
  summarise(estabs_entry_bg = sum(estabs_entry_bg))

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
    rename(GISJOIN = DestinationID) %>%
    left_join(wages_df, by = "GISJOIN") %>%
    mutate(
      cost_adj = exp(Total_TravelTime * 0.01),
      RCMA_contrib = (bg_avg_weekly_wage / cost_adj)^3.4
    )
}

RCMA_list <- lapply(names(wages_destination), function(y) {
  
  od_mat <- matrix_list[[paste0(y,"CostMatrix")]]
  wages <- wages_destination[[y]]
  
  compute_RCMA(od_mat, wages) %>%
    group_by(OriginID) %>%
    summarise(RCMA = sum(RCMA_contrib, na.rm = TRUE)) %>%
    mutate(year = y) %>%
    rename(GISJOIN = OriginID)
})

names(RCMA_list) <- names(wages_destination)

#-------- Load ACS data and assign spatially --------

acs_codebooks <- list.files(
  "nhgis0014_csv\\nhgis0014_csv",
  pattern = "*_blck_grp_codebook\\.txt$",
  full.names = TRUE
)

parse_nhgis_codebook <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  tibble(line = lines) %>%
    mutate(
      # grab the census table ID, e.g., B08301
      census_table = ifelse(
        str_detect(line, "Source code:"),
        str_extract(line, "B\\d{5}"),
        NA_character_
      )
    ) %>%
    tidyr::fill(census_table) %>%
    # grab variable lines: anything that starts with letters/numbers and ends with colon
    filter(str_detect(line, "\\b[A-Z0-9]{3,5}[E]\\d{3}:")) %>%
    mutate(
      variable = str_extract(line, "\\b[A-Z0-9]{3,5}E\\d{3}"),
      prefix   = str_extract(variable, "^[A-Z0-9]{3,5}")
    ) %>%
    select(census_table, prefix, variable)
}

nhgis_prefix_lookup <- purrr::map_dfr(
  acs_codebooks,
  parse_nhgis_codebook
) %>%
  mutate(
    prefix = sub("E\\d+$", "", variable)
  ) %>%
  distinct(census_table, prefix)

tables_needed <- c(
  "B01003",  # total population
  "B08301",  # means of transportation
  "B23025",   # employment status
  "B25002",   # housing units by occupancy
  "B25077"   # median value of housing units
)

prefix_to_table <- nhgis_prefix_lookup %>%
  filter(census_table %in% tables_needed) %>%
  distinct(prefix, census_table) %>%
  group_by(prefix) %>%
  slice(1) %>%          # keep first match only
  ungroup() %>%
  deframe()

prefixes_needed <- names(prefix_to_table)

acs_2010 <- list.files(
  "nhgis0014_csv\\nhgis0014_csv",
  pattern = "^nhgis0014_ds2.*_201.*_blck_grp\\.csv$",
  full.names = TRUE
)

acs_2020 <- list.files(
  "nhgis0014_csv\\nhgis0014_csv",
  pattern = "^nhgis0014_ds2.*_202.*_blck_grp\\.csv$",
  full.names = TRUE
)

acs_2010_list <- setNames(
  lapply(acs_2010, read.csv, stringsAsFactors = FALSE),
  tools::file_path_sans_ext(basename(acs_2010))
)

acs_2020_list <- setNames(
  lapply(acs_2020, read.csv, stringsAsFactors = FALSE),
  tools::file_path_sans_ext(basename(acs_2020))
)

denver_acs_2010_list <- lapply(acs_2010_list,function(df){
  cols_keep <- names(df)[
    grepl("E\\d{3}$", names(df)) &
      sub("E\\d{3}$", "", names(df)) %in% names(prefix_to_table)
  ]
  df %>%
    filter(
      STATEA == 8,
      COUNTYA %in% c(1,5,13,14,31,35,59)
    ) %>%
    select(GISJOIN, any_of(cols_keep)) %>%
    rename_with(
      ~ purrr::map_chr(.x, function(col) {
        
        if (col == "GISJOIN") return(col)
        
        prefix <- sub("E\\d+$", "", col)
        cell   <- sub("^.*E", "E", col)
        
        if (!prefix %in% names(prefix_to_table)) {
          col
        } else {
          paste0(prefix_to_table[[prefix]], "_", cell)
        }
      })
    )
})

denver_acs_2020_list <- lapply(acs_2020_list, function(df) {
  cols_keep <- names(df)[
    grepl("E\\d{3}$", names(df)) &
      sub("E\\d{3}$", "", names(df)) %in% names(prefix_to_table)
  ]
  df %>%
    filter(
      STATEA == 8,
      COUNTYA %in% c(1,5,13,14,31,35,59)
    ) %>%
    select(GISJOIN, any_of(cols_keep)) %>%
    rename_with(
      ~ purrr::map_chr(.x, function(col) {
        
        if (col == "GISJOIN") return(col)
        
        prefix <- sub("E\\d+$", "", col)
        cell   <- sub("^.*E", "E", col)
        
        if (!prefix %in% names(prefix_to_table)) {
          col
        } else {
          paste0(prefix_to_table[[prefix]], "_", cell)
        }
      })
    )
})

#-------- Crosswalk 2010 data to 2020 geos --------

bg_xwalk <- st_intersection(
  denver_bgs_2010 %>% select(GISJOIN_2010 = GISJOIN),
  denver_bgs_2020 %>% select(GISJOIN_2020 = GISJOIN)
)

bg_xwalk <- bg_xwalk %>%
  mutate(overlap_area = st_area(.)) %>%
  group_by(GISJOIN_2010) %>%
  mutate(
    bg2010_area = sum(overlap_area),
    weight_2010_to_2020 = overlap_area / bg2010_area
  ) %>%
  ungroup()

acs2010_to_2020 <- lapply(denver_acs_2010_list,function(df){
  left_join(df, bg_xwalk, by = c("GISJOIN" = "GISJOIN_2010")) %>%
    mutate(across(
      where(is.numeric),
      ~ .x * weight_2010_to_2020
    ))
  })

denver_acs_2010_list_est <- lapply(acs2010_to_2020, function(df) {
  
  df %>%
    st_drop_geometry() %>%
    select(GISJOIN_2020, where(is.numeric)) %>%
    group_by(GISJOIN_2020) %>%
    summarise(
      across(everything(), sum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(GISJOIN = GISJOIN_2020)
  
})

#-------- Join and calculate ACS vars --------

get_year <- function(x) {
  as.integer(stringr::str_extract(x, "20\\d{2}"))
}

acs_2010_all <- purrr::imap_dfr(
  denver_acs_2010_list_est,
  ~ dplyr::mutate(
    .x,
    dplyr::across(where(~ inherits(.x, "units")), units::drop_units),
    year = get_year(.y)
  )
)

acs_2020_all <- purrr::imap_dfr(
  denver_acs_2020_list,
  ~ dplyr::mutate(
    .x,
    dplyr::across(where(~ inherits(.x, "units")), units::drop_units),
    year = get_year(.y)
  )
)

denver_acs <- bind_rows(acs_2010_all, acs_2020_all)

lodes <- lodes_list %>%
  bind_rows(.id = "source") %>%
  mutate(year = str_extract(source, "\\d{4}"))

employment <- denver_acs %>%
  mutate(year = as.character(year)) %>%
  select(year, GISJOIN, B23025_E001, B23025_E004) %>%
  inner_join(
    lodes %>% select(year, GISJOIN, c000), 
    by = c("GISJOIN","year")) %>%
  mutate(
      pct_empl = B23025_E004/B23025_E001, 
      jobs_per_100 = c000/B23025_E001 * 100,
      ljobs_per_100 = log1p(jobs_per_100)
    ) %>%
  rename(
    total_jobs = c000,
    total_working_age_pop = B23025_E001,
    employed_pop = B23025_E004
  )

housing <- denver_acs %>%
  select(year,GISJOIN,B01003_E001,B25002_E001,B25002_E002,B25002_E003) %>%
  mutate(
    year = as.character(year),
    lhousing_per_100 = log1p(B25002_E001 / B01003_E001 * 100),
    occupancy_rate = B25002_E002 / B25002_E001,
    vacancy_rate = B25002_E003 / B25002_E001,
  ) %>%
  rename(
    population = B01003_E001,
    housing_supply = B25002_E001,
    occupied_housing = B25002_E002,
    vacant_housing = B25002_E003,
    )

# -------- Calculate FCMA --------

employment_list <- employment %>%
  split(.$year) %>%
  lapply(\(df) select(df, -year))

compute_FCMA <- function(od_df, employment_df, rcma_df, epsilon = 1e-3) {
  od_df %>%
    rename(GISJOIN = OriginID) %>%   # <-- join on ORIGIN
    left_join(employment_df, by = "GISJOIN") %>%
    left_join(rcma_df, by = "GISJOIN") %>%
    mutate(
      cost_adj = exp(Total_TravelTime * 0.01),
      FCMA_contrib = cost_adj^(-3.4) * total_working_age_pop / RCMA
    )
}

FCMA_list <- lapply(names(employment_list), function(y) {
  
  od_mat <- matrix_list[[paste0(y, "CostMatrix")]]
  employ <- employment_list[[y]]
  rcma   <- RCMA_list[[y]]
  
  compute_FCMA(od_mat, employ, rcma) %>%
    group_by(DestinationID) %>%
    summarise(FCMA = sum(FCMA_contrib, na.rm = TRUE)) %>%
    mutate(year = y) %>%
    rename(GISJOIN = DestinationID)
})

names(FCMA_list) <- names(employment_list)

FCMA <- bind_rows(FCMA_list)

# -------- Point-Pattern Analysis using Centroids --------

control.treatment <- function(
    blckgrp,
    rail_stations,
    highways = NULL,
    phase = c("phase3_flow", "phase2_stock"),
    treat_radius = 1000,
    donut_inner = 1000,
    donut_outer = 2000,
    road_buffer = 1000,
    use_continuous = FALSE,
    use_poly_intersect = FALSE,
    target_crs = 26913
) {
  
  phase <- match.arg(phase)
  
  # ---- CRS ----
  blckgrp <- st_transform(blckgrp, target_crs)
  rail_stations <- st_transform(rail_stations, target_crs)
  if (!is.null(highways)) highways <- st_transform(highways, target_crs)
  
  # ---- Centroids ----
  center <- st_centroid(blckgrp)
  
  # ---- Nearest station distances ----
  dist_mat <- st_distance(center, rail_stations)
  nearest_idx <- apply(dist_mat, 1, which.min)
  nearest_dist <- apply(dist_mat, 1, min)
  rail_attrs <- st_drop_geometry(rail_stations)
  
  blckgrp <- blckgrp %>%
    mutate(
      center_dist_to_rail = as.numeric(nearest_dist),
      nearest_station = rail_attrs$NAME[nearest_idx],
      cohort_year = rail_attrs$year[nearest_idx],
      rail_exposure = -log1p(center_dist_to_rail)
    )
  
  # ============================================================
  # PHASE 3: FIRM ENTRY
  # ============================================================
  if (phase == "phase3_flow") {
    
    if (is.null(highways)) stop("Highways object required for phase3_flow.")
    
    # ---- CRS check
    highways <- st_make_valid(highways)
    highway_union <- st_union(highways)
    
    # ---- Distance to nearest highway
    nearest_highway_dist <- as.numeric(st_distance(center, highway_union))
    near_highway <- nearest_highway_dist <= road_buffer  # restrict sample
    
    # ---- Rail distance already computed
    rail_dist <- blckgrp$center_dist_to_rail
    
    # ---- Treatment: within radius of station
    treated <- rail_dist <= treat_radius
    
    # ---- Donut: 1–2 km from rail
    rail_donut <- rail_dist > donut_inner & rail_dist <= donut_outer
    
    # ---- Control: within highway sample, NOT treated, NOT in donut
    control <- near_highway & !treated & !rail_donut
    
    # ---- Mark in-sample: only block groups near highways
    in_sample <- near_highway
    
    # ---- Assign to dataset
    blckgrp <- blckgrp %>%
      mutate(
        treated = treated,
        control = control,
        in_sample = in_sample,
        cohort_year = if_else(treated, cohort_year, 0)
      )
  }
  
  # ============================================================
  # PHASE 2: EMPLOYMENT / FIRM STOCKS
  # ============================================================
  if (phase == "phase2_stock") {
    blckgrp <- blckgrp %>%
      mutate(
        treated = center_dist_to_rail <= treat_radius,
        donut = center_dist_to_rail > donut_inner & center_dist_to_rail <= donut_outer,
        control_general = center_dist_to_rail > donut_outer,
        in_sample = !donut
      )
  }
  
  # ============================================================
  # Continuous-only option
  # ============================================================
  if (use_continuous) {
    blckgrp <- blckgrp %>%
      mutate(
        treated = NULL,
        control = NULL
      )
  }
  
  blckgrp <- blckgrp %>%
    filter(in_sample == TRUE)
  
  return(blckgrp)
}

denver_controlled <- control.treatment(
  denver_bgs_2020,
  RTD_stations_commuter,
  highways = highways
  )

#-------- Join spatially-assigned data to regression vars --------

denver_regression_set <- denver_controlled %>%
  left_join(births, by = "GISJOIN") %>%
  left_join(births_by_naics, by = c("GISJOIN","year")) %>%
  left_join(employment, by = c("GISJOIN","year")) %>%
  left_join(housing, by = c("GISJOIN","year")) %>%
  arrange(GISJOIN, year) %>%
  mutate(
    year = as.numeric(year),
    event_time = year - cohort_year,      # 0 = first treated year
    treated_post = if_else(year >= cohort_year & treated == TRUE, 1, 0)
  )


#-------- Run sunab regressions --------

birth_outcomes <- c("estabs_entry_bg","estabs_entry_11","estabs_entry_21",
              "estabs_entry_22","estabs_entry_23","estabs_entry_31-33","estabs_entry_42",
              "estabs_entry_44-45","estabs_entry_48-49","estabs_entry_51",
              "estabs_entry_52","estabs_entry_53","estabs_entry_54","estabs_entry_55",
              "estabs_entry_56","estabs_entry_61","estabs_entry_62","estabs_entry_71",
              "estabs_entry_72","estabs_entry_81")

birth_models <- map(
  birth_outcomes,
  ~ feols(
    as.formula(paste0("log1p(`", .x, "`) ~ sunab(cohort_year,year) | GISJOIN + year")),
    data = denver_regression_set,
    cluster = ~ GISJOIN
  )
)

names(birth_models) <- birth_outcomes

employment_model <- feols(
  ljobs_per_100 ~ sunab(cohort_year,year) | GISJOIN + year,
  data = denver_regression_set,
  cluster = ~ GISJOIN
)

# -------- Plot birth regression models --------

get_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

# Extract ATT for all outcomes
att_table <- map_df(names(birth_models), function(x) {
  
  # Get summary aggregated ATT
  att_sum <- summary(birth_models[[x]], agg = "ATT")
  
  # Extract estimate and SE
  tibble(
    outcome = x,
    ATT = att_sum$coeftable["ATT", "Estimate"],
    SE = att_sum$coeftable["ATT", "Std. Error"],
    t_value = att_sum$coeftable["ATT", "t value"],
    p_value = att_sum$coeftable["ATT", "Pr(>|t|)"]
  )
}) %>%
  arrange(p_value)

extract_sunab_dynamic <- function(model, outcome_name) {
  
  sum_model <- summary(model)
  coefs <- sum_model$coeftable
  
  # Keep dynamic terms
  keep <- grepl("^year::", rownames(coefs))
  coefs <- coefs[keep, , drop = FALSE]
  
  if (nrow(coefs) == 0) return(NULL)
  
  tibble(
    event_time = as.numeric(sub("year::", "", rownames(coefs))),
    estimate   = coefs[, "Estimate"],
    se         = coefs[, "Std. Error"],
    outcome    = outcome_name
  )
}

# Extract across all models
all_coefs <- map_df(names(birth_models), function(x) {
  extract_sunab_dynamic(birth_models[[x]], x)
})

# Add confidence intervals
all_coefs <- all_coefs %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )

# Select sectors for plotting
sector_labels <- c(
  "estabs_entry_bg" = "All Sectors",
  "estabs_entry_23" = "Construction",
  "estabs_entry_42" = "Wholesale Trade",
  "estabs_entry_52" = "Finance & Insurance",
  "estabs_entry_62" = "Health Care & Social Assistance"
)

att_table_plot <- att_table %>%
  filter(outcome %in% names(sector_labels)) %>%
  mutate(
    sector_name = sector_labels[outcome],
    legend_label = paste0(
      sector_name,
      "\nATT = ", round(ATT, 3),
      ", p = ", format.pval(p_value, digits = 2, eps = .001)
    )
  )

# Make sure plot_data contains the right sectors
plot_data <- all_coefs %>%
  filter(outcome %in% names(sector_labels)) %>%
  left_join(att_table_plot %>% select(outcome, legend_label),
            by = "outcome")

agg_info <- att_table_plot %>%
  filter(outcome == "estabs_entry_bg")

stars <- get_stars(agg_info$p_value)

agg_label <- paste0(
  "All Sectors",
  "\nATT = ", round(agg_info$ATT, 3),
  stars,
  "  (SE = ", round(agg_info$SE, 3), ")"
)

# Data for aggregate only
plot_all <- plot_data %>%
  filter(outcome == "estabs_entry_bg") %>%
  mutate(series = agg_label)

ggplot(plot_all,
       aes(x = event_time,
           y = estimate,
           color = series)) +
  
  geom_rect(aes(xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.2,
            inherit.aes = FALSE) +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = series,
                  group = 1),
              alpha = 0.2,
              color = NA) +
  
  geom_line(size = 2) +
  geom_point(size = 4) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_color_manual(
    values = setNames("steelblue", agg_label)
  ) +
  
  scale_fill_manual(
    values = setNames("#ae0c00", agg_label),
    guide = "none"
  ) +
  
  labs(
    x = "Years Relative to Treatment",
    y = "Dynamic ATT",
    title = "Event Study: Log Establishment Births (All Sectors)",
    color = NULL
  ) +
  
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(size = 48, face = "bold"),
    axis.title = element_text(size = 36),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 26),
    panel.grid.minor = element_blank()
  )

plot_sectors <- plot_data %>%
  filter(outcome != "estabs_entry_bg")

sector_att_info <- att_table %>%
  filter(outcome %in% c(
    "estabs_entry_23",
    "estabs_entry_42",
    "estabs_entry_52",
    "estabs_entry_62"
  )) %>%
  mutate(
    sector_name = sector_labels[outcome],
    stars = get_stars(p_value),
    legend_label = paste0(
      sector_name,
      "\nATT = ", round(ATT, 3),
      stars,
      "  (SE = ", round(SE, 3), ")"
    )
  ) %>%
  select(outcome, legend_label)

plot_sectors <- all_coefs %>%
  filter(outcome %in% sector_att_info$outcome) %>%
  mutate(outcome = as.character(outcome)) %>%
  left_join(sector_att_info, by = "outcome")

sector_colors <- c(
  "estabs_entry_23" = "#1f77b4",
  "estabs_entry_42" = "#ff7f0e",
  "estabs_entry_52" = "#2ca02c",
  "estabs_entry_62" = "#d62728"
)

# Map colors to full legend labels
color_map <- sector_colors[sector_att_info$outcome]
names(color_map) <- sector_att_info$legend_label

# Dashed confidence intervals
ggplot(plot_sectors,
       aes(x = event_time,
           y = estimate,
           color = legend_label)) +
  
  geom_rect(aes(xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.2,
            inherit.aes = FALSE) +
  
  geom_line(size = 2) +
  geom_point(size = 4) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = color_map) +
  
  labs(
    x = "Years Relative to Treatment",
    y = "Dynamic ATT",
    title = "Event Study: Log Establishment Births by Sector",
    color = NULL
  ) +
  
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(size = 48, face = "bold"),
    axis.title = element_text(size = 36),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 24),
    panel.grid.minor = element_blank()
  )

# -------- Plot employment regression model --------

employment_coefs <- broom::tidy(employment_model) %>%
  filter(grepl("year::", term)) %>%
  mutate(
    event_time = as.numeric(sub(".*::", "", term)),
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

# Extract overall ATT from the feols model
att_employment <- summary(employment_model, agg = "ATT")$coeftable

# Pull estimate, SE, and p-value
att_est <- att_employment["ATT", "Estimate"]
att_se  <- att_employment["ATT", "Std. Error"]
att_p   <- att_employment["ATT", "Pr(>|t|)"]

# Create annotation label with stars and SE
att_label <- paste0(
  "ATT = ", round(att_est, 3),
  get_stars(att_p),
  "  (SE = ", round(att_se, 3), ")"
)

# Plot
ggplot(employment_coefs, aes(x = event_time, y = estimate)) +
  
  # Shaded post-treatment window
  geom_rect(aes(xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf),
            fill = "lightgray", alpha = 0.1, inherit.aes = FALSE) +
  
  # Confidence band
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#ae0c00", alpha = 0.2) +
  
  # Line for point estimates
  geom_line(color = "steelblue", size = 5) +
  
  # Points
  geom_point(color = "steelblue", size = 7) +
  
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # Add dynamic annotation
  annotate("text",
           x = max(employment_coefs$event_time) - 1,  
           y = max(employment_coefs$upper) + 0.05,
           label = att_label,
           size = 10,
           color = "black") +
  
  # Labels
  labs(
    x = "Years Relative to Treatment",
    y = "Dynamic ATT: log(Jobs per 100 residents)",
    title = "Event Study: Local Employment Density"
  ) +
  
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(size = 48, face = "bold"),
    axis.title = element_text(size = 36),
    axis.text = element_text(size = 28),
    panel.grid.minor = element_blank()
  )

# -------- Plot study area + Treatment/Control Groups --------

counties <- st_read("nhgis0016_shape\\nhgis0016_shapefile_tl2023_us_county_2023\\US_county_2023.shp")
counties <- counties %>%
  filter(STATEFP == "08",
         COUNTYFP %in% c("001","005","013","014","031","035","059"))

# Prepare treatments as a plain data frame
treatments <- denver_controlled %>%
  st_drop_geometry() %>%
  select(GISJOIN, treated, control, cohort_year)

# 1️⃣ Prepare plot_year_factor
blckgrp_plot <- denver_bgs_2020 %>%
  left_join(
    denver_controlled %>% st_drop_geometry() %>% select(GISJOIN, treated, control, cohort_year),
    by = "GISJOIN"
  ) %>%
  mutate(
    group = case_when(
      treated ~ "Treated",
      control ~ "Control",
      TRUE ~ "Excluded"
    ),
    plot_year_factor = case_when(
      group == "Control"  ~ "Control",
      group == "Excluded" ~ "Excluded",
      TRUE ~ as.character(cohort_year)
    )
  )

# 2️⃣ Factor levels for consistent ordering
treated_years <- sort(unique(blckgrp_plot$cohort_year[blckgrp_plot$group == "Treated"]))
blckgrp_plot$plot_year_factor <- factor(
  blckgrp_plot$plot_year_factor,
  levels = c("Control", as.character(treated_years), "Excluded")
)

# 3️⃣ Fill colors matching factor levels
fill_colors <- c(
  "Control" = "steelblue",             # deep blue
  setNames(viridis::viridis(length(treated_years), option = "C"), as.character(treated_years)), # treated
  "Excluded" = "grey90"               # light grey
)

# 4️⃣ Legend labels with N
legend_labels <- blckgrp_plot %>%
  st_drop_geometry() %>%
  group_by(plot_year_factor) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(factor(plot_year_factor, levels = levels(blckgrp_plot$plot_year_factor))) %>%
  mutate(label = paste0(plot_year_factor, " (N=", N, ")")) %>%
  { setNames(.$label, .$plot_year_factor) }

# 5️⃣ Plot

RTD_commuter$plot_year_factor <- factor(RTD_commuter$year)

ggplot() +
  geom_sf(data = blckgrp_plot, aes(fill = plot_year_factor)) +
  geom_sf(data = counties, fill = NA, color = "grey40", linetype = "dashed", size = 1) +
  geom_sf(
    data = RTD_commuter,
    aes(color = plot_year_factor),
    size = 1
  ) +
  geom_sf(data = RTD_stations_commuter, color = "#ae0c00", size = 1, shape = 21, fill = "#ae0c00") +
  scale_fill_manual(
    name = "Group / Cohort (N)",
    values = fill_colors,
    labels = legend_labels
  ) +
  scale_color_manual(
    name = "Group / Cohort (N)",
    values = fill_colors,
    labels = legend_labels
  )


# -------- Commuter rail context map --------

county_centroids <- st_centroid(counties)

ggplot() +
  # County background with dashed borders
  geom_sf(data = counties, fill = "grey90", color = "grey40", linetype = "dashed", size = 1) +
  
  # Commuter rail lines
  geom_sf(data = RTD_commuter, aes(color = paste0(route_shor," (",year,")")), size = 2, show.legend = "line") +
  
  # Rail stations
  geom_sf(data = RTD_stations_commuter, color = "#ae0c00", size = 4, shape = 21, fill = "#ae0c00") +
  
  geom_sf_text(
    data = county_centroids,
    aes(label = NAME),  # replace NAME with your county name column
    size = 10,           # adjust for poster
    color = "black",
    nudge_x = -3000,
    nudge_y = -2000
  ) +
  
  # Poster-friendly theme
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(size = 48, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 28)
  ) +
  
  # Labels
  labs(
    title = "Map of Study Area with Commuter Rail Lines",
    color = "Route"  # legend title for lines
  ) +
  
  # Optional: nice color palette for routes
  scale_color_viridis_d(option = "C")   # Plasma palette for line colors

# -------- Plot Housing Model ----------

# Dynamic coefficients
housing_coefs <- broom::tidy(housing_supply_model) %>%
  filter(grepl("year::", term)) %>%
  mutate(
    event_time = as.numeric(sub(".*::", "", term)),
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error
  )

# Extract ATT table
att_housing <- summary(housing_supply_model, agg = "ATT")$coeftable

att_est <- att_housing["ATT", "Estimate"]
att_se  <- att_housing["ATT", "Std. Error"]
att_p   <- att_housing["ATT", "Pr(>|t|)"]

att_label <- paste0(
  "ATT = ", round(att_est, 3),
  get_stars(att_p),
  "  (SE = ", round(att_se, 3), ")"
)

ggplot(housing_coefs, aes(x = event_time, y = estimate)) +
  
  # Post-treatment shading
  geom_rect(aes(xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf),
            fill = "lightgray", alpha = 0.1,
            inherit.aes = FALSE) +
  
  # 95% CI
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#ae0c00", alpha = 0.2) +
  
  # Line + points
  geom_line(color = "steelblue", size = 5) +
  geom_point(color = "steelblue", size = 7) +
  
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # ATT annotation
  annotate("text",
           x = max(housing_coefs$event_time) - 1,
           y = max(housing_coefs$upper) + 0.05,
           label = att_label,
           size = 10,
           fontface = "bold") +
  
  labs(
    x = "Years Relative to Treatment",
    y = "Dynamic ATT: log(Housing Units)",
    title = "Event Study: Housing Supply"
  ) +
  
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(size = 48, face = "bold"),
    axis.title = element_text(size = 36),
    axis.text = element_text(size = 28),
    panel.grid.minor = element_blank()
  )

# RTD context map

ggplot() +
  geom_sf(data = counties,
          fill = "grey90",
          color = "black",
          linetype = "dashed",
          size = 1) +
  
  geom_sf(data = RTD_boundaries,
          fill = "blue",
          alpha = 0.1,
          color = "blue",
          size = 1.5) +
  
  geom_sf_text(
    data = county_centroids,
    aes(label = NAME),  # replace NAME with your county name column
    size = 10,           # adjust for poster
    color = "black",
    nudge_x = -3000,
    nudge_y = -2000
  ) +
  
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(size = 48, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Study Counties with RTD Service Area Overlay")
