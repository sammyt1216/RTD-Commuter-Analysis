library(sf)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)

parse_nhgis_codebook <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  tibble(line = lines) %>%
    mutate(
      # grab the census table ID
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

clean_acs_df <- function(df, prefix_to_table, valid_gisjoin) {
  
  prefix_names <- names(prefix_to_table)
  
  cols_keep <- names(df)[
    grepl("E\\d{3}$", names(df)) &
      sub("E\\d{3}$", "", names(df)) %in% prefix_names
  ]
  
  # Validation at the top, before the pipe
  found_prefixes <- unique(sub("E\\d{3}$", "", cols_keep))
  
  if (length(found_prefixes) == 0) {
    stop(
      "clean_acs_df: no expected prefixes found in df.\nLooking for: ",
      paste(prefix_names, collapse = ", "),
      "\nActual E-variable prefixes in df: ",
      paste(unique(sub("E\\d{3}$", "", names(df)[grepl("E\\d{3}$", names(df))])), collapse = ", ")
    )
  }
  
  # Return the pipe result explicitly
  df |>
    dplyr::filter(GISJOIN %in% valid_gisjoin) |>
    dplyr::select(GISJOIN, dplyr::any_of(cols_keep)) |>
    dplyr::rename_with(function(cols) {
      ifelse(
        cols == "GISJOIN",
        cols,
        {
          prefix <- sub("E\\d+$", "", cols)
          cell   <- sub("^.*E", "E", cols)
          ifelse(
            prefix %in% prefix_names,
            paste0(prefix_to_table[prefix], "_", cell),
            cols
          )
        }
      )
    }) |>
    dplyr::mutate(dplyr::across(
      dplyr::starts_with("B"),
      ~ dplyr::if_else(is.na(.x) | as.numeric(.x) < 0, NA_real_, as.numeric(.x))
    ))
}

NAICS2_CODES <- c("11","21","22","23","31-33","42","44-45",
                  "48-49","51","52","53","54","55","56",
                  "61","62","71","72","81","92")

geoid_to_gisjoin <- function(id) {
  id <- as.character(id)
  
  # Truncate block to BG if 15 chars
  id <- dplyr::case_when(
    nchar(id) == 15 ~ substr(id, 1, 12),
    TRUE            ~ id
  )
  
  # Now handle the two GEOID formats separately
  dplyr::case_when(
    # 11-char LODES format: leading zero dropped (8 CCC TTTT BBB)
    nchar(id) == 11 ~ id |>
      stringr::str_replace("^8", "G080") |>
      (\(x) stringr::str_c(stringr::str_sub(x, 1, 7), "0", stringr::str_sub(x, 8)))(),
    
    # 12-char standard GEOID format: (08 CCC TTTTTT BBB)
    nchar(id) == 12 ~ id |>
      stringr::str_replace("^08", "G080") |>
      (\(x) stringr::str_c(stringr::str_sub(x, 1, 7), "0", stringr::str_sub(x, 8)))(),
    
    TRUE ~ NA_character_
  )
}

geoid_to_gisjoin2 <- function(id) {
  id <- as.character(id)
  id <- stringr::str_pad(id, width = 15, side = "left", pad = "0")  # restore leading zero
  id <- substr(id, 1, 12)  # truncate block to BG
  
  id |>
    stringr::str_replace("^08", "G080") |>
    (\(x) stringr::str_c(
      stringr::str_sub(x, 1, 7), "0",
      stringr::str_sub(x, 8)
    ))()
}

process_lodes_file <- function(df, file_name, geos) {
  # 1. Fix GISJOIN for NHGIS-style
  df <- df %>%
    mutate(
      GISJOIN = geoid_to_gisjoin2(w_geocode)
    )
  
  # 2. Rename c000 to total_jobs
  if ("C000" %in% colnames(df)) df <- df %>% rename(total_jobs = C000)
  
  # 3. Rename NAICS columns
  NAICS <- paste0("CNS", str_pad(1:20, 2, pad = "0"))
  matched <- match(NAICS, colnames(df))
  colnames(df)[matched[!is.na(matched)]] <- paste0("NAICS", NAICS2_CODES)
  
  # 4. Add year column from filename and drop IDs except for GISJOIN
  df <- df  %>%
    mutate(year = as.integer(str_extract(file_name, "20\\d{2}"))) %>%
    filter(GISJOIN %in% geos) %>% 
    dplyr::select(-w_geocode, -createdate) %>%
    dplyr::group_by(GISJOIN, year) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE),
                     .groups = "drop")
  
  df
}

process_county_lodes_file <- function(df, file_name) {
  # 1. Format id
  df <- df %>%
    mutate(
      county_fips = str_replace(id, "^(.{1})", "\\10")
    )
  
  # 2. Rename c000 to total_jobs
  if ("c000" %in% colnames(df)) df <- df %>% 
      rename(total_jobs_county = c000)
  
  # 3. Rename NAICS columns
  NAICS <- paste0("cns", str_pad(1:20, 2, pad = "0"))
  matched <- match(NAICS, colnames(df))
  colnames(df)[matched[!is.na(matched)]] <- paste0("NAICS", NAICS2_CODES, "_county")
  
  # 4. Add year column from filename
  df <- df %>% 
    mutate(year_county = as.integer(str_extract(file_name, "20\\d{2}"))) %>%
    select(year_county,
           county_fips,
           total_jobs_county,
           all_of(paste0("NAICS", NAICS2_CODES, "_county")))
  df
}

clean_bds <- function(df){
  df <- df %>%
    filter(
      year >= 2013,
      st == 8,
      if("cty" %in% colnames(df)) cty %in% c(1,5,13,14,31,35,59,123) else TRUE
    ) %>%
    mutate(
      firms = as.numeric(firms),
      estabs = as.numeric(estabs),
      estabs_entry = as.numeric(estabs_entry),
      estabs = ifelse(is.na(estabs), 1e-8, estabs),
      county_fips = if("cty" %in% colnames(df)) {
        str_pad(cty, width = 5, side = "left", pad = "0") %>%
          { paste0("8", str_sub(., -4)) }
      } else NA
    )
  
  df
}

process_qcew <- function(qcew_files) {
  
  fips_keep  <- c("08001","08005","08013","08014",
                  "08031","08035","08059","08123")
  naics_keep <- NAICS2_CODES
  
  purrr::map(
    qcew_files,
    function(f) {
      data.table::fread(
        f,
        select = c("year", "area_fips", "industry_code", "qtr",
                   "month1_emplvl", "month2_emplvl", "month3_emplvl",
                   "total_qtrly_wages"),
        colClasses = list(character = c("area_fips", "industry_code"))
      ) |>
        dplyr::filter(
          area_fips     %in% fips_keep,
          industry_code %in% naics_keep,
          qtr           == 2
        ) |>
        dplyr::group_by(year, area_fips, industry_code) |>
        dplyr::summarise(
          avg_empl        = sum(month1_emplvl, month2_emplvl, month3_emplvl) / 3,
          total_wages     = sum(total_qtrly_wages),
          avg_weekly_wage = (total_wages / avg_empl) / 13,
          .groups         = "drop"
        ) |>
        dplyr::mutate(
          area_fips = paste0("80", stringr::str_sub(area_fips, 3, 5))
        )
    }
  ) |>
    purrr::set_names(
      tools::file_path_sans_ext(basename(qcew_files))
    )
}

# ------- OD Matrix Filtering for Inclusive Value of TT -------

parse_mode_scenario <- function(od_long, label = NULL) {
  
  # Use the label argument if provided, otherwise parse from existing scenario col
  x <- if (!is.null(label)) label else od_long$scenario
  
  mode_lookup <- c(
    "drive"        = "drive",
    "walk_transit" = "transit",
    "walk_alone"   = "walk",
    "cycle_alone"  = "bike"
  )
  
  # Match the longest mode key present in the label string
  matched_mode <- map_chr(x, function(s) {
    hits <- names(mode_lookup)[str_detect(s, names(mode_lookup))]
    if (length(hits) == 0) NA_character_
    else mode_lookup[hits[which.max(nchar(hits))]]  # longest match wins
  })
  
  # Everything after the mode prefix (and trailing underscore) is the scenario
  matched_scenario <- str_remove(x, paste0(names(mode_lookup), collapse = "|")) |>
    str_remove("^_|_$") |>
    na_if("")
  # e.g. "walk_transit_base" -> mode = "transit", scenario = "base"
  # e.g. "drive_2013"        -> mode = "drive",   scenario = "2013"
  # e.g. "walk_alone"        -> mode = "walk",     scenario = NA
  
  od_long |>
    mutate(
      mode     = matched_mode,
      scenario = coalesce(matched_scenario, "base")
    )
}

# Mode-specific ceilings reflecting realistic maximum commute times (minutes).
# Rationale: beyond these thresholds, the mode is effectively unavailable for
# that OD pair — retaining them would bias t_bar downward for other modes
# that serve only the short-haul pairs.
default_time_ceilings <- function() {
  c(
    walk    =  60,   # ~3 miles at avg pace; beyond this walk is not a real option
    bike    =  90,   # ~15 miles; aggressive but captures e-bike users
    transit = 120,   # accounts for multi-transfer trips in sparse networks
    drive   =  90
  )
}

apply_time_ceiling <- function(od_long, ceilings = default_time_ceilings()) {
  ceiling_df <- enframe(ceilings, name = "mode", value = "ceiling")
  
  od_long |>
    left_join(ceiling_df, by = "mode") |>
    filter(TotalTravelTime <= ceiling) |>
    select(-ceiling)
}

process_lodes_OD <- function(lodes_od, study_gisjoin) {
  lodes_od %>%
      mutate(
        from_GISJOIN = geoid_to_gisjoin2(h_geocode),
        to_GISJOIN   = geoid_to_gisjoin2(w_geocode)
      ) %>%
      filter(from_GISJOIN %in% study_gisjoin,
             to_GISJOIN %in% study_gisjoin) %>%
      group_by(from_GISJOIN,to_GISJOIN) %>%
      summarize(SE01 = sum(SE01),
                SE02 = sum(SE01),
                SE03 = sum(SE03),
                total_flow = sum(S000))
}
