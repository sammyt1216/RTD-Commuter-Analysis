library(targets)
library(tarchetypes)
library(here)
library(data.table)
library(alpaca)
library(cppRouting)
library(dplyr)
library(furrr)
library(ggplot2)
library(sf)
library(lwgeom)
library(gravity)
library(purrr)
library(here)
library(stringr)
library(modelsummary)
library(tibble)
library(tidyr)
library(units)
library(fixest)

tar_source()

# Set target-specific options such as packages:
tar_option_set(
  packages = c(
    "data.table", "cppRouting", "dplyr", "ggplot2",
    "purrr", "stringr", "sf", "tidyr",
    "modelsummary", "lwgeom", "fixest", "dplyr", "tidyr")
)

# Network parameters
CRS       <- 26954

BOARD_PENALTIES <- list(
  rail_board  = 2.0,
  rail_alight = 1.0,
  bus_board   = 0.5,
  bus_alight  = 0.25
)

GTFS_DIRS <- c(
  "data\\GTFS\\RTD_Purchased_Bus",
  "data\\GTFS\\RTD_Purchased_CR",
  "data\\GTFS\\RTD_Operated_Bus",
  "data\\GTFS\\RTD_Operated_CR",
  "data\\GTFS\\RTD_Operated_LR"
)

# B and N lines for counterfactual speed reference
CF_RAIL_LINES <- c("113B", "117N")

SCENARIOS <- tibble::tribble(
  ~scenario,        ~max_cohort_year,
  "pre_2016",       0,
  "post_2016",      2016,
  "post_2017",      2017,
  "post_2019",      2019,
  "post_2020",      2020,
  "counterfactual", 2045    # includes everything including cf lines
)

OD_SCENARIOS <- tibble::tribble(
  ~scenario,        ~transit_od,
  "pre_2016",       "od_walk_transit_pre_2016",
  "post_2016",      "od_walk_transit_post_2016",
  "post_2017",      "od_walk_transit_post_2017",
  "post_2019",      "od_walk_transit_post_2019",
  "post_2020",      "od_walk_transit_post_2020",
  "counterfactual", "od_walk_transit_counterfactual"
)

IV_SCENARIOS <- tibble::tribble(
  ~scenario,        ~p_car_target,
  "pre_2016",       "p_car_2013",
  "post_2016",      "p_car_2013",
  "post_2017",      "p_car_2013",
  "post_2019",      "p_car_2023",
  "post_2020",      "p_car_2023",
  "counterfactual", "p_car_2023"
)

# End this file with a list of target objects.
list(
  # -------- Load shapefiles --------
  
  # Load 2020 bgs
  tar_target(map2020_raw,st_read("data\\nhgis0014_shapefile_tl2020_us_blck_grp_2020\\US_blck_grp_2020.shp")),
  
  # Load 2010 bgs
  tar_target(map2010_raw,st_read("data\\nhgis0014_shapefile_tl2010_us_blck_grp_2010\\US_blck_grp_2010.shp")),
  
  # Load highways
  tar_target(highways_raw,st_read("data\\Highways\\Highways.shp")),
  
  # Load RTD boundary
  tar_target(RTD_raw,st_read("data\\RTDBoundary\\RTDBoundary.shp")),
  
  # Load RTD rail stations
  tar_target(RTD_stations_raw,st_read("data\\LightrailStations\\LightrailStations.shp")),
  
  # Load counterfactuals
  tar_target(counterfactual_sf,st_read("data\\Network_Shapefiles\\Counterfactual_CR.shp")),
  
  # Filter RTD stations to commuter
  tar_target(RTD_commuter_raw, RTD_stations_raw %>%
               filter(str_detect(RAIL_LINE,"A|B|G|N")) %>%
               mutate(
                 year = case_when(
                   str_detect(RAIL_LINE, "A") ~ 2016,
                   str_detect(RAIL_LINE, "N") & !str_detect(RAIL_LINE, "A") ~ 2020,
                   CITY == "Westminster" ~ 2016,
                   TRUE ~ 2019
                 )
               )),
  
  # ------- Filter to sample -------
  
  # Project all to EPSG 26954
  tar_target(map2020_sf, project_to_stateplane(map2020_raw)),
  tar_target(map2010_sf, project_to_stateplane(map2010_raw)),
  tar_target(highways_sf, project_to_stateplane(highways_raw)),
  tar_target(RTD_sf, project_to_stateplane(RTD_raw)),
  tar_target(RTD_commuter_sf, project_to_stateplane(RTD_commuter_raw)),
  
  # Filter bgs to Colorado
  tar_target(
    bgs_co_2020,
    dplyr::filter(map2020_sf, STATEFP == "08")
  ),
  
  tar_target(
    bgs_co_2010,
    dplyr::filter(map2010_sf, STATEFP10 == "08")
  ),
  
  # Pre-filter BGs that touch RTD
  tar_target(
    bgs_touch_rtd_2020,
    {
      bgs <- bgs_co_2020
      rtd_valid <- sf::st_make_valid(RTD_sf)        
      bgs[ lengths(sf::st_intersects(bgs, rtd_valid)) > 0, ]
      }
    ),
  
  # Intersect BGs with RTD to compute share
  tar_target(
    bgs_intersect_rtd,
    {
      bgs <- bgs_touch_rtd_2020
      rtd_valid <- sf::st_make_valid(RTD_sf)
      bgs <- dplyr::mutate(bgs, bg_area = sf::st_area(geometry))
      bgs_int <- sf::st_intersection(bgs, rtd_valid)
      bgs_share <- bgs_int %>%
        dplyr::mutate(int_area = sf::st_area(geometry)) %>%
        dplyr::group_by(GISJOIN) %>%
        dplyr::summarise(int_area = sum(int_area), .groups = "drop") %>%
        sf::st_drop_geometry()
      bgs <- bgs %>%
        dplyr::left_join(bgs_share, by = "GISJOIN") %>%
        dplyr::mutate(
          int_area = tidyr::replace_na(int_area, units::set_units(0, "m^2")),
          share = as.numeric(int_area / bg_area)
        ) %>%
        dplyr::filter(share >= 0.5) %>%
        dplyr::select(-bg_area, -int_area, -share)
      bgs
      }
    ),
  
  # ------ Load ACS data ------
  
  # Return list of ACS codebook files
  tar_target(
    acs_codebooks,
    list.files(
      "data/nhgis0019_csv",
      pattern = "_blck_grp_codebook\\.txt$",
      full.names = TRUE
    )
  ),
  
  # Load codebooks
  tar_target(
    nhgis_prefix_lookup,
    purrr::map_dfr(acs_codebooks, parse_nhgis_codebook) |>
      dplyr::mutate(prefix = sub("E\\d+$", "", variable)) |>
      dplyr::distinct(census_table, prefix)
  ),
  
  # Look up table prefixes in codebooks
  tar_target(
    prefix_to_table,
    {
      tables_needed <- c("B01003","B23025","B25003","B19013","B25064") # Total population, employment, housing by occupancy, MHI, rent
      
      out <- nhgis_prefix_lookup |>
        dplyr::filter(census_table %in% tables_needed) |>
        dplyr::distinct(prefix, census_table) |>
        dplyr::arrange(prefix) |>
        dplyr::group_by(prefix) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
      
      stats::setNames(out$census_table, out$prefix)
    }
  ),
  
  # Load ACS files
  tar_target(
    acs_2010_files,
    list.files(
      "data/nhgis0019_csv",
      pattern = "^nhgis0019_ds2.*_201.*_blck_grp\\.csv$",
      full.names = TRUE
    )
  ),
  
  tar_target(
    acs_2020_files,
    list.files(
      "data/nhgis0019_csv",
      pattern = "^nhgis0019_ds2.*_202.*_blck_grp\\.csv$",
      full.names = TRUE
    )
  ),
  
  tar_target(
    acs_2010_list,
    {
      files <- acs_2010_files
      
      purrr::map(
        files,
        ~ read.csv(.x, stringsAsFactors = FALSE)
      ) |>
        purrr::set_names(
          tools::file_path_sans_ext(basename(files))
        )
    }
  ),
  
  tar_target(
    acs_2020_list,
    {
      files <- acs_2020_files
      
      purrr::map(
        files,
        ~ read.csv(.x, stringsAsFactors = FALSE)
      ) |>
        purrr::set_names(
          tools::file_path_sans_ext(basename(files))
        )
    }
  ),
  
  # Create vectors for BG identifiers
  tar_target(
    study_gisjoin,
    bgs_intersect_rtd$GISJOIN
  ),
  
  tar_target(
    gisjoin2010,
    bgs_co_2010$GISJOIN
  ),
  
  # Clean datasets
  tar_target(
    denver_acs_2020_list,
    lapply(
      acs_2020_list,
      clean_acs_df,
      prefix_to_table = prefix_to_table,
      valid_gisjoin = study_gisjoin
    )
  ),
  
  tar_target(
    denver_acs_2010_list,
    lapply(
      acs_2010_list,
      clean_acs_df,
      prefix_to_table = prefix_to_table,
      valid_gisjoin = gisjoin2010
    )
  ),
  
  # Crosswalk 2010 boundaries to 2020 by area
  tar_target(
    bg_xwalk,
    {
      bgs_2010 <- bgs_co_2010 |>
        dplyr::select(GISJOIN_2010 = GISJOIN) |>
        dplyr::mutate(bg2010_area = units::drop_units(sf::st_area(geometry)))
      
      overlapping <- sf::st_intersection(
        bgs_2010,
        bgs_intersect_rtd |> dplyr::select(GISJOIN_2020 = GISJOIN)  # was study_bgs
      ) |>
        dplyr::mutate(
          overlap_area        = units::drop_units(sf::st_area(geometry)),
          weight_2010_to_2020 = overlap_area / bg2010_area
        ) |>
        sf::st_drop_geometry() |>
        dplyr::select(GISJOIN_2010, GISJOIN_2020, overlap_area, bg2010_area, weight_2010_to_2020)
      
      outside <- bgs_2010 |>
        dplyr::filter(!GISJOIN_2010 %in% overlapping$GISJOIN_2010) |>
        sf::st_drop_geometry() |>
        dplyr::mutate(
          GISJOIN_2020        = GISJOIN_2010,
          overlap_area        = bg2010_area,
          weight_2010_to_2020 = 1
        ) |>
        dplyr::select(GISJOIN_2010, GISJOIN_2020, overlap_area, bg2010_area, weight_2010_to_2020)
      
      dplyr::bind_rows(overlapping, outside)
    }
  ),
  
  # Crosswalk 2010-2019 data to 2020 boundaries
  tar_target(
    denver_acs_2010_list_est,
    lapply(denver_acs_2010_list, function(df) {
      
      df |>
        dplyr::inner_join(bg_xwalk, by = c("GISJOIN" = "GISJOIN_2010")) %>%
        dplyr::mutate(across(where(is.numeric), ~ .x * weight_2010_to_2020)) %>%
        sf::st_drop_geometry() %>%
        dplyr::select(GISJOIN_2020, where(is.numeric)) %>%
        dplyr::group_by(GISJOIN_2020) %>%
        dplyr::summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
        dplyr::rename(GISJOIN = GISJOIN_2020) %>%
        dplyr::filter(GISJOIN %in% study_gisjoin) %>%
        units::drop_units()
    })
  ),
  
  # Bind datasets
  tar_target(
    acs_2010_all,
    purrr::imap_dfr(
      denver_acs_2010_list_est,
      ~ dplyr::mutate(.x, year = as.integer(stringr::str_extract(.y, "20\\d{2}")))
    )
  ),
  
  tar_target(
    acs_2020_all,
    purrr::imap_dfr(
      denver_acs_2020_list,
      ~ dplyr::mutate(.x, year = as.integer(stringr::str_extract(.y, "20\\d{2}")))
    )
  ),
  
  # Combine to form full ACS dataset
  tar_target(
    denver_acs,
    dplyr::bind_rows(acs_2010_all, acs_2020_all)
  ),
  
  # ------ Load LODES data ------
  
  # List LODES files
  tar_target(
    lodes_files, list.files(
    "data/LODES_WAC",
    pattern = "^co_wac_S000_JT00_20.*\\.csv$",
    full.names = TRUE
  )),
  
  # Load LODES files
  tar_target(
    lodes_list_raw,
    {
      files <- lodes_files
      purrr::map(
        files,
        ~ read.csv(.x, stringsAsFactors = FALSE)
      ) |>
        purrr::set_names(
          tools::file_path_sans_ext(basename(files))
        )
    }),
  
  # Process LODES data files
  tar_target(
    lodes_list,
    purrr::imap_dfr(
      lodes_list_raw,   # already read CSVs
      ~ process_lodes_file(.x, .y, geos = study_gisjoin)
    )
  ),
  
  # Bind rows
  tar_target(
    denver_lodes,
    dplyr::bind_rows(lodes_list)
  ),
  
  # Join with ACS
  tar_target(
    acs_lodes_data,
    dplyr::left_join(
      denver_acs,
      denver_lodes,
      by = c("GISJOIN","year")
    )
  ),
  
  # List LODES files
  tar_target(
    lodes_files_private, list.files(
      "data/LODES_WAC_JT02",
      pattern = "^co_wac_S000_JT02_20.*\\.csv$",
      full.names = TRUE
    )),
  
  # Load LODES files
  tar_target(
    lodes_list_private_raw,
    {
      files <- lodes_files_private
      purrr::map(
        files,
        ~ read.csv(.x, stringsAsFactors = FALSE)
      ) |>
        purrr::set_names(
          tools::file_path_sans_ext(basename(files))
        )
    }),
  
  # Process LODES data files
  tar_target(
    lodes_list_private,
    purrr::imap_dfr(
      lodes_list_private_raw,   # already read CSVs
      ~ process_lodes_file(.x, .y, geos = study_gisjoin)
    )
  ),
  
  # Bind rows
  tar_target(
    denver_lodes_private,
    dplyr::bind_rows(lodes_list_private)
  ),
  
  # ------- Load BDS data and estimate size by county by sector ------
  
  # Load BDS tables
  tar_target(
    bds_county_sec_raw,
    read.csv("data/BDS/bds2023_st_cty_sec.csv")
  ),
  
  tar_target(
    bds_county_size_raw,
    read.csv("data/BDS/bds2023_st_cty_fzc.csv")
  ),
  
  tar_target(
    bds_st_sec_size_raw,
    read.csv("data/BDS/bds2023_st_sec_fzc.csv")
  ),
  
  # Cleaned BDS tables
  tar_target(
    bds_sec_clean,
    clean_bds(bds_county_sec_raw)
  ),
  
  tar_target(
    bds_size_clean,
    clean_bds(bds_county_size_raw)
  ),
  
  tar_target(
    bds_state_clean,
    clean_bds(bds_st_sec_size_raw)
  ),
  
  # Generate list of prior matrices by year
  tar_target(
    prior_matrix_year,
    {
      tmp <- bds_state_clean %>%
        group_by(year, sector, fsizecoarse) %>%
        summarize(estabs = sum(estabs, na.rm = TRUE), .groups = "drop") %>%
        group_by(year, sector) %>%
        mutate(p_estabs = estabs / sum(estabs)) %>%
        ungroup()
      
      yrs <- sort(unique(tmp$year))
      
      mats <- lapply(yrs, function(y){
        xtabs(
          p_estabs ~ sector + fsizecoarse,
          data = tmp[tmp$year == y, ]
        ) |>
          as.matrix() |>
          pmax(1e-8)
      })
      names(mats) <- yrs   # ← THIS IS THE IMPORTANT PART
      
      mats
    }
  ),
  
  # Iterative proportional fitting to estimate county x sector x size
  tar_target(
    county_sec_size,
    {
      combos <- expand_grid(
        cty = unique(bds_sec_clean$county_fips),
        year = unique(bds_sec_clean$year)
      )
      
      purrr::pmap_dfr(combos, function(cty_i, year_i){
        
        prior <- prior_matrix_year[[as.character(year_i)]]
        
        # sector vector with fallback
        sec_vec <- setNames(rep(1e-6, nrow(prior)), rownames(prior))
        sec_tab <- xtabs(estabs ~ sector, data = bds_sec_clean %>% filter(county_fips == cty_i, year == year_i))
        sec_vec[names(sec_tab)] <- pmax(sec_tab, 1e-6)  # ensure non-zero
        
        # size vector with fallback
        size_vec <- setNames(rep(1e-6, ncol(prior)), colnames(prior))
        size_tab <- xtabs(estabs ~ fsizecoarse, data = bds_size_clean %>% filter(county_fips == cty_i, year == year_i))
        size_vec[names(size_tab)] <- pmax(size_tab, 1e-6)
        
        # scale sector totals to match size totals
        sec_vec <- sec_vec * (sum(size_vec) / sum(sec_vec))
        
        # now run IPF
        mat <- ipf_county(sec_vec, size_vec, prior)
        
        as.data.frame(as.table(mat)) %>%
          dplyr::rename(
            sector = 1,
            fsizecoarse = 2,
            estabs = 3
          ) %>%
          dplyr::mutate(county_fips = cty_i,
                        year = year_i)
      })
    }
  ),
  
  # ------ Estimate bg-level firm size bins ------
  
  # Load county-level LODES for private primary jobs
  tar_target(
    county_lodes_files,
    list.files(
      "data/OnTheMap_Counties",
      pattern = "^polygon_.*\\.csv$",
      full.names = TRUE
    )
  ),
  
  tar_target(
    county_lodes_list_raw,
    {
      files <- county_lodes_files
      purrr::map(
        files,
        ~ read.csv(.x, stringsAsFactors = FALSE)
      ) |>
        purrr::set_names(
          tools::file_path_sans_ext(basename(files))
        )
    }
  ),
  
  # Format county-level LODES data
  tar_target(
    county_lodes_list,
    purrr::imap_dfr(
      county_lodes_list_raw,   # already read CSVs
      ~ process_county_lodes_file(.x, .y)
    )
  ),
  
  # Bind county-level LODES data into panel set
  tar_target(
    county_lodes,
    dplyr::bind_rows(county_lodes_list)
  ),
  
  # Join to bg-level private primary RTD set keying county FIPS to GISJOIN
  tar_target(
    denver_lodes_county,
    county_lodes %>%
      rowwise() %>%
      mutate(
        match_rows = list(
          denver_lodes_private %>%
            filter(year == year_county, str_detect(GISJOIN, county_fips))
        )
      ) %>%
      unnest(match_rows)
  ),
  
  # Pivot LODES and estimate bg-county job count weights
  tar_target(
    lodes_long,
    denver_lodes_county %>%
      pivot_longer(
        cols = starts_with("NAICS"),
        names_to = "col_name",
        values_to = "emp"
      ) %>%
      mutate(
        naics2 = str_extract(col_name, "\\d{2}(-\\d{2})?"),
        type   = dplyr::if_else(str_detect(col_name, "_county$"), "county", "bg")
      ) %>%
      select(-col_name) %>%
      pivot_wider(
        names_from = type,
        values_from = emp,
        values_fill = 0
      ) %>%
      rename(
        emp_bg = bg,
        emp_county = county
      ) %>%
      filter(!is.na(emp_bg)) %>%
      mutate(
        sector_share = ifelse(emp_county > 0, emp_bg / emp_county, 0)
      )
  ),
  
  # Estimate firms by size using BDS county x sector x size estimates
  tar_target(
    bg_size_estimates,
    {
      # Join county-level size-by-sector counts
      lodes_with_size <- lodes_long %>%
        left_join(
          county_sec_size,
          by = c("county_fips", "year", "naics2" = "sector")
        ) %>%
        mutate(
          # scale county sector-size counts to BG using sector share
          estabs_bg = estabs * sector_share
        )
      
      # Aggregate across sectors to get BG-level size bins
      bg_size_final <- lodes_with_size %>%
        filter(!is.na(fsizecoarse)) %>%
        group_by(year, county_fips, GISJOIN, fsizecoarse) %>%
        summarise(
          estabs_bg = sum(estabs_bg, na.rm = TRUE),
          .groups = "drop"
        ) %>%
      
      return(bg_size_final)
    }
  ),
  
  # Pivot wider for more regression-friendly formatting
  tar_target(
    lodes_acs_size_data,
    {
      bg_size_estimates %>%
        mutate(
          size_bin = case_when(
            fsizecoarse == "a) 1 to 19"   ~ "small",
            fsizecoarse == "b) 20 to 499" ~ "mid",
            fsizecoarse == "c) 500+"      ~ "large",
            TRUE ~ NA_character_
          )
        ) %>%
        group_by(year, GISJOIN, county_fips, size_bin) %>%
        summarise(estabs_sum = sum(estabs_bg, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(
          names_from = size_bin,
          values_from = estabs_sum,
          names_prefix = "estabs_",
          values_fill = 0
        ) %>%
        # join to regression data
        left_join(acs_lodes_data, by = c("year", "GISJOIN"))
    }
  ),
  
  # ------ Estimate bg-level firms by sector ------
  
  # Join LODES to BDS and estimate firms, estabs, and entry by sector
  tar_target(
    bg_sec_estimates,
    lodes_long %>%
      left_join(
        bds_sec_clean,
        by = c("county_fips", "year", "naics2" = "sector")
      ) %>%
      mutate(
        firms_bg = firms * sector_share,
        estabs_bg = estabs * sector_share,
        estabs_entry_bg = estabs_entry * sector_share
      )
  ),
  
  # Pivot wider and join to dataset
  tar_target(
    lodes_acs_bds_data,
    bg_sec_estimates %>%
      group_by(year, GISJOIN, naics2) %>%
      summarise(
        firms = sum(firms_bg, na.rm = TRUE),
        estabs = sum(estabs_bg, na.rm = TRUE), 
        estabs_entry = sum(estabs_entry_bg, na.rm = TRUE),
        .groups = "drop") %>%
      pivot_wider(
        names_from = naics2,
        values_from = c(firms, estabs, estabs_entry),
        names_sep = "_",
        values_fill = 0
      ) %>%
      mutate(
        firms_total = rowSums(across(starts_with("firms_")), na.rm = TRUE),
        estabs_total = rowSums(across(starts_with("estabs_") & !starts_with("estabs_entry")), na.rm = TRUE),
        estabs_entry_total = rowSums(across(starts_with("estabs_entry_")), na.rm = TRUE)
      ) %>%
      left_join(lodes_acs_size_data, by=c("year","GISJOIN"))
  ),
  
  # ------ Estimate wages -------
  
  # Load county-level LODES for private primary jobs
  tar_target(
    qcew_files,
    list.files(
      "data/QCEW",
      pattern = "\\.csv$",
      full.names = TRUE
    )
  ),
  
  tar_target(qcew_list, process_qcew(qcew_files)),
  
  tar_target(
    qcew,
    bind_rows(qcew_list)
  ),
  
  tar_target(wages_bg, {
    lodes_long |>
      dplyr::left_join(
        qcew |>
          dplyr::select(
            year,
            county_fips     = area_fips,
            naics2          = industry_code,
            avg_weekly_wage
          ),
        by = c("year", "county_fips", "naics2")
      ) |>
      dplyr::group_by(GISJOIN, year) |>
      dplyr::mutate(
        # Weight by share of BG employment in each sector
        # emp_bg is jobs in this BG-sector, sum gives total BG jobs
        bg_total_jobs = sum(emp_bg, na.rm = TRUE),
        bg_sector_share = emp_bg / bg_total_jobs
      ) |>
      dplyr::summarise(
        bg_avg_weekly_wage = sum(bg_sector_share * avg_weekly_wage,
                                 na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        bg_avg_annual_wage = bg_avg_weekly_wage * 52,
        lwage              = log(pmax(bg_avg_annual_wage, 1))
      )
  }),
  
  # ------ Assign treatment and controls by cohort -------
  
  # Buffer highways
  tar_target(
    highway_buffer,
    sf::st_buffer(highways_sf, 1000)
  ),
  
  tar_target(
    cf_buffer,
    sf::st_buffer(counterfactual_sf, 1000)
  ),
  
  # Final BG filter near highways
  tar_target(
    study_bgs,
    {
      bgs <- bgs_intersect_rtd
      bgs[ lengths(sf::st_intersects(bgs, highway_buffer)) > 0 | lengths(sf::st_intersects(bgs, cf_buffer)) > 0, ]
    }
  ),
  
  # Assign groups by center dist from rail and join data
  tar_target(
    regression_set_raw,
    {
      bg_treatments <- lapply(
        c("treatment","spillover","continuous"),
        function(x) {assign_groups(lodes_acs_bds_data,
                                  study_bgs,
                                  RTD_commuter_sf,
                                  type = x)}
      )
      names(bg_treatments) <- c("treatment","spillover","continuous")
      
      bg_treatments
    }
    ),
  
  # ------ FEOLS regression analysis ------
  
  # Clean regression set
  tar_target(
    regression_set,
    lapply(regression_set_raw, function (df) var_estimation(df, wages_bg))
  ),
  
  tar_target(
    regression_analysis,
    lapply(names(regression_set), function(type) {
      run_regression(regression_set[[type]],
                     outcomes = c("lpop", "lmhi", "lwage", "ljobs",
                                  "lhousing", "ljobs_per_100", "lmedian_rent",
                                  "lestabs_entry", "estabs_small_share",
                                  "estabs_mid_share", "estabs_large_share"),
                     controls = type)
    }) |> setNames(names(regression_set))
  ),
  
  tar_target(var_dict, c(
    "lhs: lpop"               = "Log Population",
    "lhs: lhousing"           = "Log Housing Units",
    "lhs: lmedian_rent"       = "Log Median Rent",
    "lhs: lmhi"               = "Log Median Household Income",
    "lhs: ljobs"              = "Log Jobs",
    "lhs: ljobs_per_100"      = "Log Jobs per 100 Residents",
    "lhs: lwage"              = "Log Avg Annual Wage",
    "lhs: lestabs_entry"      = "Log Establishment Entries",
    "lhs: estabs_small_share" = "Share of Small Establishments",
    "lhs: estabs_mid_share"   = "Share of Mid-Size Establishments",
    "lhs: estabs_large_share" = "Share of Large Establishments"
  )),
  
  tar_target(att_raw, 
             extract_all_atts(regression_analysis[c("treatment", "spillover")])),
  
  tar_target(att_raw_continuous, 
             extract_all_atts(regression_analysis["continuous"])),
  
  tar_target(att_formatted,
             format_att_table(att_raw, clean_names=var_dict)),
  
  tar_target(att_export,
             export_att_table(att_formatted, "att_table.tex")),
  
  tar_target(continuous_coefs,
              format_continuous_table(att_raw_continuous, clean_names=var_dict)),
  
  tar_target(continuous_export,
              export_continuous_table(continuous_coefs, "continuous_table.tex")),
  
  # ------ Extract coefficients ------
  
  tar_target(es_coefs_treatment,
             extract_sunab_coefs(regression_analysis, design = "treatment")
  ),
  
  tar_target(es_coefs_spillover,
             extract_sunab_coefs(regression_analysis, design = "spillover")
  ),
  
  tar_target(es_coefs_continuous,
             extract_event_study_coefs(regression_analysis, design = "continuous")
  ),
  
  # ------ Single design plots ------
  
  tar_target(plot_sunab_residential, {
    p <- make_event_study_panel(
      es_coefs_treatment,
      outcomes = OUTCOME_GROUPS$residential,
      title    = "Residential Outcomes: Treatment Group (<1km)",
      ncol     = 2
    )
    save_plot(p, "sunab_residential_treatment.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_sunab_labor, {
    p <- make_event_study_panel(
      es_coefs_treatment,
      outcomes = OUTCOME_GROUPS$labor,
      title    = "Labor Market Outcomes: Treatment Group (<1km)",
      ncol     = 2
    )
    save_plot(p, "sunab_labor_treatment.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_sunab_composition, {
    p <- make_event_study_panel(
      es_coefs_treatment,
      outcomes  = OUTCOME_GROUPS$composition,
      title     = "Firm Size Composition: Treatment Group (<1km)",
      ncol      = 3,
      ci_limit  = 0.05
    )
    save_plot(p, "sunab_composition_treatment.png", width = 15, height = 5)
    p
  }),
  
  tar_target(plot_sunab_residential_spillover, {
    p <- make_event_study_panel(
      es_coefs_spillover,
      outcomes = OUTCOME_GROUPS$residential,
      title    = "Residential Outcomes: Spillover Group (2-5km)",
      ncol     = 2,
      color    = "#d7191c"
    )
    save_plot(p, "sunab_residential_spillover.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_sunab_labor_spillover, {
    p <- make_event_study_panel(
      es_coefs_spillover,
      outcomes = OUTCOME_GROUPS$labor,
      title    = "Labor Market Outcomes: Spillover Group (2-5km)",
      ncol     = 2,
      color    = "#d7191c"
    )
    save_plot(p, "sunab_labor_spillover.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_sunab_composition_spillover, {
    p <- make_event_study_panel(
      es_coefs_spillover,
      outcomes  = OUTCOME_GROUPS$composition,
      title     = "Firm Size Composition: Spillover Group (2-5km)",
      ncol      = 3,
      color     = "#d7191c",
      ci_limit  = 0.05
    )
    save_plot(p, "sunab_composition_spillover.png", width = 15, height = 5)
    p
  }),
  
  tar_target(plot_continuous_residential, {
    p <- make_event_study_panel(
      es_coefs_continuous,
      outcomes = OUTCOME_GROUPS$residential,
      title    = "Residential Outcomes: Continuous DiD",
      ncol     = 2,
      color    = "#1a9641",
      caption  = "Continuous DiD: i(ref_year, log_dist, ref = -1) | GISJOIN + year."
    )
    save_plot(p, "continuous_residential.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_continuous_labor, {
    p <- make_event_study_panel(
      es_coefs_continuous,
      outcomes = OUTCOME_GROUPS$labor,
      title    = "Labor Market Outcomes: Continuous DiD",
      ncol     = 2,
      color    = "#1a9641",
      caption  = "Continuous DiD: i(ref_year, log_dist, ref = -1) | GISJOIN + year."
    )
    save_plot(p, "continuous_labor.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_continuous_composition, {
    p <- make_event_study_panel(
      es_coefs_continuous,
      outcomes  = OUTCOME_GROUPS$composition,
      title     = "Firm Size Composition: Continuous DiD",
      ncol      = 3,
      color     = "#1a9641",
      ci_limit  = 0.05,
      caption   = "Continuous DiD: i(ref_year, log_dist, ref = -1) | GISJOIN + year."
    )
    save_plot(p, "continuous_composition.png", width = 15, height = 5)
    p
  }),
  
  # ------ Overlay plots (all three designs) ------
  
  tar_target(plot_overlay_residential, {
    p <- make_event_study_overlay(
      coef_list = list(
        "Treatment (<1km)"  = es_coefs_treatment,
        "Spillover (2-5km)" = es_coefs_spillover,
        "Continuous DiD"    = es_coefs_continuous
      ),
      outcomes = OUTCOME_GROUPS$residential,
      title    = "Residential Outcomes: Treatment vs Spillover vs Continuous",
      ncol     = 2
    )
    save_plot(p, "overlay_residential.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_overlay_labor, {
    p <- make_event_study_overlay(
      coef_list = list(
        "Treatment (<1km)"  = es_coefs_treatment,
        "Spillover (2-5km)" = es_coefs_spillover,
        "Continuous DiD"    = es_coefs_continuous
      ),
      outcomes = OUTCOME_GROUPS$labor,
      title    = "Labor Market Outcomes: Treatment vs Spillover vs Continuous",
      ncol     = 2
    )
    save_plot(p, "overlay_labor.png", width = 12, height = 9)
    p
  }),
  
  tar_target(plot_overlay_composition, {
    p <- make_event_study_overlay(
      coef_list = list(
        "Treatment (<1km)"  = es_coefs_treatment,
        "Spillover (2-5km)" = es_coefs_spillover,
        "Continuous DiD"    = es_coefs_continuous
      ),
      outcomes  = OUTCOME_GROUPS$composition,
      title     = "Firm Size Composition: Treatment vs Spillover vs Continuous",
      ncol      = 3,
      ci_limit  = 0.05
    )
    save_plot(p, "overlay_composition.png", width = 15, height = 5)
    p
  }),
  
  # ------ Load network shapefiles ------
  
  tar_target(origins_raw,
             st_centroid(bgs_intersect_rtd) %>%
               select(GISJOIN, geometry)
  ),
  
  # Counterfactual commuter rail pre-split at proposed stations
  # No speed field — will derive from GTFS B/N line averages
  tar_target(transit_counterfactual,
             st_read("data\\Network_Shapefiles\\Counterfactual_CR.shp")
  ),
  
  # Pedestrian network — used for both walk and cycle weights
  tar_target(pedestrian_routes,
             st_read("data\\Network_Shapefiles\\RTD_Walking_Network.shp")
  ),
  
  # Local driving network — free flow speeds
  tar_target(road_routes,
             st_read("data\\Network_Shapefiles\\RTD_Driving_Network.shp")
  ),
  
  # Highway networks with congestion adjusted travel times
  tar_target(highway_routes_2013,
             st_read("data\\Network_Shapefiles\\RTD_Traffic_2013.shp")
  ),
  
  tar_target(highway_routes_2023,
             st_read("data\\Network_Shapefiles\\RTD_Traffic_2023.shp")
  ),
  
  # ── GTFS Processing ───────────────────────────────────────────────────────────
  
  # Build all transit segments from GTFS
  tar_target(gtfs_segments,
             build_transit_edges_from_gtfs(
               gtfs_dirs       = GTFS_DIRS,
               crs             = CRS,
               peak_hour_start = 6.5,
               peak_hour_end   = 9.5
             )
  ),
  
  # Extract average speed for B and N lines
  # Used to calculate counterfactual travel times
  tar_target(cf_rail_speed,
             gtfs_segments %>%
               st_drop_geometry() %>%
               filter(route_id %in% CF_RAIL_LINES) %>%
               mutate(
                 length_m     = as.numeric(st_length(gtfs_segments %>%
                                                       filter(route_id %in% CF_RAIL_LINES))),
                 speed_kph    = (length_m / 1000) / (TravelTime / 60)
               ) %>%
               summarise(avg_speed_kph = mean(speed_kph, na.rm = TRUE)) %>%
               pull(avg_speed_kph)
  ),
  
  # Add travel times to counterfactual lines using B/N average speed
  tar_target(transit_counterfactual_tt,
             transit_counterfactual %>%
               mutate(
                 length_m    = as.numeric(st_length(.)),
                 TravelTime  = length_m / (cf_rail_speed * 1000 / 60),
                 cohort_year = 2045,    # all counterfactual lines are new
                 route_id    = case_when(SUBDIV == "B" ~ "113B",
                                        SUBDIV == "N" ~ "117N"),
                 mode = "commuter_rail"
               )
  ),
  
  # ── Edge Construction ─────────────────────────────────────────────────────────
  
  # Walk edges
  tar_target(walk_edges,
             make_road_edges(
               lines          = pedestrian_routes,
               weight_field   = "TravelTime",
               mode_tag       = "walk"
             )
  ),
  
  # Cycle edges
  tar_target(cycle_edges,
             make_road_edges(
               lines          = pedestrian_routes,
               weight_field   = "Cycle_TT",
               mode_tag       = "cycle"
             )
  ),
  
  # Drive edges — 2013 and 2023 highway + local roads
  tar_target(drive_edges_2013,
             get_drive_edges(highway_routes_2013, road_routes)
  ),
  
  # All transit edges from GTFS
  # Split into bus and rail for separate connector handling
  tar_target(bus_edges_gtfs,
             make_transit_edges(
               lines          = gtfs_segments %>% filter(mode == "bus"),
               weight_field   = "TravelTime",
               mode_tag       = "bus",
               bidirectional  = TRUE
             )
  ),
  
  tar_target(rail_edges_existing,
             make_transit_edges(
               lines          = gtfs_segments %>%
                 filter(mode %in% c("commuter_rail", "light_rail")),
               weight_field   = "TravelTime",
               mode_tag       = "rail",
               bidirectional  = TRUE
             )
  ),
  
  tar_target(rail_edges_counterfactual,
             make_transit_edges(
               lines          = transit_counterfactual_tt,
               weight_field   = "TravelTime",
               mode_tag       = "rail",
               bidirectional  = TRUE
             )
  ),
  
  # All rail edges combined
  tar_target(rail_edges_all,
             bind_rows(
               rail_edges_existing,
               rail_edges_counterfactual
             )
  ),
  
  # ── Connector Edges ───────────────────────────────────────────────────────────
  
  tar_target(lr_connector_edges,
             make_connector_edges_from_nodes(
               transit_edges  = rail_edges_all,
               road_edges     = walk_edges,
               board_penalty  = BOARD_PENALTIES$rail_board,
               alight_penalty = BOARD_PENALTIES$rail_alight,
               mode_tag       = "lr",
               max_dist       = 500
             )
  ),
  
  tar_target(bus_connector_edges,
             make_connector_edges_from_nodes(
               transit_edges  = bus_edges_gtfs,
               road_edges     = walk_edges,
               board_penalty  = BOARD_PENALTIES$bus_board,
               alight_penalty = BOARD_PENALTIES$bus_alight,
               mode_tag       = "bus",
               max_dist       = 200
             )
  ),
  
  # ------ Network Analysis ------
  
  # Driving graphs
  
  tar_target(graph_drive_alone,
             assemble_graph(road_edges = drive_edges_2013)
  ),
  
  tar_target(origin_nodes_drive_alone,
             snap_points_to_graph(
               points   = origins_raw,
               graph    = graph_drive_alone,
               max_dist = 2000,
               id_field = "GISJOIN"
             )
  ),
  
  tar_target(od_drive_alone,
             run_od_matrix(
               graph        = graph_drive_alone,
               origin_nodes = origin_nodes_drive_alone,
               dest_nodes   = origin_nodes_drive_alone,
               label        = "drive_alone"
             )
  ),
  
  # ── Walk Alone ────────────────────────────────────────────────────────────────
  
  tar_target(graph_walk_alone,
             assemble_graph(road_edges = walk_edges)
  ),
  
  tar_target(origin_nodes_walk,
             snap_points_to_graph(
               points   = origins_raw,
               graph    = graph_walk_alone,
               max_dist = 2000,
               id_field = "GISJOIN"
             )
  ),
  
  tar_target(od_walk_alone,
             run_od_matrix(
               graph        = graph_walk_alone,
               origin_nodes = origin_nodes_walk,
               dest_nodes   = origin_nodes_walk,
               label        = "walk_alone"
             )
  ),
  
  # ── Cycle Alone ───────────────────────────────────────────────────────────────
  
  tar_target(graph_cycle_alone,
             assemble_graph(road_edges = cycle_edges)
  ),
  
  tar_target(origin_nodes_cycle,
             snap_points_to_graph(
               points   = origins_raw,
               graph    = graph_cycle_alone,
               max_dist = 2000,
               id_field = "GISJOIN"
             )
  ),
  
  tar_target(od_cycle_alone,
             run_od_matrix(
               graph        = graph_cycle_alone,
               origin_nodes = origin_nodes_cycle,
               dest_nodes   = origin_nodes_cycle,
               label        = "cycle_alone"
             )
  ),
  
  # ── Base Transit Graph For Snapping ──────────────────────────────────────────
  
  tar_target(graph_base_transit,
             assemble_graph(
               road_edges      = walk_edges,
               transit_edges   = bind_rows(
                 rail_edges_existing %>% filter(cohort_year == 0),
                 bus_edges_gtfs
               ),
               connector_edges = bind_rows(
                 lr_connector_edges,
                 bus_connector_edges
               )
             )
  ),
  
  tar_target(origin_nodes_transit,
             snap_points_to_graph(
               points   = origins_raw,
               graph    = graph_base_transit,
               max_dist = 1000,
               id_field = "GISJOIN"
             )
  ),
  
  # ── Transit Graphs + OD Matrices (one per network state) ─────────────────────
  
  tar_map(
    values = SCENARIOS,
    names  = "scenario",
    
    tar_target(rail_filtered,
               filter_scenario(
                 edges    = rail_edges_all,
                 scenario = scenario
               )
    ),
    
    tar_target(graph_walk_transit,
               assemble_graph(
                 road_edges      = walk_edges,
                 transit_edges   = bind_rows(rail_filtered, bus_edges_gtfs),
                 connector_edges = bind_rows(lr_connector_edges, bus_connector_edges)
               )
    ),
    
    tar_target(od_walk_transit,
               run_od_matrix(
                 graph        = graph_walk_transit,
                 origin_nodes = origin_nodes_transit,
                 dest_nodes   = origin_nodes_transit,
                 label        = paste0("walk_transit_", scenario)
               )
    )
  ),
  
  # ── Scenario OD Bundles ───────────────────────────────────────────────────────
  # Each scenario binds walk, cycle, scenario-specific transit, and base drive
  # Transit referenced explicitly by name from tar_map above
  # Drive held at 2013 baseline for all scenarios (isolates transit effect)
  
  tar_target(od_scenario_pre_2016,
             parse_mode_scenario(bind_rows(
               od_walk_alone,
               od_cycle_alone,
               od_walk_transit_pre_2016,
               od_drive_alone
             ))
  ),
  
  tar_target(od_scenario_post_2016,
             parse_mode_scenario(bind_rows(
               od_walk_alone,
               od_cycle_alone,
               od_walk_transit_post_2016,
               od_drive_alone
             ))
  ),
  
  tar_target(od_scenario_post_2017,
             parse_mode_scenario(bind_rows(
               od_walk_alone,
               od_cycle_alone,
               od_walk_transit_post_2017,
               od_drive_alone
             ))
  ),
  
  tar_target(od_scenario_post_2019,
             parse_mode_scenario(bind_rows(
               od_walk_alone,
               od_cycle_alone,
               od_walk_transit_post_2019,
               od_drive_alone
             ))
  ),
  
  tar_target(od_scenario_post_2020,
             parse_mode_scenario(bind_rows(
               od_walk_alone,
               od_cycle_alone,
               od_walk_transit_post_2020,
               od_drive_alone
             ))
  ),
  
  tar_target(od_scenario_counterfactual,
             parse_mode_scenario(bind_rows(
               od_walk_alone,
               od_cycle_alone,
               od_walk_transit_counterfactual,
               od_drive_alone
             ))
  ),
  
  # ── OD Cleaning (ceiling only — no common pair filter for aggregate t_bar) ───
  
  tar_target(required_modes,
             c("walk", "bike", "transit", "drive")
  ),
  
  tar_target(od_clean_pre_2016,
             od_scenario_pre_2016 |> apply_time_ceiling()
  ),
  
  tar_target(od_clean_post_2016,
             od_scenario_post_2016 |> apply_time_ceiling()
  ),
  
  tar_target(od_clean_post_2017,
             od_scenario_post_2017 |> apply_time_ceiling()
  ),
  
  tar_target(od_clean_post_2019,
             od_scenario_post_2019 |> apply_time_ceiling()
  ),
  
  tar_target(od_clean_post_2020,
             od_scenario_post_2020 |> apply_time_ceiling()
  ),
  
  tar_target(od_clean_counterfactual,
             od_scenario_counterfactual |> apply_time_ceiling()
  ),
  
  # ── Export ────────────────────────────────────────────────────────────────────
  
  tar_target(export_results, {
    
    all_od <- dplyr::bind_rows(
      od_scenario_pre_2016,
      od_scenario_post_2016,
      od_scenario_post_2017,
      od_scenario_post_2019,
      od_scenario_post_2020,
      od_scenario_counterfactual
    )
    
    dir.create("outputs", showWarnings = FALSE)
    
    data.table::fwrite(all_od, "outputs/OD_All_Scenarios.csv")
    
    all_od |>
      dplyr::group_by(scenario) |>
      dplyr::group_walk(~ data.table::fwrite(
        .x,
        paste0("outputs/OD_", .y$scenario, ".csv")
      ))
    
    "outputs/OD_All_Scenarios.csv"
  },
  format = "file"
  ),
  
  # ── Collect transit ODs into named list for downstream use ───────────────────
  
  tar_target(od_transit_all, {
    list(
      pre_2016       = od_walk_transit_pre_2016,
      post_2016      = od_walk_transit_post_2016,
      post_2017      = od_walk_transit_post_2017,
      post_2019      = od_walk_transit_post_2019,
      post_2020      = od_walk_transit_post_2020,
      counterfactual = od_walk_transit_counterfactual
    )
  }),
  
  # ------ Inclusive value preference-shifter ------
  
  # Return list of ACS codebook files
  tar_target(
    tt_codebooks,
    list.files(
      "data/nhgis0017_csv",
      pattern = "_codebook\\.txt$",
      full.names = TRUE
    )
  ),
  
  # Load codebooks
  tar_target(
    nhgis_tt_prefix_lookup,
    purrr::map_dfr(tt_codebooks, parse_nhgis_codebook) |>
      dplyr::mutate(prefix = sub("E\\d+$", "", variable)) |>
      dplyr::distinct(census_table, prefix)
  ),
  
  # Look up table prefixes in codebooks
  tar_target(
    bg_prefix_to_table,
    {
      out <- nhgis_tt_prefix_lookup |>
        dplyr::filter(census_table %in% "B08301") |>
        dplyr::distinct(prefix, census_table) |>
        dplyr::group_by(prefix) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
      stats::setNames(out$census_table, out$prefix)
    }
  ),
  
  tar_target(
    tract_prefix_to_table,
    {
      out <- nhgis_tt_prefix_lookup |>
        dplyr::filter(census_table %in% c("B08141","B08119")) |>
        dplyr::distinct(prefix, census_table) |>
        dplyr::group_by(prefix) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
      stats::setNames(out$census_table, out$prefix)
    }
  ),
  
  # Create list of RTD tracts from earlier valid geos
  tar_target(study_tracts,
             substr(study_gisjoin, 1, nchar(study_gisjoin) - 1)
  ),
  
  tar_target(tracts2010,
             substr(gisjoin2010, 1, nchar(gisjoin2010) - 1)
  ),
  
  # Load ACS files
  tar_target(
    tt_2013_bg,
    read.csv("data/nhgis0017_csv/nhgis0017_ds201_20135_blck_grp.csv")
  ),
  
  tar_target(
    tt_2013_tract,
    read.csv("data/nhgis0017_csv/nhgis0017_ds202_20135_tract.csv")
  ),
  
  tar_target(
    tt_2023_bg,
    read.csv("data/nhgis0017_csv/nhgis0017_ds267_20235_blck_grp.csv")
  ),
  
  tar_target(
    tt_2023_tract,
    read.csv("data/nhgis0017_csv/nhgis0017_ds268_20235_tract.csv")
  ),
  
  # Clean datasets
  tar_target(
    denver_tt_2023_bg,
    clean_acs_df(
      tt_2023_bg,
      prefix_to_table = bg_prefix_to_table,
      valid_gisjoin = study_gisjoin
    )
  ),
  
  tar_target(
    denver_tt_2013_bg,
    clean_acs_df(
      tt_2013_bg,
      prefix_to_table = bg_prefix_to_table,
      valid_gisjoin = gisjoin2010
    )
  ),
  
  tar_target(tract_2010, st_read("data\\nhgis0018_shapefile_tl2010_us_tract_2010\\US_tract_2010.shp")),
  tar_target(tract_2020, st_read("data\\nhgis0018_shapefile_tl2020_us_tract_2020\\US_tract_2020.shp")),
  tar_target(tract_2010_sf, project_to_stateplane(tract_2010)),
  tar_target(tract_2020_sf, project_to_stateplane(tract_2020)),
  
  # Filter bgs to Colorado
  tar_target(
    tract_co_2020,
    dplyr::filter(tract_2020_sf, STATEFP == "08")
  ),
  
  tar_target(
    tract_co_2010,
    dplyr::filter(tract_2010_sf, STATEFP10 == "08")
  ),
  
  tar_target(
    tracts_touch_rtd_2020,
    {
      tracts <- tract_co_2020
      rtd_valid <- sf::st_make_valid(RTD_sf)        
      tracts[ lengths(sf::st_intersects(tracts, rtd_valid)) > 0, ]
    }
  ),
  
  tar_target(
    tracts_touch_rtd_2010,
    {
      tracts <- tract_co_2010
      rtd_valid <- sf::st_make_valid(RTD_sf)        
      tracts[ lengths(sf::st_intersects(tracts, rtd_valid)) > 0, ]
    }
  ),
  
  tar_target(tract_gisjoin_2020, tracts_touch_rtd_2020$GISJOIN),
  
  tar_target(tract_gisjoin_2010, tracts_touch_rtd_2010$GISJOIN),
  
  tar_target(
    denver_tt_2023_tract,
    clean_acs_df(
      tt_2023_tract,
      prefix_to_table = tract_prefix_to_table,
      valid_gisjoin = tract_gisjoin_2020
    )
  ),
  
  tar_target(
    denver_tt_2013_tract,
    clean_acs_df(
      tt_2013_tract,
      prefix_to_table = tract_prefix_to_table,
      valid_gisjoin = tract_gisjoin_2010
    )
  ),
  
  # Crosswalk 2010-2019 data to 2020 boundaries
  tar_target(
    tract_xwalk,
    {
      xwalk <- sf::st_intersection(
        tracts_touch_rtd_2010 |> dplyr::select(GISJOIN_2010 = GISJOIN),
        tracts_touch_rtd_2020 |> dplyr::select(GISJOIN_2020 = GISJOIN)
      )
      
      xwalk |>
        dplyr::mutate(
          overlap_area = sf::st_area(geometry)
        ) |>
        dplyr::group_by(GISJOIN_2010) |>
        dplyr::mutate(
          tract2010_area = sum(overlap_area),
          weight_2010_to_2020 = overlap_area / tract2010_area
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          overlap_area = units::drop_units(sf::st_area(geometry))
        )
    }
  ),
  
  tar_target(
    denver_tt_2013_bg_est,
    denver_tt_2013_bg %>%
      dplyr::inner_join(bg_xwalk, by = c("GISJOIN" = "GISJOIN_2010")) %>%
      dplyr::mutate(across(where(is.numeric), ~ .x * weight_2010_to_2020)) %>%
      sf::st_drop_geometry() %>%
      dplyr::select(GISJOIN_2020, where(is.numeric)) %>%
      dplyr::group_by(GISJOIN_2020) %>%
      dplyr::summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
      dplyr::rename(GISJOIN = GISJOIN_2020) %>%
      dplyr::filter(GISJOIN %in% study_gisjoin) %>%
      units::drop_units()
  ),
  
  tar_target(
    denver_tt_2013_tract_est,
    denver_tt_2013_tract %>%
      dplyr::left_join(tract_xwalk, by = c("GISJOIN" = "GISJOIN_2010"))  %>%
      dplyr::mutate(across(where(is.numeric),
                           ~ .x * weight_2010_to_2020))  %>%
      sf::st_drop_geometry()  %>%
      dplyr::select(GISJOIN_2020, where(is.numeric))  %>%
      dplyr::group_by(GISJOIN_2020)  %>%
      dplyr::summarise(across(everything(), sum, na.rm = TRUE),
                       .groups = "drop")  %>%
      dplyr::rename(GISJOIN = GISJOIN_2020)  %>%
      dplyr::filter(!is.na(GISJOIN))  %>%
      units::drop_units()
  ),
  
  # ------ Aggregate modal shares ------
  
  tar_target(agg_shares_bg_2013,
             compute_aggregate_shares_bg(denver_tt_2013_bg_est, bg_mode_cols)
  ),
  
  tar_target(agg_shares_bg_2023,
             compute_aggregate_shares_bg(denver_tt_2023_bg, bg_mode_cols)
  ),
  
  tar_target(agg_shares_tract_2013,
             compute_aggregate_shares_tract(denver_tt_2013_tract_est, tract_mode_cols)
  ),
  
  tar_target(agg_shares_tract_2023,
             compute_aggregate_shares_tract(denver_tt_2023_tract, tract_mode_cols)
  ),
  
  # Combine: total stratum from BG, no_car from tract, bike fallback
  tar_target(agg_shares_2013, {
    dplyr::bind_rows(
      agg_shares_bg_2013,
      dplyr::filter(agg_shares_tract_2013, car_stratum == "no_car"),
      dplyr::filter(agg_shares_bg_2013, mode == "bike") |>
        dplyr::mutate(car_stratum = "no_car")
    ) |>
      dplyr::distinct(mode, car_stratum, .keep_all = TRUE)
  }),
  
  tar_target(agg_shares_2023, {
    dplyr::bind_rows(
      agg_shares_bg_2023,
      dplyr::filter(agg_shares_tract_2023, car_stratum == "no_car"),
      dplyr::filter(agg_shares_bg_2023, mode == "bike") |>
        dplyr::mutate(car_stratum = "no_car")
    ) |>
      dplyr::distinct(mode, car_stratum, .keep_all = TRUE)
  }),
  
  # Skill shares from B08119 (tract-level)
  tar_target(agg_shares_skill_2013,
             compute_aggregate_shares_skill(denver_tt_2013_tract_est, skill_mode_cols)
  ),
  
  tar_target(agg_shares_skill_2023,
             compute_aggregate_shares_skill(denver_tt_2023_tract, skill_mode_cols)
  ),
  
  tar_target(agg_shares_skill_pooled, {
    dplyr::bind_rows(
      agg_shares_skill_car_2013 |> dplyr::mutate(weight = 6),
      agg_shares_skill_car_2023 |> dplyr::mutate(weight = 5)
    ) |>
      dplyr::group_by(mode, skill_stratum, car_stratum) |>
      dplyr::summarise(
        S = weighted.mean(S, w = weight, na.rm = TRUE),
        .groups = "drop"
      )
  }),
  
  # IPF reconciliation of skill × car joint distribution
  tar_target(agg_shares_skill_car_2013,
             compute_ipf_aggregate(agg_shares_skill_2013, agg_shares_tract_2013)
  ),
  
  tar_target(agg_shares_skill_car_2023,
             compute_ipf_aggregate(agg_shares_skill_2023, agg_shares_tract_2023)
  ),
  
  # t_bar now needs to be computed for each network state
  tar_target(agg_t_bar_pre_2016,
             compute_aggregate_t_bar(od_clean_pre_2016)),
  tar_target(agg_t_bar_post_2016,
             compute_aggregate_t_bar(od_clean_post_2016)),
  tar_target(agg_t_bar_post_2017,
             compute_aggregate_t_bar(od_clean_post_2017)),
  tar_target(agg_t_bar_post_2019,
             compute_aggregate_t_bar(od_clean_post_2019)),
  tar_target(agg_t_bar_post_2020,
             compute_aggregate_t_bar(od_clean_post_2020)),
  
  # ------ b_m scalars ------
  
  # General (car-stratified only)
  tar_target(agg_t_bar_pooled, {
    dplyr::bind_rows(
      compute_aggregate_t_bar(od_clean_pre_2016)  |> dplyr::mutate(weight = 3),
      compute_aggregate_t_bar(od_clean_post_2016) |> dplyr::mutate(weight = 1),
      compute_aggregate_t_bar(od_clean_post_2017) |> dplyr::mutate(weight = 2),
      compute_aggregate_t_bar(od_clean_post_2019) |> dplyr::mutate(weight = 1),
      compute_aggregate_t_bar(od_clean_post_2020) |> dplyr::mutate(weight = 4)
    ) |>
      dplyr::group_by(mode) |>
      dplyr::summarise(
        t_bar = weighted.mean(t_bar, w = weight, na.rm = TRUE),
        .groups = "drop"
      )
  }),
  
  # Single b_m estimated from pooled t_bar and pooled shares
  tar_target(agg_shares_pooled, {
    dplyr::bind_rows(
      agg_shares_2013 |> dplyr::mutate(weight = 6),
      agg_shares_2023 |> dplyr::mutate(weight = 5)
    ) |>
      dplyr::group_by(mode, car_stratum) |>
      dplyr::summarise(
        S = weighted.mean(S, w = weight, na.rm = TRUE),
        .groups = "drop"
      )
  }),
  
  tar_target(agg_bm_pooled,
             compute_aggregate_bm(
               agg_t_bar_pooled,
               agg_shares_pooled,
               scenario_label = "pooled"
             )
  ),
  
  tar_target(agg_bm_skill_pooled,
             compute_bm_stratified(
               agg_t_bar_pooled,
               agg_shares_skill_pooled,
               scenario_label = "pooled"
             )
  ),
  
  # ------ P(car) and P(low skill) weights ------
  
  tar_target(p_car_2013,
             compute_car_access_prob(denver_tt_2013_bg_est, denver_tt_2013_tract_est)
  ),
  
  tar_target(p_car_2023,
             compute_car_access_prob(denver_tt_2023_bg, denver_tt_2023_tract)
  ),
  
  tar_target(p_low_skill_2013,
             compute_skill_prob(denver_tt_2013_bg_est, denver_tt_2013_tract_est, skill_stratum = "low")),
  
  tar_target(p_low_skill_2023, 
             compute_skill_prob(denver_tt_2023_bg, denver_tt_2023_tract, skill_stratum = "low")),
  
  # ------ Inclusive value of travel time ------
  
  # Then in tar_map, pass all OD targets and select inside:
  tar_map(
    values = IV_SCENARIOS,
    names  = "scenario",
    
    tar_target(iv_public, {
      od <- list(
        pre_2016       = od_clean_pre_2016,
        post_2016      = od_clean_post_2016,
        post_2017      = od_clean_post_2017,
        post_2019      = od_clean_post_2019,
        post_2020      = od_clean_post_2020,
        counterfactual = od_clean_counterfactual
      )[[scenario]]
      compute_inclusive_tt(od, agg_bm_pooled, nest = "public")
    }),
    
    tar_target(iv_private, {
      od <- list(
        pre_2016       = od_clean_pre_2016,
        post_2016      = od_clean_post_2016,
        post_2017      = od_clean_post_2017,
        post_2019      = od_clean_post_2019,
        post_2020      = od_clean_post_2020,
        counterfactual = od_clean_counterfactual
      )[[scenario]]
      compute_inclusive_tt(od, agg_bm_pooled, nest = "private")
    }),
    
    tar_target(iv_reconciled, {
      p_car <- if (p_car_target == "p_car_2013") p_car_2013 else p_car_2023
      inclusive_tt_final(iv_public, iv_private, p_car)
    }),
    
    tar_target(iv_skill_low_public, {
      od <- list(
        pre_2016       = od_clean_pre_2016,
        post_2016      = od_clean_post_2016,
        post_2017      = od_clean_post_2017,
        post_2019      = od_clean_post_2019,
        post_2020      = od_clean_post_2020,
        counterfactual = od_clean_counterfactual
      )[[scenario]]
      compute_inclusive_tt(od, agg_bm_skill_pooled,
                           nest = "public", skill = "low_skill")
    }),
    
    tar_target(iv_skill_low_private, {
      od <- list(
        pre_2016       = od_clean_pre_2016,
        post_2016      = od_clean_post_2016,
        post_2017      = od_clean_post_2017,
        post_2019      = od_clean_post_2019,
        post_2020      = od_clean_post_2020,
        counterfactual = od_clean_counterfactual
      )[[scenario]]
      compute_inclusive_tt(od, agg_bm_skill_pooled,
                           nest = "private", skill = "low_skill")
    }),
    
    tar_target(iv_skill_high_public, {
      od <- list(
        pre_2016       = od_clean_pre_2016,
        post_2016      = od_clean_post_2016,
        post_2017      = od_clean_post_2017,
        post_2019      = od_clean_post_2019,
        post_2020      = od_clean_post_2020,
        counterfactual = od_clean_counterfactual
      )[[scenario]]
      compute_inclusive_tt(od, agg_bm_skill_pooled,
                           nest = "public", skill = "high_skill")
    }),
    
    tar_target(iv_skill_high_private, {
      od <- list(
        pre_2016       = od_clean_pre_2016,
        post_2016      = od_clean_post_2016,
        post_2017      = od_clean_post_2017,
        post_2019      = od_clean_post_2019,
        post_2020      = od_clean_post_2020,
        counterfactual = od_clean_counterfactual
      )[[scenario]]
      compute_inclusive_tt(od, agg_bm_skill_pooled,
                           nest = "private", skill = "high_skill")
    }),
    
    tar_target(iv_skill_low_reconciled, {
      p_car <- if (p_car_target == "p_car_2013") p_car_2013 else p_car_2023
      inclusive_tt_final(
        public_iv  = iv_skill_low_public,
        private_iv = iv_skill_low_private,
        p_car      = p_car
        # no p_low_skill — keeps low skill separate
      )
    }),
    
    tar_target(iv_skill_high_reconciled, {
      p_car <- if (p_car_target == "p_car_2013") p_car_2013 else p_car_2023
      inclusive_tt_final(
        public_iv  = iv_skill_high_public,
        private_iv = iv_skill_high_private,
        p_car      = p_car
        # no p_low_skill — keeps high skill separate
      )
    })
  ),
  
  # ------ LODES OD data ------
  
  tar_target(lodes_od_files,
             list.files(
               "data/LODES_OD",
               pattern = "co_od_main_JT00_20[0-9]{2}\\.csv$",
               full.names = TRUE
             )
  ),
  
  tar_target(lodes_od_list_raw, {
               files <- lodes_od_files
               purrr::map(
                 files,
                 ~ read.csv(.x, stringsAsFactors = FALSE)
               ) |>
                 purrr::set_names(
                   tools::file_path_sans_ext(basename(files))
                 )}
  ),
  
  tar_target(lodes_od_list,
             lapply(lodes_od_list_raw, function (x) process_lodes_OD(x, study_gisjoin))
  ),
  
  # ------ PPML estimation ------
  
  # General
  tar_target(lodes_ppml_set, {
    
    year_to_scenario <- c(
      "2013" = "pre_2016",  "2014" = "pre_2016",  "2015" = "pre_2016",
      "2016" = "post_2016", "2017" = "post_2017",  "2018" = "post_2017",
      "2019" = "post_2019", "2020" = "post_2020",  "2021" = "post_2020",
      "2022" = "post_2020", "2023" = "post_2020"
    )
    
    iv_list <- list(
      pre_2016  = iv_reconciled_pre_2016,
      post_2016 = iv_reconciled_post_2016,
      post_2017 = iv_reconciled_post_2017,
      post_2019 = iv_reconciled_post_2019,
      post_2020 = iv_reconciled_post_2020
    )
    
    purrr::imap_dfr(lodes_od_list, function(lodes, name) {
      yr       <- stringr::str_extract(name, "20[0-9]{2}")
      scenario <- year_to_scenario[yr]
      iv       <- iv_list[[scenario]]
      build_ppml_set(lodes, iv) |>
        dplyr::mutate(year = as.integer(yr)) |>
        dplyr::mutate(
          orig_year = factor(paste0(from_GISJOIN, "_", year)),
          dest_year = factor(paste0(to_GISJOIN,   "_", year))
        )
    })
  }),
  
  tar_target(lodes_ppml_set_low, {
    
    year_to_scenario <- c(
      "2013" = "pre_2016",  "2014" = "pre_2016",  "2015" = "pre_2016",
      "2016" = "post_2016", "2017" = "post_2017",  "2018" = "post_2017",
      "2019" = "post_2019", "2020" = "post_2020",  "2021" = "post_2020",
      "2022" = "post_2020", "2023" = "post_2020"
    )
    
    iv_list <- list(
      pre_2016 = iv_skill_low_reconciled_pre_2016,
      post_2016 = iv_skill_low_reconciled_post_2016,
      post_2017 = iv_skill_low_reconciled_post_2017,
      post_2019 = iv_skill_low_reconciled_post_2019,
      post_2020 = iv_skill_low_reconciled_post_2020
    )
    
    purrr::imap_dfr(lodes_od_list, function(lodes, name) {
      yr       <- stringr::str_extract(name, "20[0-9]{2}")
      scenario <- year_to_scenario[yr]
      iv       <- iv_list[[scenario]]
      build_ppml_set(lodes, iv) |>
        dplyr::mutate(year = as.integer(yr)) |>
        dplyr::mutate(
          orig_year = factor(paste0(from_GISJOIN, "_", year)),
          dest_year = factor(paste0(to_GISJOIN,   "_", year))
        )
    })
  }),
  
  tar_target(lodes_ppml_set_high, {
    
    year_to_scenario <- c(
      "2013" = "pre_2016",  "2014" = "pre_2016",  "2015" = "pre_2016",
      "2016" = "post_2016", "2017" = "post_2017",  "2018" = "post_2017",
      "2019" = "post_2019", "2020" = "post_2020",  "2021" = "post_2020",
      "2022" = "post_2020", "2023" = "post_2020"
    )
    
    iv_list <- list(
      pre_2016 = iv_skill_high_reconciled_pre_2016,
      post_2016 = iv_skill_high_reconciled_post_2016,
      post_2017 = iv_skill_high_reconciled_post_2017,
      post_2019 = iv_skill_high_reconciled_post_2019,
      post_2020 = iv_skill_high_reconciled_post_2020
    )
    
    purrr::imap_dfr(lodes_od_list, function(lodes, name) {
      yr       <- stringr::str_extract(name, "20[0-9]{2}")
      scenario <- year_to_scenario[yr]
      iv       <- iv_list[[scenario]]
      build_ppml_set(lodes, iv) |>
        dplyr::mutate(year = as.integer(yr)) |>
        dplyr::mutate(
          orig_year = factor(paste0(from_GISJOIN, "_", year)),
          dest_year = factor(paste0(to_GISJOIN,   "_", year))
        )
    })
  }),
  
  tar_target(ppml_fit,
             alpaca::feglm(
               total_flow ~ inclusive_traveltime | orig_year + dest_year,
               data   = lodes_ppml_set,
               family = poisson()
             )
  ),
  
  tar_target(ppml_fit_low,
             alpaca::feglm(
               total_flow ~ inclusive_traveltime | orig_year + dest_year,
               data   = lodes_ppml_set_low,
               family = poisson()
             )
  ),
  
  tar_target(ppml_fit_high,
             alpaca::feglm(
               total_flow ~ inclusive_traveltime | orig_year + dest_year,
               data   = lodes_ppml_set_high,
               family = poisson()
             )
  ),
  
  # ------- CMA Estimation ------
  
  tar_target(lodes_wac_2013_raw,
             read.csv("data/LODES_WAC/co_wac_S000_JT00_2013.csv")),
  
  tar_target(L_jobs_2013, {
    lodes_wac_2013_raw |>
      dplyr::mutate(GISJOIN = geoid_to_gisjoin2(w_geocode)) |>
      dplyr::filter(GISJOIN %in% study_gisjoin) |>
      dplyr::group_by(GISJOIN) |>
      dplyr::summarise(Lf = sum(C000, na.rm = TRUE), .groups = "drop") |>
      tibble::deframe()   # named vector
  }),
  
  tar_target(L_residence_2013, {
    denver_acs_2010_list_est[["nhgis0019_ds201_20135_blck_grp"]] |>
      dplyr::select(GISJOIN, B23025_E004) |>
      tibble::deframe()
  }),
  
  # Extract theta from PPML regressions
  tar_target(theta_general,
             extract_theta(ppml_fit)
  ),
  
  tar_target(theta_low_skill,
             extract_theta(ppml_fit_low)
  ),
  
  tar_target(theta_high_skill,
             extract_theta(ppml_fit_high)
  ),
  
  # Compute FCMA
  tar_target(cma_pre_2016,
             compute_fcma(
               od_long     = iv_reconciled_pre_2016,
               L_residence = L_residence_2013,
               L_jobs      = L_jobs_2013,
               theta       = theta_general
             )
  ),
  
  tar_target(cma_post_2016,
             compute_fcma(
               od_long     = iv_reconciled_post_2016,
               L_residence = L_residence_2013,
               L_jobs      = L_jobs_2013,
               theta       = theta_general
             )
  ),
  
  tar_target(cma_post_2017,
             compute_fcma(
               od_long     = iv_reconciled_post_2017,
               L_residence = L_residence_2013,
               L_jobs      = L_jobs_2013,
               theta       = theta_general
             )
  ),
  
  tar_target(cma_post_2019,
             compute_fcma(
               od_long     = iv_reconciled_post_2019,
               L_residence = L_residence_2013,
               L_jobs      = L_jobs_2013,
               theta       = theta_general
             )
  ),
  
  tar_target(cma_post_2020,
             compute_fcma(
               od_long     = iv_reconciled_post_2020,
               L_residence = L_residence_2013,
               L_jobs      = L_jobs_2013,
               theta       = theta_general
             )
  ),
  
  tar_target(cma_counterfactual,
             compute_fcma(
               od_long     = iv_reconciled_counterfactual,
               L_residence = L_residence_2013,
               L_jobs      = L_jobs_2013,
               theta       = theta_general
             )
  ),
  
  # ── Delta ln(CMA) relative to pre-2016 baseline ──────────────────────────────
  
  tar_target(delta_cma_panel, {
    
    # Year to network state mapping
    year_to_scenario <- c(
      "2013" = "pre_2016",  "2014" = "pre_2016",  "2015" = "pre_2016",
      "2016" = "post_2016", "2017" = "post_2017",  "2018" = "post_2017",
      "2019" = "post_2019", "2020" = "post_2020",  "2021" = "post_2020",
      "2022" = "post_2020", "2023" = "post_2020"
    )
    
    # CMA by network state
    cma_list <- list(
      pre_2016  = cma_pre_2016,
      post_2016 = cma_post_2016,
      post_2017 = cma_post_2017,
      post_2019 = cma_post_2019,
      post_2020 = cma_post_2020
    )
    
    # Pre-2016 baseline for log change computation
    baseline <- cma_pre_2016 |>
      dplyr::select(
        GISJOIN,
        RCMA_base = RCMA,
        FCMA_base = FCMA
      )
    
    # Build one row per GISJOIN per year
    purrr::map_dfr(names(year_to_scenario), function(yr) {
      scenario <- year_to_scenario[yr]
      cma      <- cma_list[[scenario]]
      
      cma |>
        dplyr::left_join(baseline, by = "GISJOIN") |>
        dplyr::mutate(
          year             = as.integer(yr),
          scenario         = scenario,
          delta_ln_FCMA    = log(FCMA) - log(FCMA_base),
          delta_ln_RCMA    = log(RCMA) - log(RCMA_base)
        ) |>
        dplyr::select(
          GISJOIN, year, scenario,
          RCMA, FCMA,
          delta_ln_FCMA, delta_ln_RCMA
        )
    })
  }),
  
  # ── Counterfactual delta (relative to post-2020) ──────────────────────────────
  # Kept separate since it's a forward-looking comparison not a panel year
  
  tar_target(delta_cma_cf, {
    post2020_base <- cma_post_2020 |>
      dplyr::select(GISJOIN, RCMA_2020 = RCMA, FCMA_2020 = FCMA)
    
    cma_counterfactual |>
      dplyr::left_join(post2020_base, by = "GISJOIN") |>
      dplyr::mutate(
        delta_ln_FCMA_cf = log(FCMA) - log(FCMA_2020),
        delta_ln_RCMA_cf = log(RCMA) - log(RCMA_2020)
      ) |>
      dplyr::select(GISJOIN, delta_ln_FCMA_cf, delta_ln_RCMA_cf)
  }),
  
  # ------ FCMA Regression Analysis ------
  
  tar_target(cma_regression_sets, {
    
    base <- regression_set[["continuous"]] |>
      dplyr::filter(!is.na(distance_m)) |>
      dplyr::left_join(
        delta_cma_panel |> dplyr::mutate(GISJOIN = trimws(GISJOIN)),
        by = c("GISJOIN", "year")
      ) |>
      dplyr::left_join(
        delta_cma_cf |> dplyr::mutate(GISJOIN = trimws(GISJOIN)),
        by = "GISJOIN",
        relationship = "many-to-one"
      ) |>
      dplyr::filter(!is.na(delta_ln_FCMA), !is.na(delta_ln_RCMA)) |>
      dplyr::mutate(
        post      = as.integer(year >= cohort_year),
        treat     = as.integer(group == "treatment"),
        post_2020 = as.integer(year >= 2020)
      )
    
    list(
      fcma = list(
        jobs        = base %>%
          dplyr::filter(
            !is.infinite(ljobs), !is.na(ljobs)
          ) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup(),
        jobs_per_100 = base %>%
          dplyr::filter(
            !is.infinite(ljobs_per_100), !is.na(ljobs_per_100)
          ) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup(),
        wage        = base %>%
          dplyr::filter(
            !is.infinite(lwage), !is.na(lwage)
          ) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup(),
        estab       = base |>
          dplyr::filter(
            !is.infinite(lestabs_entry), !is.na(lestabs_entry)
          ) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup()
      ),
      rcma = list(
        pop     = base %>%
          dplyr::filter(!is.infinite(lpop), !is.na(lpop)) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup(),
        income = base %>%
          dplyr::filter(!is.infinite(lmhi), !is.na(lmhi)) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup(),
        housing = base %>%
          dplyr::filter(!is.infinite(lhousing), !is.na(lhousing)) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup(),
        rent    = base %>%
          dplyr::filter(!is.infinite(lmedian_rent), !is.na(lmedian_rent)) %>%
          dplyr::group_by(GISJOIN) %>%
          dplyr::filter(dplyr::n_distinct(year) == 11) %>%
          dplyr::ungroup()
      )
    )
  }),
  
  tar_target(fcma_regs_general,
             purrr::imap(
               cma_regression_sets[["fcma"]],
               function(data, outcome_name) {
                 run_cma_regressions(
                   data,
                   delta_var    = "delta_ln_FCMA",
                   delta_cf_var = "delta_ln_FCMA_cf",
                   outcome      = outcome_name
                 )
               }
             )
  ),
  
  tar_target(rcma_regs_general,
             purrr::imap(
               cma_regression_sets[["rcma"]],
               function(data, outcome_name) {
                 run_cma_regressions(
                   data,
                   delta_var    = "delta_ln_RCMA",
                   delta_cf_var = "delta_ln_RCMA_cf",
                   outcome      = outcome_name
                 )
               }
             )
  ),
  
  tar_target(fcma_boot_general,
             purrr::imap(
               cma_regression_sets[["fcma"]],
               function(data, outcome_name) {
                 run_cma_bootstrap(
                   data,
                   delta_var    = "delta_ln_FCMA",
                   delta_cf_var = "delta_ln_FCMA_cf",
                   outcome      = outcome_name,
                   n_boot       = 200
                 )
               }
             )
  ),
  
  tar_target(rcma_boot_general,
             purrr::imap(
               cma_regression_sets[["rcma"]],
               function(data, outcome_name) {
                 run_cma_bootstrap(
                   data,
                   delta_var    = "delta_ln_RCMA",
                   delta_cf_var = "delta_ln_RCMA_cf",
                   outcome      = outcome_name,
                   n_boot       = 200
                 )
               }
             )
  ),
  
  # ------ Export CMA Result Tables ------
  
  # ------ FCMA tables ------
  
  tar_target(table_fcma_base, {
    make_cma_table(fcma_boot_general, spec = "base", measure = "fcma") |>
      export_cma_table("table_fcma_base.tex", spec = "base", measure = "fcma")
  }),
  
  tar_target(table_fcma_treatment, {
    make_cma_table(fcma_boot_general, spec = "treatment", measure = "fcma") |>
      export_cma_table("table_fcma_treatment.tex", spec = "treatment", measure = "fcma")
  }),
  
  tar_target(table_fcma_high_rent, {
    make_cma_table(fcma_boot_general, spec = "high_rent", measure = "fcma") |>
      export_cma_table("table_fcma_high_rent.tex", spec = "high_rent", measure = "fcma")
  }),
  
  # ------ RCMA tables ------
  
  tar_target(table_rcma_base, {
    make_cma_table(rcma_boot_general, spec = "base", measure = "rcma") |>
      export_cma_table("table_rcma_base.tex", spec = "base", measure = "rcma")
  }),
  
  tar_target(table_rcma_treatment, {
    make_cma_table(rcma_boot_general, spec = "treatment", measure = "rcma") |>
      export_cma_table("table_rcma_treatment.tex", spec = "treatment", measure = "rcma")
  }),
  
  tar_target(table_rcma_high_rent, {
    make_cma_table(rcma_boot_general, spec = "high_rent", measure = "rcma") |>
      export_cma_table("table_rcma_high_rent.tex", spec = "high_rent", measure = "rcma")
  }),
  
  # ------ FCMA coef plots ------
  
  tar_target(plot_fcma_base, {
    p <- extract_coef_df(fcma_boot_general, spec = "base", measure = "fcma") |>
      make_coef_plot(
        title    = "FCMA Effects on Labor Market Outcomes",
        ci_limit = 2
      )
    save_plot(p, "coef_fcma_base.png", width = 12, height = 7)
    p
  }),
  
  tar_target(plot_fcma_treatment, {
    p <- extract_coef_df(fcma_boot_general, spec = "treatment", measure = "fcma") |>
      make_coef_plot(
        title    = "FCMA Treatment Interaction: Labor Market Outcomes",
        ci_limit = 2
      )
    save_plot(p, "coef_fcma_treatment.png", width = 12, height = 7)
    p
  }),
  
  tar_target(plot_fcma_high_rent, {
    p <- extract_coef_df(fcma_boot_general, spec = "high_rent", measure = "fcma") |>
      make_coef_plot(
        title    = "FCMA High-Rent Heterogeneity: Labor Market Outcomes",
        ci_limit = 2
      )
    save_plot(p, "coef_fcma_high_rent.png", width = 12, height = 7)
    p
  }),
  
  # ------ RCMA coef plots ------
  
  tar_target(plot_rcma_base, {
    p <- extract_coef_df(rcma_boot_general, spec = "base", measure = "rcma") |>
      make_coef_plot(
        title    = "RCMA Effects on Residential Outcomes",
        ci_limit = 2
      )
    save_plot(p, "coef_rcma_base.png", width = 12, height = 7)
    p
  }),
  
  tar_target(plot_rcma_treatment, {
    p <- extract_coef_df(rcma_boot_general, spec = "treatment", measure = "rcma") |>
      make_coef_plot(
        title    = "RCMA Treatment Interaction: Residential Outcomes",
        ci_limit = 2
      )
    save_plot(p, "coef_rcma_treatment.png", width = 12, height = 7)
    p
  }),
  
  tar_target(plot_rcma_high_rent, {
    p <- extract_coef_df(rcma_boot_general, spec = "high_rent", measure = "rcma") |>
      make_coef_plot(
        title    = "RCMA High-Rent Heterogeneity: Residential Outcomes",
        ci_limit = 2
      )
    save_plot(p, "coef_rcma_high_rent.png", width = 12, height = 7)
    p
  })
)

