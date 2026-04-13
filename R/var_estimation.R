library(sf)
library(tibble)
library(alpaca)
library(dplyr)
library(tidyr)
library(stringr)

ipf_county <- function(row_totals, col_totals, prior, max_iter = 50, tol = 1e-8){
  
  m <- prior * sum(row_totals)
  
  for(i in seq_len(max_iter)){
    
    m_old <- m
    
    m <- sweep(m, 1, rowSums(m), "/")
    m <- sweep(m, 1, row_totals, "*")
    
    m <- sweep(m, 2, colSums(m), "/")
    m <- sweep(m, 2, col_totals, "*")
    
    if(max(abs(m - m_old), na.rm = TRUE) < tol) break
  }
  
  m
}

compute_weighted_median <- function(x, w) {
  # Remove NA pairs
  keep <- !is.na(x) & !is.na(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  
  # Sort by value
  ord <- order(x)
  x   <- x[ord]
  w   <- w[ord]
  
  # Find weighted median
  cumw <- cumsum(w) / sum(w)
  x[which(cumw >= 0.5)[1]]
}

var_estimation <- function(df, wage) {
  
  baseline_2013 <- df |>
    dplyr::filter(year == 2013) |>
    dplyr::select(
      GISJOIN,
      baseline_rent   = B25064_E001,
      renter_units    = B25003_E003
    ) |>
    dplyr::mutate(
      baseline_rent = dplyr::if_else(
        is.na(baseline_rent) | baseline_rent <= 0,
        NA_real_,
        as.numeric(baseline_rent)
      ),
      renter_units = dplyr::if_else(
        is.na(renter_units) | renter_units <= 0,
        NA_real_,
        as.numeric(renter_units)
      )
    )
  
  rent_threshold <- compute_weighted_median(
    baseline_2013$baseline_rent,
    baseline_2013$renter_units
  )
  
  df <- df %>%
    left_join(wage, by = c("year","GISJOIN")) %>%
    left_join(baseline_2013, by = "GISJOIN") %>%
    mutate(
      ljobs    = log1p(total_jobs),
      lpop     = log1p(B01003_E001),
      lmhi     = log1p(B19013_E001),
      lhousing = log1p(B25003_E001),
      lmedian_rent = log1p(B25064_E001),
      ljobs_per_100 = log1p(total_jobs / B23025_E004),
      lestabs_entry = log1p(estabs_entry_total),
      estabs_small_share = estabs_small / estabs_total,
      estabs_mid_share = estabs_mid / estabs_total,
      estabs_large_share = estabs_large / estabs_total,
      high_rent_int = as.integer(baseline_rent > rent_threshold)
    ) %>%
    select(GISJOIN, year, starts_with("l"), starts_with("estabs"),
           total_jobs, high_rent_int, nearest_station_id, distance_m, cohort_year,
           group, any_of(c("log_dist","ref_year")))
  
  df
}

# ------ Modal share lookup tables ------

bg_mode_cols <- list(
  drive   = "B08301_E003",
  transit = "B08301_E010",
  bike    = "B08301_E018",
  walk    = "B08301_E019"
)

tract_mode_cols <- list(
  total = list(
    drive   = "B08141_E006",
    transit = "B08141_E016",
    walk    = "B08141_E021"
  ),
  no_car = list(
    transit = "B08141_E017",
    bike    = "B08141_E027",
    walk    = "B08141_E022"
  )
)

skill_mode_cols <- list(
  low_skill = list(
    denom   = c("B08119_E002","B08119_E003","B08119_E004"),
    drive   = c("B08119_E011","B08119_E012","B08119_E013"),
    transit = c("B08119_E029","B08119_E030","B08119_E031"),
    walk    = c("B08119_E038","B08119_E039","B08119_E040"),
    bike    = c("B08119_E047","B08119_E048","B08119_E049")
  ),
  high_skill = list(
    denom   = c("B08119_E007","B08119_E008","B08119_E009"),
    drive   = c("B08119_E016","B08119_E017","B08119_E018"),
    transit = c("B08119_E034","B08119_E035","B08119_E036"),
    walk    = c("B08119_E043","B08119_E044","B08119_E045"),
    bike    = c("B08119_E052","B08119_E053","B08119_E054")
  )
)

# ------ Aggregate share functions (study-area scalars) ------

compute_aggregate_shares_bg <- function(df, mode_cols) {
  total_workers <- sum(df[["B08301_E001"]] - df[["B08301_E021"]], na.rm = TRUE)
  
  purrr::map_dfr(names(mode_cols), function(m) {
    tibble::tibble(
      mode        = m,
      car_stratum = "total",
      S           = sum(df[[mode_cols[[m]]]], na.rm = TRUE) / total_workers
    )
  })
}

compute_aggregate_shares_tract <- function(tract_df, mode_cols) {
  purrr::map_dfr(names(mode_cols), function(stratum) {
    cols  <- mode_cols[[stratum]]
    denom <- if (stratum == "total") {
      sum(tract_df[["B08141_E001"]] - tract_df[["B08141_E031"]], na.rm = TRUE)
    } else {
      sum(tract_df[["B08141_E002"]] - tract_df[["B08141_E032"]], na.rm = TRUE)
    }
    purrr::map_dfr(names(cols), function(m) {
      tibble::tibble(
        mode        = m,
        car_stratum = stratum,
        S           = sum(tract_df[[cols[[m]]]], na.rm = TRUE) / denom
      )
    })
  })
}

# Aggregate skill shares from B08119 — returns mode × skill_stratum scalars
compute_aggregate_shares_skill <- function(tract_df, skill_mode_cols) {
  purrr::map_dfr(names(skill_mode_cols), function(skill) {
    cols  <- skill_mode_cols[[skill]]
    
    # Exclude WFH from denominator
    wfh_cols <- paste0("B08119_E0",
                       switch(skill,
                              low_skill  = c("56","57","58"),
                              high_skill = c("61","62","63")))
    denom <- sum(rowSums(tract_df[, cols$denom], na.rm = TRUE)) -
      sum(rowSums(tract_df[, wfh_cols],   na.rm = TRUE))
    
    modes <- setdiff(names(cols), "denom")
    purrr::map_dfr(modes, function(m) {
      tibble::tibble(
        mode          = m,
        skill_stratum = skill,
        S             = sum(rowSums(tract_df[, cols[[m]]], na.rm = TRUE)) / denom
      )
    })
  })
}

# IPF reconciliation — aggregate version (scalar inputs, scalar outputs)
compute_ipf_aggregate <- function(skill_shares, car_shares,
                                  max_iter = 100, tol = 1e-8) {
  modes        <- c("drive", "transit", "walk", "bike")
  skill_strata <- c("low_skill", "high_skill")
  car_strata   <- c("total", "no_car")
  
  seed <- matrix(0, nrow = length(modes), ncol = length(skill_strata),
                 dimnames = list(modes, skill_strata))
  for (m in modes) for (s in skill_strata) {
    val <- skill_shares$S[skill_shares$mode == m &
                            skill_shares$skill_stratum == s]
    seed[m, s] <- if (length(val) > 0) val else 0
  }
  
  margin <- matrix(0, nrow = length(modes), ncol = length(car_strata),
                   dimnames = list(modes, car_strata))
  for (m in modes) for (c in car_strata) {
    val <- car_shares$S[car_shares$mode == m & car_shares$car_stratum == c]
    margin[m, c] <- if (length(val) > 0) val else 0
  }
  
  fitted <- ipf_shares(seed, margin, max_iter, tol)
  
  purrr::map_dfr(modes, function(m)
    purrr::map_dfr(skill_strata, function(s)
      purrr::map_dfr(car_strata, function(c)
        tibble::tibble(mode = m, skill_stratum = s,
                       car_stratum = c, S = fitted[m, s, c])
      )
    )
  )
}

# ------ IPF core (shared by aggregate and BG-level) ------

ipf_shares <- function(seed_skill, margin_car, max_iter = 100, tol = 1e-8) {
  modes  <- rownames(seed_skill)
  skills <- colnames(seed_skill)
  cars   <- colnames(margin_car)
  
  arr <- array(0,
               dim      = c(length(modes), length(skills), length(cars)),
               dimnames = list(mode = modes, skill = skills, car = cars)
  )
  for (m in modes) for (s in skills) for (c in cars)
    arr[m, s, c] <- seed_skill[m, s] * margin_car[m, c]
  
  for (iter in seq_len(max_iter)) {
    arr_old <- arr
    for (m in modes) for (c in cars) {
      cur <- sum(arr[m, , c])
      if (cur > 0) arr[m, , c] <- arr[m, , c] * (margin_car[m, c] / cur)
    }
    for (m in modes) for (s in skills) {
      cur <- sum(arr[m, s, ])
      if (cur > 0) arr[m, s, ] <- arr[m, s, ] * (seed_skill[m, s] / cur)
    }
    if (max(abs(arr - arr_old)) < tol) {
      message("IPF converged at iteration ", iter); break
    }
  }
  arr
}

# BG-level IPF shares (spatial, broadcast from tract)
compute_shares_skill_car <- function(bg_df, tract_skill_df, tract_car_df,
                                     skill_mode_cols, car_mode_cols) {
  bg_with_tract <- bg_df |>
    dplyr::mutate(tract_GISJOIN = substr(GISJOIN, 1, nchar(GISJOIN) - 1))
  
  purrr::map_dfr(unique(bg_with_tract$tract_GISJOIN), function(tr) {
    tract_skill_row <- tract_skill_df |> dplyr::filter(GISJOIN == tr)
    tract_car_row   <- tract_car_df   |> dplyr::filter(GISJOIN == tr)
    if (nrow(tract_skill_row) == 0 | nrow(tract_car_row) == 0) return(NULL)
    
    modes  <- c("drive", "transit", "walk", "bike")
    skills <- c("low_skill", "high_skill")
    cars   <- c("total", "no_car")
    
    seed_skill <- matrix(0, nrow = length(modes), ncol = length(skills),
                         dimnames = list(modes, skills))
    for (skill in skills) {
      cols  <- skill_mode_cols[[skill]]
      denom <- sum(as.numeric(tract_skill_row[, cols$denom]), na.rm = TRUE)
      for (m in modes) {
        num <- sum(as.numeric(tract_skill_row[, cols[[m]]]), na.rm = TRUE)
        seed_skill[m, skill] <- if (denom > 0) num / denom else 0
      }
    }
    
    margin_car <- matrix(0, nrow = length(modes), ncol = length(cars),
                         dimnames = list(modes, cars))
    denom_total  <- as.numeric(tract_car_row[["B08141_E001"]]) -
      as.numeric(tract_car_row[["B08141_E031"]])
    denom_no_car <- as.numeric(tract_car_row[["B08141_E002"]]) -
      as.numeric(tract_car_row[["B08141_E032"]])
    car_cols <- list(
      total  = list(drive = "B08141_E006", transit = "B08141_E016",
                    walk  = "B08141_E021", bike    = NA),
      no_car = list(drive = "B08141_E007", transit = "B08141_E017",
                    walk  = "B08141_E022", bike    = "B08141_E027")
    )
    for (m in modes) for (car in cars) {
      col   <- car_cols[[car]][[m]]
      denom <- if (car == "total") denom_total else denom_no_car
      if (!is.na(col) && denom > 0)
        margin_car[m, car] <- as.numeric(tract_car_row[[col]]) / denom
    }
    
    fitted       <- ipf_shares(seed_skill, margin_car)
    bgs_in_tract <- bg_with_tract |>
      dplyr::filter(tract_GISJOIN == tr) |> dplyr::pull(GISJOIN)
    
    purrr::map_dfr(bgs_in_tract, function(bg)
      purrr::map_dfr(modes,  function(m)
        purrr::map_dfr(skills, function(s)
          purrr::map_dfr(cars,   function(c)
            tibble::tibble(GISJOIN = bg, mode = m,
                           skill_stratum = s, car_stratum = c,
                           S = fitted[m, s, c])
          )
        )
      )
    )
  })
}

# ------ b_m computation ------

# Aggregate b_m — study-area scalars, grouped by car_stratum only
compute_aggregate_bm <- function(t_bar, shares, scenario_label, kappa = 0.011) {
  modes  <- c("drive", "transit", "bike")
  t_walk <- t_bar$t_bar[t_bar$mode == "walk"]
  
  bm_df <- shares |>
    dplyr::filter(mode %in% c(modes, "walk")) |>
    dplyr::group_by(car_stratum) |>
    dplyr::group_modify(function(stratum_shares, key) {
      S_walk <- stratum_shares$S[stratum_shares$mode == "walk"]
      purrr::map_dfr(modes, function(m) {
        tibble::tibble(
          mode      = m,
          scenario  = scenario_label,
          log_share = log(max(stratum_shares$S[stratum_shares$mode == m], 1e-8) /
                            max(S_walk, 1e-8)),
          t_diff    = t_bar$t_bar[t_bar$mode == m] - t_walk,
          b_m       = log_share + kappa * t_diff
        )
      })
    }) |>
    dplyr::ungroup() |>
    dplyr::filter(!(mode == "drive" & car_stratum == "no_car")) |>
    dplyr::group_by(car_stratum) |>
    dplyr::mutate(C = -log(sum(exp(b_m)))) |>
    dplyr::ungroup() |>
    dplyr::mutate(b_m = b_m + C)
  
  walk_rows <- bm_df |>
    dplyr::distinct(car_stratum, C) |>
    dplyr::mutate(mode = "walk", scenario = scenario_label,
                  log_share = 0, t_diff = 0, b_m = C) |>
    dplyr::select(-C)
  
  bm_df |> dplyr::select(-C) |>
    dplyr::bind_rows(walk_rows) |>
    dplyr::arrange(car_stratum, mode)
}

# Skill-stratified b_m — study-area scalars, grouped by skill × car_stratum
compute_bm_stratified <- function(t_bar, shares, scenario_label, kappa = 0.011) {
  modes        <- c("drive", "transit", "bike")
  skill_strata <- c("low_skill", "high_skill")
  car_strata   <- c("total", "no_car")
  t_walk       <- t_bar$t_bar[t_bar$mode == "walk"]
  
  bm_df <- shares |>
    dplyr::filter(mode %in% c(modes, "walk")) |>
    dplyr::group_by(skill_stratum, car_stratum) |>
    dplyr::group_modify(function(stratum_shares, key) {
      S_walk <- stratum_shares$S[stratum_shares$mode == "walk"]
      purrr::map_dfr(modes, function(m) {
        if (m == "drive" && key$car_stratum == "no_car") return(NULL)
        tibble::tibble(
          mode      = m,
          scenario  = scenario_label,
          log_share = log(max(stratum_shares$S[stratum_shares$mode == m], 1e-8) /
                            max(S_walk, 1e-8)),
          t_diff    = t_bar$t_bar[t_bar$mode == m] - t_walk,
          b_m       = log_share + kappa * t_diff
        )
      })
    }) |>
    dplyr::ungroup() |>
    dplyr::group_by(skill_stratum, car_stratum) |>
    dplyr::mutate(C = -log(sum(exp(b_m)))) |>
    dplyr::ungroup() |>
    dplyr::mutate(b_m = b_m + C)
  
  walk_rows <- bm_df |>
    dplyr::group_by(skill_stratum, car_stratum) |>
    dplyr::summarise(C = dplyr::first(C), .groups = "drop") |>
    dplyr::mutate(
      mode      = "walk",
      scenario  = scenario_label,
      log_share = 0,
      t_diff    = 0,
      b_m       = C
    ) |>
    dplyr::select(-C)
  
  bm_df |> dplyr::select(-C) |>
    dplyr::bind_rows(walk_rows) |>
    dplyr::arrange(skill_stratum, car_stratum, mode)
}

# ------ Inclusive value computation ------

# Now handles both aggregate (no skill arg) and stratified (with skill arg)
compute_inclusive_tt <- function(tt_od, pref_shifter, nest = "private",
                                 skill = NULL, kappa = 0.011) {
  car_filter <- if (nest == "public") "no_car" else "total"
  
  prefs <- pref_shifter |>
    dplyr::filter(car_stratum == car_filter)
  
  if (!is.null(skill)) {
    prefs <- prefs |> dplyr::filter(skill_stratum == skill)
  }
  
  prefs <- prefs |> dplyr::select(mode, b_m)
  
  tt_od |>
    dplyr::left_join(prefs, by = "mode") |>
    dplyr::group_by(from_GISJOIN, to_GISJOIN) |>
    dplyr::summarize(
      iv_tt   = -1/kappa * log(sum(exp(b_m - kappa * TotalTravelTime),
                                   na.rm = TRUE)),
      .groups = "drop"
    )
}

# ------ IV reconciliation ------

# Unified final IV — handles both aggregate (no skill weighting) and stratified
inclusive_tt_final <- function(public_iv, private_iv, p_car,
                               public_iv_high = NULL, private_iv_high = NULL,
                               p_low_skill = NULL) {
  
  base <- dplyr::left_join(
    public_iv  |> dplyr::rename(tt_public  = iv_tt),
    private_iv |> dplyr::rename(tt_private = iv_tt),
    by = c("from_GISJOIN", "to_GISJOIN")
  )
  
  # Aggregate (no skill dimension)
  if (is.null(p_low_skill)) {
    return(
      base |>
        dplyr::left_join(p_car, by = c("from_GISJOIN" = "GISJOIN")) |>
        dplyr::mutate(
          inclusive_traveltime = p_car * tt_private + (1 - p_car) * tt_public
        ) |>
        dplyr::select(from_GISJOIN, to_GISJOIN, inclusive_traveltime)
    )
  }
  
  # Skill-stratified — weight across four strata
  high_base <- dplyr::left_join(
    public_iv_high  |> dplyr::rename(tt_public_high  = iv_tt),
    private_iv_high |> dplyr::rename(tt_private_high = iv_tt),
    by = c("from_GISJOIN", "to_GISJOIN")
  )
  
  base |>
    dplyr::left_join(high_base, by = c("from_GISJOIN", "to_GISJOIN")) |>
    dplyr::left_join(p_car,       by = c("from_GISJOIN" = "GISJOIN")) |>
    dplyr::left_join(p_low_skill, by = c("from_GISJOIN" = "GISJOIN")) |>
    dplyr::mutate(
      inclusive_traveltime =
        p_low_skill       * p_car       * tt_private      +   # low, private
        p_low_skill       * (1-p_car)   * tt_public       +   # low, public
        (1-p_low_skill)   * p_car       * tt_private_high +   # high, private
        (1-p_low_skill)   * (1-p_car)   * tt_public_high      # high, public
    ) |>
    dplyr::select(from_GISJOIN, to_GISJOIN, inclusive_traveltime)
}

# ------ PPML helpers (unchanged) ------

compute_aggregate_t_bar <- function(od_clean) {
  od_clean |>
    dplyr::group_by(mode) |>
    dplyr::summarise(t_bar = mean(TotalTravelTime, na.rm = TRUE), .groups = "drop")
}

compute_car_access_prob <- function(bg_df, tract_df) {
  tract_p_car <- tract_df |>
    dplyr::mutate(p_car = (B08141_E001 - B08141_E002) / B08141_E001) |>
    dplyr::select(tract_GISJOIN = GISJOIN, p_car)
  
  bg_df |>
    dplyr::select(GISJOIN) |>
    dplyr::mutate(tract_GISJOIN = substr(GISJOIN, 1, nchar(GISJOIN) - 1)) |>
    dplyr::left_join(tract_p_car, by = "tract_GISJOIN") |>
    dplyr::select(GISJOIN, p_car)
}

compute_skill_prob <- function(bg_df, tract_df, skill_stratum) {
  
  if(skill_stratum == "low") {
    tract_p_skill <- tract_df |>
      dplyr::mutate(p_low_skill = (B08119_E002 + B08119_E003 + B08119_E004 + B08119_E005 + B08119_E006) / B08119_E001) |>
      dplyr::select(tract_GISJOIN = GISJOIN, p_low_skill)
  } else {
    tract_p_skill <- tract_df |>
      dplyr::mutate(p_high_skill = (B08119_E007 + B08119_E008 + B08119_E009) / B08119_E001) |>
      dplyr::select(tract_GISJOIN = GISJOIN, p_high_skill)
  }
  
  bg_df |>
    dplyr::select(GISJOIN) |>
    dplyr::mutate(tract_GISJOIN = substr(GISJOIN, 1, nchar(GISJOIN) - 1)) |>
    dplyr::left_join(tract_p_skill, by = "tract_GISJOIN") |>
    dplyr::select(GISJOIN, any_of(c("p_low_skill","p_high_skill")))
}

# Define which OD clean target each scenario uses
# by creating a helper that switches on scenario name
get_od_clean <- function(scenario) {
  switch(scenario,
         pre_2016      = od_clean_pre_2016,
         post_2016     = od_clean_post_2016,
         post_2017     = od_clean_post_2017,
         post_2019     = od_clean_post_2019,
         post_2020     = od_clean_post_2020,
         counterfactual = od_clean_counterfactual
  )
}

build_ppml_set <- function(flows, tt, id = c("from_GISJOIN","to_GISJOIN")) {
  dplyr::full_join(flows, tt, by = id) |>
    dplyr::filter(from_GISJOIN != to_GISJOIN) |>
    dplyr::mutate(
      inclusive_traveltime = dplyr::case_when(
        is.infinite(inclusive_traveltime) ~ 120,
        is.na(inclusive_traveltime)       ~ 120,
        inclusive_traveltime > 120        ~ 120,
        inclusive_traveltime >= 0         ~ inclusive_traveltime,
        inclusive_traveltime < 0          ~ 0.5,
        TRUE                              ~ NA_real_
      ),
      total_flow = dplyr::coalesce(total_flow, 0L)
    )
}

elasticity_ppml <- function(data) {
  alpaca::feglm(
    total_flow ~ inclusive_traveltime | from_GISJOIN + to_GISJOIN,
    data   = data,
    family = poisson()
  )
}

# helpers/fcma.R

extract_theta <- function(ppml_fit, kappa = 0.011) {
  beta  <- ppml_fit$coefficients[["inclusive_traveltime"]]
  theta <- -beta / kappa
  message("Extracted theta: ", round(theta, 4))
  theta
}

build_od_matrix <- function(od_long, id_col = "from_GISJOIN",
                            jd_col = "to_GISJOIN",
                            val_col = "inclusive_traveltime") {
  all_ids <- sort(union(od_long[[id_col]], od_long[[jd_col]]))
  
  mat <- od_long |>
    tidyr::pivot_wider(
      id_cols     = dplyr::all_of(id_col),
      names_from  = dplyr::all_of(jd_col),
      values_from = dplyr::all_of(val_col),
      values_fill = 120
    ) |>
    tibble::column_to_rownames(id_col) |>
    as.matrix()
  
  # Add any missing rows (origins not in data)
  missing_rows <- setdiff(all_ids, rownames(mat))
  if (length(missing_rows) > 0) {
    empty <- matrix(NA_real_, nrow = length(missing_rows), ncol = ncol(mat),
                    dimnames = list(missing_rows, colnames(mat)))
    mat <- rbind(mat, empty)
  }
  
  # Add any missing cols (destinations not in data)
  missing_cols <- setdiff(all_ids, colnames(mat))
  if (length(missing_cols) > 0) {
    empty <- matrix(NA_real_, nrow = nrow(mat), ncol = length(missing_cols),
                    dimnames = list(rownames(mat), missing_cols))
    mat <- cbind(mat, empty)
  }
  
  # Now safe to reindex
  mat <- mat[all_ids, all_ids]
  mat[is.na(mat)] <- 0
  diag(mat) <- 0
  mat
}

compute_cost_matrix <- function(t_matrix, theta, kappa = 0.011) {
  # C_ij = exp(-theta * kappa * t_ij)
  # NA travel times -> 0 cost (excluded from sum)
  cost <- exp(-theta * kappa * t_matrix)
  cost[is.na(cost)] <- 0
  diag(cost) <- 0
  cost
}

solve_cma <- function(cost_matrix, L_residence, L_jobs,
                      max_iter = 1000, tol = 1e-8) {
  n   <- nrow(cost_matrix)
  ids <- rownames(cost_matrix)
  
  Lr <- as.numeric(L_residence[ids])
  Lf <- as.numeric(L_jobs[ids])
  Lr[is.na(Lr)] <- 0
  Lf[is.na(Lf)] <- 0
  
  # Validate inputs before iterating
  if (all(Lr == 0)) stop("All residential population values are zero or NA")
  if (all(Lf == 0)) stop("All job counts are zero or NA")
  if (any(is.na(cost_matrix))) {
    warning("NA values in cost matrix — setting to 0")
    cost_matrix[is.na(cost_matrix)] <- 0
  }
  
  RCMA <- rep(1, n)
  FCMA <- rep(1, n)
  
  for (iter in seq_len(max_iter)) {
    RCMA_old <- RCMA
    FCMA_old <- FCMA
    
    RCMA_new <- as.vector(cost_matrix %*% (Lf / pmax(FCMA, 1e-10)))
    FCMA_new <- as.vector(t(cost_matrix) %*% (Lr / pmax(RCMA_new, 1e-10)))
    
    # Guard against NA/NaN/Inf before updating
    if (any(!is.finite(RCMA_new)) || any(!is.finite(FCMA_new))) {
      warning("Non-finite values at iteration ", iter, " — stopping early")
      break
    }
    
    RCMA <- RCMA_new
    FCMA <- FCMA_new
    
    # Normalize
    rcma_mean <- mean(RCMA[RCMA > 0])
    fcma_mean <- mean(FCMA[FCMA > 0])
    if (is.finite(rcma_mean) && rcma_mean > 0) RCMA <- RCMA / rcma_mean
    if (is.finite(fcma_mean) && fcma_mean > 0) FCMA <- FCMA / fcma_mean
    
    delta <- max(abs(FCMA - FCMA_old), abs(RCMA - RCMA_old), na.rm = TRUE)
    
    if (iter %% 50 == 0) message("Iter ", iter, " | delta: ", round(delta, 8))
    if (is.finite(delta) && delta < tol) {
      message("Converged at iteration ", iter)
      break
    }
  }
  
  tibble::tibble(GISJOIN = ids, RCMA = RCMA, FCMA = FCMA)
}

compute_fcma <- function(od_long, L_residence, L_jobs, theta, kappa = 0.011,
                         max_iter = 1000, tol = 1e-8) {
  
  message("Building cost matrix...")
  t_mat    <- build_od_matrix(od_long)
  cost_mat <- compute_cost_matrix(t_mat, theta, kappa)
  
  message("Solving CMA system...")
  solve_cma(cost_mat, L_residence, L_jobs, max_iter, tol)
}
