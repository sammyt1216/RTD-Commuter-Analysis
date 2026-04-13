library(fixest)
library(furrr)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(modelsummary)
library(kableExtra)
library(ggplot2)
library(viridis)

run_regression <- function(data,
                           outcomes = c(
                             "lpop", "lmhi", "lwage", "lhousing",
                             "ljobs", "ljobs_per_100", "lmedian_rent",
                             "lestabs", "lestabs_entry", "lestabs_small",
                             "lestabs_mid", "lestabs_large"
                           ),
                           controls){
  
  if(controls == "continuous") {
    fml <- as.formula(
      paste0(
        "c(", paste(outcomes, collapse = ", "), ")",
        " ~ i(ref_year, log_dist, ref = -1) | GISJOIN + year"
      )
    )
  } else {
    fml <- as.formula(
      paste0(
        "c(", paste(outcomes, collapse = ", "), ")",
        " ~ sunab(cohort_year, year, ref.p = -1) | GISJOIN + year"
      )
    )
  }
  
  fixest::feols(fml, data = data)
}

extract_model_atts <- function(model, design, outcome) {
  coefs <- if (design == "continuous") {
    summary(model)$coeftable
  } else {
    summary(model, agg = "ATT")$coeftable
  }
  
  result <- tibble::tibble(
    term      = rownames(coefs),
    estimate  = coefs[, "Estimate"],
    std.error = coefs[, "Std. Error"],
    p.value   = coefs[, "Pr(>|t|)"],
    design    = design,
    outcome   = outcome
  )
  
  # Add reference period (k = -1) as explicit zero row for sunab
  if (design %in% c("treatment", "spillover")) {
    ref_row <- tibble::tibble(
      term      = "year::-1",   # matches fixest's term naming convention
      estimate  = 0,
      std.error = 0,
      p.value   = NA_real_,
      design    = design,
      outcome   = outcome
    )
    result <- dplyr::bind_rows(result, ref_row)
  }
  
  result
}

extract_design_atts <- function(design_models, design_name) {
  purrr::map_dfr(
    seq_along(design_models),
    function(i) extract_model_atts(
      design_models[[i]],
      design_name,
      names(design_models)[i]
    )
  ) |> purrr::compact()
}

extract_all_atts <- function(regression_analysis) {
  purrr::map_dfr(
    names(regression_analysis),
    function(type) extract_design_atts(regression_analysis[[type]], type)
  )
}

format_att_table <- function(att_raw, clean_names) {
  att_raw |>
    mutate(
      outcome = clean_names[outcome],
      stars = case_when(
        p.value < 0.01 ~ "***",
        p.value < 0.05 ~ "**",
        p.value < 0.1  ~ "*",
        TRUE           ~ ""
      ),
      estimate_fmt = paste0(round(estimate, 3), stars),
      se_fmt       = paste0("(", round(std.error, 3), ")")
    ) |>
    select(outcome, design, estimate_fmt, se_fmt) |>
    pivot_longer(c(estimate_fmt, se_fmt),
                 names_to  = "stat",
                 values_to = "value") |>
    mutate(stat = factor(stat, levels = c("estimate_fmt", "se_fmt"))) |>
    arrange(outcome, stat) |>
    pivot_wider(names_from = design, values_from = value) |>
    mutate(outcome = ifelse(as.integer(stat) == 2, "", outcome)) |>
    select(-stat)
}

export_att_table <- function(att_table_wide, file) {
  att_table_wide |>
    kableExtra::kbl(
      format     = "latex",
      booktabs   = TRUE,
      escape     = FALSE,
      align      = "lcc",
      col.names  = c("Outcome", "Treatment ($<$1km)", "Spillover (2-5km)"),
      caption    = "Average Treatment Effects by Design",
      label      = "tab:att_results"
    ) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position")
    ) |>
    kableExtra::save_kable(paste0("outputs\\tables\\",file))
  
  file
}

format_continuous_table <- function(att_raw_continuous, clean_names = var_dict) {
  att_raw_continuous |>
    filter(design == "continuous") |>
    mutate(post = !grepl("-", term)) |>
    group_by(outcome) |>
    summarise(
      avg_post_estimate = mean(estimate[post]),
      avg_post_se       = mean(std.error[post]),
      avg_pre_estimate  = mean(estimate[!post]),
      avg_pre_se        = mean(std.error[!post])
    ) |>
    mutate(
      outcome   = clean_names[outcome],
      post_fmt  = paste0(round(avg_post_estimate, 3),
                         case_when(avg_post_estimate/avg_post_se > 2.576 ~ "***",
                                   avg_post_estimate/avg_post_se > 1.960 ~ "**",
                                   avg_post_estimate/avg_post_se > 1.645 ~ "*",
                                   TRUE ~ "")),
      post_se   = paste0("(", round(avg_post_se, 3), ")"),
      pre_fmt   = paste0(round(avg_pre_estimate, 3),
                         case_when(avg_pre_estimate/avg_pre_se > 2.576 ~ "***",
                                   avg_pre_estimate/avg_pre_se > 1.960 ~ "**",
                                   avg_pre_estimate/avg_pre_se > 1.645 ~ "*",
                                   TRUE ~ ""))  ,
      pre_se    = paste0("(", round(avg_pre_se, 3), ")")
    ) |>
    select(outcome, pre_fmt, pre_se, post_fmt, post_se) |>
    pivot_longer(c(pre_fmt, pre_se, post_fmt, post_se),
                 names_to  = "stat",
                 values_to = "value") |>
    mutate(
      period = ifelse(grepl("pre", stat), "Pre-Period", "Post-Period"),
      type   = ifelse(grepl("fmt", stat), "estimate", "se")
    ) |>
    select(-stat) |>
    pivot_wider(names_from = period, values_from = value) |>
    arrange(outcome, type) |>
    mutate(outcome = ifelse(type == "se", "", outcome)) |>
    select(-type)
}

export_continuous_table <- function(att_continuous_wide, file) {
  att_continuous_wide |>
    kableExtra::kbl(
      format     = "latex",
      booktabs   = TRUE,
      escape     = FALSE,
      align      = "lcc",
      col.names  = c("Outcome", "Pre-Period Avg", "Post-Period Avg"),
      caption    = "Continuous Treatment: Avg Log Distance Coefficients",
      label      = "tab:continuous_results"
    ) |>
    kableExtra::kable_styling(latex_options = "hold_position") |>
    kableExtra::save_kable(paste0("outputs\\", file))
  
  file
}

# Null coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------ Outcome groupings ------

OUTCOME_GROUPS <- list(
  residential = c(
    "Log Population",
    "Log Median Household Income",
    "Log Housing Units",
    "Log Median Rent"
  ),
  labor = c(
    "Log Total Jobs",
    "Log Jobs per 100 Residents",
    "Log Avg Annual Wage",
    "Log Establishment Entry"
  ),
  composition = c(
    "Share Small Establishments",
    "Share Mid Establishments",
    "Share Large Establishments"
  )
)

OUTCOME_LABELS <- c(
  "lhs: lpop"               = "Log Population",
  "lhs: lmhi"               = "Log Median Household Income",
  "lhs: lhousing"           = "Log Housing Units",
  "lhs: lmedian_rent"       = "Log Median Rent",
  "lhs: ljobs"              = "Log Total Jobs",
  "lhs: ljobs_per_100"      = "Log Jobs per 100 Residents",
  "lhs: lwage"              = "Log Avg Annual Wage",
  "lhs: lestabs_entry"      = "Log Establishment Entry",
  "lhs: estabs_small_share" = "Share Small Establishments",
  "lhs: estabs_mid_share"   = "Share Mid Establishments",
  "lhs: estabs_large_share" = "Share Large Establishments"
)

# ------ Extraction functions ------

extract_event_study_coefs <- function(regression_analysis, design = "continuous") {
  model_list <- as.list(regression_analysis[[design]])
  
  result <- purrr::map_dfr(seq_along(model_list), function(i) {
    model   <- model_list[[i]]
    outcome <- names(model_list)[i]
    coefs   <- summary(model)$coeftable
    
    tibble::tibble(
      term      = rownames(coefs),
      estimate  = coefs[, "Estimate"],
      std.error = coefs[, "Std. Error"],
      p.value   = coefs[, "Pr(>|t|)"],
      outcome   = outcome
    )
  }) |>
    dplyr::mutate(
      year_rel      = as.numeric(stringr::str_extract(term, "-?\\d+")),
      post          = year_rel >= 0,
      ci_low        = estimate - 1.96 * std.error,
      ci_high       = estimate + 1.96 * std.error,
      outcome_label = dplyr::coalesce(OUTCOME_LABELS[outcome], outcome) |>
        stringr::str_remove("^lhs: ")
    ) |>
    dplyr::filter(!is.na(year_rel))
  
  # Add k=-1 reference rows separately using the result object
  ref_rows <- result |>
    dplyr::distinct(outcome, outcome_label) |>
    dplyr::mutate(
      term      = "ref",
      estimate  = 0,
      std.error = 0,
      p.value   = NA_real_,
      year_rel  = -1,
      post      = FALSE,
      ci_low    = 0,
      ci_high   = 0
    )
  
  dplyr::bind_rows(result, ref_rows)
}

extract_sunab_coefs <- function(regression_analysis, design = "treatment") {
  # Same as event study but uses sunab summary
  model_list <- as.list(regression_analysis[[design]])
  
  result <- purrr::map_dfr(seq_along(model_list), function(i) {
    model   <- model_list[[i]]
    outcome <- names(model_list)[i]
    coefs   <- summary(model)$coeftable
    
    tibble::tibble(
      term      = rownames(coefs),
      estimate  = coefs[, "Estimate"],
      std.error = coefs[, "Std. Error"],
      p.value   = coefs[, "Pr(>|t|)"],
      outcome   = outcome
    )
  }) |>
    dplyr::mutate(
      year_rel      = as.numeric(stringr::str_extract(term, "-?\\d+")),
      post          = year_rel >= 0,
      ci_low        = estimate - 1.96 * std.error,
      ci_high       = estimate + 1.96 * std.error,
      outcome_label = dplyr::coalesce(OUTCOME_LABELS[outcome], outcome) |>
        stringr::str_remove("^lhs: ")
    ) |>
    dplyr::filter(!is.na(year_rel))
  
  # Add k=-1 reference rows separately using the result object
  ref_rows <- result |>
    dplyr::distinct(outcome, outcome_label) |>
    dplyr::mutate(
      term      = "ref",
      estimate  = 0,
      std.error = 0,
      p.value   = NA_real_,
      year_rel  = -1,
      post      = FALSE,
      ci_low    = 0,
      ci_high   = 0
    )
  
  dplyr::bind_rows(result, ref_rows)
}

# ------ Single-design plot ------

make_event_study_panel <- function(coef_df, outcomes, title,
                                   ncol = 2, ci_limit = NULL,
                                   color = "#2c7bb6",
                                   caption = "Sun-Abraham estimator: sunab(cohort_year, year) | GISJOIN + year.") {
  
  plot_df <- coef_df |>
    dplyr::filter(outcome_label %in% outcomes) |>
    dplyr::mutate(
      outcome_label = factor(outcome_label, levels = outcomes)
    )
  
  if (!is.null(ci_limit)) {
    plot_df <- plot_df |>
      dplyr::mutate(
        ci_low  = pmax(ci_low,  -ci_limit),
        ci_high = pmin(ci_high,  ci_limit),
        clamped = ci_low == -ci_limit | ci_high == ci_limit
      )
  } else {
    plot_df <- plot_df |> dplyr::mutate(clamped = FALSE)
  }
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = year_rel, y = estimate)
  ) +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed",
                        color = "grey40", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey60", linewidth = 0.4) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.15, fill = color
    ) +
    ggplot2::geom_line(color = color, linewidth = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(shape = post),
      size = 2, color = color
    ) +
    ggplot2::scale_shape_manual(
      values = c("FALSE" = 1, "TRUE" = 16),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(-7, 7, by = 2)) +
    ggplot2::facet_wrap(~outcome_label, ncol = ncol, scales = "free_y") +
    ggplot2::labs(
      x        = "Years Relative to Station Opening",
      y        = "ATT Estimate",
      title    = title,
      subtitle = if (!is.null(ci_limit))
        paste0("95% CIs clamped to \u00b1", ci_limit,
               ". Open points = pre-period.")
      else
        "Shaded area = 95% CI. Open points = pre-period.",
      caption  = caption
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(size = 8)
    )
}

# ------ Overlay plot (treatment vs spillover vs continuous) ------

make_event_study_overlay <- function(coef_list, outcomes, title,
                                     ncol = 2, ci_limit = NULL) {
  # coef_list: named list of coef_dfs e.g.
  # list("Treatment (<1km)" = coefs_treatment,
  #      "Spillover (2-5km)" = coefs_spillover)
  # or add continuous too
  
  group_colors <- c(
    "Treatment (<1km)"   = "#2c7bb6",
    "Spillover (2-5km)"  = "#d7191c",
    "Continuous DiD"     = "#1a9641"
  )
  
  combined <- purrr::imap_dfr(coef_list, function(df, grp) {
    df |> dplyr::mutate(group = grp)
  }) |>
    dplyr::filter(outcome_label %in% outcomes) |>
    dplyr::mutate(
      outcome_label = factor(outcome_label, levels = outcomes),
      group         = factor(group, levels = names(coef_list))
    )
  
  if (!is.null(ci_limit)) {
    combined <- combined |>
      dplyr::mutate(
        ci_low  = pmax(ci_low,  -ci_limit),
        ci_high = pmin(ci_high,  ci_limit)
      )
  }
  
  ggplot2::ggplot(
    combined,
    ggplot2::aes(x = year_rel, y = estimate,
                 color = group, fill = group)
  ) +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed",
                        color = "grey40", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey60", linewidth = 0.4) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.08, color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(shape = post), size = 1.5
    ) +
    ggplot2::scale_color_manual(
      values = group_colors, name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors, name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c("FALSE" = 1, "TRUE" = 16),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(-7, 7, by = 2)) +
    ggplot2::facet_wrap(~outcome_label, ncol = ncol, scales = "free_y") +
    ggplot2::labs(
      x        = "Years Relative to Station Opening",
      y        = "Estimate",
      title    = title,
      subtitle = "Shaded area = 95% CI. Open points = pre-period.",
      caption  = "Blue = Treatment (<1km). Red = Spillover (2-5km). Green = Continuous DiD."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(size = 8)
    )
}

# Legacy wrappers for backward compatibility
make_event_study_plot <- function(coef_df, outcomes_to_plot = NULL,
                                  ncol = 2, title = NULL) {
  make_event_study_panel(
    coef_df,
    outcomes  = outcomes_to_plot %||% unique(coef_df$outcome_label),
    title     = title %||% "Event Study: RTD Station Opening Effects",
    ncol      = ncol,
    color     = "#2c7bb6",
    caption   = "Continuous DiD: i(ref_year, log_dist, ref = -1) | GISJOIN + year."
  )
}

make_sunab_plot <- function(coef_df, outcomes_to_plot = NULL,
                            ncol = 2, title = NULL) {
  make_event_study_panel(
    coef_df,
    outcomes = outcomes_to_plot %||% unique(coef_df$outcome_label),
    title    = title %||% "Event Study: Sun-Abraham ATT Estimates",
    ncol     = ncol,
    color    = "#d7191c"
  )
}

# ------ FCMA Regression Setup -------

# ------ FCMA specs ------

fcma_spec_single <- function(outcome, delta_var = "delta_ln_FCMA",
                             delta_cf_var = "delta_ln_FCMA_cf") {
  
  dep_var <- switch(outcome,
                    jobs         = "ljobs",
                    jobs_per_100 = "ljobs_per_100",
                    wage         = "lwage",
                    estab        = "lestabs_entry"
  )
  
  list(
    base = list(
      as.formula(paste(
        dep_var, "~", delta_var, ":post +",
        delta_cf_var, ":post_2020 | GISJOIN + year"
      ))
    ) |> setNames(outcome),
    treatment = list(
      as.formula(paste(
        dep_var, "~", delta_var, ":post +",
        delta_cf_var, ":post_2020 +",
        delta_var, ":post:treat +",
        delta_cf_var, ":post_2020:treat",
        "| GISJOIN + year"
      ))
    ) |> setNames(outcome),
    high_rent = list(
      as.formula(paste(
        dep_var, "~", delta_var, ":post +",
        delta_cf_var, ":post_2020 +",
        delta_var, ":post:high_rent_int +",
        delta_cf_var, ":post_2020:high_rent_int",
        "| GISJOIN + year"
      ))
    ) |> setNames(outcome)
  )
}

rcma_spec_single <- function(outcome, delta_var = "delta_ln_RCMA",
                             delta_cf_var = "delta_ln_RCMA_cf") {
  
  dep_var <- switch(outcome,
                    pop     = "lpop",
                    income  = "lmhi",
                    housing = "lhousing",
                    rent    = "lmedian_rent"
  )
  
  list(
    base = list(
      as.formula(paste(
        dep_var, "~", delta_var, ":post +",
        delta_cf_var, ":post_2020 | GISJOIN + year"
      ))
    ) |> setNames(outcome),
    treatment = list(
      as.formula(paste(
        dep_var, "~", delta_var, ":post +",
        delta_cf_var, ":post_2020 +",
        delta_var, ":post:treat +",
        delta_cf_var, ":post_2020:treat",
        "| GISJOIN + year"
      ))
    ) |> setNames(outcome),
    high_rent = list(
      as.formula(paste(
        dep_var, "~", delta_var, ":post +",
        delta_cf_var, ":post_2020 +",
        delta_var, ":post:high_rent_int +",
        delta_cf_var, ":post_2020:high_rent_int",
        "| GISJOIN + year"
      ))
    ) |> setNames(outcome)
  )
}

select_specs <- function(data, delta_var, delta_cf_var) {
  # Detect which single outcome this dataset is for
  outcome <- dplyr::case_when(
    "ljobs"         %in% names(data) & 
      !"ljobs_per_100" %in% names(data)    ~ "jobs",
    "ljobs_per_100" %in% names(data) &
      !"lwage"         %in% names(data)    ~ "jobs_per_100",
    "lwage"         %in% names(data) &
      !"ljobs_per_100" %in% names(data)    ~ "wage",
    "lestabs_entry" %in% names(data)       ~ "estab",
    "lmedian_rent"  %in% names(data) &&
      sum(!is.na(data$lmedian_rent)) > 0   ~ "rent",
    "lpop"          %in% names(data)       ~ "pop",
    "lmhi"          %in% names(data)       ~ "income",
    "lhousing"      %in% names(data)       ~ "housing",
    TRUE                                   ~ stop("Unknown outcome dataset")
  )
  
  if (delta_var == "delta_ln_FCMA") {
    fcma_spec_single(outcome, delta_var, delta_cf_var)
  } else {
    rcma_spec_single(outcome, delta_var, delta_cf_var)
  }
}

run_cma_regressions <- function(data, delta_var = "delta_ln_FCMA",
                                delta_cf_var = "delta_ln_FCMA_cf",
                                outcome = NULL) {
  
  specs <- if (delta_var == "delta_ln_FCMA") {
    fcma_spec_single(outcome, delta_var, delta_cf_var)
  } else {
    rcma_spec_single(outcome, delta_var, delta_cf_var)
  }
  
  purrr::map(specs, function(spec) {
    purrr::map(spec, function(formula) {
      fixest::feols(formula, data = data, cluster = ~GISJOIN)
    })
  })
}

run_cma_bootstrap <- function(data, delta_var = "delta_ln_FCMA",
                              delta_cf_var  = "delta_ln_FCMA_cf",
                              outcome       = NULL,
                              n_boot        = 200,
                              seed          = 42,
                              parallel      = TRUE) {
  set.seed(seed)
  
  specs <- if (delta_var == "delta_ln_FCMA") {
    fcma_spec_single(outcome, delta_var, delta_cf_var)
  } else {
    rcma_spec_single(outcome, delta_var, delta_cf_var)
  }
  
  dt <- data |>
    dplyr::filter(!is.na(year)) |>
    dplyr::mutate(tract_id = substr(GISJOIN, 1, nchar(GISJOIN) - 1)) |>
    data.table::as.data.table()
  data.table::setkey(dt, tract_id)
  
  all_formula_vars <- unique(unlist(
    lapply(specs, function(spec)
      lapply(spec, all.vars)
    )
  ))
  req_cols <- unique(c(
    all_formula_vars,
    "GISJOIN", "year", "treat", "post", "post_2020", "high_rent_int"
  ))
  req_cols <- intersect(req_cols, names(dt))
  dt       <- dt[, .SD, .SDcols = c(req_cols, "tract_id")]
  
  tracts   <- unique(dt$tract_id)
  n_tracts <- length(tracts)
  
  # Set up parallel if requested
  if (parallel) {
    n_workers <- max(1, parallel::detectCores() - 1)
    future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential))
    map_fn <- function(n, f, ...) {
      furrr::future_map_dfr(
        n, f, ...,
        .options = furrr::furrr_options(
          seed     = TRUE,
          globals  = list(dt = dt, tracts = tracts, n_tracts = n_tracts),
          packages = c("data.table", "dplyr", "fixest")
        ),
        .progress = TRUE
      )
    }
  } else {
    map_fn <- function(n, f, ...) purrr::map_dfr(n, f, ...)
  }
  
  purrr::map(specs, function(spec) {
    purrr::map(spec, function(fml) {
      
      fit <- fixest::feols(fml, data = data, cluster = ~GISJOIN)
      
      boot_coefs <- map_fn(seq_len(n_boot), function(b) {
        on.exit(gc())
        boot_tracts <- sample(tracts, n_tracts, replace = TRUE)
        
        boot_data <- data.table::rbindlist(
          lapply(seq_along(boot_tracts), function(k) {
            rows <- dt[tract_id == boot_tracts[k]]
            rows[, GISJOIN := paste0(GISJOIN, "_", k)]
            rows
          })
        ) |>
          dplyr::select(-dplyr::any_of("tract_id")) |>
          dplyr::filter(dplyr::if_all(
            dplyr::any_of(c("ljobs", "ljobs_per_100", "lwage",
                            "lestabs_entry", "lhousing",
                            "lmedian_rent", "lpop", "lmhi")),
            ~ !is.infinite(.x)
          ))
        
        tryCatch({
          boot_fit <- fixest::feols(
            fml, data = boot_data,
            cluster = ~GISJOIN,
            warn = FALSE, notes = FALSE
          )
          result <- as.data.frame(t(coef(boot_fit)))
          rm(boot_data, boot_fit)
          result
        }, error = function(e) {
          rm(boot_data)
          NULL
        })
      })
      
      if (nrow(boot_coefs) == 0) {
        warning("All bootstrap iterations failed for: ", deparse(fml))
        return(list(
          fit        = fit,
          boot_se    = sqrt(diag(vcov(fit))),
          boot_ci    = t(confint(fit)),
          boot_coefs = boot_coefs
        ))
      }
      
      boot_se <- apply(boot_coefs, 2, sd,       na.rm = TRUE)
      boot_ci <- apply(boot_coefs, 2, quantile,
                       probs = c(0.025, 0.975), na.rm = TRUE)
      
      list(fit        = fit,
           boot_se    = boot_se,
           boot_ci    = boot_ci,
           boot_coefs = boot_coefs)
    })
  })
}

# Table output with bootstrap SEs substituted in
fcma_etable <- function(..., boot_results = NULL, title = NULL) {
  
  fits <- list(...)
  
  if (is.null(boot_results)) {
    # No bootstrap — use cluster SEs from fit directly
    do.call(fixest::etable, c(fits, list(title = title)))
    
  } else {
    # Substitute bootstrap SEs and recompute z-stats/p-values
    # boot_results should be a list of the same length as fits,
    # each element being the corresponding boot result list
    # e.g. list(fcma_boot_general$base$jobs, fcma_boot_low$base$jobs)
    
    stopifnot(length(fits) == length(boot_results))
    
    # Build se and ci lists for etable
    se_list <- purrr::map(boot_results, ~ .x$boot_se)
    ci_list <- purrr::map(boot_results, function(br) {
      # etable expects ci as a 2-row matrix: lower, upper
      br$boot_ci
    })
    
    do.call(fixest::etable, c(
      fits,
      list(
        se      = se_list,
        ci      = ci_list,
        title   = title,
        notes   = "Bootstrap standard errors clustered at tract level (500 replications)."
      )
    ))
  }
}

# ------ Coefficient plot helpers ------

save_plot <- function(plot_obj, filename, width = 10, height = 7) {
  dir.create("outputs/plots", recursive = TRUE, showWarnings = FALSE)
  path <- file.path("outputs/plots", filename)
  ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = 300)
  invisible(path)
}

# Updated make_coef_table_latex to handle single boot object
# ------ Coefficient labels ------

FCMA_COEF_LABELS <- c(
  "delta_ln_FCMA:post"               = "$\\Delta\\ln(\\text{FCMA}) \\times \\text{Post}$",
  "delta_ln_FCMA_cf:post_2020"  = "$\\Delta\\ln(\\text{FCMA})_{CF} \\times \\text{Post}_{2020}$",
  "delta_ln_FCMA:post:treat"         = "$\\Delta\\ln(\\text{FCMA}) \\times \\text{Post} \\times \\text{Treat}$",
  "delta_ln_FCMA_cf:post_2020:treat" = "$\\Delta\\ln(\\text{FCMA})_{CF} \\times \\text{Post}_{2020} \\times \\text{Treat}$",
  "delta_ln_FCMA:post:high_rent_int" = "$\\Delta\\ln(\\text{FCMA}) \\times \\text{Post} \\times \\text{High Rent}$",
  "delta_ln_FCMA_cf:post_2020:high_rent_int" = "$\\Delta\\ln(\\text{FCMA})_{CF} \\times \\text{Post}_{2020} \\times \\text{High Rent}$"
)

RCMA_COEF_LABELS <- c(
  "delta_ln_RCMA:post"               = "$\\Delta\\ln(\\text{RCMA}) \\times \\text{Post}$",
  "delta_ln_RCMA_cf:post_2020"  = "$\\Delta\\ln(\\text{RCMA})_{CF} \\times \\text{Post}_{2020}$",
  "delta_ln_RCMA:post:treat"         = "$\\Delta\\ln(\\text{RCMA}) \\times \\text{Post} \\times \\text{Treat}$",
  "delta_ln_RCMA_cf:post_2020:treat" = "$\\Delta\\ln(\\text{RCMA})_{CF} \\times \\text{Post}_{2020} \\times \\text{Treat}$",
  "delta_ln_RCMA:post:high_rent_int" = "$\\Delta\\ln(\\text{RCMA}) \\times \\text{Post} \\times \\text{High Rent}$",
  "delta_ln_RCMA_cf:post_2020:high_rent_int" = "$\\Delta\\ln(\\text{RCMA})_{CF} \\times \\text{Post}_{2020} \\times \\text{High Rent}$"
)

FCMA_OUTCOME_LABELS <- c(
  jobs         = "Log Total Jobs",
  jobs_per_100 = "Log Jobs per 100 Residents",
  wage         = "Log Avg Annual Wage",
  estab        = "Log Establishment Entry"
)

RCMA_OUTCOME_LABELS <- c(
  pop     = "Log Population",
  income  = "Log Med Household Income",
  housing = "Log Housing Units",
  rent    = "Log Median Rent"
)

COEF_LABELS_PLOT <- c(
  "delta_ln_FCMA:post"               = "\u0394ln(FCMA) \u00d7 Post",
  "delta_ln_FCMA_cf:post_2020"  = "\u0394ln(FCMA) CF \u00d7 Post\u2082\u2080\u2082\u2080",
  "delta_ln_FCMA:post:treat"         = "\u0394ln(FCMA) \u00d7 Post \u00d7 Treat",
  "delta_ln_FCMA_cf:post_2020:treat" = "\u0394ln(FCMA) CF \u00d7 Post\u2082\u2080\u2082\u2080 \u00d7 Treat",
  "delta_ln_FCMA:post:high_rent_int" = "\u0394ln(FCMA) \u00d7 Post \u00d7 High Rent",
  "delta_ln_FCMA_cf:post_2020:high_rent_int" = "\u0394ln(FCMA) CF \u00d7 Post\u2082\u2080\u2082\u2080 \u00d7 High Rent",
  "delta_ln_RCMA:post"               = "\u0394ln(RCMA) \u00d7 Post",
  "delta_ln_RCMA_cf:post_2020"  = "\u0394ln(RCMA) CF \u00d7 Post\u2082\u2080\u2082\u2080",
  "delta_ln_RCMA:post:treat"         = "\u0394ln(RCMA) \u00d7 Post \u00d7 Treat",
  "delta_ln_RCMA_cf:post_2020:treat" = "\u0394ln(RCMA) CF \u00d7 Post\u2082\u2080\u2082\u2080 \u00d7 Treat",
  "delta_ln_RCMA:post:high_rent_int" = "\u0394ln(RCMA) \u00d7 Post \u00d7 High Rent",
  "delta_ln_RCMA_cf:post_2020:high_rent_int" = "\u0394ln(RCMA) CF \u00d7 Post\u2082\u2080\u2082\u2080 \u00d7 High Rent"
)

# ------ Cell extractor ------

extract_cell <- function(fit_obj, coef_nm) {
  if (is.null(fit_obj)) return(list(est = "---", se = ""))
  
  est <- coef(fit_obj$fit)[coef_nm]
  if (is.na(est)) return(list(est = "---", se = ""))
  
  if (length(fit_obj$boot_se) > 0 && coef_nm %in% names(fit_obj$boot_se)) {
    se       <- fit_obj$boot_se[coef_nm]
    se_flag  <- ""
  } else {
    se       <- sqrt(diag(vcov(fit_obj$fit)))[coef_nm]
    se_flag  <- "*"   # flags clustered SE fallback
  }
  
  if (is.na(se) || se == 0) return(list(est = round(est, 4), se = ""))
  
  p_approx <- 2 * pnorm(-abs(est / se))
  stars <- dplyr::case_when(
    p_approx < 0.01 ~ "***",
    p_approx < 0.05 ~ "**",
    p_approx < 0.1  ~ "*",
    TRUE            ~ ""
  )
  
  list(
    est = paste0(round(est, 4), stars),
    se  = paste0("(", round(se, 4), ")", se_flag)
  )
}

# ------ Table builder ------
# boot_list: named list of boot objects e.g.
# list(jobs = fcma_boot_general[["jobs"]], estab = fcma_boot_general[["estab"]])

make_cma_table <- function(boot_list, spec = "base", measure = "fcma") {
  
  coef_labels    <- if (measure == "fcma") FCMA_COEF_LABELS else RCMA_COEF_LABELS
  outcome_labels <- if (measure == "fcma") FCMA_OUTCOME_LABELS else RCMA_OUTCOME_LABELS
  
  # Filter coef labels to spec
  if (spec == "base") {
    coef_labels <- coef_labels[!grepl("Treat|High Rent", coef_labels)]
  } else if (spec == "treatment") {
    coef_labels <- coef_labels[!grepl("High Rent", coef_labels)]
  } else if (spec == "high_rent") {
    coef_labels <- coef_labels[!grepl("Treat", coef_labels)]
  }
  
  # Flatten boot_list — each element may contain multiple outcomes
  # e.g. boot_list[["jobs"]][[spec]][["jobs_per_100"]]
  all_outcomes <- unlist(
    lapply(names(boot_list), function(grp) {
      outcomes <- names(boot_list[[grp]][[spec]])
      setNames(
        lapply(outcomes, function(o) boot_list[[grp]][[spec]][[o]]),
        outcomes
      )
    }),
    recursive = FALSE
  )
  
  purrr::map_dfr(
    intersect(names(outcome_labels), names(all_outcomes)),
    function(outcome) {
      fit_obj <- all_outcomes[[outcome]]
      purrr::map_dfr(names(coef_labels), function(coef_nm) {
        cell <- extract_cell(fit_obj, coef_nm)
        dplyr::bind_rows(
          tibble::tibble(
            outcome    = outcome_labels[outcome],
            coef_label = coef_labels[coef_nm],
            estimate   = cell$est,
            row_type   = "estimate"
          ),
          tibble::tibble(
            outcome    = "",
            coef_label = "",
            estimate   = cell$se,
            row_type   = "se"
          )
        )
      })
    }
  ) |>
    dplyr::select(outcome, coef_label, estimate)
}

# ------ Table export ------

export_cma_table <- function(table_df, file, spec = "base", measure = "fcma") {
  
  dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
  
  caption <- dplyr::case_when(
    spec == "base"      & measure == "fcma" ~
      "Effects of Firm Commuter Market Access on Labor Market Outcomes",
    spec == "treatment" & measure == "fcma" ~
      "FCMA Treatment Interaction Effects on Labor Market Outcomes",
    spec == "high_rent" & measure == "fcma" ~
      "FCMA High-Rent Heterogeneity: Labor Market Outcomes",
    spec == "base"      & measure == "rcma" ~
      "Effects of Resident Commuter Market Access on Residential Outcomes",
    spec == "treatment" & measure == "rcma" ~
      "RCMA Treatment Interaction Effects on Residential Outcomes",
    spec == "high_rent" & measure == "rcma" ~
      "RCMA High-Rent Heterogeneity: Residential Outcomes"
  )
  
  table_df |>
    kableExtra::kbl(
      format     = "latex",
      booktabs   = TRUE,
      escape     = FALSE,
      align      = "llc",
      col.names  = c("Outcome", "Coefficient", "Estimate"),
      caption    = caption,
      label      = paste0(measure, "_", spec)
    ) |>
    kableExtra::kable_styling(
      latex_options = c("hold_position", "scale_down")
    ) |>
    kableExtra::footnote(
      general = paste0(
        "Bootstrap SEs clustered at tract level in parentheses. ",
        "$^{*}$p$<$0.1, $^{**}$p$<$0.05, $^{***}$p$<$0.01. ",
        "Fixed effects: block group + year. ",
        "$^{(*)}$ indicates clustered SE used in place of bootstrap SE."
      ),
      escape = FALSE
    ) |>
    kableExtra::save_kable(file.path("outputs", "tables", file))
  
  file
}

# ------ Coefficient plot extractor ------

extract_coef_df <- function(boot_list, spec = "base", measure = "fcma") {
  
  outcome_labels <- if (measure == "fcma") FCMA_OUTCOME_LABELS else RCMA_OUTCOME_LABELS
  
  # Flatten boot_list
  all_outcomes <- unlist(
    lapply(names(boot_list), function(grp) {
      outcomes <- names(boot_list[[grp]][[spec]])
      setNames(
        lapply(outcomes, function(o) boot_list[[grp]][[spec]][[o]]),
        outcomes
      )
    }),
    recursive = FALSE
  )
  
  purrr::map_dfr(
    intersect(names(outcome_labels), names(all_outcomes)),
    function(outcome) {
      fit_obj <- all_outcomes[[outcome]]
      if (is.null(fit_obj)) return(NULL)
      
      if (length(fit_obj$boot_se) > 0) {
        se_vec  <- fit_obj$boot_se
        ci_low  <- fit_obj$boot_ci[1, ]
        ci_high <- fit_obj$boot_ci[2, ]
      } else {
        se_vec  <- sqrt(diag(vcov(fit_obj$fit)))
        ci_low  <- coef(fit_obj$fit) - 1.96 * se_vec
        ci_high <- coef(fit_obj$fit) + 1.96 * se_vec
      }
      
      tibble::tibble(
        coef_name     = names(coef(fit_obj$fit)),
        estimate      = as.numeric(coef(fit_obj$fit)),
        boot_se       = as.numeric(se_vec),
        ci_low        = as.numeric(ci_low),
        ci_high       = as.numeric(ci_high),
        outcome       = outcome,
        outcome_label = outcome_labels[outcome]
      )
    }
  ) |>
    dplyr::mutate(
      coef_label = dplyr::coalesce(
        COEF_LABELS_PLOT[coef_name],
        coef_name
      )
    )
}

# ------ Coefficient plot ------

make_coef_plot <- function(coef_df, title, coefs_to_plot = NULL,
                           ci_limit = NULL) {
  plot_df <- coef_df
  if (!is.null(coefs_to_plot)) {
    plot_df <- plot_df |> dplyr::filter(coef_label %in% coefs_to_plot)
  }
  if (!is.null(ci_limit)) {
    plot_df <- plot_df |>
      dplyr::mutate(
        ci_low  = pmax(ci_low,  -ci_limit),
        ci_high = pmin(ci_high,  ci_limit),
        clamped = ci_low == -ci_limit | ci_high == ci_limit
      )
  } else {
    plot_df <- plot_df |> dplyr::mutate(clamped = FALSE)
  }
  
  # Clean up any remaining raw coef names in coef_label
  plot_df <- plot_df |>
    dplyr::mutate(
      coef_label = dplyr::coalesce(
        COEF_LABELS_PLOT[coef_label],
        coef_label
      ),
      # Shorten outcome labels for x axis
      outcome_short = dplyr::case_when(
        outcome_label == "Log Total Jobs"            ~ "Total Jobs",
        outcome_label == "Log Jobs per 100 Residents"~ "Jobs/100",
        outcome_label == "Log Avg Weekly Wage"       ~ "Wage",
        outcome_label == "Log Avg Annual Wage"       ~ "Wage",
        outcome_label == "Log Establishment Entry"   ~ "Estab Entry",
        outcome_label == "Log Housing Units"         ~ "Housing",
        outcome_label == "Log Median Rent"           ~ "Rent",
        outcome_label == "Log Population"            ~ "Population",
        outcome_label == "Log Median HH Income"      ~ "MHI",
        TRUE ~ outcome_label
      )
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = outcome_short, y = estimate)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey50", linewidth = 0.4) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = ci_low, ymax = ci_high,
                   linetype = clamped),
      size = 0.5, linewidth = 0.7, color = "#2c7bb6"
    ) +
    ggplot2::scale_linetype_manual(
      values = c("FALSE" = "solid", "TRUE" = "dotted"),
      guide  = "none"
    ) +
    ggplot2::facet_wrap(
      ~coef_label,
      scales = "free_y",   # free y since scales differ substantially
      ncol   = 2
    ) +
    ggplot2::labs(
      x        = NULL,
      y        = "Coefficient Estimate",
      title    = title,
      subtitle = if (!is.null(ci_limit))
        paste0("95% CIs clamped to \u00b1", ci_limit,
               ". Dotted = clamped.")
      else
        "Point estimates with 95% bootstrap confidence intervals.",
      caption  = "Bootstrap SEs clustered at tract level. FEs: block group + year."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(size = 9, angle = 30,
                                               hjust = 1),
      panel.spacing    = ggplot2::unit(1.5, "lines")
    )
}

make_fcma_map <- function(bg_fcma, study_bgs_sf, var, title,
                          midpoint = 0, limit = NULL) {
  
  map_df <- study_bgs_sf |>
    dplyr::left_join(bg_fcma, by = "GISJOIN")
  
  if (is.null(limit)) {
    limit <- max(abs(map_df[[var]]), na.rm = TRUE)
  }
  
  ggplot2::ggplot(map_df) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), color = NA) +
    ggplot2::scale_fill_gradient2(
      low      = "red",
      mid      = "white",
      high     = "blue",
      midpoint = midpoint,
      limits   = c(-limit, limit),
      name     = "Δln(FCMA)",
      na.value = "grey80"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
}

make_treatment_map <- function(study_bgs_sf, regression_set) {
  
  map_df <- study_bgs_sf |>
    dplyr::left_join(
      regression_set |> dplyr::distinct(GISJOIN, group),
      by = "GISJOIN"
    )
  
  ggplot2::ggplot(map_df) +
    ggplot2::geom_sf(ggplot2::aes(fill = group), color = NA) +
    ggplot2::scale_fill_manual(
      values = c(
        "treatment" = "#d7191c",
        "exclude"   = "#fdae61",
        "spillover" = "#abd9e9",
        "control"   = "#2c7bb6"
      ),
      na.value = "grey80",
      name     = "Group",
      labels   = c(
        "treatment" = "Treatment (<1km)",
        "exclude"   = "Exclusion (1-2km)",
        "spillover" = "Spillover (2-5km)",
        "control"   = "Control (>5km)"
      )
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::labs(title = "Treatment Classification") +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
}
