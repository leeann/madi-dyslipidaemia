# Compare sensitivity estimates with Model 2 on matched samples.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("broom", "dplyr", "gt", "survey", "tibble"))
create_output_directories()

sensitivity_models <- readRDS(
  file.path("results", "models", "sensitivity-models.rds")
)

comparison_specifications <- c(
  age_base = "Main Model 2",
  age_spline = "Age spline",
  locality = "+ Urban/rural locality",
  lifestyle_base = "Model 2 on lifestyle sample",
  lifestyle = "+ Lifestyle block"
)

comparison_rows <- lapply(names(outcome_labels), function(outcome) {
  outcome_key <- outcome
  dplyr::bind_rows(lapply(names(comparison_specifications), function(model_name) {
    model <- sensitivity_models[[outcome_key]]$models[[model_name]]
    sample_name <- if (model_name %in% c("lifestyle_base", "lifestyle")) {
      "lifestyle"
    } else if (model_name == "locality") {
      "locality"
    } else {
      "age_spline"
    }

    tidy_madi_terms(
      model,
      outcome,
      comparison_specifications[[model_name]]
    ) |>
      dplyr::filter(grepl("Q5", term)) |>
      dplyr::mutate(
        n = unname(sensitivity_models[[outcome_key]]$n[[sample_name]]),
        global_madi_p = global_madi_p(model)
      )
  }))
})

comparison_results <- dplyr::bind_rows(comparison_rows) |>
  dplyr::mutate(
    `Q5 vs Q1 PR (95% CI)` = sprintf(
      "%.2f (%.2f, %.2f)", estimate, conf.low, conf.high
    )
  )

comparison_csv <- comparison_results |>
  dplyr::transmute(
    outcome,
    specification,
    n,
    estimate,
    conf_low = conf.low,
    conf_high = conf.high,
    q5_p_value = p.value,
    global_madi_p
  )
utils::write.csv(
  comparison_csv,
  file.path("results", "data", "sensitivity-q5-comparison.csv"),
  row.names = FALSE
)

comparison_gt <- comparison_results |>
  dplyr::select(
    Outcome = outcome,
    Specification = specification,
    n,
    `Q5 vs Q1 PR (95% CI)`,
    `Q5 p-value` = p.value,
    `Global MADI p-value` = global_madi_p
  ) |>
  gt::gt(groupname_col = "Outcome") |>
  gt::fmt_integer(columns = n) |>
  gt::fmt_number(
    columns = c(`Q5 p-value`, `Global MADI p-value`),
    decimals = 3
  ) |>
  gt::tab_header(
    title = gt::md(
      "**Table S4. Q5-versus-Q1 estimates across sensitivity specifications**"
    )
  ) |>
  gt::tab_source_note(
    gt::md(
      paste0(
        "Q1 (least deprived) is the reference. The lifestyle-sample Model 2 ",
        "row separates complete-case selection from lifestyle adjustment."
      )
    )
  ) |>
  gt::opt_row_striping()
save_gt_html(comparison_gt, "sensitivity-q5-comparison.html")

calculate_percent_change <- function(base_model_name, sensitivity_model_name) {
  dplyr::bind_rows(lapply(names(outcome_labels), function(outcome) {
    base_results <- tidy_madi_terms(
      sensitivity_models[[outcome]]$models[[base_model_name]],
      outcome,
      "base"
    ) |>
      dplyr::select(outcome, term, base_estimate = estimate)

    sensitivity_results <- tidy_madi_terms(
      sensitivity_models[[outcome]]$models[[sensitivity_model_name]],
      outcome,
      "sensitivity"
    ) |>
      dplyr::select(outcome, term, sensitivity_estimate = estimate)

    dplyr::left_join(
      base_results,
      sensitivity_results,
      by = c("outcome", "term")
    ) |>
      dplyr::mutate(
        absolute_percent_change =
          100 * abs(sensitivity_estimate / base_estimate - 1)
      )
  }))
}

age_change <- calculate_percent_change("age_base", "age_spline")
locality_change <- calculate_percent_change("locality_base", "locality")
lifestyle_change <- calculate_percent_change("lifestyle_base", "lifestyle")

change_summary <- tibble::tibble(
  sensitivity_analysis = c("Age spline", "Urban/rural locality", "Lifestyle block"),
  median_absolute_percent_change = c(
    stats::median(age_change$absolute_percent_change),
    stats::median(locality_change$absolute_percent_change),
    stats::median(lifestyle_change$absolute_percent_change)
  ),
  maximum_absolute_percent_change = c(
    max(age_change$absolute_percent_change),
    max(locality_change$absolute_percent_change),
    max(lifestyle_change$absolute_percent_change)
  )
)
utils::write.csv(
  change_summary,
  file.path("results", "data", "sensitivity-change-summary.csv"),
  row.names = FALSE
)

message("Saved sensitivity comparison table and change metrics.")
