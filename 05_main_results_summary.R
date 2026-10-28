# Create the main fully adjusted summary table and multiplicity-adjusted results.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("broom", "dplyr", "gt", "survey", "tibble"))
create_output_directories()
options(survey.lonely.psu = "adjust")

analysis_data <- readRDS(
  file.path("data", "nhms2023_dyslipidaemia_analysis.rds")
)
primary_models <- readRDS(
  file.path("results", "models", "primary-models.rds")
)
full_design <- make_survey_design(analysis_data)

summary_rows <- lapply(names(outcome_labels), function(outcome) {
  outcome_key <- outcome
  complete <- complete_case_mask(
    full_design,
    c(outcome_key, "madi_q", primary_covariates)
  )
  outcome_design <- full_design[complete, ]
  model <- primary_models[[outcome_key]]$models$model2

  prevalence_formula <- stats::as.formula(paste0("~", outcome_key))
  weighted_prevalence <- survey::svyby(
    prevalence_formula,
    ~madi_q,
    outcome_design,
    survey::svymean,
    na.rm = TRUE
  )

  q5_result <- tidy_madi_terms(model, outcome_key, "Model 2") |>
    dplyr::filter(grepl("Q5", term))

  tibble::tibble(
    outcome = outcome_labels[[outcome_key]],
    n = primary_models[[outcome_key]]$n,
    events = primary_models[[outcome_key]]$events,
    q1_prevalence = weighted_prevalence[[outcome_key]][1],
    q5_prevalence = weighted_prevalence[[outcome_key]][5],
    q5_pr = q5_result$estimate,
    q5_conf_low = q5_result$conf.low,
    q5_conf_high = q5_result$conf.high,
    q5_p_value = q5_result$p.value,
    global_madi_p = global_madi_p(model)
  )
})

summary_results <- dplyr::bind_rows(summary_rows) |>
  dplyr::mutate(
    global_fdr_q = stats::p.adjust(global_madi_p, method = "BH"),
    `Adjusted PR (95% CI)` = sprintf(
      "%.2f (%.2f, %.2f)", q5_pr, q5_conf_low, q5_conf_high
    )
  )

utils::write.csv(
  summary_results,
  file.path("results", "data", "main-results-summary.csv"),
  row.names = FALSE
)

summary_gt <- summary_results |>
  dplyr::select(
    Outcome = outcome,
    n,
    `Q1 prevalence` = q1_prevalence,
    `Q5 prevalence` = q5_prevalence,
    `Adjusted PR (95% CI)`,
    `Q5 p-value` = q5_p_value,
    `Global MADI p-value` = global_madi_p,
    `Global FDR q-value` = global_fdr_q
  ) |>
  gt::gt() |>
  gt::fmt_integer(columns = n) |>
  gt::fmt_percent(
    columns = c(`Q1 prevalence`, `Q5 prevalence`),
    decimals = 1
  ) |>
  gt::fmt_number(
    columns = c(`Q5 p-value`, `Global MADI p-value`, `Global FDR q-value`),
    decimals = 3
  ) |>
  gt::tab_header(
    title = gt::md("**Summary of fully adjusted findings**"),
    subtitle = "Model 2; Q1 (least deprived) is the reference"
  ) |>
  gt::tab_source_note(
    gt::md(
      paste0(
        "FDR q-values use Benjamini–Hochberg adjustment across the five ",
        "global MADI tests."
      )
    )
  ) |>
  gt::opt_row_striping()

save_gt_html(summary_gt, "main-results-summary.html")
message("Saved the fully adjusted summary of findings.")
