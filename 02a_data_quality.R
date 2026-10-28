# Produce structural checks, missingness, and sample-retention outputs.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("dplyr", "gt", "tibble"))
create_output_directories()

analysis_data <- readRDS(
  file.path("data", "nhms2023_dyslipidaemia_analysis.rds")
)

module_data <- analysis_data |>
  dplyr::filter(adult %in% TRUE, hyperchol_module_respondent %in% TRUE)

audit_vars <- c(
  "cluster_id", "design_stratum", "final_wt", "madi_q", "age", "age_band",
  "sex", "ethnicity_f", "marital_f", "education_f", "occupation_f",
  "income_f", "locality", "smoking_f", "physical_activity_f",
  "fruit_veg_servings", "alcohol_f", names(outcome_labels)
)

variable_labels <- c(
  cluster_id = "Survey cluster",
  design_stratum = "Survey stratum",
  final_wt = "Final weight",
  madi_q = "MADI quintile",
  age = "Age",
  age_band = "Age band",
  sex = "Sex",
  ethnicity_f = "Ethnicity",
  marital_f = "Marital status",
  education_f = "Highest education level",
  occupation_f = "Occupation",
  income_f = "Household-income quintile",
  locality = "Urban/rural locality",
  smoking_f = "Smoking status",
  physical_activity_f = "Physical activity status",
  fruit_veg_servings = "Fruit and vegetable servings per day",
  alcohol_f = "Alcohol-use status",
  outcome_labels
)

missingness <- tibble::tibble(
  variable = unname(variable_labels[audit_vars]),
  missing_n = vapply(
    module_data[audit_vars],
    function(variable) sum(is.na(variable)),
    FUN.VALUE = integer(1)
  )
) |>
  dplyr::mutate(missing_percent = 100 * missing_n / nrow(module_data))

structural_checks <- tibble::tibble(
  check = c(
    "Source records retained for survey design",
    "Adults aged 18 years or older",
    "Hypercholesterolaemia-module respondents",
    "Unique survey clusters",
    "Design strata",
    "MADI quintiles represented among adults",
    "Non-positive final weights",
    "Missing final weights",
    "Missing cluster IDs",
    "Missing design strata"
  ),
  result = c(
    nrow(analysis_data),
    sum(analysis_data$adult %in% TRUE),
    nrow(module_data),
    dplyr::n_distinct(analysis_data$cluster_id, na.rm = TRUE),
    dplyr::n_distinct(analysis_data$design_stratum, na.rm = TRUE),
    dplyr::n_distinct(module_data$madi_q, na.rm = TRUE),
    sum(analysis_data$final_wt <= 0, na.rm = TRUE),
    sum(is.na(analysis_data$final_wt)),
    sum(is.na(analysis_data$cluster_id)),
    sum(is.na(analysis_data$design_stratum))
  )
)

retention_rows <- lapply(names(outcome_labels), function(outcome) {
  required_variables <- c(outcome, "madi_q", primary_covariates)
  complete <- stats::complete.cases(
    module_data[, required_variables, drop = FALSE]
  )

  tibble::tibble(
    outcome = outcome_labels[[outcome]],
    outcome_observed_n = sum(!is.na(module_data[[outcome]])),
    complete_case_n = sum(complete),
    excluded_n = sum(!complete),
    retention_percent = 100 * mean(complete),
    events_n = sum(module_data[[outcome]][complete] == 1L)
  )
})

primary_retention <- dplyr::bind_rows(retention_rows)

primary_complete <- stats::complete.cases(
  module_data[, c("madi_q", primary_covariates), drop = FALSE]
)
lifestyle_complete <- primary_complete & stats::complete.cases(
  module_data[, lifestyle_covariates, drop = FALSE]
)

sensitivity_retention <- tibble::tibble(
  analysis = c("Primary Model 2 / age spline", "Urban/rural locality", "Lifestyle block"),
  complete_case_n = c(
    sum(primary_complete),
    sum(primary_complete & !is.na(module_data$locality)),
    sum(lifestyle_complete)
  )
) |>
  dplyr::mutate(
    excluded_n = nrow(module_data) - complete_case_n,
    retention_percent = 100 * complete_case_n / nrow(module_data)
  )

utils::write.csv(
  structural_checks,
  file.path("results", "data", "structural-checks.csv"),
  row.names = FALSE
)
utils::write.csv(
  missingness,
  file.path("results", "data", "missingness.csv"),
  row.names = FALSE
)
utils::write.csv(
  primary_retention,
  file.path("results", "data", "primary-sample-retention.csv"),
  row.names = FALSE
)
utils::write.csv(
  sensitivity_retention,
  file.path("results", "data", "sensitivity-sample-retention.csv"),
  row.names = FALSE
)

structural_gt <- structural_checks |>
  gt::gt() |>
  gt::fmt_integer(columns = result) |>
  gt::tab_header(title = gt::md("**Routine structural and range checks**"))

missingness_gt <- missingness |>
  gt::gt() |>
  gt::fmt_integer(columns = missing_n) |>
  gt::fmt_number(columns = missing_percent, decimals = 2) |>
  gt::cols_label(
    variable = "Variable",
    missing_n = "Missing n",
    missing_percent = "Missing %"
  ) |>
  gt::tab_header(
    title = gt::md("**Missingness among module respondents**")
  )

retention_gt <- primary_retention |>
  gt::gt() |>
  gt::fmt_integer(
    columns = c(outcome_observed_n, complete_case_n, excluded_n, events_n)
  ) |>
  gt::fmt_number(columns = retention_percent, decimals = 1) |>
  gt::cols_label(
    outcome = "Outcome",
    outcome_observed_n = "Outcome-observed n",
    complete_case_n = "Complete-case n",
    excluded_n = "Excluded n",
    retention_percent = "Retention %",
    events_n = "Events n"
  ) |>
  gt::tab_header(title = gt::md("**Primary analytic-sample retention**"))

save_gt_html(structural_gt, "data-quality-structural-checks.html")
save_gt_html(missingness_gt, "data-quality-missingness.html")
save_gt_html(retention_gt, "data-quality-primary-retention.html")

message("Saved data-quality outputs.")
