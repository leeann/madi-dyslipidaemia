# Fit the three prespecified sensitivity analyses for each outcome.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("survey"))
create_output_directories()
options(survey.lonely.psu = "adjust")

analysis_data <- readRDS(
  file.path("data", "nhms2023_dyslipidaemia_analysis.rds")
)
full_design <- make_survey_design(analysis_data)

age_spline_terms <- c(
  "madi_q", "splines::ns(age, df = 4)", "sex", "ethnicity_f", "marital_f",
  "education_f", "occupation_f", "income_f"
)
locality_terms <- c(model2_terms, "locality")
lifestyle_terms <- c(
  model2_terms,
  "smoking_f",
  "physical_activity_f",
  "splines::ns(fruit_veg_servings, df = 3)",
  "alcohol_f"
)

sensitivity_models <- lapply(names(outcome_labels), function(outcome) {
  age_complete <- complete_case_mask(
    full_design,
    c(outcome, "madi_q", "age", primary_covariates)
  )
  locality_complete <- complete_case_mask(
    full_design,
    c(outcome, "madi_q", primary_covariates, "locality")
  )
  lifestyle_complete <- complete_case_mask(
    full_design,
    c(outcome, "madi_q", primary_covariates, lifestyle_covariates)
  )

  age_design <- full_design[age_complete, ]
  locality_design <- full_design[locality_complete, ]
  lifestyle_design <- full_design[lifestyle_complete, ]

  list(
    n = c(
      age_spline = nrow(age_design$variables),
      locality = nrow(locality_design$variables),
      lifestyle = nrow(lifestyle_design$variables)
    ),
    models = list(
      age_base = fit_pr_model(age_design, outcome, model2_terms),
      age_spline = fit_pr_model(age_design, outcome, age_spline_terms),
      locality_base = fit_pr_model(locality_design, outcome, model2_terms),
      locality = fit_pr_model(locality_design, outcome, locality_terms),
      lifestyle_base = fit_pr_model(lifestyle_design, outcome, model2_terms),
      lifestyle = fit_pr_model(lifestyle_design, outcome, lifestyle_terms)
    )
  )
})
names(sensitivity_models) <- names(outcome_labels)

saveRDS(
  sensitivity_models,
  file.path("results", "models", "sensitivity-models.rds"),
  compress = "xz"
)

sensitivity_counts <- do.call(
  rbind,
  lapply(names(sensitivity_models), function(outcome) {
    data.frame(
      outcome = outcome_labels[[outcome]],
      analysis = names(sensitivity_models[[outcome]]$n),
      n = as.numeric(sensitivity_models[[outcome]]$n),
      row.names = NULL
    )
  })
)
utils::write.csv(
  sensitivity_counts,
  file.path("results", "data", "sensitivity-model-counts.csv"),
  row.names = FALSE
)

message("Saved sensitivity model objects for five outcomes.")
