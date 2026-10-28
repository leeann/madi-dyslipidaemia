# Fit the three primary survey-weighted prevalence-ratio models per outcome.

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

primary_models <- lapply(names(outcome_labels), function(outcome) {
  complete <- complete_case_mask(
    full_design,
    c(outcome, "madi_q", primary_covariates)
  )
  outcome_design <- full_design[complete, ]

  list(
    n = nrow(outcome_design$variables),
    events = sum(outcome_design$variables[[outcome]] == 1L),
    models = fit_primary_model_set(outcome_design, outcome)
  )
})
names(primary_models) <- names(outcome_labels)

saveRDS(
  primary_models,
  file.path("results", "models", "primary-models.rds"),
  compress = "xz"
)

model_counts <- data.frame(
  outcome = unname(outcome_labels),
  n = vapply(primary_models, `[[`, numeric(1), "n"),
  events = vapply(primary_models, `[[`, numeric(1), "events")
)
utils::write.csv(
  model_counts,
  file.path("results", "data", "primary-model-counts.csv"),
  row.names = FALSE
)

message("Saved primary model objects for five outcomes.")
