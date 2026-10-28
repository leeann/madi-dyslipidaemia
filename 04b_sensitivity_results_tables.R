# Create paired tables for the age-spline, locality, and lifestyle analyses.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("gt", "gtsummary"))
create_output_directories()

sensitivity_models <- readRDS(
  file.path("results", "models", "sensitivity-models.rds")
)

age_tables <- lapply(names(outcome_labels), function(outcome) {
  merge_regression_tables(
    models = sensitivity_models[[outcome]]$models[c("age_base", "age_spline")],
    spanners = c("**Model 2: age bands**", "**Age-spline sensitivity**")
  )
})
names(age_tables) <- names(outcome_labels)

locality_tables <- lapply(names(outcome_labels), function(outcome) {
  merge_regression_tables(
    models = sensitivity_models[[outcome]]$models[c("locality_base", "locality")],
    spanners = c("**Model 2 on matched sample**", "**Model 2 + locality**")
  )
})
names(locality_tables) <- names(outcome_labels)

lifestyle_tables <- lapply(names(outcome_labels), function(outcome) {
  merge_regression_tables(
    models = sensitivity_models[[outcome]]$models[c("lifestyle_base", "lifestyle")],
    spanners = c("**Model 2 on matched sample**", "**Model 2 + lifestyle**")
  )
})
names(lifestyle_tables) <- names(outcome_labels)

age_table <- stack_outcome_tables(
  age_tables,
  "**Table S1. Sensitivity to non-linear continuous-age adjustment**"
)
locality_table <- stack_outcome_tables(
  locality_tables,
  "**Table S2. Sensitivity to additional urban/rural locality adjustment**"
)
lifestyle_table <- stack_outcome_tables(
  lifestyle_tables,
  "**Table S3. Sensitivity to additional lifestyle adjustment**"
)

save_gtsummary_html(age_table, "sensitivity-age-spline.html")
save_gtsummary_html(locality_table, "sensitivity-locality.html")
save_gtsummary_html(lifestyle_table, "sensitivity-lifestyle.html")

message("Saved three sensitivity-analysis tables.")
