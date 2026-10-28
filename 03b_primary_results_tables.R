# Create individual and combined tables for the primary models.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("dplyr", "gt", "gtsummary"))
create_output_directories()

primary_models <- readRDS(
  file.path("results", "models", "primary-models.rds")
)

model_spanners <- c(
  "**Unadjusted**",
  "**Model 1: demographic adjusted**",
  "**Model 2: fully adjusted**"
)

primary_tables <- lapply(names(outcome_labels), function(outcome) {
  table <- merge_regression_tables(
    models = primary_models[[outcome]]$models,
    spanners = model_spanners,
    caption = paste0("**", outcome_labels[[outcome]], "**")
  )

  save_gtsummary_html(
    table,
    paste0("primary-", outcome_slugs[[outcome]], ".html")
  )
  table
})
names(primary_tables) <- names(outcome_labels)

joint_primary_table <- stack_outcome_tables(
  primary_tables,
  "**Combined prevalence-ratio results for dyslipidaemia outcomes**"
)
save_gtsummary_html(joint_primary_table, "primary-combined.html")

joint_stars_table <- add_significance_stars(joint_primary_table) |>
  gtsummary::modify_caption(
    "**Combined prevalence-ratio results with significance stars**"
  )

joint_stars_gt <- gtsummary::as_gt(joint_stars_table) |>
  gt::tab_source_note(
    "Significance: * p < 0.1; ** p < 0.05; *** p < 0.01."
  )
save_gt_html(joint_stars_gt, "primary-combined-significance-stars.html")

message("Saved individual and combined primary regression tables.")
