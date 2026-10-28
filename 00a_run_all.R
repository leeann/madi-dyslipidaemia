# Run the complete dyslipidaemia analysis pipeline.
# Execute from the repository root:
#   Rscript R/00a_run_all.R

options(stringsAsFactors = FALSE)

source(file.path("R", "00b_functions.R"))

analysis_scripts <- c(
  "01_prepare_analysis_data.R",
  "02a_data_quality.R",
  "02b_table1.R",
  "03a_fit_primary_models.R",
  "03b_primary_results_tables.R",
  "03c_primary_forest_plots.R",
  "04a_fit_sensitivity_models.R",
  "04b_sensitivity_results_tables.R",
  "04c_sensitivity_summary.R",
  "05_main_results_summary.R",
  "99_session_info.R"
)

for (script in analysis_scripts) {
  message("\n--- Running R/", script, " ---")
  source(file.path("R", script), local = new.env(parent = globalenv()))
}

message("\nAnalysis pipeline completed successfully.")
