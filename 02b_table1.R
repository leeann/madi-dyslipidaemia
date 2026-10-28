# Create the survey-weighted participant characteristics table.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("dplyr", "gt", "gtsummary", "survey"))
create_output_directories()
options(survey.lonely.psu = "adjust")

analysis_data <- readRDS(
  file.path("data", "nhms2023_dyslipidaemia_analysis.rds")
)
full_design <- make_survey_design(analysis_data)
module_design <- full_design[study_domain(full_design), ]

module_counts <- analysis_data |>
  dplyr::filter(adult %in% TRUE, hyperchol_module_respondent %in% TRUE) |>
  dplyr::count(madi_q, name = "n")

overall_n <- sum(module_counts$n)
madi_n <- stats::setNames(module_counts$n, as.character(module_counts$madi_q))

table_1 <- module_design |>
  gtsummary::tbl_svysummary(
    by = madi_q,
    include = c(
      age, sex, ethnicity_f, marital_f, education_f, occupation_f, income_f,
      hyperchol, undiagnosed_hyperchol, high_ldl, low_hdl, high_tg
    ),
    label = list(
      age = "Age, years",
      sex = "Sex",
      ethnicity_f = "Ethnicity",
      marital_f = "Marital status",
      education_f = "Highest education level",
      occupation_f = "Occupation",
      income_f = "Household-income quintile",
      hyperchol = "Total hypercholesterolaemia",
      undiagnosed_hyperchol = "Undiagnosed hypercholesterolaemia",
      high_ldl = "High LDL cholesterol",
      low_hdl = "Low HDL cholesterol",
      high_tg = "High triglycerides"
    ),
    statistic = list(
      age ~ "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{p}"
    ),
    type = list(
      c(hyperchol, undiagnosed_hyperchol, high_ldl, low_hdl, high_tg) ~
        "dichotomous"
    ),
    digits = list(age ~ 1, gtsummary::all_categorical() ~ 1),
    missing = "no"
  ) |>
  gtsummary::add_overall(last = FALSE) |>
  gtsummary::add_ci(
    include = gtsummary::all_categorical(),
    method = gtsummary::all_categorical() ~ "svyprop.logit",
    statistic = gtsummary::all_categorical() ~ "{conf.low}, {conf.high}",
    style_fun = gtsummary::all_categorical() ~
      gtsummary::label_style_number(digits = 1, scale = 100),
    pattern = "{stat} ({ci})"
  ) |>
  gtsummary::add_p(
    test = list(
      age ~ "svy.kruskal.test",
      gtsummary::all_categorical() ~ "svy.chisq.test"
    ),
    pvalue_fun = gtsummary::label_style_pvalue(digits = 3)
  ) |>
  gtsummary::bold_labels() |>
  gtsummary::modify_header(
    stat_0 = paste0(
      "**Overall**<br>n = ", format(overall_n, big.mark = ","),
      "<br>% (95% CI)"
    ),
    stat_1 = paste0(
      "**Q1 (least deprived)**<br>n = ",
      format(madi_n[["Q1 (least deprived)"]], big.mark = ","),
      "<br>% (95% CI)"
    ),
    stat_2 = paste0(
      "**Q2**<br>n = ", format(madi_n[["Q2"]], big.mark = ","),
      "<br>% (95% CI)"
    ),
    stat_3 = paste0(
      "**Q3**<br>n = ", format(madi_n[["Q3"]], big.mark = ","),
      "<br>% (95% CI)"
    ),
    stat_4 = paste0(
      "**Q4**<br>n = ", format(madi_n[["Q4"]], big.mark = ","),
      "<br>% (95% CI)"
    ),
    stat_5 = paste0(
      "**Q5 (most deprived)**<br>n = ",
      format(madi_n[["Q5 (most deprived)"]], big.mark = ","),
      "<br>% (95% CI)"
    ),
    p.value = "**p-value**"
  ) |>
  gtsummary::modify_caption(
    "**Table 1. Survey-weighted characteristics and dyslipidaemia outcomes by MADI quintile**"
  )

table_1_gt <- gtsummary::as_gt(table_1) |>
  gt::tab_source_note(
    gt::md(
      paste0(
        "Column sample sizes are unweighted; percentages and 95% CIs are ",
        "survey weighted. Age is reported as median (IQR)."
      )
    )
  ) |>
  gt::tab_options(
    table.font.size = gt::px(12),
    data_row.padding = gt::px(3),
    column_labels.padding = gt::px(4)
  )

save_gt_html(table_1_gt, "table-1.html")
message("Saved Table 1.")
