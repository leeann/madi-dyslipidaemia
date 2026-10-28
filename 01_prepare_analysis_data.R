# Validate the linked NHMS source file and create a condensed analysis dataset.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("dplyr", "tibble"))
create_output_directories()

source_path <- locate_source_data()
message("Reading source data: ", source_path)
nhms_raw <- readRDS(source_path)

required_source_vars <- c(
  "ebid10", "state_st", "final_wt", "madi_quintile",
  "a22", "a25", "ethnicity", "marital_status", "education_level",
  "occupation", "hhincome_quintile", "strata", "pa_status", "csmk",
  "f04", "fv_dy", "n01", "n02_new2", "overall_hchol", "undx_hchol",
  "overall_hi_ldl", "overall_low_hdl", "overall_hi_tg"
)

missing_source_vars <- setdiff(required_source_vars, names(nhms_raw))
if (length(missing_source_vars) > 0L) {
  stop(
    "Required variables missing from source RDS: ",
    paste(missing_source_vars, collapse = ", "),
    call. = FALSE
  )
}

analysis_data <- nhms_raw |>
  dplyr::transmute(
    cluster_id = ebid10,
    design_stratum = state_st,
    final_wt = final_wt,
    locality = factor(
      strata,
      levels = c("1", "2"),
      labels = c("Urban", "Rural")
    ),
    adult = !is.na(a25) & a25 >= 18,
    madi_q = factor(
      as.character(madi_quintile),
      levels = paste0("Q", 1:5),
      labels = c(
        "Q1 (least deprived)", "Q2", "Q3", "Q4",
        "Q5 (most deprived)"
      )
    ),
    age = a25,
    age_band = cut(
      a25,
      breaks = c(18, 30, 40, 50, 60, 70, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = c("18–29", "30–39", "40–49", "50–59", "60–69", "≥70")
    ),
    sex = factor(a22, levels = 1:2, labels = c("Male", "Female")),
    ethnicity_f = factor(
      dplyr::case_when(
        ethnicity %in% 1:3 ~ ethnicity,
        ethnicity %in% c(4, 5) ~ 4,
        ethnicity == 6 ~ 5,
        TRUE ~ NA_real_
      ),
      levels = 1:5,
      labels = c(
        "Malay", "Chinese", "Indian", "Bumiputera Sabah/Sarawak", "Other"
      )
    ),
    marital_f = factor(
      marital_status,
      levels = 1:3,
      labels = c(
        "Never married", "Married/living with partner",
        "Separated/divorced/widowed"
      )
    ),
    education_f = factor(
      education_level,
      levels = 1:4,
      labels = c("No formal education", "Primary", "Secondary", "Tertiary")
    ),
    occupation_f = factor(
      occupation,
      levels = 1:7,
      labels = c(
        "Government employee", "Private employee", "Self-employed",
        "Homemaker/caregiver/unpaid", "Retired", "Student", "Not working"
      )
    ),
    income_f = factor(
      hhincome_quintile,
      levels = 1:5,
      labels = paste("Household income quintile", 1:5)
    ),
    smoking_f = factor(
      dplyr::case_when(
        csmk == 1 ~ "Current",
        csmk == 0 & f04 %in% c(1, 2) ~ "Former",
        csmk == 0 & f04 == 3 ~ "Never",
        TRUE ~ NA_character_
      ),
      levels = c("Never", "Former", "Current")
    ),
    physical_activity_f = factor(
      pa_status,
      levels = c(1, 2),
      labels = c("Active", "Inactive")
    ),
    fruit_veg_servings = as.numeric(fv_dy),
    alcohol_f = factor(
      dplyr::case_when(
        n02_new2 == 1 ~ "Current",
        n01 == 1 & n02_new2 == 2 ~ "Former",
        n01 == 2 ~ "Never",
        TRUE ~ NA_character_
      ),
      levels = c("Never", "Former", "Current")
    ),
    hyperchol_module_respondent =
      !is.na(overall_hchol) &
      !is.na(undx_hchol) &
      !is.na(overall_hi_ldl) &
      !is.na(overall_low_hdl) &
      !is.na(overall_hi_tg),
    hyperchol = dplyr::case_when(
      overall_hchol == 0 ~ 0L,
      overall_hchol == 1 ~ 1L,
      TRUE ~ NA_integer_
    ),
    undiagnosed_hyperchol = dplyr::case_when(
      undx_hchol == 0 ~ 0L,
      undx_hchol == 1 ~ 1L,
      TRUE ~ NA_integer_
    ),
    high_ldl = dplyr::case_when(
      overall_hi_ldl == 0 ~ 0L,
      overall_hi_ldl == 1 ~ 1L,
      TRUE ~ NA_integer_
    ),
    low_hdl = dplyr::case_when(
      overall_low_hdl == 0 ~ 0L,
      overall_low_hdl == 1 ~ 1L,
      TRUE ~ NA_integer_
    ),
    high_tg = dplyr::case_when(
      overall_hi_tg == 0 ~ 0L,
      overall_hi_tg == 1 ~ 1L,
      TRUE ~ NA_integer_
    )
  )

binary_vars <- names(outcome_labels)
binary_valid <- vapply(
  analysis_data[binary_vars],
  function(variable) all(stats::na.omit(variable) %in% c(0L, 1L)),
  FUN.VALUE = logical(1)
)

stopifnot(
  all(analysis_data$final_wt > 0, na.rm = TRUE),
  all(is.finite(analysis_data$final_wt[!is.na(analysis_data$final_wt)])),
  nlevels(analysis_data$madi_q) == 5L,
  all(binary_valid)
)

saveRDS(
  analysis_data,
  file.path("data", "nhms2023_dyslipidaemia_analysis.rds"),
  compress = "xz"
)

cohort_counts <- tibble::tibble(
  population = c(
    "All source records",
    "Adults aged 18 years or older",
    "Adult hypercholesterolaemia-module respondents"
  ),
  n = c(
    nrow(analysis_data),
    sum(analysis_data$adult %in% TRUE),
    sum(
      analysis_data$adult %in% TRUE &
        analysis_data$hyperchol_module_respondent %in% TRUE
    )
  )
)

utils::write.csv(
  cohort_counts,
  file.path("results", "data", "cohort-counts.csv"),
  row.names = FALSE
)

message(
  "Saved condensed analysis data with ",
  nrow(analysis_data), " records and ", ncol(analysis_data), " variables."
)
