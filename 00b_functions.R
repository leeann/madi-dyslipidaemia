# Shared constants and helper functions for the dyslipidaemia study.

outcome_labels <- c(
  hyperchol = "Total hypercholesterolaemia",
  undiagnosed_hyperchol = "Undiagnosed hypercholesterolaemia",
  high_ldl = "High LDL cholesterol",
  low_hdl = "Low HDL cholesterol",
  high_tg = "High triglycerides"
)

outcome_slugs <- c(
  hyperchol = "total-hypercholesterolaemia",
  undiagnosed_hyperchol = "undiagnosed-hypercholesterolaemia",
  high_ldl = "high-ldl",
  low_hdl = "low-hdl",
  high_tg = "high-triglycerides"
)

primary_covariates <- c(
  "age_band", "sex", "ethnicity_f", "marital_f",
  "education_f", "occupation_f", "income_f"
)

model1_terms <- c("madi_q", "age_band", "sex", "ethnicity_f")

model2_terms <- c(
  "madi_q", "age_band", "sex", "ethnicity_f", "marital_f",
  "education_f", "occupation_f", "income_f"
)

lifestyle_covariates <- c(
  "smoking_f", "physical_activity_f", "fruit_veg_servings", "alcohol_f"
)

required_packages <- c(
  "broom", "dplyr", "ggplot2", "gt", "gtsummary", "survey", "tibble"
)

load_analysis_packages <- function(packages = required_packages) {
  missing_packages <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]

  if (length(missing_packages) > 0L) {
    stop(
      "Install the following packages before running the analysis: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(
    packages,
    function(package) {
      suppressPackageStartupMessages(
        library(package, character.only = TRUE)
      )
    }
  ))
}

create_output_directories <- function() {
  directories <- c(
    file.path("data"),
    file.path("results", ".sass-cache"),
    file.path("results", "data"),
    file.path("results", "figures"),
    file.path("results", "models"),
    file.path("results", "tables")
  )

  invisible(lapply(
    directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  options(
    sass.cache = normalizePath(
      file.path("results", ".sass-cache"),
      mustWork = TRUE
    )
  )
}

locate_source_data <- function() {
  environment_path <- Sys.getenv("NHMS2023_RDS", unset = "")
  candidate_paths <- c(
    environment_path,
    file.path("data-raw", "nhms2023_master_linked.rds")
  )
  candidate_paths <- candidate_paths[nzchar(candidate_paths)]
  available_paths <- candidate_paths[file.exists(candidate_paths)]

  if (length(available_paths) == 0L) {
    stop(
      paste0(
        "The restricted NHMS source RDS was not found. Set NHMS2023_RDS to ",
        "its full path or place nhms2023_master_linked.rds in data-raw/."
      ),
      call. = FALSE
    )
  }

  normalizePath(available_paths[[1]], mustWork = TRUE)
}

make_survey_design <- function(data) {
  survey::svydesign(
    ids = ~cluster_id,
    strata = ~design_stratum,
    weights = ~final_wt,
    data = data,
    nest = TRUE
  )
}

study_domain <- function(design) {
  design$variables$adult %in% TRUE &
    design$variables$hyperchol_module_respondent %in% TRUE
}

complete_case_mask <- function(design, variables) {
  study_domain(design) &
    stats::complete.cases(design$variables[, variables, drop = FALSE])
}

fit_pr_model <- function(design, outcome, terms) {
  survey::svyglm(
    stats::reformulate(termlabels = terms, response = outcome),
    design = design,
    family = stats::quasipoisson(link = "log")
  )
}

fit_primary_model_set <- function(design, outcome) {
  list(
    unadjusted = fit_pr_model(design, outcome, "madi_q"),
    model1 = fit_pr_model(design, outcome, model1_terms),
    model2 = fit_pr_model(design, outcome, model2_terms)
  )
}

regression_table <- function(model) {
  gtsummary::tbl_regression(
    model,
    exponentiate = TRUE,
    include = "madi_q"
  )
}

merge_regression_tables <- function(models, spanners, caption = NULL) {
  table <- gtsummary::tbl_merge(
    tbls = lapply(models, regression_table),
    tab_spanner = spanners
  ) |>
    gtsummary::modify_header(
      tidyselect::starts_with("estimate") ~ "**PR**",
      tidyselect::starts_with("conf.low") ~ "**95% CI**",
      tidyselect::starts_with("p.value") ~ "**p-value**"
    ) |>
    gtsummary::modify_abbreviation("PR = Prevalence Ratio")

  if (!is.null(caption)) {
    table <- gtsummary::modify_caption(table, caption)
  }

  table
}

stack_outcome_tables <- function(tables, caption) {
  gtsummary::tbl_stack(
    tbls = tables,
    group_header = unname(outcome_labels[names(tables)])
  ) |>
    gtsummary::modify_caption(caption)
}

add_significance_stars <- function(table) {
  estimate_columns <- grep("^estimate_", names(table$table_body), value = TRUE)

  for (estimate_column in estimate_columns) {
    suffix <- sub("^estimate_", "", estimate_column)
    p_column <- paste0("p.value_", suffix)

    if (!p_column %in% names(table$table_body)) {
      next
    }

    estimate <- table$table_body[[estimate_column]]
    p_value <- table$table_body[[p_column]]
    stars <- ifelse(
      is.na(p_value), "",
      ifelse(
        p_value < 0.01, "***",
        ifelse(p_value < 0.05, "**", ifelse(p_value < 0.10, "*", ""))
      )
    )

    table$table_body[[estimate_column]] <- ifelse(
      is.na(estimate),
      NA_character_,
      paste0(sprintf("%.2f", estimate), stars)
    )
  }

  table |>
    gtsummary::modify_table_styling(
      columns = tidyselect::starts_with("estimate"),
      fmt_fun = identity
    ) |>
    gtsummary::modify_table_styling(
      columns = tidyselect::starts_with("p.value"),
      hide = TRUE
    ) |>
    gtsummary::modify_header(
      tidyselect::starts_with("estimate") ~ "**PR**",
      tidyselect::starts_with("conf.low") ~ "**95% CI**"
    )
}

save_gtsummary_html <- function(table, filename) {
  gt::gtsave(
    data = gtsummary::as_gt(table),
    filename = file.path("results", "tables", filename)
  )
}

save_gt_html <- function(table, filename) {
  gt::gtsave(
    data = table,
    filename = file.path("results", "tables", filename)
  )
}

tidy_madi_terms <- function(model, outcome, specification) {
  outcome_variable <- outcome
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::filter(grepl("^madi_q", term)) |>
    dplyr::mutate(
      outcome = outcome_labels[[outcome_variable]],
      outcome_variable = outcome_variable,
      specification = specification
    )
}

global_madi_p <- function(model) {
  as.numeric(survey::regTermTest(model, ~madi_q)$p)
}

save_plot_outputs <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = file.path("results", "figures", paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path("results", "figures", paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}
